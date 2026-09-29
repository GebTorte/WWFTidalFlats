// ============================================================
// 0. Epoch — change dates for each epoch
// ============================================================
var epochStart = '2017-01-01';
var epochEnd   = '2019-12-31';
var epochLabel = '2017_2019';

// ============================================================
// 1. Region — coast buffer, simplified to cut vertex count for export
// ============================================================
var region = coastBuffer.geometry().simplify(30); // 30m ≈ pixel scale
Map.centerObject(region, 8);

// ============================================================
// 2. Options + Collection 2 cloud mask
// ============================================================
var globOptions = {
  maskAltitude: 100,
  maskDepth: -100,
  parallelScale: 8,
  nTrees: 10,
  conPixels: 100
};

function maskAndScaleC2(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(62).eq(0);
  var optical = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  return image.addBands(optical, null, true).updateMask(qaMask);
}

// ============================================================
// 3. Build the collection for this epoch, scoped to the region
// ============================================================
var bandSelect = ['green', 'swir1', 'swir2', 'nir', 'red'];
var bandsTM  = ['SR_B2', 'SR_B5', 'SR_B7', 'SR_B4', 'SR_B3']; // L5, L7
var bandsOLI = ['SR_B3', 'SR_B6', 'SR_B7', 'SR_B5', 'SR_B4']; // L8, L9

function loadCollection(id, bandMap) {
  return ee.ImageCollection(id)
    .filterBounds(region)
    .filterDate(epochStart, epochEnd)
    .map(maskAndScaleC2)
    .select(bandMap, bandSelect);
}

var L7 = loadCollection('LANDSAT/LE07/C02/T1_L2', bandsTM);
var L8 = loadCollection('LANDSAT/LC08/C02/T1_L2', bandsOLI);
var L9 = loadCollection('LANDSAT/LC09/C02/T1_L2', bandsOLI);
var collection = ee.ImageCollection(L7.merge(L8).merge(L9));
print('Image count for ' + epochLabel + ':', collection.size());

// ============================================================
// 4. Indices, reducer, covariates
// ============================================================
function applyNDWI(img)  { return img.normalizedDifference(['green','nir']).rename('ndwi'); }
function applyMNDWI(img) { return img.normalizedDifference(['green','swir1']).rename('mndwi'); }
function applyNDVI(img)  { return img.normalizedDifference(['nir','red']).rename('ndvi'); }
function applyAWEI(img)  { return img.expression("4*(b('green')-b('swir1'))-(0.25*b('nir')+2.75*b('swir2'))").rename('awei'); }

var reducer = ee.Reducer.min()
  .combine(ee.Reducer.max(), '', true)
  .combine(ee.Reducer.stdDev(), '', true)
  .combine(ee.Reducer.median(), '', true)
  .combine(ee.Reducer.percentile([10, 25, 50, 75, 90]), '', true)
  .combine(ee.Reducer.intervalMean(0, 10).setOutputs(['intMn0010']), '', true)
  .combine(ee.Reducer.intervalMean(10, 25).setOutputs(['intMn1025']), '', true)
  .combine(ee.Reducer.intervalMean(25, 50).setOutputs(['intMn2550']), '', true)
  .combine(ee.Reducer.intervalMean(50, 75).setOutputs(['intMn5075']), '', true)
  .combine(ee.Reducer.intervalMean(75, 90).setOutputs(['intMn7590']), '', true)
  .combine(ee.Reducer.intervalMean(90, 100).setOutputs(['intMn90100']), '', true)
  .combine(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090']), '', true)
  .combine(ee.Reducer.intervalMean(25, 75).setOutputs(['intMn2575']), '', true);

var covariates = {
  aweiReduced:  collection.map(applyAWEI).reduce(reducer, globOptions.parallelScale),
  ndwiReduced:  collection.map(applyNDWI).reduce(reducer, globOptions.parallelScale),
  mndwiReduced: collection.map(applyMNDWI).reduce(reducer, globOptions.parallelScale),
  ndvi:      collection.map(applyNDVI).reduce(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090'])),
  nirBand:   collection.select(['nir']).reduce(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090'])),
  swir1Band: collection.select(['swir1']).reduce(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090'])),
  etopo: ee.Image('NOAA/NGDC/ETOPO1').select(['bedrock'], ['etopo']).resample('bicubic'),
  swOccurrence: ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select(['occurrence'], ['surfaceWater']).unmask()
};

var trainComposite = covariates.aweiReduced
  .addBands(covariates.ndwiReduced)
  .addBands(covariates.mndwiReduced)
  .addBands(covariates.ndvi)
  .addBands(covariates.nirBand)
  .addBands(covariates.swir1Band)
  .addBands(covariates.etopo)
  .addBands(covariates.swOccurrence)
  .clip(region);

var bands = trainComposite.bandNames();

// ============================================================
// 5. Topo mask
// ============================================================
var topoMask = covariates.etopo.gte(globOptions.maskDepth).and(covariates.etopo.lte(globOptions.maskAltitude));
var finalMask = topoMask.updateMask(topoMask).rename('datamask').byte();

// ============================================================
// 6. Spatially blocked train/validation split
// ============================================================
var blockSizeDeg = 0.1; // should exceed the spatial autocorrelation range

var withBlock = trainingtable.map(function(f) {
  var coords = f.geometry().coordinates();
  var lon = ee.Number(coords.get(0));
  var lat = ee.Number(coords.get(1));
  var blockLon = lon.divide(blockSizeDeg).floor();
  var blockLat = lat.divide(blockSizeDeg).floor();
  var blockId = blockLon.multiply(100000).add(blockLat);
  return f.set('blockId', blockId);
});

var blockIds = ee.List(withBlock.distinct('blockId').aggregate_array('blockId'));
var blockFC = ee.FeatureCollection(blockIds.map(function(id) {
  return ee.Feature(null, {blockId: id});
})).randomColumn('random', 42);

var joined = ee.Join.saveFirst('blockMatch').apply(
  withBlock, blockFC, ee.Filter.equals({leftField: 'blockId', rightField: 'blockId'})
);
var withBlockRandom = joined.map(function(f) {
  return f.set('random', ee.Feature(f.get('blockMatch')).get('random'));
});

var trainSet = withBlockRandom.filter(ee.Filter.lt('random', 0.7));
var validationSet = withBlockRandom.filter(ee.Filter.gte('random', 0.7));

// ============================================================
// 7. Train + validate
// ============================================================
var classifier = ee.Classifier.smileRandomForest({
    numberOfTrees: globOptions.nTrees,
    variablesPerSplit: null,
    bagFraction: 0.5,
    seed: 0
  })
  .train(trainSet, 'CLASS', bands)
  .setOutputMode('CLASSIFICATION');

var validated = validationSet.classify(classifier);
var errorMatrix = validated.errorMatrix('CLASS', 'classification');
print('Validation confusion matrix:', errorMatrix);
print('Validation overall accuracy:', errorMatrix.accuracy());
print('Validation kappa:', errorMatrix.kappa());

// ============================================================
// 8. Save the trained classifier as an asset (run once per epoch, reusable after)
// ============================================================
Export.classifier.toAsset({
  classifier: classifier,
  description: 'export_classifier_' + epochLabel,
  assetId: 'users/rjones/yellowsea_rf_classifier_' + epochLabel // match your existing asset path style
});

// ============================================================
// 9. Classify + denoise
// ============================================================
var classified = trainComposite.select(bands).classify(classifier);
var finalOut = classified.byte().mask(finalMask);
finalOut = finalOut.mask(finalOut.connectedPixelCount(globOptions.conPixels).gte(globOptions.conPixels));
finalOut = finalOut.updateMask(finalOut.eq(2));

// ============================================================
// 10. Tile the region and export each tile separately
// ============================================================
var boundsCoords = region.bounds().coordinates().get(0).getInfo();
var lons = boundsCoords.map(function(c){ return c[0]; });
var lats = boundsCoords.map(function(c){ return c[1]; });
var minLon = Math.min.apply(null, lons), maxLon = Math.max.apply(null, lons);
var minLat = Math.min.apply(null, lats), maxLat = Math.max.apply(null, lats);

var tileSize = 1.5; // in degrees - decrease/increase if needed
var tiles = [];
for (var lon = minLon; lon < maxLon; lon += tileSize) {
  for (var lat = minLat; lat < maxLat; lat += tileSize) {
    tiles.push(ee.Geometry.Rectangle([lon, lat, Math.min(lon + tileSize, maxLon), Math.min(lat + tileSize, maxLat)]));
  }
}
print('Number of tiles:', tiles.length);

var validTiles = [];
tiles.forEach(function(tileGeom, i) {
  var tileRegion = tileGeom.intersection(region, ee.ErrorMargin(30));
  var area = tileRegion.area(ee.ErrorMargin(30)).getInfo();
  if (area > 0) {
    validTiles.push({region: tileRegion, index: i});
  } else {
    print('Skipping empty tile ' + i);
  }
});
print('Valid (non-empty) tiles:', validTiles.length);

validTiles.forEach(function(t) {
  var finalOutTile = finalOut.clip(t.region);
  Export.image.toDrive({
    image: finalOutTile,
    description: 'yellowsea_tidalflats_' + epochLabel + '_tile' + t.index,
    region: t.region,
    scale: 30,
    crs: 'EPSG:32652',
    maxPixels: 1e13
  });
});