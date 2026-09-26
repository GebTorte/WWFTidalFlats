// ============================================================
// Yellow Sea tidal flats — classify + export
// Notable choice: training labels are from 2014-2016 but re-sampled
// against this epoch's imagery in STEP 4.5. Fine for stable points,
// risky near shifting tidal edges — worth spot-checking.
// ============================================================

// ---- STEP 0: epoch + inputs ----
var epochStart = '2017-01-01';
var epochEnd   = '2019-12-31';
var epochLabel = '2017_2019';

var coastBuffer    = ee.FeatureCollection('projects/lstcalculation/assets/yellowsea_coasts_1km');
var trainingtable   = ee.FeatureCollection('projects/lstcalculation/assets/TrainingTable_YellowSea_2014_2016');

// ---- STEP 1: study region ----
var region = coastBuffer.geometry().simplify(30);
var regionBounds = region.bounds();
Map.centerObject(region, 8);

// ---- STEP 2: settings ----
var opt = {
  maskAltitude: 100,
  maskDepth: -100,
  parallelScale: 2,   // was 8; not needed once the reducer dropped percentiles
  nTrees: 10,
  conPixels: 100,
  maxCloud: 70
};

// Cloud mask + scale reflectance (Collection 2)
function maskAndScaleC2(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(62).eq(0);
  var optical = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  return image.addBands(optical, null, true).updateMask(qaMask);
}

// ---- STEP 3: load Landsat, filtered early ----
var bandSelect = ['green', 'swir1', 'swir2', 'nir', 'red'];
var bandsTM  = ['SR_B2', 'SR_B5', 'SR_B7', 'SR_B4', 'SR_B3']; // L7
var bandsOLI = ['SR_B3', 'SR_B6', 'SR_B7', 'SR_B5', 'SR_B4']; // L8/L9

function loadCollection(id, bandMap) {
  return ee.ImageCollection(id)
    .filterDate(epochStart, epochEnd)
    .filterBounds(regionBounds)
    .filter(ee.Filter.lte('CLOUD_COVER', opt.maxCloud))
    .select(bandMap.concat(['QA_PIXEL']))
    .map(maskAndScaleC2)
    .select(bandMap, bandSelect);
}

var collections = [
  loadCollection('LANDSAT/LE07/C02/T1_L2', bandsTM),
  loadCollection('LANDSAT/LC08/C02/T1_L2', bandsOLI)
];
if (new Date(epochEnd) > new Date('2021-10-31')) {
  collections.push(loadCollection('LANDSAT/LC09/C02/T1_L2', bandsOLI)); // L9 launched Oct 2021
}
var collection = collections.reduce(function(a, b) { return a.merge(b); });
print('Image count:', collection.size());

// ---- STEP 4: indices + covariates ----
function applyNDWI(img)  { return img.normalizedDifference(['green','nir']).rename('ndwi'); }
function applyMNDWI(img) { return img.normalizedDifference(['green','swir1']).rename('mndwi'); }
function applyNDVI(img)  { return img.normalizedDifference(['nir','red']).rename('ndvi'); }
function applyAWEI(img)  { return img.expression("4*(b('green')-b('swir1'))-(0.25*b('nir')+2.75*b('swir2'))").rename('awei'); }

// Detailed stats for the water indices (min/max/mean/median/stdDev — no percentiles, no variance)
var reducer = ee.Reducer.min()
  .combine(ee.Reducer.max(), '', true)
  .combine(ee.Reducer.mean(), '', true)
  .combine(ee.Reducer.median(), '', true)
  .combine(ee.Reducer.stdDev(), '', true);

var waterIndices = collection.map(function(img) {
  return ee.Image.cat([applyNDWI(img), applyMNDWI(img), applyAWEI(img)]);
});
var vegAndBands = collection.map(function(img) {
  return ee.Image.cat([applyNDVI(img), img.select(['nir', 'swir1'])]);
});

var covariates = {
  waterStats: waterIndices.reduce(reducer, opt.parallelScale),
  otherStats: vegAndBands.reduce(ee.Reducer.mean()), // just the mean for these — cheap, sufficient
  etopo: ee.Image('NOAA/NGDC/ETOPO1').select(['bedrock'], ['etopo']).resample('bicubic'),
  swOccurrence: ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select(['occurrence'], ['surfaceWater']).unmask()
};

var trainComposite = covariates.waterStats
  .addBands(covariates.otherStats)
  .addBands(covariates.etopo)
  .addBands(covariates.swOccurrence);

var bands = trainComposite.bandNames();
print('Covariate bands:', bands);

// ---- STEP 4.5: re-sample training points against this epoch's covariates ----
// Keeps CLASS labels, replaces old columns with fresh values at the same points.

// Shift class codes by 1 so 0 is never a real class — it stays reserved
// for masked/background pixels, which otherwise look identical to class 0.
var shiftedClasses = trainingtable.map(function(f) {
  return f.set('CLASS', ee.Number(f.get('CLASS')).add(1));
});

var updatedTrainingTable = trainComposite.sampleRegions({
  collection: shiftedClasses.select(['CLASS']),
  properties: ['CLASS'],
  scale: 30,
  tileScale: 4,
  geometries: true
});

// ---- STEP 5: elevation mask ----
var topoMask = covariates.etopo.gte(opt.maskDepth).and(covariates.etopo.lte(opt.maskAltitude));
var finalMask = topoMask.updateMask(topoMask).rename('datamask').byte();

// ---- STEP 6: spatially blocked train/validation split ----
var blockSizeDeg = 0.1;

var withBlock = updatedTrainingTable.map(function(f) {
  var coords = f.geometry().coordinates();
  var blockLon = ee.Number(coords.get(0)).divide(blockSizeDeg).floor();
  var blockLat = ee.Number(coords.get(1)).divide(blockSizeDeg).floor();
  return f.set('blockId', blockLon.multiply(100000).add(blockLat));
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

// ---- STEP 7: train + validate ----
var classifier = ee.Classifier.smileRandomForest({
    numberOfTrees: opt.nTrees,
    variablesPerSplit: null,
    bagFraction: 0.5,
    seed: 0
  })
  .train(trainSet, 'CLASS', bands)
  .setOutputMode('CLASSIFICATION');

var validated = validationSet.classify(classifier);
var errorMatrix = validated.errorMatrix('CLASS', 'classification');
print('Confusion matrix:', errorMatrix);
print('Overall accuracy:', errorMatrix.accuracy());
print('Kappa:', errorMatrix.kappa());

// Feature importance table
var impDict  = ee.Dictionary(classifier.explain().get('importance'));
var impTotal = ee.Number(impDict.values().reduce(ee.Reducer.sum()));
var impFC = ee.FeatureCollection(impDict.keys().map(function(k) {
  var v = ee.Number(impDict.get(k));
  return ee.Feature(null, {band: k, importance: v, share_pct: v.divide(impTotal).multiply(100)});
})).sort('importance', false);
print('Feature importance:', ui.Chart.feature.byFeature(impFC, 'band', ['importance', 'share_pct']).setChartType('Table'));

// ---- STEP 8: save classifier ----
Export.classifier.toAsset({
  classifier: classifier,
  description: 'export_classifier_' + epochLabel,
  assetId: 'projects/lstcalculation/assets/yellowsea_rf_classifier_' + epochLabel
});

// ---- STEP 9: classify + clean up ----
var classified = trainComposite.select(bands).classify(classifier);
var finalOut = classified.byte().mask(finalMask);
finalOut = finalOut.mask(finalOut.connectedPixelCount(opt.conPixels).gte(opt.conPixels));
// No class filter here: all 3 classes (now coded 1, 2, 3) stay in the output.
// Anything masked out above (elevation mask, denoising) exports as 0.
var finalClipped = finalOut.clip(region);

// ---- STEP 10: tile + export ----
var boundsCoords = regionBounds.coordinates().get(0).getInfo(); // one round trip, for the loop limits
var lons = boundsCoords.map(function(c){ return c[0]; });
var lats = boundsCoords.map(function(c){ return c[1]; });
var minLon = Math.min.apply(null, lons), maxLon = Math.max.apply(null, lons);
var minLat = Math.min.apply(null, lats), maxLat = Math.max.apply(null, lats);

var tileSize = 1.5;
var tileGeoms = [];
var tileFeatures = [];
for (var lon = minLon; lon < maxLon; lon += tileSize) {
  for (var lat = minLat; lat < maxLat; lat += tileSize) {
    var idx = tileGeoms.length;
    var geom = ee.Geometry.Rectangle(
      [lon, lat, Math.min(lon + tileSize, maxLon), Math.min(lat + tileSize, maxLat)], null, false);
    tileGeoms.push(geom);
    tileFeatures.push(ee.Feature(geom, {index: idx}));
  }
}

// One server call for all tiles, instead of one per tile
var validIdx = ee.FeatureCollection(tileFeatures).filterBounds(region).aggregate_array('index').getInfo();
print('Exporting tiles:', validIdx.length, 'of', tileGeoms.length);

validIdx.forEach(function(i) {
  Export.image.toDrive({
    image: finalClipped,
    description: 'yellowsea_tidalflats_' + epochLabel + '_tile' + i,
    region: tileGeoms[i],
    scale: 30,
    crs: 'EPSG:32652',
    maxPixels: 1e13
  });
});
