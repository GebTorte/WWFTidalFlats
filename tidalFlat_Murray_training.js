//////////////////////////////////////////// 
// See Murray et al. (Nature, 2019) for further information.
// Developed in Google Earth Engine:
// https://code.earthengine.google.com
// Cite:
// Murray N. J., Phinn S. R., DeWitt M., Ferrari R., Johnston R., Lyons M. B., Clinton N., Thau D. & 
// Fuller R. A. (2019) The global distribution and trajectory of tidal flats. *Nature*. 565, 222-225.
//////////////////////////////////////////// 

// ANGEPASST

// 0. Global Variables
var aoi = ee.FeatureCollection('projects/eagle-exercises/assets/Yellow_Sea_coastal_ii');
var trainingPoints = ee.FeatureCollection('projects/eagle-exercises/assets/Murray_tidalflat_2019_trainingset');

var bandSelect = ['green', 'swir1', 'swir2', 'nir', 'red'];
var bandsOLI   = ['SR_B3', 'SR_B6', 'SR_B7', 'SR_B5', 'SR_B4']; // Landsat 8 (OLI), Collection 2
var bandsTM    = ['SR_B2', 'SR_B5', 'SR_B7', 'SR_B4', 'SR_B3']; // Landsat 4/5/7 (TM/ETM+), Collection 2

var epochs = [
  // {start: '1984-01-01', end: '1986-12-31', label: '1984_1986'},
  // {start: '1987-01-01', end: '1989-12-31', label: '1987_1989'},
  // {start: '1990-01-01', end: '1992-12-31', label: '1990_1992'},
  // {start: '1993-01-01', end: '1995-12-31', label: '1993_1995'},
  // {start: '1996-01-01', end: '1998-12-31', label: '1996_1998'},
  // {start: '1999-01-01', end: '2001-12-31', label: '1999_2001'},
  // {start: '2002-01-01', end: '2004-12-31', label: '2002_2004'},
  // {start: '2005-01-01', end: '2007-12-31', label: '2005_2007'},
  // {start: '2008-01-01', end: '2010-12-31', label: '2008_2010'},
  // {start: '2011-01-01', end: '2013-12-31', label: '2011_2013'},
  {start: '2014-01-01', end: '2016-12-31', label: '2014_2016'}
];
var trainingEpochLabel = '2014_2016';

Map.centerObject(aoi, 8);
Map.setOptions('SATELLITE');
Map.addLayer(aoi.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'AOI');
Map.addLayer(trainingPoints, {color: 'yellow'}, 'Trainingspunkte (roh)', false);
 

// 1. Functions
function maskQA(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(parseInt('11111', 2)).eq(0); // Fill/DilatedCloud/Cirrus/Cloud/Shadow
  return image.updateMask(mask);
}

// Scaling Factors for Collection 2 Level 2 Surface Reflectance

function scaleSR(image) {
  var optical = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  return image.addBands(optical, null, true);
}

var landsatFunctions = {
  applyNDWI: function(image) {
    return image.normalizedDifference(['green', 'nir']).rename('ndwi');
  },
  applyMNDWI: function(image) {
    return image.normalizedDifference(['green', 'swir1']).rename('mndwi');
  },
  applyAWEI: function(image) {
    return image.expression("4*(b('green')-b('swir1'))-(0.25*b('nir')+2.75*b('swir2'))").rename('awei');
  },
  applyNDVI: function(image) {
    return image.normalizedDifference(['nir', 'red']).rename('ndvi');
  }
};

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

// 2. Data Imports & Processing

// images
function generateLandsatCollection(startDate, endDate) {
  var l4 = ee.ImageCollection('LANDSAT/LT04/C02/T1_L2')
      .filterDate(startDate, endDate).filterBounds(aoi)
      .map(maskQA).map(scaleSR).select(bandsTM, bandSelect);
  var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
      .filterDate(startDate, endDate).filterBounds(aoi)
      .map(maskQA).map(scaleSR).select(bandsTM, bandSelect);
  var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
      .filterDate(startDate, endDate).filterBounds(aoi)
      .map(maskQA).map(scaleSR).select(bandsTM, bandSelect);
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
      .filterDate(startDate, endDate).filterBounds(aoi)
      .map(maskQA).map(scaleSR).select(bandsOLI, bandSelect);
  return ee.ImageCollection(l4.merge(l5).merge(l7).merge(l8));
}
 
function buildFeatureStack(startDate, endDate) {
  var collection = generateLandsatCollection(startDate, endDate);
  
  var covariates = {
    aweiReduced: collection.map(landsatFunctions.applyAWEI).reduce(reducer),
    ndwiReduced: collection.map(landsatFunctions.applyNDWI).reduce(reducer),
    mndwiReduced: collection.map(landsatFunctions.applyMNDWI).reduce(reducer),
    ndvi: collection.map(landsatFunctions.applyNDVI)
        .reduce(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090'])),
    nirBand: collection.select(['nir'])
        .reduce(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090'])),
    swir1Band: collection.select(['swir1'])
        .reduce(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090'])),
    etopo: ee.Image('NOAA/NGDC/ETOPO1').select(['bedrock'], ['etopo']).resample('bicubic'),
    swOccurrence: ee.Image('JRC/GSW1_0/GlobalSurfaceWater')
        .select(['occurrence'], ['surfaceWater']).unmask()
  };
 
  return covariates.aweiReduced
      .addBands(covariates.ndwiReduced)
      .addBands(covariates.mndwiReduced)
      .addBands(covariates.ndvi)
      .addBands(covariates.nirBand)
      .addBands(covariates.swir1Band)
      .addBands(covariates.etopo)
      .addBands(covariates.swOccurrence)
      .clip(aoi);
}

// building Feature Stacks
var trainComposite;
epochs.forEach(function(epoch) {
  var stack = buildFeatureStack(epoch.start, epoch.end);
 
  if (epoch.label === trainingEpochLabel) {
    trainComposite = stack; 
  }
 
  Export.image.toDrive({
    image: stack,
    description: 'FeatureStack_YellowSea_' + epoch.label,
    folder: 'GEE_exports',
    fileNamePrefix: 'FeatureStack_YellowSea_' + epoch.label,
    region: aoi.geometry(),
    scale: 30,
    crs: 'EPSG:32651', 
    maxPixels: 1e13
  });
});

// Check
Map.addLayer(trainComposite.select('ndwi_median'),
    {min: -0.5, max: 0.5, palette: ['brown', 'white', 'blue']},
    'NDWI Median 2014-2016 (Check)', false);
print('Anzahl Feature-Baender im Stack:', trainComposite.bandNames());

// Training Table & Export
var sampledPoints = trainComposite.sampleRegions({
  collection: trainingPoints,
  properties: ['CLASS'],
  scale: 30,
  tileScale: 4,
  geometries: true
});
 
var validPoints = sampledPoints.filter(ee.Filter.notNull(trainComposite.bandNames()));
print('Trainingspunkte gesamt:', trainingPoints.size());
print('Punkte mit gueltigen Feature-Werten:', validPoints.size());
 
Map.addLayer(validPoints, {color: 'lime'}, 'Punkte mit extrahierten Features (Check)', false);
 
Export.table.toDrive({
  collection: validPoints,
  description: 'TrainingTable_YellowSea_2014_2016',
  folder: 'GEE_exports',
  fileFormat: 'CSV'
});
 