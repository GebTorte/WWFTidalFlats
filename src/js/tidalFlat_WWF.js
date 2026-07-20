//////////////////////////////////////////// 
// Script based on:
// Murray N. J., Phinn S. R., DeWitt M., Ferrari R., Johnston R., Lyons M. B., Clinton N., Thau D. & 
// Fuller R. A. (2019) The global distribution and trajectory of tidal flats. *Nature*. 565, 222-225.

// Adapted for the Yellow Sea (2017-2025) to monitor tidal flats as a habitat for migratory birds.

// Developed in Google Earth Engine:
// https://code.earthengine.google.com
//////////////////////////////////////////// 

// 0. AOI and GLOBAL VARIABLES

var site = ee.Geometry.Polygon([117.5, 32, 127, 41], null, false); // rectangular bounding box of the yellow sea
var globOptions = { 
  versionID: '_SR',
  outFolder: 'SR',
  startDate: '2014-01-01',
  endDate: '2016-12-31',
  bandSelect: ['green', 'swir1', 'swir2', 'nir', 'red'],
  bands8: ['B3', 'B6', 'B7', 'B5', 'B4'],
  bands7: ['B2', 'B5', 'B7', 'B4', 'B3'], 
  maskAltitude: 100,  
  maskDepth: -100, 
  maskDistance: 50000,
  maskApplySRTM: false,
  parallelScale: 8,
  trainingValidationRatio: 0.0001,
  nTrees: 10, 
  outScale: 30, 
  conPixels: 100
};


// 1. FUNCTIONS
var landsatFunctions = {
  applyFMask: function(image) {
    // Mask out SHADOW, SNOW, and CLOUD classes. SR data.
    return image
      .updateMask(image.select('cfmask')
      .lt(2)); 
  },

  applyNDWI: function(image) {
    // apply NDWI to an image
    var ndwi = image.normalizedDifference(['green','nir']);
    return ndwi.select([0], ['ndwi']);
  },

  applyMNDWI: function(image) {
    // apply MNDWI to an image
    var mndwi = image.normalizedDifference(["green","swir1"]);
    return image.select([0], ['mndwi']);
  },

  applyAWEI: function(image) {
    // apply AWEI to an image
    var awei = image.expression("4*(b('green')-b('swir1'))-(0.25*b('nir')+2.75*b('swir2'))");
    return awei.select([0], ['awei']);
  },

  applyNDVI: function(image) {
    // apply NDVI to an image
    var ndvi = image.normalizedDifference(['nir','red']);
    return ndvi.select([0], ['ndvi']);
  },
};

var reducer = ee.Reducer.min()
    .combine(ee.Reducer.max(), '', true)
    .combine(ee.Reducer.stdDev(), '', true)
    .combine(ee.Reducer.median(), '', true)
    .combine(ee.Reducer.percentile([10, 25, 50, 75,90]), '', true)
    .combine(ee.Reducer.intervalMean(0, 10).setOutputs(['intMn0010']), '', true)
    .combine(ee.Reducer.intervalMean(10, 25).setOutputs(['intMn1025']), '', true)
    .combine(ee.Reducer.intervalMean(25, 50).setOutputs(['intMn2550']), '', true)
    .combine(ee.Reducer.intervalMean(50, 75).setOutputs(['intMn5075']), '', true)
    .combine(ee.Reducer.intervalMean(75, 90).setOutputs(['intMn7590']), '', true)
    .combine(ee.Reducer.intervalMean(90, 100).setOutputs(['intMn90100']), '', true)
    .combine(ee.Reducer.intervalMean(10, 90).setOutputs(['intMn1090']), '', true)
    .combine(ee.Reducer.intervalMean(25, 75).setOutputs(['intMn2575']), '', true);

// 2. DATA IMPORTS

// vectors
var globCoast = ee.FeatureCollection('ft:1Hsoe_WwULJ23Nuj1wikGzfH_WQMtpDWOR3XpWkHk');
var randomPointsPreComputed = ee.FeatureCollection('ft:1hVC5uIlWZQxtsapNsr5AzcLm1Vzo4I_DqjfNgNmN'); // Precomputed Training

// landsat
function generateLandsatCollection(){ //starting with landsat 7 because of the observation period
  var L7collection = ee.ImageCollection('LANDSAT/LE7_SR')
      .filterDate(globOptions.startDate,globOptions.endDate)
      .map(landsatFunctions.applyFMask)
      .select(globOptions.bands7, globOptions.bandSelect);
  var L8collection = ee.ImageCollection('LANDSAT/LC8_SR')
      .filterDate(globOptions.startDate, globOptions.endDate)
      .map(landsatFunctions.applyFMask)
      .select(globOptions.bands8, globOptions.bandSelect);
  var collectionFull = ee.ImageCollection(L4collection
      .merge(L5collection)
      .merge(L7collection)
      .merge(L8collection))
      .filterBounds(site)
      .filter(ee.Filter.intersects('.geo', globCoast.geometry(), null, null, 1000))
      .filterMetadata('WRS_ROW', 'less_than', 120); 
  return collectionFull;
}
var collection = generateLandsatCollection();
