// define AOI of classification
var coastBuffer = ee.FeatureCollection('projects/earthobserver-cc-476010/assets/yellowsea_coasts_1km');

// TODO: clip to coast buffer, so background is removed for analysis
var loadedClassifiedTif = ee.Image('projects/earthobserver-cc-476010/assets/yellowsea_tidalflats_2017_2019_most');

var nativeProj = loadedClassifiedTif.projection();
var nativeScale = nativeProj.nominalScale();
var imageExtent = loadedClassifiedTif.geometry().bounds();

// Clip image extent to coast buffer to isolate tidal /analysis  region
var analysisArea = imageExtent.intersection(coastBuffer.geometry(), 30);

// Your Classified Image (0: Background/Land, 1: Tidal Flat, 2: Water)
//var myClassified2017 = loadedClassifiedTif.select(0).rename('class_2017');
var myClassified2017 = loadedClassifiedTif.select(0)
  .eq(1) // our tidal-flat classification code
  .rename('tidal_flat') // name same as murray
  .set({
    'system:time_start': ee.Date('2018-01-01').millis(),
    //'system:time_end': ee.Date('2019-12-31').millis()
    'source': 'User Classified'
  });

var myCollection = ee.ImageCollection([myClassified2017]);

// -------------------------------------------------------------------------
// combine ImageCollections (murray + ours)
// -------------------------------------------------------------------------
var murrayCollection = ee.ImageCollection('UQ/murray/Intertidal/v1_1/global_intertidal');
var murrayCollection = murrayCollection.filterBounds(imageExtent);
// 3. Map over the collection to manually reduce the region using a specific CRS
var murrayCollection = murrayCollection.map(function(img) {
  var tfBinary = img.clip(coastBuffer.geometry()).select('classification').eq(1).rename('tidal_flat');
  return tfBinary.copyProperties(img, ['system:time_start', 'system:time_end']);
});


var mergedCollection = murrayCollection.merge(myCollection)
  .sort('system:time_start');
  
  
//var murrayCollection = murrayCollection.map(function(image) {
//  var stats = image.reduceRegion({
//    reducer: ee.Reducer.mean(),
//    geometry: imageExtent,
//    scale: 250,
//    crs: nativeProj, 
//    bestEffort: false
//  });
  
  // Set the extracted classification code and system time as properties on a feature
//  return ee.Feature(null, {
//    'classification': stats.get('classification'),
//    'system:time_start': image.get('system:time_start')
//  });
//});

// -------------------------------------------------------------------------
// 4. Compute Quantitative Tidal Flat Area (Sq Km) per Timestep
// -------------------------------------------------------------------------

var pixelAreaKm2 = ee.Image.pixelArea().divide(1e6);

var timeSeriesFeatures = mergedCollection.map(function(img) {
  // Multiply binary tidal flat mask by pixel area
  var tfAreaImg = img.select('tidal_flat').multiply(pixelAreaKm2);
  
  var stats = tfAreaImg.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: analysisArea,
    scale: 30, // Reasonable scale for regional aggregate
    crs: nativeProj,
    maxPixels: 1e13,
    bestEffort: true
  });
  
  return ee.Feature(null, {
    'system:time_start': img.get('system:time_start'),
    'area_sqkm': stats.get('tidal_flat'),
    'source': img.get('source')
  });
});

// -------------------------------------------------------------------------
// 5. Build and Display Chart
// -------------------------------------------------------------------------
var chart = ui.Chart.feature.byFeature({
  features: timeSeriesFeatures,
  xProperty: 'system:time_start',
  yProperties: ['area_sqkm']
}).setOptions({
  title: 'Yellow Sea Intertidal Flat Extent Over Time',
  vAxis: {title: 'Tidal Flat Area (sq km)'},
  hAxis: {title: 'Year', format: 'YYYY'},
  lineWidth: 2,
  pointSize: 5,
  series: {
    0: {color: '#1b9e77'}
  }
});

// Center map & print chart
Map.centerObject(analysisArea, 8);
Map.addLayer(analysisArea, {color: 'red'}, 'AOI Analysis Area', false);
print(chart);
