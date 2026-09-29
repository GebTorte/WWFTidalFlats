// define AOI of classification
var coastBuffer = ee.FeatureCollection('projects/earthobserver-cc-476010/assets/korea_studyarea');

// TODO: clip to coast buffer, so background is removed for analysis
var loadedClassifiedTif_2017 = ee.Image('projects/earthobserver-cc-476010/assets/classified/2017_2019_koreaSA');
var loadedClassifiedTif_2020 = ee.Image('projects/earthobserver-cc-476010/assets/classified/2020_2022_koreaSA');
var loadedClassifiedTif_2023 = ee.Image('projects/earthobserver-cc-476010/assets/classified/2023_2025_koreaSA');

var nativeProj = loadedClassifiedTif_2017.projection();
var nativeScale = nativeProj.nominalScale();
var imageExtent = loadedClassifiedTif_2017.geometry().bounds();

// Clip image extent to coast buffer to isolate tidal /analysis  region
var analysisArea = imageExtent.intersection(coastBuffer.geometry(), 30);

// Your Classified Image (0: Background/Land, 1: Tidal Flat, 2: Water)
var prepareClassified = function(img) {
  return img.select(0)
  .eq(2) // our tidal-flat classification code
  .rename('tidal_flat').copyProperties(img, ['system:time_start', 'system:time_end']) // name same as murray
};


var myCollection = ee.ImageCollection([loadedClassifiedTif_2017, loadedClassifiedTif_2020, loadedClassifiedTif_2023]);
var myCollection = myCollection.map(function(img) {
  return prepareClassified(img);
});

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
  title: 'Korea AOI Intertidal Flat Extent Over Time',
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
