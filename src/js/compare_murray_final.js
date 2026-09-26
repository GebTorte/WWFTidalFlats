// =========================================================================
// 1. INPUT ASSETS & GEOMETRY
// =========================================================================
var loadedClassifiedImage = ee.Image('projects/earthobserver-cc-476010/assets/classified/yellowsea_tidalflats_2017_2019_tile3');

// Clean 4-corner bounding box to prevent geometry clipping errors
var imageExtent = loadedClassifiedImage.geometry().bounds();

// =========================================================================
// 2. DATASET PREPARATION & ALIGNMENT
// =========================================================================
// Murray et al. Intertidal Dataset (2014-2016 timestep)
// Band 'classification': 1 = Tidal Flat, 0 = Non-intertidal
var murray2014 = ee.ImageCollection('UQ/murray/Intertidal/v1_1/global_intertidal')
  .filterDate('2014-01-01', '2016-12-31')
  .first()
  .select('classification')
  .rename('class_2014');

// Your classified image (0: background, 1: tidal flat, 2: water, 3: land)
var myClassified2026 = loadedClassifiedImage.select(0).rename('class_2026');

// Combine both rasters clipped to image extent
var combined = ee.Image.cat([myClassified2026, murray2014]).clip(imageExtent);

// =========================================================================
// 3. TRANSITION ANALYSIS & VISUALIZATION LAYERS
// =========================================================================
var myTF = 1;     // Your Tidal Flat Class ID
var murrayTF = 1; // Murray Tidal Flat Class ID

// A. Stable Tidal Flat (Intertidal in both 2014 & 2026)
var tfStable = murray2014.eq(murrayTF).and(myClassified2026.eq(myTF));

// B. Tidal Flat Loss (2014 Intertidal -> 2026 Non-Tidal Flat)
var tfLoss = murray2014.eq(murrayTF).and(myClassified2026.neq(myTF));

// C. Tidal Flat Gain (2014 Non-Intertidal -> 2026 Tidal Flat)
var tfGain = murray2014.neq(murrayTF).and(myClassified2026.eq(myTF));

// Map Display
Map.centerObject(imageExtent, 10);

// Base Classified Image
Map.addLayer(myClassified2026, {min: 0, max: 3, palette: ['black', 'cyan', 'blue', 'green']}, 'My Classified Raster');

// Change Dynamics
Map.addLayer(tfStable.selfMask(), {palette: ['#808080']}, 'Stable Tidal Flat');
Map.addLayer(tfLoss.selfMask(), {palette: ['#FF0000']}, 'Tidal Flat Loss');
Map.addLayer(tfGain.selfMask(), {palette: ['#00FF00']}, 'Tidal Flat Gain');

// =========================================================================
// 4. TRANSITION MATRIX & QUANTITATIVE AREA CALCULATION
// =========================================================================
// Encode transitions into 2-digit codes: (Murray2014 * 10) + MyClass2026
// e.g., 11 = Stable TF | 10, 12, 13 = Loss | 01 = Gain
var transitionImage = combined.select('class_2014').multiply(10)
  .add(combined.select('class_2026'))
  .rename('transition');

// Add area in square kilometers
var areaImage = ee.Image.pixelArea().divide(1e6).addBands(transitionImage);

// Compute transition statistics
var transitionStats = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName: 'code',
  }),
  geometry: imageExtent,
  scale: 30,
  crs: 'EPSG:32652' , // hardcode instead // loadedClassifiedImage.projection(),
  maxPixels: 1e13,
  tileScale: 4
});

print("Raw Transition Area Output (Sq Km):", transitionStats);

// =========================================================================
// 5. PARSE MATRIX & UI CHART GENERATION
// =========================================================================
var processStats = function(stats) {
  var groups = ee.List(stats.get('groups'));
  
  var fc = ee.FeatureCollection(groups.map(function(item) {
    var dict = ee.Dictionary(item);
    var code = ee.Number(dict.get('code'));
    var area = ee.Number(dict.get('sum'));
    
    var mClass = code.divide(10).floor();
    var myClass = code.mod(10);
    
    return ee.Feature(null, {
      'Transition_Code': code,
      'Murray_2014_Class': mClass,
      'My_2026_Class': myClass,
      'Area_SqKm': area
    });
  }));
  
  return fc;
};

var transitionFC = processStats(transitionStats);

// Console Summary Chart
var chart = ui.Chart.feature.byFeature({
  features: transitionFC,
  xProperty: 'Transition_Code',
  yProperties: ['Area_SqKm']
}).setChartType('ColumnChart')
  .setOptions({
    title: 'Intertidal Area Transition Breakdown (Sq Km)',
    hAxis: {title: 'Transition Code (11=Stable, 10/12/13=Loss, 01=Gain)'},
    vAxis: {title: 'Area (sq km)'},
    colors: ['#2b5c8f']
  });

print(chart);

// Export transition table to Google Drive
Export.table.toDrive({
  collection: transitionFC,
  description: 'Intertidal_Change_Matrix_2014_2026',
  fileFormat: 'CSV'
});