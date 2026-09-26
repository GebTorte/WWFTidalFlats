
var coastBuffer    = ee.FeatureCollection('projects/earthobserver-cc-476010/assets/yellowsea_coasts_1km');

// TODO define all tiles (load combined TIF)
var loadedClassifiedTif = ee.Image('projects/earthobserver-cc-476010/assets/classified/yellowsea_tidalflats_2017_2019_tile2');
var imageExtent = loadedClassifiedTif.geometry().bounds();
// 1. Load Murray et al. 2014 dataset

var murrayCollection = ee.ImageCollection('UQ/murray/Intertidal/v1_1/global_intertidal');
var murray2014 = murrayCollection
  //.filterBounds(imageExtent)
  .filterDate('2014-01-01', '2016-12-31')
  .first()
  //.select('classification')
  //.rename('class_2014')
  //.clip(imageExtent);

// Reclassify Murray 2014 to match your class indices:
// [0: Background/Land, 1: Tidal Flat, 2: Water]
var murray2014 = murray2014.remap(
  [0, 1, 2], // Original Murray classes
  [1, 2, 3] // our mapping, shifted by 1
).rename('class_2014');


// 2. Prepare your classified image (ensure it is renamed)
var myClassified2026 = loadedClassifiedTif.select(0).rename('class_2026');

// 3. Combine both rasters aligned to your classified image's projection
var combined = ee.Image.cat([myClassified2026, murray2014])
  .clip(loadedClassifiedTif.geometry());
  //.reproject({
  //  crs: loadedClassifiedTif.projection(),
  //  scale: 30
  //});
  
// 4. Stable Tidal Flat calculation
var tfStable = combined.select('class_2014').eq(1)
  .and(combined.select('class_2026').eq(1));
  
  
/// trainsition matrix
// Combine classes into single 2-digit transition codes
// e.g., 12 = Tidal Flat (2014) -> Water (2026) [Loss/Erosion]
// e.g., 21 = Water (2014) -> Tidal Flat (2026) [Gain/Accretion]
var transitionImage = combined.select('class_2014').multiply(10)
  .add(combined.select('class_2026'))
  .rename('transition');

// Add pixel area in square kilometers
var areaImage = ee.Image.pixelArea().divide(1e6).addBands(transitionImage);

// Calculate total area per transition type
var transitionStats = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName: 'code',
  }),
  geometry: imageExtent,
  scale: 30,
  maxPixels: 1e13,
  tileScale: 4,
  crs: loadedClassifiedTif.projection()
});

// TODO: show
//print("Transition Breakdown (Sq Km):", transitionStats);


Map.centerObject(coastBuffer, 8);

// Display your raw image to verify it renders
Map.addLayer(myClassified2026.clip(imageExtent), {min: 0, max: 3, palette: ['black', 'cyan', 'blue', 'green']}, 'My Classified 2026');

// Display Murray 2014 to verify overlap
Map.addLayer(murray2014.clip(imageExtent).selfMask(), {palette: ['orange']}, 'Murray 2014 Tidal Flat');  
  
/////////////////////////////////
// THE below part does not work yet
/////////////////////////////////

var tfStable = murray2014.eq(1)
  .and(myClassified2026.eq(1))
  .clip(imageExtent);

// Display Stable Tidal Flat
Map.addLayer(tfStable.selfMask(), {palette: ['red']}, 'Stable Tidal Flat (Overlap)');  

// gain loss viz
// 1. Unchanged Tidal Flat (Persistence)
var tfStable = combined.select('class_2014').eq(1)
  .and(combined.select('class_2026').eq(1));

// 2. Tidal Flat Loss (2014 Tidal Flat -> 2026 Land/Water)
var tfLoss = combined.select('class_2014').eq(1)
  .and(combined.select('class_2026').neq(1));

// 3. Tidal Flat Gain (2014 Non-Tidal Flat -> 2026 Tidal Flat)
var tfGain = combined.select('class_2014').neq(1)
  .and(combined.select('class_2026').eq(1));

// Vis Layers
Map.addLayer(tfStable.selfMask(), {palette: ['gray']}, 'Stable Tidal Flat');
Map.addLayer(tfLoss.selfMask(), {palette: ['red']}, 'Tidal Flat Loss');
Map.addLayer(tfGain.selfMask(), {palette: ['green']}, 'Tidal Flat Gain');