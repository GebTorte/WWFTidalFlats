/**
 * Improved Yellow Sea Intertidal Change & Quantitative Analysis Code
 */

// -------------------------------------------------------------------------
// 0. Configuration & Inputs
// -------------------------------------------------------------------------
var coastBuffer = ee.FeatureCollection('projects/earthobserver-cc-476010/assets/yellowsea_coasts_1km');
var loadedClassifiedTif = ee.Image('projects/earthobserver-cc-476010/assets/classified/yellowsea_tidalflats_2017_2019_tile3');

var nativeProj = loadedClassifiedTif.projection();
var nativeScale = nativeProj.nominalScale();
var imageExtent = loadedClassifiedTif.geometry().bounds();

// Your Classified Image (0: Background/Land, 1: Tidal Flat, 2: Water)
var myClassified2026 = loadedClassifiedTif.select(0).rename('class_2026');

// -------------------------------------------------------------------------
// 1. Murray Intertidal Collection Processing
// -------------------------------------------------------------------------
var murrayCollection = ee.ImageCollection('UQ/murray/Intertidal/v1_1/global_intertidal');

// Helper function to extract and reclass Murray epoch
function prepareMurrayEpoch(startDate, endDate, newName) {
  var img = murrayCollection
    .filterDate(startDate, endDate)
    .first()
    .clip(imageExtent);
  
  // Remap Murray: Keep class values consistent:
  // Original: 0 = Non-intertidal, 1 = Intertidal (Tidal Flat)
  // We map directly to: 0 = Background/Land, 1 = Tidal Flat
  return img.select('classification').eq(1).rename(newName);
}

//var murray2011 = prepareMurrayEpoch('2011-01-01', '2013-12-31', 'class_2011');
// latest murray epoch
var murray2014 = prepareMurrayEpoch('2014-01-01', '2016-12-31', 'class_2014');

// -------------------------------------------------------------------------
// 2. Multi-Timestep Alignment & Combinations
// -------------------------------------------------------------------------
// Stack rasters using your classified image's projection and scale
var combined = ee.Image.cat([murray2011, murray2014, myClassified2026])
  .clip(imageExtent);

// -------------------------------------------------------------------------
// 3. Trajectory & Dynamics Analysis (2014 vs 2026)
// -------------------------------------------------------------------------
var m14 = combined.select('class_2014');
var my26 = combined.select('class_2026');

// Binary Change Mask Logic (1 = Tidal Flat, 0 = Other)
var tfStable = m14.eq(1).and(my26.eq(1)).rename('stable');
var tfLoss   = m14.eq(1).and(my26.neq(1)).rename('loss');
var tfGain   = m14.neq(1).and(my26.eq(1)).rename('gain');

// Dynamic Transition Map
// Coding Scheme: 10 * T1 + T2
// 11: Stable Flat, 10: Flat -> Non-Flat (Loss), 01: Non-Flat -> Flat (Gain)
var transitionImage = m14.multiply(10).add(my26).rename('transition');

// -------------------------------------------------------------------------
// 4. Quantitative Statistics (Area Breakdown)
// -------------------------------------------------------------------------
var pixelAreaSqKm = ee.Image.pixelArea().divide(1e6);
var areaByTransition = pixelAreaSqKm.addBands(transitionImage);

var transitionStats = areaByTransition.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName: 'code',
  }),
  geometry: imageExtent,
  crs: nativeProj,
  scale: nativeScale,
  maxPixels: 1e13,
  tileScale: 8
});

// Post-process & Print Area Breakdown
print('--- QUANTITATIVE AREA ANALYSIS ---');
var statsList = ee.List(transitionStats.get('groups'));
var formattedStats = statsList.map(function(item) {
  var dict = ee.Dictionary(item);
  var code = dict.get('code');
  var area = dict.get('sum');
  
  // Readable transition dictionary decoding
  var label = ee.Algorithms.If(ee.Number(code).eq(11), 'Stable Tidal Flat (11)',
              ee.Algorithms.If(ee.Number(code).eq(10), 'Tidal Flat Loss (10)',
              ee.Algorithms.If(ee.Number(code).eq(1),  'Tidal Flat Gain (01)', 'Other Transition')));
              
  return ee.Dictionary({'Transition': label, 'Code': code, 'Area_sqkm': area});
});

print('Transition Matrix Breakdown (Sq Km):', formattedStats);

// -------------------------------------------------------------------------
// 5. Spatial Metrics: Tidal Flat Patch Analysis
// -------------------------------------------------------------------------
// Connected components analysis to evaluate fragmentation in stable areas
var patchSize = tfStable.selfMask().connectedPixelCount(100, true);
var patchArea = patchSize.multiply(pixelAreaSqKm).rename('patch_area_sqkm');

var maxPatchStats = patchArea.reduceRegion({
  reducer: ee.Reducer.max().combine({
    reducer2: ee.Reducer.mean(),
    sharedInputs: true
  }),
  geometry: imageExtent,
  crs: nativeProj,
  scale: nativeScale,
  maxPixels: 1e13
});

print('--- SPATIAL PATTERN METRICS ---');
print('Stable Patch Metrics (Sq Km):', maxPatchStats);

// -------------------------------------------------------------------------
// 6. Map Visualization
// -------------------------------------------------------------------------
Map.centerObject(coastBuffer, 9);

// Vis Palettes
var classVis = {min: 0, max: 2, palette: ['#222222', '#00ffff', '#0000ff']};

// Background Layers
Map.addLayer(m14.selfMask(), {palette: ['#ffaa00']}, 'Murray 2014 Intertidal', false);
Map.addLayer(my26, classVis, 'My Classified 2026', false);

// Change Map Layers
Map.addLayer(tfStable.selfMask(), {palette: ['#888888']}, 'Stable Tidal Flat (Persistence)');
Map.addLayer(tfLoss.selfMask(),   {palette: ['#ff0000']}, 'Tidal Flat Loss');
Map.addLayer(tfGain.selfMask(),   {palette: ['#00ff00']}, 'Tidal Flat Gain');

// Patch Area Layer (Color ramp by size)
Map.addLayer(patchArea, {min: 0, max: 10, palette: ['blue', 'yellow', 'red']}, 'Patch Sizes (Sq Km)', false);