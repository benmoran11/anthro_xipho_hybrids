

var point = ee.Geometry.Point(-98.58305741232668, 20.896058641104556);

var start = ee.Date('2022-05-01');
var finish = ee.Date('2022-05-30');

var filteredCollection = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
  .filterBounds(point)
  .filterDate(start, finish)
  .filter(ee.Filter.lt('CLOUD_COVER', 9))
  .sort('CLOUD_COVER', true)
  .map(function(image){return image.clip(HGO_study_area)});
  
var focalimage1 = ee.Image("LANDSAT/LC09/C02/T1_L2/LC09_026045_20220521").clip(HGO_study_area).select('SR_B.*|ST.*').toInt16();
var focalimage2 = ee.Image("LANDSAT/LC09/C02/T1_L2/LC09_026046_20220521").clip(HGO_study_area).select('SR_B.*|ST.*').toInt16();
var May_study_area = filteredCollection.median().select('SR_B.*|ST.*');
  
var visParamsTrue = {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 6000, max: 15000}
  

Map.addLayer(May_study_area, visParamsTrue, 'Landsat May 2022');


//Export.image.toDrive({
//  image: May_study_area,
//  description: 'LANDSAT_May_2022_StudyArea',
//  folder: 'ee_exports',
//  scale: 30,
//  region: HGO_study_area,
//  maxPixels: 1e13
//});

var filteredCollection2 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
  .filterBounds(point)
  .filterDate(ee.Date('1986-05-09'), ee.Date('1986-05-11'))
  //.filter(ee.Filter.lt('CLOUD_COVER', 9))
  .sort('CLOUD_COVER', true);
  //.map(function(image){return image.clip(HGO_study_area)});

print(filteredCollection2);  
var early_study_area = filteredCollection2.median().clip(HGO_study_area).select('SR_B.*|ST.*');
var early_image1 =ee.Image("LANDSAT/LT05/C02/T1_L2/LT05_026045_19850131").clip(HGO_study_area).select('SR_B.*|ST.*').toInt16();
var early_image2 =ee.Image("LANDSAT/LT05/C02/T1_L2/LT05_026046_19860510").clip(HGO_study_area).select('SR_B.*|ST.*').toInt16();
var visParamsTrue = {bands: ['SR_B3', 'SR_B2', 'SR_B1'], min: 7000, max: 15000}

print(early_image1);  
//Map.addLayer(early_image1, visParamsTrue, 'Landsat Jan 1985');
//Map.addLayer(early_study_area, visParamsTrue, 'Landsat May 1986');

/**
 * Function to mask clouds using the Sentinel-2 QA band
 * @param {ee.Image} image Sentinel-2 image
 * @return {ee.Image} cloud masked Sentinel-2 image
 */
function maskS2clouds(image) {
  var qa = image.select('SCL');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.neq(3)
      .and(qa.neq(8))
      .and(qa.neq(9));

  return image.updateMask(mask);
}

var SentinelCollection = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(HGO_study_area)
  .filterDate(ee.Date('2020-05-01'), ee.Date('2023-05-31'))
  .filter(ee.Filter.calendarRange(5, 5,'month'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .sort('CLOUD_COVER', true)
  .map(maskS2clouds)
  .map(function(image){return image.clip(HGO_study_area)});
 
 print(SentinelCollection)
 var Sentinel_study_area = SentinelCollection.median();
 var Sentinel_Atlapexco = Sentinel_study_area.clip(table);
 var visParamsSentinel = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000}
 Map.addLayer(Sentinel_study_area, visParamsSentinel, 'Sentinel May 2022');
 Map.addLayer(Sentinel_Atlapexco, visParamsSentinel, 'Atlapexco May 2022');

Map.addLayer(table, undefined, "Atlapexco");

Export.image.toDrive({
  image: Sentinel_Atlapexco,
  description: 'Sentinel2_May_2020-23_AtlapexcoDrainage',
  folder: 'ee_exports',
  scale: 10,
  region: table,
  maxPixels: 1e13
});
