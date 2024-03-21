
///This Land Productivity Dynamics (LPD) algorithm developed by   
// Dr. Cesar Luis Garcia and Dr. Ingrid Teich .- cesarnon@gmail.com - ingridteich@gmail.com 
// Under a WOCAT-FAO letter of agreement for the CACILM2-Project

// Below there are some explanation of what can be chosen when the script... 
//  and the particular choices (parametrization) made for this general version:
//1.- Making correction to Bad Quality Pixels when doing the annual mean
//2.- Using MannKendall p <0.5 for the positive or negative trends
//3.- MTID adapted to 3 year mean to reduce outlier impacts (more biological/climate meaning)
//4.- Steadiness index is based on the a combination of Mann-Kendall Trend and MTDI helps when there is no significance trend
//5.- Biomass Low is NDVI <0.4 and High is NDVI >0.7.
//6.- Water is masked to NoData 
//7.- Hyper-Arid areas are defined as never going over NDVI 0.25 in all the period, and are set as "Stable"

// Given the special conditions and particularities of countries in CA users might prefer to have more control on thresholds, masks and methods to calculate the different LPD components. So in this script there are sections where you need to choose options for making the final LPD -  By Commenting or un-Commenting you can activate the preferred method for each component.
//Some of the questions you may ask yourself or the team and decisions to make:

// Which is the: Period of time (Start and End dates)
// Should a correction of the Bad quality pixels be implemented in the Modis series (cloud/ice/etc.)?
// How to calculate Steadiness (trend + MDTI: Trend is LinearFit or MannKendall (parametric or non-Parametric)? use only significance Trends and set a threshold?
// MTDI should be made with the Last year integral? Or use a Last Period should be used to avoid that one atypical last year govern the destiny of the analysis?
// Initial biomass thresholds.... using user defined values or should it be computed using a Mean and Standard Deviation over an area (country, sampling area, etc.)
// Parametrization of the State calculation: which will be the Length of periods (Baseline, initial, end) and how many percentiles means that we should consider a change? 
// Should water be masked out? Should desert and bare areas also be masked out considering the biological meaning in the results?

//Palletes:
var mk_palette = {"opacity":1,"bands":["NDVI"],"min":-125,"max":125,"palette":["d21603","ffee75","16ce1b"]}
var traject_palette = {"opacity":1,"bands":["constant"],"palette":["080808","da310b","fff79c","33a117"]}

var SubfixName = 'MK_MTDI3_Des025'
// Select Country
/*
var LSIB = ee.FeatureCollection("USDOS/LSIB/2017");
var roiName = 'YourCountry'
var country = LSIB.filterMetadata('COUNTRY_NA', 'equals', 'YourCountry')
*/

// Choose All World

var roiName = 'World'
var globe = ee.Geometry.Polygon([-180, 88, 0, 88, 180, 88, 180, -88, 0, -88, -180, -88], null, false);
var country = globe


// definition of the initial and final year - This script is for standar calendar year and not hydrological year
var yrStart = 2005;
var yrEnd = 2019; 

// for calendar year 2000 is not complete!!!
var years = ee.List.sequence(yrStart, yrEnd);
print(years, 'years for calendar analysis')

var initDate = yrStart + '-01-01'
var endDate = yrEnd + '-12-31'

//Load NDVI collection
var modis = ee.ImageCollection('MODIS/006/MOD13Q1')
.filterDate(initDate, endDate);
print(modis, 'collection to process');

//Load water cover from same sensor
//make a mask for more permanent water bodies
var modiswater = ee.ImageCollection('MODIS/006/MOD44W').select('water_mask')
var waterSum = modiswater.sum()//.reproject('EPSG:4326', null, 250)
var waterMask = waterSum.gte(12)// limit set on 12 years

//Map.addLayer(waterSum,{},'waterMask',false)

//create a variable for band name
var nn = ee.String('ndvi_')

// Get the NDVI Annual mean for every calendar Year 
// and replace bad quality pixels with anual mean


var byYearCal = ee.ImageCollection.fromImages(
      years.map(function (y) {
        //get the subset for the target year
        var modisYear = modis.filter(ee.Filter.calendarRange(y, y, 'year'));
        //get the mean for NDVI
        var modisYearMean = modisYear.select('NDVI').mean();
        // Make a funtion to replace bad pixels with the mean
        var maskQAYear = function(image) {
          var image2 = image.select('NDVI');
          var image3 = image2.where(image.select("SummaryQA").gte(2),modisYearMean);
        return image3.rename('NDVI');
        }
        var ModisYearCorrected = modisYear.map(maskQAYear)
        var ModisYearCorrmean = ModisYearCorrected.mean()
                        .set('year', y)
                        .set('system:time_start',ee.Date.fromYMD(y,1,1))
                        
        return ModisYearCorrmean
}));
print(byYearCal,'Colection of calendar year means');

/*
Map.addLayer(byYearCal)
var stacked = byYearCal.toBands()
print(stacked,'stacked');
Map.addLayer(stacked,{},'stacked')
*/


//------- // This is the Default for making a annual mean collection----
// Get the NDVI Annual mean for every calendar Year 
// do not apply correction to bad quality pixels with anual mean
/*
var byYearCal = ee.ImageCollection.fromImages(
      years.map(function (y) {
        //get the subset for the target year
        var modisYear = modis.filter(ee.Filter.calendarRange(y, y, 'year'));
        var ModisYearCorrmean = modisYear.select('NDVI').mean()
                        .set('year', y)
                        .set('system:time_start',ee.Date.fromYMD(y,1,1))
                        
        return ModisYearCorrmean
}));
var stackedMean = byYearCal.toBands()
print(stackedMean,'stackedMean');
Map.addLayer(stackedMean,{},'stackedMean',false)
*/


//---------------------------------------------------------------------
//******************** Calculation starts *****************************
//---------------------------------------------------------------------


//// Calculate trend with MannKendall
//-------------------------------------------------------------------------
//       -----Calculate the significance with MannKendall no parametric------
// **
// ** Code to generate temporal ndvi analysis for the project 'Enabling the use of global
// ** data sources to assess and monitor land degradation at multiple scales' 
// ** http://www.conservation.org/about/gef/Pages/NDVI.aspx
// ** by Mariano Gonzalez-Roglich (mgonzalez-roglich@conservation.org)

// define kendal parameter values for a significance of 0.1, 0.05 and  0.01

var period = yrEnd - yrStart + 1
var coefficients90 = ee.Array([4,6,7,9,10,12,15,17,18,22,23,27,
                                28,32,35,37,40,42,
                               45,49,52,56,59,61,66,68,73,75,80,84,87,91,94,98,103,
                               107,110,114,119,123,128,132,135,141,144,150,153,159,
                               162,168,173,177,182,186,191,197,202]); 
var coefficients95 = ee.Array([4,6,9,11,14,16,19,21,24,26,31,33,
                              36,40,43,47,50,54,
                               59,63,66,70,75,79,84,88,93,97,102,106,111,115,120,
                               126,131,137,142,146,151,157,162,168,173,179,186,190,
                               197,203,208,214,221,227,232,240,245,251,258]);
var coefficients99 = ee.Array([6,8,11,18,22,25,29,34,38,41,47,50,
                              56,61,65,70,76,81,
                               87,92,98,105,111,116,124,129,135,142,150,155,163,170,
                               176,183,191,198,206,213,221,228,236,245,253,260,268,
                               277,285,294,302,311,319,328,336,345,355,364]); 
var kendal90 = coefficients90.get([period - 4]);
var kendal95 = coefficients95.get([period - 4]);
var kendal99 = coefficients99.get([period - 4]);

// Choose your level of significance:

var kendalS = kendal95

// ===Mann Kendalls S statistic=========================================================================
//    This function returns the Mann Kendalls S statistic, assuming that n is less than 40,             
//  the significance of a calculated S statistic is found in table A.30 of Nonparametric                
//  Statistical Methods, second edition by Hollander & Wolfe.                                           
//  reproduced for conveniance:                                                                         
//  https://www.dropbox.com/s/zc5u3oc1e3ou2me/v2appendixc.pdf?dl=0   _ Thanks Cesar! :-)
// =====================================================================================================
var f_MannKendallStat = function(imageCollection) {
  var TimeSeriesList = imageCollection.toList(50);
  var NumberOfItems = TimeSeriesList.length().getInfo();
  var ConcordantArray = [];
  var DiscordantArray = [];
  for (var k = 0; k <= NumberOfItems-2; k += 1) {
    var CurrentImage = ee.Image(TimeSeriesList.get(k));
    var l = k + 1;
    for (l; l <= NumberOfItems-1; l += 1) {
      var nextImage = ee.Image(TimeSeriesList.get(l));
      var Concordant = CurrentImage.lt(nextImage);
      ConcordantArray.push(Concordant);
      var Discordant = CurrentImage.gt(nextImage);
      DiscordantArray.push(Discordant);}}
  var ConcordantSum = ee.ImageCollection(ConcordantArray).sum();
  var DiscordantSum = ee.ImageCollection(DiscordantArray).sum();
  var MKSstat = ConcordantSum.subtract(DiscordantSum);   
  return MKSstat;};
  
// Compute Kendall statistics

var mk_trend_CalNDVI  = f_MannKendallStat(byYearCal.select('NDVI'))

Map.addLayer(mk_trend_CalNDVI, mk_palette,'mk_trend_cal',false)

//------ make the Coded Maps---------------

//Make coded maps  for 3 categories where:
// 1: Negative trend - Significative
// 2: No significative Trend
// 3: Possitive trend - Significative

var TrajectCalNDVI_MK = ee.Image(0)
  .where(mk_trend_CalNDVI.lt(0).and(mk_trend_CalNDVI.abs().gte(kendalS)),1) //Significant negative
  .where(mk_trend_CalNDVI.abs().lt(kendalS),2)                              // No Significant
  .where(mk_trend_CalNDVI.gt(0).and(mk_trend_CalNDVI.abs().gte(kendalS)),3) // Significant possitive
  .where(waterMask.eq(1),0)                                                 // No da
  //.clip(zona)


Map.addLayer(TrajectCalNDVI_MK,traject_palette,'TrajectCalNDVI_MK',false)

var finalTrend = TrajectCalNDVI_MK

/*
//------- // This is the Default for Trend----
//-------------------------------------------------------------------------
//// Calculate trend with LinearFit
//-------------------------------------------------------------------------
// Calculate Slopes for NDVI using parametric linearFit...
// Alternative it can be done with ThielSen non-parametric, which also is less sensitive to outliers.

var createTimeBand = function(image) {
    return image.addBands(image.metadata('year').rename('time'));
  };
var byYearCalTimed = byYearCal.map(createTimeBand)

// Calculate the slope and intersect for NDVI (linearfit and thielSen)
var LiFitCalNDVI = byYearCalTimed.select(['time', 'NDVI'])
    .reduce(ee.Reducer.linearFit());

var LiFitScale = LiFitCalNDVI.select(["scale"]);

var TrajectCalNDVI_LF = ee.Image(0)
  .where(LiFitScale.lt(0),1) // negative
  .where(LiFitScale.eq(0),1)// not many values here,so it may go to either way in steadiness
  .where(LiFitScale.gt(0),3) // possitive
  .where(waterMask.eq(1),0) // No data
  //.clip(zona)

var finalTrend = TrajectCalNDVI_LF
*/

//-------------
//// Calculate MTID 
//-------------------------------------------------------------------------

// The Multi Temporal Image Differencing (MTID).


// Traditional way based in (Guo et al. 2008) see below to avoid default
var MTID = function(imageCollection) {
  var TimeSeriesList = imageCollection.toList(50);
//  var TimeSeriesList = byMean.toList(50);
  var NumberOfItems = TimeSeriesList.length().getInfo();
  print (NumberOfItems,'NumberOfItems')
  var sumArray = [];
// ------ Choose the Last image Default or not!!!!
  //var lastImage = ee.Image(TimeSeriesList.get(NumberOfItems-1)); // This is the default considering just last year 
      // Alternative with more biological meaning based on the computation of anomalies
      // the logic is to detect the direction of net change of last period to reduce outlier impacts (more biological/climate meaning)
      var lastImage = imageCollection.select('NDVI')
                  .filterMetadata('year', 'greater_than', yrEnd -3)//last 3 years
                .mean()
  //print (lastImage, 'lastImage')
  for (var k = 0; k <= NumberOfItems-2; k += 1) {
    var currentImage = ee.Image(TimeSeriesList.get(k));
    sumArray.push(lastImage.subtract(currentImage))
    //print (sumArray,'sumArray')
    //print (currentImage, 'currentImage'+k)
    }
  var finalSum = ee.ImageCollection(sumArray).sum();
  //print (finalSum,'finalSum')
  return finalSum;};


//Compute MTDI
var MTDINDVI  = MTID(byYearCal.select('NDVI'))


Map.addLayer(MTDINDVI,{palette:["d21603","ffee75","16ce1b"]},'MTDINDVI',false)
var MTDIcode= MTDINDVI.lte(0)
              .where(MTDINDVI.gt(0),2)
              .where(waterMask.eq(1),0)

Map.addLayer(MTDIcode,{min:1, palette: ["#d21603","16ce1b"]}, "MTDIcode",false);

////-----------------------------------------------------------------------
//// Calculate steadiness 
//-------------------------------------------------------------------------
// Calculate the 4 value steadiness index  based on the a combination of Mann-Kendall Trend and MTDI
//where MTDI helps when there is no significance

var steadiness = ee.Image(0)
  .where(finalTrend.eq(1).and(MTDIcode.eq(1)),1) // T- MTDI-
  .where(finalTrend.eq(1).and(MTDIcode.eq(2)),2) // T- MTDI+
  .where(finalTrend.eq(2).and(MTDIcode.eq(1)),2) // T 0 MTDI-
  .where(finalTrend.eq(2).and(MTDIcode.eq(2)),3) // T 0 MTDI+
  .where(finalTrend.eq(3).and(MTDIcode.eq(1)),3) // T+ MTDI-
  .where(finalTrend.eq(3).and(MTDIcode.eq(2)),4) // T+ MTDI+
  //.where(waterMask.eq(1),0);   
  

// Calculate the 4 value steadiness index  based on only Mann-Kendall
/*
var steadiness = ee.Image(0)
  .where(finalTrend.eq(1),1) // T- significant
  .where(finalTrend.eq(2).and(mk_trend_CalNDVI.lt(0)),2) // slope - No significant
  .where(finalTrend.eq(2).and(mk_trend_CalNDVI.gt(0)),3) // slope + No significant
  .where(finalTrend.eq(3).and(MTDIcode.eq(2)),4) // T+ significant
  //.where(waterMask.eq(1),0);   
*/

// Calculate the 4 value steadiness index  based on LinearFit and MTDI  ------ This is the default
//there is no significance involved here as proposed also in https://drive.google.com/file/d/14XvPHvIV0_xh5ujz5ETawE0SOKw-maJi/view?usp=sharing
/*
var steadiness = ee.Image(0)
  .where(finalTrend.eq(1).and(MTDIcode.eq(1)),1) // T- MTDI-
  .where(finalTrend.eq(1).and(MTDIcode.eq(2)),2) // T - MTDI+
  .where(finalTrend.eq(3).and(MTDIcode.eq(1)),3) // T + MTDI-
  .where(finalTrend.eq(3).and(MTDIcode.eq(2)),4) // T+ MTDI+
  //.where(waterMask.eq(1),0);   
*/

Map.addLayer(steadiness,{min:1, palette: ["d21603","ffee75","16ce1b"]}, "steadiness",false);  
//*************************************************************************
//*************************************************************************
//*************************************************************************
// Calculation of initial Biomass

var NDVIinitialMean = byYearCal.select('NDVI')
              .filterMetadata('year', 'less_than', yrStart+3)// first 3 years
              .mean()
              .divide(10000);
              
var lowBio = ee.Number(0.4);
var highBio = ee.Number(0.7);

var initialBiomass = ee.Image(0)
      .where(NDVIinitialMean.lte(lowBio),1) // Low
      .where(NDVIinitialMean.gt(lowBio).and(NDVIinitialMean.lte(highBio)),2) // medium
      .where(NDVIinitialMean.gte(highBio),3)  // high

Map.addLayer(initialBiomass,{min:1, palette: ["d21603","ffee75","16ce1b"]}, "initialBiomass",false);

//*************************************************************************
//*************************************************************************
//*************************************************************************
// Calculation of State o % change performance

// State Baseline for percentiles
var NDVImeanTbase = byYearCal.select('NDVI')
              .filterMetadata('year', 'less_than', yrStart+15)//baseline first 15 years
              //.mean()
print (NDVImeanTbase,'NDVImeanTbase')
var bl_ndvi_perc = NDVImeanTbase.reduce(ee.Reducer.percentile([0,10,20,30,40,50,60,70,80,90,100]));
print(bl_ndvi_perc,'bl_ndvi_perc')
Map.addLayer(bl_ndvi_perc,{},'bl_ndvi_perc',false)

// State Time 1
var NDVImeanT1 = byYearCal.select('NDVI')
              .filterMetadata('year', 'less_than', yrStart+4)//last 4 years
              //.mean()
print (NDVImeanT1,'NDVImeanT1')

// State Time 2
var NDVImeanT2 = byYearCal.select('NDVI')
              .filterMetadata('year', 'greater_than', yrEnd -4)//last 4 years
              //.mean()
print (NDVImeanT2,'NDVImeanT2')

//get the mean
var t1_ndvi_mean = NDVImeanT1.mean().rename('ndvi');
var t2_ndvi_mean = NDVImeanT2.mean().rename('ndvi');

Map.addLayer(t1_ndvi_mean,{},'t1_ndvi_mean',false)
Map.addLayer(t2_ndvi_mean,{},'t2_ndvi_mean',false)

// emerging degradation
// reclassify mean ndvi for baseline period based on the percentiles
var t1_classes = ee.Image(-32768).where(t1_ndvi_mean.lte(bl_ndvi_perc.select('NDVI_p10')), 1)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p10')), 2)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p20')), 3)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p30')), 4)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p40')), 5)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p50')), 6)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p60')), 7)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p70')), 8)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p80')), 9)
                                 .where(t1_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p90')),10)
Map.addLayer(t1_classes,{},'t1_classes',false)
// reclassify mean ndvi for target period based on the percentiles
var t2_classes = ee.Image(-32768).where(t2_ndvi_mean.lte(bl_ndvi_perc.select('NDVI_p10')), 1)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p10')), 2)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p20')), 3)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p30')), 4)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p40')), 5)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p50')), 6)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p60')), 7)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p70')), 8)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p80')), 9)
                                 .where(t2_ndvi_mean.gt(bl_ndvi_perc.select('NDVI_p90')),10)
Map.addLayer(t2_classes,{},'t2_classes',false)

// emerging state: difference between start and end clusters >= 2 percentile
var classes_chg = t2_classes.subtract(t1_classes).updateMask(t1_ndvi_mean.neq(-32768));

var emeState = ee.Image(0).where(classes_chg.lte(-2),1)                            // negative
                              .where(classes_chg.gte(-1).and(classes_chg.lte( 1)), 2)   // stable
                              .where(classes_chg.gte( 2), 3)                           //  possitive
                              .where(t1_ndvi_mean.subtract(t2_ndvi_mean).abs().lte(500),2);//Deal with small differences
Map.addLayer(emeState,{min:0,max:3, palette:['#000000',"d21603","ffee75","16ce1b"]},'emeState',false)


//*************************************************************************
//*************************************************************************
//*************************************************************************
// Calculation of Final expression
/*
logical tabulation/combination of 4 steadiness classes, 
with 3 classes (low, medium, high) of initial standing biomass 
and 3 classes of change of GPP proxy between beginning and end of time series - State
to 36 combinations and assignment to the final 5 LPD classes
*/

var semifinal_lpd = ee.Image().expression(
  '(a == 1 && b==1 && c==1 )?1:'+ //1
  '(a == 1 && b==2 && c==1 )?2:'+ //2
  '(a == 1 && b==3 && c==1 )?3:'+ //3
  '(a == 1 && b==1 && c==2 )?4:'+  //4
  '(a == 1 && b==2 && c==2 )?5:'+  //5
  '(a == 1 && b==3 && c==2 )?6:'+  //6
  '(a == 1 && b==1 && c==3 )?7:'+  //7
  '(a == 1 && b==2 && c==3 )?8:'+  //8
  '(a == 1 && b==3 && c==3 )?9:'+  //9
  '(a == 2 && b==1 && c==1 )?10:'+ //10
  '(a == 2 && b==2 && c==1 )?11:'+ //11
  '(a == 2 && b==3 && c==1 )?12:'+ //12
  '(a == 2 && b==1 && c==2 )?13:'+  //13
  '(a == 2 && b==2 && c==2 )?14:'+  //14
  '(a == 2 && b==3 && c==2 )?15:'+  //15
  '(a == 2 && b==1 && c==3 )?16:'+  //16
  '(a == 2 && b==2 && c==3 )?17:'+  //17
  '(a == 2 && b==3 && c==3 )?18:'+  //18
  '(a == 3 && b==1 && c==1 )?19:'+ //19
  '(a == 3 && b==2 && c==1 )?20:'+ //20
  '(a == 3 && b==3 && c==1 )?21:'+ //21
  '(a == 3 && b==1 && c==2 )?22:'+  //22
  '(a == 3 && b==2 && c==2 )?23:'+  //23
  '(a == 3 && b==3 && c==2 )?24:'+  //24
  '(a == 3 && b==1 && c==3 )?25:'+  //25
  '(a == 3 && b==2 && c==3 )?26:'+  //26
  '(a == 3 && b==3 && c==3 )?27:'+  //27
  '(a == 4 && b==1 && c==1 )?28:'+ //28
  '(a == 4 && b==2 && c==1 )?29:'+ //29
  '(a == 4 && b==3 && c==1 )?30:'+ //30
  '(a == 4 && b==1 && c==2 )?31:'+  //31
  '(a == 4 && b==2 && c==2 )?32:'+   //32
  '(a == 4 && b==3 && c==2 )?33:'+  //33
  '(a == 4 && b==1 && c==3 )?34:'+  //34
  '(a == 4 && b==2 && c==3 )?35:'+  //35
  '(a == 4 && b==3 && c==3 )?36:99'   //36
  ,
    {
        a: steadiness, //Steadiness (1,2,3,4)
        b: initialBiomass, // (1,2,3)
        c: emeState, // GPP change (1,2,3)
    }
  );
  
print('Semi-Final', semifinal_lpd);
Map.addLayer(semifinal_lpd,{},'Semi-Final',false);

var final_lpd = semifinal_lpd.expression(
  '(a > 0 && a<9)?1:'+ //1
  '(a > 8 && a<15)?2:'+ //2
  '(a > 14 && a<23)?3:'+ //3
  '(a > 22 && a<33)?4:'+ //4
  '(a > 32 && a<37)?5:0' //5
      
  ,
    {
        a: semifinal_lpd,
    }
  );

print('Land Productivity Country Image', final_lpd);


// Make a mask to remove desert areas also need to un-comment the line in the final step that creates LPD
// alternatively you can make all desert areas also Stable category
var modisMax = ee.ImageCollection('MODIS/006/MOD13Q1')
.filterDate(initDate, endDate).select('NDVI').max();
var desertMask = modisMax.lte(2500) // This will identify all areas that in the period never went over 0.25 NDVI;
Map.addLayer(desertMask,{},'desertmask',false)

//Final LPD vis param
var lpd_vispar = {max: 5,
min: 1,
opacity: 1,
palette: ['#f23c46', '#e9a358', '#e5e6b3', '#a9afae', '#267300'],
};



var LPD = final_lpd
          .where(desertMask.eq(1),4).updateMask(1) //Desert mask (make sure desertMask line is not commented)
          .where(waterMask.eq(1),0).updateMask(1) //Water mask (make sure waterMask line is not commented)
          .selfMask().rename('LPD')
          .clip(country) // commnet this line if wish to see the global layer

//Map.centerObject(country,5,false)

//Map.addLayer(lpd_unccd,lpd_vispar,'LPD UNCCD',false);
Map.addLayer(LPD,lpd_vispar,'Final LPD_'+ yrStart+'_'+yrEnd);


Export.image.toDrive({
  image:LPD.toByte(), 
  description: 'LPD_'+ yrStart+'_'+yrEnd+'_'+roiName+'_'+SubfixName, 
  region: country, 
  scale: 250, 
  crs: 'EPSG:4326',
  folder:'LPD_country',
  maxPixels: 1e13
    });
    
Export.image.toAsset({
  image: LPD.toByte(),
  description: 'LPD_'+ yrStart+'_'+yrEnd+'_'+roiName+'_'+SubfixName,
  assetId:'LPD_'+ yrStart+'_'+yrEnd+'_'+roiName+'_'+SubfixName,
  scale: 250,
  region: country,
  crs: 'EPSG:4326',
  //crs:'SR-ORG:6974',
  //dimensions: [172800,86400],
  maxPixels: 1e13
});    


//       Put the charts in a fix panel


// Create a panel to hold our widgets.
var panel = ui.Panel();
panel.style().set('width', '300px');

// Create an intro panel with labels.
var intro = ui.Panel([
  ui.Label({
    value: 'Explore the data',
    style: {fontSize: '20px', fontWeight: 'bold'}
  }),
  ui.Label('Click a point on the map to inspect.')
]);
panel.add(intro);

// panels to hold lon/lat values
var lon = ui.Label();
var lat = ui.Label();
panel.add(ui.Panel([lon, lat], ui.Panel.Layout.flow('horizontal')));


Map.onClick(function(coords) {
  // Update the lon/lat panel with values from the click event.
  lon.setValue('lon: ' + coords.lon.toFixed(2)),
  lat.setValue('lat: ' + coords.lat.toFixed(2));
  var point = ee.Geometry.Point(coords.lon, coords.lat); 
  
// Create an NDVI annual year chart.
  var calChart = ui.Chart.image.series(byYearCal, point, ee.Reducer.mean(), 250);
  calChart.setOptions({
    title: 'MODIS Annual mean NDVI calendar year',
    vAxis: {title: 'Index * 10000'},//, maxValue: 9000},
    hAxis: {title: 'Year', format: 'yyyy', gridlines: {count: 7}},
  });
panel.widgets().set(2, calChart);

// Create an NDVI  16day chart.
var monChart = ui.Chart.image.series(modis.select('NDVI'), point, ee.Reducer.mean(), 250);
  monChart.setOptions({
    title: 'MODIS 16 days NDVI',
    vAxis: {title: 'Index * 10000'},
    hAxis: {title: 'Calendar_Year', format: 'yyyy', gridlines: {count: 7}},
  });
panel.widgets().set(3, monChart);
})

Map.style().set('cursor', 'crosshair');
// Add the panel to the ui.root.
ui.root.insert(0, panel);
