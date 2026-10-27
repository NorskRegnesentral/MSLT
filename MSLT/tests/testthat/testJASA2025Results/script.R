library(mslt)

surveyData = readRDS("testJASA2025Results/data/lineTransectData.rds") #Line transect data (data frame)

strataFull = suppressWarnings(sf::read_sf('testJASA2025Results/data/survey_strata.geojson',crs="+proj=longlat"))#Warning that some areas do not form a ring, see fix below for ESS
strataESS <- strataFull[which(strataFull$stratum %in% c('ESS')),]
coordESS = sf::st_coordinates(strataESS)[,1:2]
coordESS = rbind(coordESS,coordESS[1,]) #NB, something wrong with polygon, last point is not equal the first
polygonESS <- sf::st_polygon(list(coordESS))

strataSS <- strataFull[which(strataFull$stratum %in% c('SS')),]
coordSS = sf::st_coordinates(strataSS)[,1:2]
polygonSS <- sf::st_polygon(list(coordSS))

strataAP <- strataFull[which(strataFull$stratum %in% c('AP')),]
coordAP = sf::st_coordinates(strataAP)[,1:2]
polygonAP <- sf::st_polygon(list(coordAP))

tmp = sf::st_union(sf::st_union(polygonESS,polygonSS), polygonAP)
fullPolygon <- sf::st_sf(geometry = sf::st_sfc(tmp), crs="+proj=longlat")


conf = defConf(detectionTrunc = 3.5, #Truncation limit 3.5 km in our case, see Eq. (20)
               g_function = 2, #Hazard detection function, see eq. (20)
               cutoff = 8, #Maximum distance between SPDE mesh nodes along transects, see Section 2.2 and S3.1.
               max.edge = c(20,80), #Maximum distance between SPDE mesh nodes in survey area and outer area, see Section 2.2  and S3.1.
               cellsize = 50, #Distance between spatial integration points, see Eq (19)
               space_intensity_matern = 1, #Activate SPDE, see Eq (2)
               spaceTime_intensity_maternAR1 = 0, #Activate space-time
               mmpp = 1, #Activate MMPP , see Eq (2)
               applyPodSize = 1, #Activate group size model, see section 2.4.
               covDetection = "vessel" ,#Use vessel as a covariate in Eq. (2),  
               UTMproj =sf::st_crs("+proj=utm +zone=22 +units=km"),
               LatLonProj =  sf::st_crs("+proj=longlat +datum=WGS84")
)

data = setData(surveyData,panels_ll = fullPolygon,conf)
par = setPar(data,conf)

run = fitMSLT(data,par,conf)



resultsOut = list()
resultsOut$logAbundance = run$rl$logAbundance
resultsOut$objective = run$opt$objective


load("testJASA2025Results/resultsExp.RData")
expect_equal(resultsOut$logAbundance, resultsExp$logAbundance,tolerance = 1e-6)
expect_equal(resultsOut$objective, resultsExp$objective,tolerance = 1e-6)


if(FALSE){
  resultsExp = list()
  resultsExp$logAbundance = run$rl$logAbundance
  resultsExp$objective  = run$opt$objective
  save(resultsExp,file = "testJASA2025Results/resultsExp.RData")
}
