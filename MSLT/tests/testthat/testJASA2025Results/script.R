library(mslt)

surveyData = readRDS("testJASA2025Results/data/lineTransectData.rds") #Line transect data (data frame)
predAreaUTM = readRDS("testJASA2025Results/data/surveyArea.rds") #Survey area (polygon)

conf = defConf(detectionTrunc = 3.5, #Truncation limit 3.5 km in our case, see Eq. (20)
               g_function = 2, #Hazard detection function, see eq. (20)
               cutoff = 8, #Maximum distance between SPDE mesh nodes along transects, see Section 2.2 and S3.1.
               max.edge = c(20,80), #Maximum distance between SPDE mesh nodes in survey area and outer area, see Section 2.2  and S3.1.
               cellsize = 50, #Distance between spatial integration points, see Eq (19)
               matern_intensity = 1, #Activate SPDE, see Eq (2)
               mmpp = 1, #Activate MMPP , see Eq (2)
               applyPodSize = 1, #Activate group size model, see section 2.4.
               covDetection = "vessel" #Use vessel as a covariate in Eq. (2)
)

data = setData(surveyData,predAreaUTM,conf)
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
