<!-- badges: start -->
[![R-CMD-check](https://github.com/NorskRegnesentral/MSLT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NorskRegnesentral/MSLT/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


# MSLT
Git page for the  multiple scale line transect (MSLT) R-package `mslt`.  
The `mslt` package uses line transect data to estimate the rate of a thinned spatial point process, which can be used to estimate the abundance of ,e.g., marine mammals. For a detailed description of the model, see https://doi.org/10.1080/01621459.2025.2566422.


### Installation

The model is installed by typing 

```R
remotes::install_github("NorskRegnesentral/MSLT/MSLT")
```

# A quick example

Here is a quick example of how to use the mslt package. For the full R code to run an example, see the section "A runnable example" blow.

### Data

The length data must be in the format of a data frame with the following columns: Description of each columns  are provided below

```R
code            time         lon     lat       vessel  Port Stbd   dist     size     perpDist
  0   2019-01-14 10:56:10   -64.886 -56.934    KPH    EMB  <NA>   0.000    NA       NA
  1   2019-01-14 10:59:45   -64.881 -56.938    KPH    EMB  <NA>   0.568    NA       NA
  1   2019-01-14 10:59:45   -64.876 -56.942    KPH    EMB  <NA>   0.568    NA       NA
  1   2019-01-14 11:01:04   -64.871 -56.946    KPH    EMB  <NA>   0.451    NA       NA
  2   2019-01-14 11:02:56   -64.866 -56.950    KPH    EMB  <NA>   0.581    2        0.221
  1   2019-01-14 11:06:01   -64.861 -56.954    KPH    EMB  <NA>   0.506    NA       NA
```

Here is a description of each column in the line transect survey data:

`code`           Description of what the row represents: 0) start of a transect. 1) no fin whale was detected. 2) a fin whale group was detected.

`time`           Time of recording (POSIXct).  

`lon`            Longitude of this point on the transect (numeric). 

`lat`            Latitude of this point on the transect (numeric).

`vessel`         Name of the vessel

`Port`           Observer on the port side (string); NA if no observer is present.

`Stbd`           Observer on the starboard side (string); NA if no observer is present.

`dist`           Distance in kilometers with sampling effort from the previous row (numeric); defined as zero when the row represents the start of a transect.

`size`           Number of animals in the detected group (numeric); NA if no group is detected.

`perpDist`  	   Perpendicular distance in kilometers to the detected whale group (numeric).


### Set up configurations and run model



Set up configurations, see `?defConf` for details:

```r
conf = defConf(...)

```

Set up the data and the parameters: 
```r
data = setData(surveyData,predAreaUTM,conf)
par = setPar(data,conf)
```

Estimating the model: 
```r
run = fitMSLT(data,par,conf)
```


## A runnable example 
This example reproduces the main results in the JASA 2025 paper: "Spatial Variation on Multiple Scales in Line Transect Data; the Case of Antarctic Fin Whales"

The paper can be found here: https://doi.org/10.1080/01621459.2025.2566422


## Load the R-package:
```r
library(mslt)

```


## Read Data
Read the line transect data and the survey area. The data is stored in folder `papers/jasa2025/data`

```R
setwd(tempdir())
files = c("lineTransectData.rds", "surveyArea.rds")
url <- "https://raw.githubusercontent.com/NorskRegnesentral/MSLT/main/papers/JASA2025/data/"
d <- lapply(files, function(f)download.file(paste(url,f,sep=""), f))
surveyData <- readRDS("lineTransectData.rds")
predAreaUTM <- readRDS("surveyArea.rds")
```

Set configurations for the model:
```r
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
```

Set up the data and the parameters: 
```r
data = setData(surveyData,predAreaUTM,conf)
par = setPar(data,conf)
```

Estimating the model:

```r
run = fitMSLT(data,par,conf)
```



### Inspect Results
We can inspect the estimated model parameters reported in Table 1 in https://doi.org/10.1080/01621459.2025.2566422:
```r
run$rep
```

We can inspect the derived model parameters reported in Table 2 by calling `partable_derived(run)`. However, we have not yet performed the bias-correction to obtain $\hat{N}_{\text{BC}}$.


### Bias-correction
Bias correction can be performed like this:

```r
index = list(run$obj$env$ADreportIndex()$logAbundance) #Only do bias-correction for estimated abundance
sdrep <- TMB::sdreport(run$obj,
                               bias.correct = TRUE,
                               bias.correct.control = list(split = index,sd = TRUE),
                               getReportCovariance = FALSE)
```



We can inspect the bias-corrected abundance estimate reported in Table 2 in https://doi.org/10.1080/01621459.2025.2566422:

```r
partable_derived(run,sdrep)
```


## Estimate simpler submodels
The different submodels can be estimated by modifying the `conf` object accordingly: 

SPDE: To estimate the SPDE-only model, set `conf$mmpp = 0` and `conf$matern_intensity = 1`. 

MMPP: To estimate the MMPP-only model, set `conf$mmpp = 1` and `conf$matern_intensity = 0`.

HPP: To estimate the HPP model, set both `conf$mmpp = 0` and `conf$matern_intensity = 0`.

Only number of groups: To estimate only the number of groups, set `conf$applyPodSize = 0`.


## Practical convenient details

All estimated parameters are included in the object `run$pl`. That object has exactly the same format as the initial parameters `par`. Similarly, the standard deviations of the estimated parameters are stored in `run$plSd`. Note that the `run` object includes everything needed to reproduce and inspect the results. For example, the SPDE mesh is located here: `attributes(run$data)$mesh`.

