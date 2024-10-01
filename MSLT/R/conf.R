#' defConf
#' 
#' @param matern_intensity 1: include spatial effect in group abundance intensity, 0: not 
#' @param mmpp  1: include MMPP in group abundance intensity, 0:not
#' @param detectionTrunc Distance at which detection is truncated
#' @param matern_size 1: include spatial effect in group size, 0: not
#' @param g_function 1: Half normal detection function, 2: Hazard detection function
#' @param covRate Covariates in group abundance intensity
#' @param covDetection Covariates in detection function
#' @param usePCpriors Use pc-priors (penalty) for spatial effect model parameters
#' @param pcPriorsRange_intensity pc-prior group abundance intensity range
#' @param pcPriorsSD_intensity pc-prior group abundance intensity marginal standard deviation
#' @param pcPriorsRange_size pc-prior group size intensity range
#' @param pcPriorsSD_size pc-prior group size intensity marginal standard deviation
#' @param spdeDetails Details for constructing mesh with use of INLA
#' @param applyPodSize 1: Pod size is included, 0: Not included
#' @param podSizeDist 1:poisson, 2:negative binomial group size
#' @param independentPodSize 1: Size independent of intensity
#' @param dependentPodSize First if beta0 is included and second if beta_SPDE is included
#' @param cellsize Distance between spatial quadratic integration points for total abundance 
#' @param UTMproj CRS-UTM projection
#' @param LatLonProj CRS-latitude and longitude projection
#' @param buffer A buffer (in km) outside of the area of interest. The buffer is used when constructing the boundary input to fmesher::fm_mesh_2d
#' @return Configurations 
#' @export
defConf = function(matern_intensity = 1,
                   mmpp = 1,
                   detectionTrunc = -1,
                   matern_size = 1,
                   g_function = 1,
                   covRate =  "",
                   covDetection =  "",
                   usePCpriors = 0,
                   pcPriorsRange_intensity = c(300,0.1),
                   pcPriorsSD_intensity = c(1,0.1),
                   pcPriorsRange_size = c(300,0.1),
                   pcPriorsSD_size = c(1,0.1),
                   penalizeMMPP = c(1,0.1),
                   spdeDetails = list(cutoff = 30,max.edge = c(80,150),min.angle = 14),
                   applyPodSize = 1,
                   podSizeDist = 2,
                   independentPodSize = 1 ,
                   dependentPodSize = c(0,1),
                   cellsize = 25,
                   UTMproj = sf::st_crs("+proj=utm +zone=22 +units=km"),
                   LatLonProj =  sf::st_crs("+proj=longlat +datum=WGS84"),
                   buffer = 80
                   ){
  conf = list()
  conf$matern_intensity = matern_intensity 
  conf$mmpp = mmpp 
  conf$matern_size = matern_size 
  conf$g_function = g_function 
  conf$UTMproj = UTMproj
  conf$LatLonProj = LatLonProj
  conf$covRate =  covRate
  conf$covDetection =  covDetection
  conf$pcPriorsRange_intensity = pcPriorsRange_intensity
  conf$pcPriorsSD_intensity = pcPriorsSD_intensity
  conf$pcPriorsRange_size = pcPriorsRange_size
  conf$pcPriorsSD_size = pcPriorsSD_size
  conf$usePCpriors = usePCpriors
  conf$penalizeMMPP = penalizeMMPP
  conf$spdeDetails = spdeDetails
  conf$applyPodSize = applyPodSize
  conf$podSizeDist = podSizeDist
  conf$independentPodSize = independentPodSize
  conf$dependentPodSize = dependentPodSize 
  conf$detectionTrunc = detectionTrunc 
  conf$cellsize = cellsize
  conf$buffer = buffer
  return(conf)
}



#' setMap
#' 
#' @param conf Model configurations returned by \code{\link{defConf}}
#' @param par Model parameters set by \code{\link{setPar}} 
#' @return map-variable used by TMB
#' @export
setMap = function(conf,par){
  map = list()
  map$log_sigma =c(0,1)
  map$log_kappa = c(0,1)
  if(conf$applyPodSize==1 & conf$independentPodSize ==0){
    map$beta_size = c(0,1)
    if(conf$dependentPodSize[1]==0){
      map$beta_size = map$beta_size-1
      map$beta_size[1] = as.factor(NA)
    }
    if(conf$dependentPodSize[2]==0){
      map$beta_size[2] = as.factor(NA)
    }
    map$beta_size = as.factor(map$beta_size)
  }
  if(conf$g_function==1){
    map$logB = as.factor(NA)
  }
  if(conf$matern_intensity==0){ #Remove spatial and mmpp effects.
    map$x_intensity = as.factor(rep(NA,length(par$x_intensity)))
    map$log_sigma[1] = NA
    map$log_kappa[1] = NA
  }
  if(conf$mmpp==0){
    map$log_c_mmpp = as.factor(NA)
    map$log_mu = as.factor(c(NA,NA))
  }
  if(conf$applyPodSize == 0 | conf$independentPodSize == 0 | conf$matern_size==0){
    map$x_size = as.factor(rep(NA,length(par$x_size)))
    map$log_sigma[2] = NA
    map$log_kappa[2] = NA
  }
  map$log_sigma = as.factor(map$log_sigma)
  map$log_kappa = as.factor(map$log_kappa)
  return(map)
}
