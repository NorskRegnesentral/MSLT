#' @param
#' @return Configurations
#' @export
#' @examples
defConf = function(d = NULL,detectionTrunc = -1){
  conf = list()

  conf$mmpp = 1 #apply MMPP
  conf$matern_intensity = 1 #apply spatial effect intensity
  conf$matern_size = 1 #apply spatial effect size
  conf$g_function = 1 #1: Normal, 2: Hazard
  conf$UTMproj = CRS("+proj=utm +zone=22 +units=km")
  conf$LatLonProj =  CRS("+proj=longlat +datum=WGS84")

  conf$covRate =  ""
  conf$covDetection =  "vessel"
  conf$minObserved = 0

  conf$pcPriorsRange_intensity = c(300,0.1)
  conf$pcPriorsSD_intensity = c(1,0.1)
  conf$pcPriorsRange_size = c(300,0.1)
  conf$pcPriorsSD_size = c(1,0.1)
  conf$usePCpriors = 0

  conf$penalizeMMPP = c(1,0.2)

  conf$spdeDetails = list()
  conf$spdeDetails$cutoff = 50 #TODO, use d
  conf$spdeDetails$max.edge = c(80,150) #TODO, use d

  conf$applyPodSize = 1 #1: Pod size is included, 0: Not included
  conf$podSizeDist = 2 #1:poisson, 2:negative binomial
  conf$independentPodSize = 1 #1: Size independent of intensity
  conf$dependentPodSize = c(0,1) #First if beta0 is included and second if beta_SPDE is included
  conf$detectionTrunc = detectionTrunc #No truncation if negative

  conf$cellsize = 20
  
  return(conf)
}



#' @param
#' @return map
#' @export
#' @examples
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
