#' setData
#' @description Define data used internally based on data porvided and configuration file
#' @param d Data provided by user
#' @param predAreaUTM prediction area provided by user
#' @param conf model configurations list  returned by \code{\link{defConf}}
#' @return data list used by TMB
#' @export
setData = function(d,predAreaUTM, conf){

  if(conf$detectionTrunc>0){
    ii = which(d$perpDist>conf$detectionTrunc & d$code==2)
    if(length(ii)>0){
      d$code[ii] = 1
    }
  }

  #Convert to UTM coordinates
  coords = data.frame(d[,c("lon","lat")])
  obsLatLon <- sf::st_as_sf(coords, coords = c("lon", "lat"),crs="+proj=longlat") 
  obsUTM <- sf::st_transform(obsLatLon, conf$UTMproj)
  obsUTM = sf::st_coordinates(obsUTM)
  obsLatLon = sf::st_coordinates(obsLatLon)
  d$utmy = obsUTM[,2]
  d$utmx = obsUTM[,1]
  
  #Define mid-points as coordinates for integration points
  integrationPointsUTM = as.matrix(d[,c("utmx","utmy")])
  endpoints = integrationPointsUTM
  for(i in 1:dim(integrationPointsUTM)[1]){
    if(d$code[i]!=0){
      integrationPointsUTM[i,] = colMeans(endpoints[(i-1):i,])
    }
  }

  #Accomondate for observer sides
  d$abortLeft = rep(0,dim(d)[1])
  d$abortRight = rep(0,dim(d)[1])
  d$abortLeft[is.na(d$Port)] =1
  d$abortRight[is.na(d$Stbd)] =1

  #Define aundance integration points
  points = sf::st_make_grid(predAreaUTM,cellsize=c(conf$cellsize,conf$cellsize),what="centers")
  points = sf::st_as_sf(points)
  points = sf::st_join(points,predAreaUTM,left=FALSE)
  predData = data.frame(sf::st_coordinates(points))
  names(predData) = c("utmx", "utmy")
  
  buffer = sf::st_buffer(predAreaUTM,conf$buffer)
  mesh <- fmesher::fm_mesh_2d(predData,
                              max.edge =conf$spdeDetails$max.edge,
                              boundary = buffer,
                              cutoff = conf$spdeDetails$cutoff,
                              offset = c(1,-0.1)
  )

  #Set up regression coefficients
  if(conf$covDetection== "vessel"){
    X_g_obs = gam(lon ~  -1 + vessel  , data = d[which(d$code==2),],fit =FALSE)$X
    X_g_left = gam(lon ~  -1 + vessel , data = d,fit =FALSE,)$X
    X_g_right = X_g_left
  }else{
    X_g_obs = gam(lon ~  1  , data = d[which(d$code==2),],fit =FALSE)$X
    X_g_left = gam(lon ~  1 , data = d,fit =FALSE)$X
    X_g_right = gam(lon ~  1 , data = d,fit =FALSE)$X
  }


  #Set up regression coefficients for rate
  X_z = gam(lon ~  1  , data = d,fit =FALSE)$X
  X_z_pred = gam(utmx ~  1 , data = predData,fit =FALSE)$X

  #matrices needed in the SPDE procedure
  AalongLines <- fmesher::fm_basis(mesh,integrationPointsUTM)
  AalongLinesEndpoints <- fmesher::fm_basis(mesh,endpoints)
  Apred <- fmesher::fm_basis(mesh,as.matrix(predData[,1:2]))
  AObs <- fmesher::fm_basis(mesh,obsUTM[d$code==2,])
  
  spde <- fmesher::fm_fem(mesh)
  spdeMatrices = list(M0 = spde$c0, M1 = spde$g1,M2 = spde$g2)
  
#  spdeI = inla.spde2.matern(mesh, alpha=2)
#  spdeMatrices = spdeI$param.inla[c("M0","M1","M2")]


  data = list(AalongLines = AalongLines,
              AalongLinesEndpoints = AalongLinesEndpoints,
              Apred = Apred,
              AObs = AObs,
              spdeMatrices = spdeMatrices,
              lineIntegralDelta = d$dist,
              code = d$code,
              distObs = d$perpDist[ which(d$code==2)],
              X_g_obs = as.matrix(X_g_obs),
              X_g_left = as.matrix(X_g_left),
              X_g_right = as.matrix(X_g_right),
              X_z = X_z,
              X_z_pred = X_z_pred,
              g_function = conf$g_function, #Which detection function to use
              penalize = conf$penalizeMMPP,
              area = as.numeric(sf::st_area(predAreaUTM)),
              abortLeft = d$abortLeft,
              abortRight = d$abortRight,
              matern_intensity = conf$matern_intensity,
              matern_size = conf$matern_size,
              pcPriorsRange_intensity = conf$pcPriorsRange_intensity,
              pcPriorsSD_intensity = conf$pcPriorsSD_intensity,
              pcPriorsRange_size = conf$pcPriorsRange_size,
              pcPriorsSD_size = conf$pcPriorsSD_size,
              usePCpriors = conf$usePCpriors,
              size = d$size[which(d$code==2 & !is.na(d$size))],
              applyPodSize = conf$applyPodSize,
              podSizeDist = conf$podSizeDist,
              independentPodSize = conf$independentPodSize,
              detectionTrunc = conf$detectionTrunc,
              useMMPP = conf$mmpp[1],
              spatialBiasCorFigure = 0,
              meanGroupFigure = 0
  )

  attributes(data)$mesh = mesh
  attributes(data)$obsUTM = obsUTM
  attributes(data)$obsLatLon = obsLatLon
  attributes(data)$predData = predData
  attributes(data)$detectionCov = colnames(X_g_right)
  return(data)
}

