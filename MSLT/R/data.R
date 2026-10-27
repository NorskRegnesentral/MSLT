#' setData
#' @description Define data used internally based on data porvided and configuration file
#' @param d Data provided by user
#' @param panels_ll prediction area provided by user. A polygon with possible several sub-polygons in standard lat-lon coordinates
#' @param conf model configurations list  returned by \code{\link{defConf}}
#' @param digits Number of digits in line transect integration points and in prediction area. Sometimes needed because different OS may result in slightly different mesh because of numerical issues with small numbers.
#' @return data list used by TMB
#' @export
setData = function(d,panels_ll, conf, digits =6){

  panels_utm <- sf::st_transform(
    panels_ll,
    conf$UTMproj
  )
  
  #Different OS sometimes result in slightly different SPDE mesh because of rounding issue##
  geom <- sf::st_geometry(panels_utm)
  for(i in seq_along(geom)) {
    geom[[i]][[1]] <- round(geom[[i]][[1]], digits)
  }
  sf::st_geometry(panels_utm) <- geom
  panels_utm <- sf::st_sf(geometry = geom)
  ###

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
  
  #Different OS sometimes result in slightly different mesh because of rounding issues
  obsUTM = round(obsUTM,digits) 
  
  d$utmx = obsUTM[,1]
  d$utmy = obsUTM[,2]
  
  #Define mid-points as coordinates for integration points
  integrationPointsUTM = as.matrix(d[,c("utmx","utmy")])
  endpoints = integrationPointsUTM
  for(i in 1:dim(integrationPointsUTM)[1]){
    if(d$code[i]!=0){
      integrationPointsUTM[i,] = colMeans(endpoints[(i-1):i,])
    }
  }

  #Accommodate for observer sides
  d$abortLeft = rep(0,dim(d)[1])
  d$abortRight = rep(0,dim(d)[1])
  d$abortLeft[is.na(d$Port)] =1
  d$abortRight[is.na(d$Stbd)] =1

  #Define abundance integration points
  predAreaUTM = sf::st_union(panels_utm)
  
  grid = sf::st_make_grid(predAreaUTM, cellsize = conf$cellsize)
  grid_sf = sf::st_as_sf(grid)
  discretized_areas = sf::st_intersection(grid_sf, predAreaUTM)# Clip grid to predAreaUTM
  discretized_areas$area_total = sf::st_area(discretized_areas)# Area represented by each integration point
  discretized_areas$cell_id = seq_len(nrow(discretized_areas))
  
  points=sf::st_centroid(sf::st_geometry(discretized_areas))
  predData = data.frame(sf::st_coordinates(points))
  names(predData) <- c("utmx", "utmy")

  # Divide area into panels for panel-specific abundance estimates
  cell_panel <- sf::st_intersection(
    discretized_areas,
    panels_utm
  )
  cell_panel$panel_area = sf::st_area(cell_panel)
  if(is.null(cell_panel[["panel", exact = TRUE]])){
    cell_panel$panel = rep("A", dim(cell_panel)[1])
  }
  
  cell_integration <- sf::st_drop_geometry(
    cell_panel[c("cell_id", "panel","panel_area")]
  )
  
  buffer = sf::st_buffer(predAreaUTM,conf$buffer)
  mesh <- fmesher::fm_mesh_2d(obsUTM,
                              max.edge =conf$spdeDetails$max.edge,
                              boundary = buffer,
                              cutoff = conf$spdeDetails$cutoff,
                              min.angle = conf$spdeDetails$min.angle,
                              offset = conf$spdeDetails$offset)
  
  
  #Set up regression coefficients for detection 
  if(conf$covDetection[1] =="" & length(conf$covDetection)==1) {
    f_det <- lon ~ 1
  }else {
    # Check that all requested covariates exist in the data
    missing_covs <- setdiff(conf$covDetection, names(d))
    if(length(missing_covs) > 0){
      stop(
        "The following detection covariates are not in the data: ",
        paste(missing_covs, collapse = ", ")
      )
    }
    
    # Build formula: no intercept + covariates
    f_det <- as.formula(
      paste("lon ~ -1 +", paste(conf$covDetection, collapse = " + "))
    )
  }
  X_g_obs = mgcv::gam(f_det  , data = d[which(d$code==2),],fit =FALSE)$X
  X_g_left = mgcv::gam(f_det , data = d,fit =FALSE,)$X
  X_g_right = X_g_left

  
  #Set up regression coefficients for rate
  X_z = mgcv::gam(lon ~  1  , data = d,fit =FALSE)$X
  X_z_pred = mgcv::gam(utmx ~  1 , data = predData,fit =FALSE)$X

  #matrices needed in the SPDE procedure
  AalongLines_S <- fmesher::fm_basis(mesh,integrationPointsUTM)
  AalongLinesEndpoints_S <- fmesher::fm_basis(mesh,endpoints)
  AObs_S <- fmesher::fm_basis(mesh,obsUTM[d$code==2,])
  Apred <- fmesher::fm_basis(mesh,as.matrix(predData[,1:2]))
  
  spde <- fmesher::fm_fem(mesh)
  spdeMatrices = list(M0 = spde$c0, M1 = spde$g1,M2 = spde$g2)
  
  
  #Space-time A-matrices
  AObs_ST = list()
  AalongLines_ST = list()
  AalongLinesEndpoints_ST = list()
  if(conf$spaceTime_intensity_maternAR1){
    nYear = diff(range(d$year)) + 1
    for(t in 1:nYear){
      AObs_ST[[t]] = fmesher::fm_basis(mesh,obsUTM[d$code==2 & d$year == min(d$year) + t-1,])
      AalongLines_ST[[t]] <- fmesher::fm_basis(mesh,integrationPointsUTM[d$year == min(d$year) + t-1,])
      AalongLinesEndpoints_ST[[t]] <- fmesher::fm_basis(mesh,endpoints[d$year == min(d$year) + t-1,])
    }
  }else{
    nYear = 1
  }

  

  data = list(AalongLines_S = AalongLines_S,
              AalongLinesEndpoints_S = AalongLinesEndpoints_S,
              AalongLines_ST = AalongLines_ST,
              AalongLinesEndpoints_ST = AalongLinesEndpoints_ST,
              AObs_S = AObs_S,
              AObs_ST = AObs_ST,
              Apred = Apred,
              year = d$year,
              nYear = nYear,
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
              areas = as.numeric(discretized_areas$area),
              cell_integration = cell_integration,
              abortLeft = d$abortLeft,
              abortRight = d$abortRight,
              space_intensity_matern = conf$space_intensity_matern,
              spaceTime_intensity_maternAR1  = conf$spaceTime_intensity_maternAR1,
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

  #Needed for ridge correct only parts of data
  data$ridgeCorrectLine = rep(0,dim(data$AalongLines_S)[1])
  data$ridgeCorrectObservations = rep(0,dim(data$AObs_S)[1])
  
  #Add attributes that is not needed in the TMB program, but neat to have
  obsLatLon = data.frame(sf::st_coordinates(obsLatLon))
  colnames(obsLatLon) = c("lon","lat")
  obsUTM = data.frame(obsUTM)
  colnames(obsUTM) = c("UTMX","UTMY")
  
  attributes(data)$mesh = mesh
  attributes(data)$obsUTM = obsUTM
  attributes(data)$obsLatLon = obsLatLon
  attributes(data)$predData = predData
  attributes(data)$detectionCov = colnames(X_g_right)
  attributes(data)$surveyData = d
  attributes(data)$panels_utm = panels_utm
  return(data)
}

