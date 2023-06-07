#' readData This provides the data in the format that the user is supposed to have
#' @description ...
#' @param maxDist maximum distance between integration points
#' @param minDist minimum distance between integration points
#' @return Processed data
#' @export
#' @examples
readData = function(maxDist = 1, minDist = 0){

  #Load data
  load("inst/extdata/GPP.RData")
  d = rbind(gpp[[1]],gpp[[2]], gpp[[3]])

  #Remove observation with NA distance
  d = d[-which(d$code==2 & is.na(d$distance) & d$species %in% c(2,66)),]

  #Remove observations of other whales species
  d[which((!d$species %in% c(2,66)) & d$code==2),]$code = 1

  #Remove incidental observations
  d$code[which(d$obs.type=="incidental" & d$code==2)] = 1

  #Distance given in km
  d$dist = d$dist/1000
  d$perpDist = d$perp.dist/1000

  #Remove integration points with very short distance.
  d = d[d$dist>minDist | d$code !=1,]

  #Fill in integration points
  done = FALSE
  counter = 1
  while(!done){
    if(d$code[counter]>0 & d$dist[counter]>maxDist){

      newLine = d[counter,]
      newLine$code = 1
      newLine$dist = newLine$dist/2
      newLine$lon = mean(d$lon[(counter-1): counter])
      newLine$lat = mean(d$lat[(counter-1): counter])
      d$dist[counter] = d$dist[counter]/2
      d = rbind(d[1:(counter-1),], newLine, d[counter:dim(d)[1],])
    }else{
      counter = counter + 1
    }
    if(counter == dim(d)[1]){
      done = TRUE
    }
  }

  return(d)
}

#' readArea
#' @description ...
#' @param file Path to file with correct format
#' @return area of interest
#' @export
#' @examples
readArea = function(conf,file = NULL){

  #TODO, how shall the user provide the area of interest?
  strata <- rgdal::readOGR('inst/extdata/survey_strata.geojson', verbose=F)
  strata <- strata[unlist(lapply(c('SS', 'AP', 'ESS'), function(x) match(x, strata$stratum))),]
  full.stratum <- rgeos::gUnion(rgeos::gUnion(strata[which(strata$stratum=='SS'),], strata[which(strata$stratum=='AP'),]), strata[which(strata$stratum=='ESS'),])
  utm.strata <- spTransform(strata, conf$UTMproj)
  utm.full.stratum <- spTransform(full.stratum, conf$UTMproj)

  predAreaUTM = utm.full.stratum

  return(predAreaUTM)
}




#' setData
#' @description Define data used internally based on data porvided and configuration file
#' @param d Raw data
#' @param predAreaUTM prediciton area
#' @param conf configuration applied
#' @return data used internally in package
#' @export
#' @examples
setData = function(d,predAreaUTM, conf){

  if(conf$detectionTrunc>0){
    ii = which(d$perpDist>conf$detectionTrunc & d$code==2)
    if(length(ii)>0){
      d$code[ii] = 1
    }
  }

  #Convert to UTM coordinates
  coords = as.matrix(d[,c("lon","lat")])
  obsLatLon = SpatialPoints(coords,proj4string = conf$LatLonProj)
  obsUTM = spTransform(obsLatLon, conf$UTMproj)
  d$utmx = obsUTM@coords[,1]
  d$utmy = obsUTM@coords[,2]

  #Define mid-points as coordinates for integration points
  integrationPointsUTM = as.matrix(d[,c("utmx","utmy")])
  for(i in 1:dim(integrationPointsUTM)[1]){
    if(d$code[i]==1){
      integrationPointsUTM[i,] = colMeans(integrationPointsUTM[(i-1):i,])
    }
  }

  #Accomondate for observer sides
  d$abortLeft = rep(0,dim(d)[1])
  d$abortRight = rep(0,dim(d)[1])
  d$abortLeft[is.na(d$Port)] =1
  d$abortRight[is.na(d$Stbd)] =1


  #Define aundance integration points
  predPoints = sp::spsample(predAreaUTM,n = conf$nPredPoints , type = "regular",offset = 0.5) #Offset is otherwise random (shifts the grid)
  meshPoints = rbind(predPoints@coords,data.frame(x1 = d$utmx,x2 = d$utmy))
  predData = data.frame(predPoints@coords)
  names(predData) = c("utmx", "utmy")


  #Define mesh
  #Make 300km buffer polygon around full stratum
  full.buffer <- rgeos::gBuffer(predAreaUTM, width=100, joinStyle='ROUND')
  boundary <- as.inla.mesh.segment(full.buffer)

  mesh <- inla.mesh.2d(meshPoints,
                       max.edge = conf$spdeDetails$max.edge ,
                       cutoff = conf$spdeDetails$cutoff,offset = c(1,-0.1),
                       boundary=boundary)

  plot(mesh)
  mesh$n



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
  AalongLines = inla.spde.make.A(mesh,integrationPointsUTM)
  Apred = inla.spde.make.A(mesh,as.matrix(predData[,1:2]))
  spde = inla.spde2.matern(mesh, alpha=2)
  spdeMatrices = spde$param.inla[c("M0","M1","M2")]

  AObs = inla.spde.make.A(mesh,obsUTM[d$code==2,])


  data = list(AalongLines = AalongLines,
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
              area = predAreaUTM@polygons[[1]]@area,
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
              useMMPP = conf$mmpp[1]
  )

  attributes(data)$mesh = mesh
  attributes(data)$obsUTM = obsUTM
  attributes(data)$obsLatLon = obsLatLon
  attributes(data)$predData = predData
  attributes(data)$detectionCov = colnames(X_g_right)

  return(data)
}

