#' plot_spatial
#' @description Plot spatial and spatio-temporal effect
#' @param fit Estimated MSLT model returned by \code{mslt} or \code{fitMSLT_fallback_ridge}
#' @param year year to plot
#' @param sd logical, if true the spatial standard deviation is plotted
#' @param spatial_stationary logical, if true the stationary spatial field is included
#' @param spatio_temporal logical, if true the spatio-temporal field is included
#' @return Plot of spatial and spatio-temporal effect
#' @export
plot_spatial = function(fit,year = NULL, sd = FALSE,spatial_stationary = TRUE,spatio_temporal = TRUE){

  panels_utm = attributes(fit$data)$panels_utm
  survey_year = attributes(fit$data)$surveyData
  survey_year = survey_year[survey_year$year==year,]
  
  #Linear Predictor
  scale = 1/((4*pi)*exp(fit$pl$log_kappa_intensity)*exp(fit$pl$log_kappa_intensity))
  if(sd){
    linPred_year =  fit$pl$x_intensity_S*0 
    if(spatial_stationary){
      linPred_year = linPred_year + fit$plSd$x_intensity_S/sqrt(scale[1])*exp(fit$plSd$log_sigma_intensity[1]) 
    }
    if(spatio_temporal){
      linPred_year = linPred_year + fit$plSd$x_intensity_ST[,year - min(fit$data$year) + 1]/sqrt(scale[2])*exp(fit$plSd$log_sigma_intensity[2])
    }
    main = "Spatial standard deviation of MAP estiamte in"
  }else{
    linPred_year =  fit$pl$beta_z  + fit$pl$x_intensity_S*0 ########TODO: Include covariates##########
    if(spatial_stationary){
      linPred_year = linPred_year + fit$pl$x_intensity_S/sqrt(scale[1])*exp(fit$pl$log_sigma_intensity[1]) 
    }
    if(spatio_temporal){
      linPred_year = linPred_year + fit$pl$x_intensity_ST[,year - min(fit$data$year) + 1]/sqrt(scale[2])*exp(fit$pl$log_sigma_intensity[2])
    }
    main = "Spatial MAP estimate in"
  }

  main = paste0(main, " year ", year)
  
  zlim = range(linPred_year)
  mesh = attributes(fit$data)$mesh
  levels = round(seq(zlim[1],zlim[2],by = 0.4),1)
  
  # akima::interp is sensitive to large coordinate values; center coordinates first
  x <- mesh$loc[,1] - mean(mesh$loc[,1])
  y <- mesh$loc[,2] - mean(mesh$loc[,2])
  interp_surface <- akima::interp(x,y,linPred_year,nx = 100,ny = 100)
  interp_surface$x <- interp_surface$x + mean(mesh$loc[,1])
  interp_surface$y <- interp_surface$y + mean(mesh$loc[,2])

  fields::image.plot(interp_surface,
                     col=viridis::turbo(500),
                     xlab = '', ylab = '',zlim = zlim,
                     main = "", useRaster=T, horizontal=T)
  contour(interp_surface ,add = TRUE,
          labcex  = 0.8,lwd = 2,levels = levels, col=rgb(1,1,1,0.7))
  title(main = main, line = 1, cex.main=1.2, xpd=NA)
  title(xlab = "Easting (km)", line = 2.5, cex.lab=1.3, xpd=NA)
  title(ylab = "Northing (km)", line = 2.5, cex.lab=1.3, xpd=NA)
  
  plot(sf::st_geometry(panels_utm), add = TRUE, border = "black",lwd = 2)
  points(survey_year$utmx,survey_year$utmy,cex = 0.3,col = 'black')
  points(survey_year$utmx[survey_year$code==2],survey_year$utmy[survey_year$code==2],cex = 0.6,col = 'purple',pch = 20)
  
  # World coastline/countries
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  norway_utm <- sf::st_transform(world, fit$conf$UTMproj)
  
  plot(
    sf::st_geometry(norway_utm),
    col = "grey90",
    border = "black",
    add = TRUE
  )
}



