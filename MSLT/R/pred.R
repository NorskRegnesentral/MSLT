#' plot_spatial
#' @description Plot spatial and spatio-temporal effect
#' @param fit Estimated MSLT model returned by \code{mslt} or \code{fitMSLT_fallback_ridge}
#' @param year_area_idx data frame that provides which year and area that we want to do bias-correction for
#' @param nSplit if year_area_idx not included. Split the biascorrect in nSplit splits. This is done to reduce memory usage. 
#' @param sd logical; Calculate bias corrected standard deviations if true
#' @return Plot of spatial and spatio-temporal effect
#' @export
pred_biascorrect = function(fit,year_area_idx = NULL,nSplit = 10,sd = FALSE){
  
  startTime <- Sys.time()
  
  if(is.null(year_area_idx)){
    idx = which(names(fit$rep$value) == "logAbundance_panels_matrix" | names(fit$rep$value) == "logAbundance_vector")
    splits <- split(idx, cut(seq_along(idx),
                           breaks = nSplit,
                           labels = FALSE))
    
    sdrep_biascorrect = RTMB::sdreport(fit$obj,bias.correct = TRUE,
                                       bias.correct.control = list(sd = sd, split = splits, nsplit = NULL))
  
  }else{
    idx = which(names(fit$rep$value) == "logAbundance_panels_matrix")
    panel_names = fit$obj$report()$panel_names
    logAbundance_panels_matrix = matrix(idx,nrow = fit$data$nYear,ncol = length(panel_names))
    colnames(logAbundance_panels_matrix) = colnames(fit$rl$logAbundance_panels_matrix)
    rownames(logAbundance_panels_matrix) = rownames(fit$rl$logAbundance_panels_matrix)
    
    split = list(logAbundance_panels_matrix[cbind(as.character(year_area_idx$year),year_area_idx$panel)])
    
    sdrep_biascorrect = RTMB::sdreport(fit$obj,bias.correct = TRUE,
                                       bias.correct.control = list(sd = sd, split = split, nsplit = NULL))
    
  }
  
  endTime = Sys.time()
  timeUsed = endTime - startTime
  print(timeUsed)
  
  rl_bias_correct = as.list(sdrep_biascorrect,what = "Est. (bias.correct)", report = TRUE)
  rlsd_bias_correct = as.list(sdrep_biascorrect,what = "Std. (bias.correct)", report = TRUE)
  
  colnames(rl_bias_correct$logAbundance_panels_matrix) = colnames(fit$rl$logAbundance_panels_matrix)
  rownames(rl_bias_correct$logAbundance_panels_matrix) = rownames(fit$rl$logAbundance_panels_matrix)
  names(rl_bias_correct$logAbundance_vector)= names(fit$rl$logAbundance_vector)
  colnames(rlsd_bias_correct$logAbundance_panels_matrix) = colnames(fit$rl$logAbundance_panels_matrix)
  rownames(rlsd_bias_correct$logAbundance_panels_matrix) = rownames(fit$rl$logAbundance_panels_matrix)
  names(rlsd_bias_correct$logAbundance_vector)= names(fit$rl$logAbundance_vector)
  
  toReturn = list(sdrep_biascorrect = sdrep_biascorrect,
                  rl_bias_correct = rl_bias_correct,
                  rlsd_bias_correct = rlsd_bias_correct)
  return(toReturn)
  
}
