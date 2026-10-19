

#' Fits the MSLT model. If the initial fit fails, the function iteratively applies
#' ridge correction to a subset of the transects and refits the model.
#'
#' Transects are selected iteratively as follows. First, the smallest eigenvalue
#' of the Hessian of the latent effects is identified. Next, the transects that
#' contribute most strongly in the direction of the corresponding eigenvector are identified.
#' 
#' @param data data list returned by \code{\link{setData}}
#' @param par parameter list  returned by \code{\link{setPar}}
#' @param conf model configurations list  returned by \code{\link{defConf}}
#' @param map map object used by TMB to couple and remove parameters
#' @param ridge.correct Activates ridge.correct in MakeADFun
#' @param nTransect_ridge Number of transects that are selected for ridge correction in each iteration
#' @param n_trials Maximum number of iterations.
#' @param ... additional parameters passed to sdreport 
#' @return fitted model
#' @export
fitMSLT_fallback_ridge = function(data,par,conf,map = setMap(conf, par),ridge.correct = FALSE,nTransect_ridge = 5,n_trials = 20,...){
  
  surveyData = attributes(data)$surveyData
  
  run = fitMSLT(data,par,conf, map = map, getJointPrecision = TRUE)
  sdrep = run$rep
  if(run$opt$convergence==0){
    return(run)
  }
  
  runs = list()
  done = FALSE
  i = 0
  par_names <- names(run$obj$env$par)
  ii_S  <- which(par_names == "x_intensity_S")- min(which(par_names == "x_intensity_S"))+1
  ii_ST <- which(par_names == "x_intensity_ST")- min(which(par_names == "x_intensity_S"))+1
  nMesh <- attributes(data)$mesh$n
  id_ridge_used = NULL
  while(!done) {
    if(i>0){
      sdrep = runs[[i]]$rep
    }
    i = i+1
    
    H <- sdrep$jointPrecision
    ii_ST_H= which(colnames(H) %in% c("x_intensity_S","x_intensity_ST"))
    H_ST = H[ii_ST_H,ii_ST_H]
    ev = eigen(H_ST, symmetric=TRUE)
    
    v <- ev$vectors[, length(ev$values)]  # minste egenverdi
    vS  <- v[ii_S]
    vST <- v[ii_ST]
    
    signal_S <- as.vector(
      data$AalongLines_S %*% vS
    )
    
    vST_mat <- matrix(
      vST,
      nrow = nMesh,
      ncol = data$nYear
    )
    
    signal_ST <- rep(0, nrow(data$AalongLines_S))
    
    for (t in 1:data$nYear) {
      idx <- which(
        data$year == min(data$year) + t - 1
      )
      signal_ST[idx] <-
        as.vector(
          data$AalongLines_ST[[t]] %*%
            vST_mat[, t]
        )
    }
    
    signal_total <- signal_S + signal_ST
    
    score <- tapply(
      signal_total^2,
      surveyData$transect_id,
      sum#mean
    )
    score = sort(score, decreasing = TRUE)
    
    if(length(id_ridge_used)>0){
      score = score[!(names(score) %in%id_ridge_used)]
    }
    print(head(score,10))
    
    problemTransects = names(score)[1:nTransect_ridge]
    id_ridge_used = c(id_ridge_used,problemTransects)
    for(id in problemTransects){
      if(!is.na(id)){
        data$ridgeCorrectLine[which(surveyData$transect==id)] = 1
        
        rawObs = surveyData[surveyData$code==2,]
        if(length(which(rawObs$transect==id))>0){
          data$ridgeCorrectObservations[which(rawObs$transect==id)] = 1
        }
      }
    }
    print(paste0("Number of ridge.correct points: " , sum(data$ridgeCorrectLine)  ))
    print("Transects with ridge.correct points:" )
    print(id_ridge_used)
    
    runs[[i]] = fitMSLT(data,par,conf,ridge.correct = TRUE,getJointPrecision = TRUE)
    
    if(runs[[i]]$opt$convergence==0 | i > n_trials){
      print(paste0("Converged in trial ",i,"."))
      done = TRUE
    }else{
      runs[[i]]$obj = NULL
      gc()
      print("Did not converge")
      print(i)
      print(runs[[i]]$rep)
    }
  }
  run_first = run
  run = runs[[i]]
  attributes(run)$run_first = run_first
  attributes(run)$n_trials = i 
  return(run)
}


