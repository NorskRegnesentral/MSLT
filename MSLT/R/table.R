#' Return table with model parameters
#' 
#' @param run run returned by \code{\link{fitMSLT}} run to extract model parameters from
#' @param CI_level Confidence level to accompany the point estimates with
#' @return Table with model parameters
#' @export
partable<-function(run,CI_level=0.95){

  CI_L <- (1-CI_level)/2
  CI_U <- 1-CI_L
  qq = c(qnorm(CI_L),qnorm(CI_U))
  pl = run$pl
  plSd = run$plSd

  sigma_int = exp(c(pl$log_sigma[1], pl$log_sigma[1] + qq*plSd$log_sigma[1]))
  log_kappa_int = c(pl$log_kappa[1], pl$log_kappa[1] + qq*plSd$log_kappa[1])
  sigma_size = exp(c(pl$log_sigma[2], pl$log_sigma[2] + qq*plSd$log_sigma[2]))
  log_kappa_size = c(pl$log_kappa[2], pl$log_kappa[2] + qq*plSd$log_kappa[2])
  log_mu1 = c(pl$log_mu[1], pl$log_mu[1] + qq*plSd$log_mu[1])
  log_mu2 = c(pl$log_mu[2], pl$log_mu[2] + qq*plSd$log_mu[2])
  c_mmpp = c(exp(pl$log_c_mmpp)+1, exp(pl$log_c_mmpp +  qq*plSd$log_c_mmpp)+1)
  betaz = c(pl$beta_z[1], pl$beta_z[1] +  qq*plSd$beta_z[1])
  
  betaG = list()
  for(i in 1:length(pl$beta_g)){
    betaG[[i]] = c(pl$beta_g[i], pl$beta_g[i] + qq*plSd$beta_g[i])
  }

  if(run$conf$g_function==2){
    b = exp(c(pl$logB, pl$logB + qq*plSd$logB ))
  }else{
    b = c(NA,NA,NA)
  }
    
  beta0_size = c(pl$beta_size[1],pl$beta_size[1] + qq*plSd$beta_size[1])
  beta1_size = c(pl$beta_size[2],pl$beta_size[2] + qq*plSd$beta_size[2])
  log_size_NB = c(pl$logSizeNB,pl$logSizeNB + qq*plSd$logSizeNB)
  
  parTabTmp = rbind(sigma_int,log_kappa_int,log_mu1,log_mu2,c_mmpp,betaz)
  parTab = rbind(parTabTmp,do.call(rbind,betaG))
  rownames(parTab) = c(rownames(parTabTmp),paste0("betaG",1:length(pl$beta_g)))
  parTab = rbind(parTab,b,sigma_size,log_kappa_size,beta0_size,beta1_size,log_size_NB)
  parTab <- as.data.frame(parTab)
  colnames(parTab) <- c("Est","CI_L","CI_U")

  # Setting estimates to NA if not estimated
  if(!is.null(run$map$x)){
    parTab[rownames(parTab)=="sigma_int",] = NA
    parTab[rownames(parTab)=="log_kappa_int",] = NA
  }
  if(!is.null(run$map$log_mu)){
    parTab[rownames(parTab)=="log_mu1",] = NA
    parTab[rownames(parTab)=="log_mu2",] = NA
    parTab[rownames(parTab)=="c_mmpp",] = NA
  }
  if(is.na(run$map$log_sigma[1])){
    parTab[rownames(parTab)=="sigma_int",] = NA
    parTab[rownames(parTab)=="log_kappa_int",] = NA
  }

  return(parTab)
}


#' Return table with derived model parameters
#' 
#' @param run run returned by \code{\link{fitMSLT}} to extract model parameters from
#' @param sdrep_bc the biascorrected sdreport object for the run to extract model parameters from
#' @param CI_level Confidence level to accompany the point estimates with
#' @return Table with derived model parameters
#' @export
partable_derived <-function(run,sdrep_bc = NULL,CI_level=0.95){
  CI_L <- (1-CI_level)/2
  CI_U <- 1-CI_L
  qq = c(qnorm(CI_L),qnorm(CI_U))
  
  range_int = c(sqrt(8)/exp(run$pl$log_kappa[1]),sqrt(8)/exp(run$pl$log_kappa[1] +  rev(qq)*run$plSd$log_kappa[1]))
  
  half_stripe_widthG = list()
  for(i in 1:length(run$pl$beta_g)){
    half_stripe_widthG[[i]] = exp(c(run$rl$logEswReport[i], run$rl$logEswReport[i] + qq*run$rlSd$logEswReport[i]))
  }

  range_psi = exp(c(run$rl$log_range_psi[1], run$rl$log_range_psi[1]+ qq*run$rlSd$log_range_psi[1]))
  k_psi = c(run$rl$k_psi, run$rl$k_psi[1] + qq*run$rlSd$k_psi[1])
  range_size = c(sqrt(8)/exp(run$pl$log_kappa[2]),sqrt(8)/exp(run$pl$log_kappa[2] + rev(qq)*run$plSd$log_kappa[2]))
  meanGroupSize = c(run$rl$meanGroupSize[1],run$rl$meanGroupSize[1] +qq*run$rlSd$meanGroupSize[1])
  abundance = exp(c(run$rl$logAbundance[1],run$rl$logAbundance[1] + qq*run$rlSd$logAbundance[1]))
  if(!is.null(sdrep_bc)){
    rlBiasCorr = as.list(sdrep_bc,"Est. (bias.correct)", report = TRUE)
    rlBiasCorrSd = as.list(sdrep_bc,"Std. (bias.correct)", report = TRUE)
    abundanceBiasCorrected = exp(c(rlBiasCorr$logAbundance[1], rlBiasCorr$logAbundance[1] + qq*rlBiasCorrSd$logAbundance[1]))
  }else{
    abundanceBiasCorrected = c(NA,NA,NA)
  }

  parTabTmp = rbind(range_int)
  parTab = rbind(parTabTmp,do.call(rbind,half_stripe_widthG))
  rownames(parTab) = c(rownames(parTabTmp),paste0("half_stripe_widthG",1:length(run$pl$beta_g)))
  parTab = rbind(parTab,range_psi,k_psi,range_size,meanGroupSize,abundance,abundanceBiasCorrected)
  parTab <- as.data.frame(parTab)
  
  colnames(parTab) <- c("Est","CI_L","CI_U")
  # Setting estimates to NA if not estimated
  if(!is.null(run$map$x)){
    parTab[rownames(parTab)=="range_int",] = NA
  }
  if(!is.null(run$map$log_mu)){
    parTab[rownames(parTab)=="range_psi",] = NA
    parTab[rownames(parTab)=="k_psi",] = NA
  }
  return(parTab)
}

