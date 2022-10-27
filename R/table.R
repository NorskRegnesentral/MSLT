#' Return table with model parameters
#' 
#' @param run run to extract model parameters from
#' @param CI_level Confidence level to accompany the point estimates with
#' @return Table with model parameters
#' @export
partable_model<-function(run,CI_level=0.95){

  CI_L <- (1-CI_level)/2
  CI_U <- 1-CI_L

  map = run$map
  pl = as.list(run$rep,"Est")
  plSd = as.list(run$rep,"Std")

  sigma_int = c(exp(pl$log_sigma[1]),
                exp(pl$log_sigma[1] + qnorm(CI_L)*plSd$log_sigma[1]),
                exp(pl$log_sigma[1] + qnorm(CI_U)*plSd$log_sigma[1]))
  kappa_int = c(exp(pl$log_kappa[1]),
                exp(pl$log_kappa[1] + qnorm(CI_L)*plSd$log_kappa[1]),
                exp(pl$log_kappa[1] + qnorm(CI_U)*plSd$log_kappa[1]))
  sigma_size = c(exp(pl$log_sigma[2]),
                 exp(pl$log_sigma[2] + qnorm(CI_L)*plSd$log_sigma[2]),
                 exp(pl$log_sigma[2] + qnorm(CI_U)*plSd$log_sigma[2]))
  kappa_size = c(exp(pl$log_kappa[2]),
                 exp(pl$log_kappa[2] + qnorm(CI_L)*plSd$log_kappa[2]),
                 exp(pl$log_kappa[2] + qnorm(CI_U)*plSd$log_kappa[2]))
  mu1 = c(exp(pl$log_mu[1]),
          exp(pl$log_mu[1] + qnorm(CI_L)*plSd$log_mu[1]),
          exp(pl$log_mu[1] + qnorm(CI_U)*plSd$log_mu[1]))
  mu2 = c(exp(pl$log_mu[2]),
          exp(pl$log_mu[2] + qnorm(CI_L)*plSd$log_mu[2]),
          exp(pl$log_mu[2] + qnorm(CI_U)*plSd$log_mu[2]))
  c_mmpp = c(exp(pl$log_c_mmpp),
             exp(pl$log_c_mmpp + qnorm(CI_L)*plSd$log_c_mmpp),
             exp(pl$log_c_mmpp + qnorm(CI_U)*plSd$log_c_mmpp))
  betaz = c(pl$beta_z[1],
            pl$beta_z[1] + qnorm(CI_L)*plSd$beta_z[1],
            pl$beta_z[1] + qnorm(CI_U)*plSd$beta_z[1])
  betaG1 = c(pl$beta_g[1],
             pl$beta_g[1] + qnorm(CI_L)*plSd$beta_g[1],
             pl$beta_g[1] + qnorm(CI_U)*plSd$beta_g[1])
  betaG2 = c(pl$beta_g[2],
             pl$beta_g[2] + qnorm(CI_L)*plSd$beta_g[2],
             pl$beta_g[2] + qnorm(CI_U)*plSd$beta_g[2])
  betaG3 = c(pl$beta_g[3],
             pl$beta_g[3] + qnorm(CI_L)*plSd$beta_g[3],
             pl$beta_g[3] + qnorm(CI_U)*plSd$beta_g[3])

  #  B = round(c(pl$logB, pl$logB - 2*plSd$logB,pl$logB + 2*plSd$logB ),2)
  beta0_size = c(pl$beta_size[1],
                 pl$beta_size[1] + qnorm(CI_L)*plSd$beta_size[1],
                 pl$beta_size[1] + qnorm(CI_U)*plSd$beta_size[1])
  beta1_size = c(pl$beta_size[2],
                 pl$beta_size[2] + qnorm(CI_L)*plSd$beta_size[2],
                 pl$beta_size[2] + qnorm(CI_U)*plSd$beta_size[2])
  SizeNB = exp(c(pl$logSizeNB,
                pl$logSizeNB + qnorm(CI_L)*plSd$logSizeNB,
                pl$logSizeNB + qnorm(CI_U)*plSd$logSizeNB))

  parTab = rbind(sigma_int,kappa_int,mu1,mu2,c_mmpp,betaz,betaG1,betaG2,betaG3,sigma_size,kappa_size,beta0_size,beta1_size,SizeNB)
  parTab <- as.data.frame(parTab)
  colnames(parTab) <- c("Est","CI_L","CI_U")

  # Setting estimates to NA if not estimated
  if(!is.null(map$x)){
    parTab[rownames(parTab)=="sigma_int",] = NA
    parTab[rownames(parTab)=="kappa_int",] = NA
  }
  if(!is.null(map$log_mu)){
    parTab[rownames(parTab)=="mu1",] = NA
    parTab[rownames(parTab)=="mu2",] = NA
    parTab[rownames(parTab)=="c_mmpp",] = NA
  }
  if(is.na(map$log_sigma[1])){
    parTab[rownames(parTab)=="sigma_int",] = NA
    parTab[rownames(parTab)=="kappa_int",] = NA
  }

  return(parTab)
}


#' Return table with derived model parameters
#' 
#' @param run run to extract model parameters from
#' @param sdreport_run the sdreport object for the run to extract model parameters from
#' @param CI_level Confidence level to accompany the point estimates with
#' @return Table with derived model parameters
#' @export
partable_model_derived <-function(run,sdreport_run,CI_level=0.95){

  CI_L <- (1-CI_level)/2
  CI_U <- 1-CI_L

  map = run$map
  pl = as.list(run$rep,"Est")
  plSd = as.list(run$rep,"Std")
  rlBiasCorr = as.list(sdreport_run,"Est. (bias.correct)", report = TRUE)
  rlBiasCorrSd = as.list(sdreport_run,"Std. (bias.correct)", report = TRUE)


  range_int = c(sqrt(8)/exp(pl$log_kappa[1]),
            sqrt(8)/exp(pl$log_kappa[1] + qnorm(CI_U)*plSd$log_kappa[1]), # CI_U here as we divide by kappa
            sqrt(8)/exp(pl$log_kappa[1] + qnorm(CI_L)*plSd$log_kappa[1])) # CI_L here as we divide by kappa

  betaG1 = c(pl$beta_g[1],
             pl$beta_g[1] + qnorm(CI_L)*plSd$beta_g[1],
             pl$beta_g[1] + qnorm(CI_U)*plSd$beta_g[1])
  betaG2 = c(pl$beta_g[2],
             pl$beta_g[2] + qnorm(CI_L)*plSd$beta_g[2],
             pl$beta_g[2] + qnorm(CI_U)*plSd$beta_g[2])
  betaG3 = c(pl$beta_g[3],
             pl$beta_g[3] + qnorm(CI_L)*plSd$beta_g[3],
             pl$beta_g[3] + qnorm(CI_U)*plSd$beta_g[3])

  half_stripe_widthG1 <- sqrt(2*pi*exp(betaG1))/2
  half_stripe_widthG2 <- sqrt(2*pi*exp(betaG2))/2
  half_stripe_widthG3 <- sqrt(2*pi*exp(betaG3))/2

  range_psi = c(rlBiasCorr$range_psi,
                rlBiasCorr$range_psi + qnorm(CI_L)*rlBiasCorrSd$range_psi,
                rlBiasCorr$range_psi + qnorm(CI_U)*rlBiasCorrSd$range_psi)

  k_psi = c(rlBiasCorr$k_psi,
                rlBiasCorr$k_psi + qnorm(CI_L)*rlBiasCorrSd$k_psi,
                rlBiasCorr$k_psi + qnorm(CI_U)*rlBiasCorrSd$k_psi)

  range_size = c(sqrt(8)/exp(pl$log_kappa[2]),
                 sqrt(8)/exp(pl$log_kappa[2] + qnorm(CI_U)*plSd$log_kappa[2]), # CI_U here as we divide by kappa
                 sqrt(8)/exp(pl$log_kappa[2] + qnorm(CI_L)*plSd$log_kappa[2])) # CI_L here as we divide by kappa

  meanGroupSize = c(rlBiasCorr$meanGroupSize,
                rlBiasCorr$meanGroupSize + qnorm(CI_L)*rlBiasCorrSd$meanGroupSize,
                rlBiasCorr$meanGroupSize + qnorm(CI_U)*rlBiasCorrSd$meanGroupSize)


  abundance = c(exp(rlBiasCorr$logAbundance),
                exp(rlBiasCorr$logAbundance + qnorm(CI_L)*rlBiasCorrSd$logAbundance),
                exp(rlBiasCorr$logAbundance + qnorm(CI_U)*rlBiasCorrSd$logAbundance))



  parTab = rbind(range_int,half_stripe_widthG1,half_stripe_widthG2,half_stripe_widthG3,range_psi,k_psi,range_size,meanGroupSize,abundance)
  parTab <- as.data.frame(parTab)
  colnames(parTab) <- c("Est","CI_L","CI_U")

  # Setting estimates to NA if not estimated
  if(!is.null(map$x)){
    parTab[rownames(parTab)=="range_int",] = NA
  }

  if(!is.null(map$log_mu)){
    parTab[rownames(parTab)=="range_psi",] = NA
    parTab[rownames(parTab)=="k_psi",] = NA
  }


  return(parTab)
}



##' Print abundance
##' @param  run
##' @details Print abundance
##' @export
logAbundance = function(run){
  rl = as.list(run$rep,"Est",report = TRUE)$logAbundance
}

