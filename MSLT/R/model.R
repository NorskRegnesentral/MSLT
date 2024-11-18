#' MSLT model
#' @description MSLT model 
#' @param par parameters in model 
#' @return Liklihood of MSLT model given parameters
#' @export
mslt = function(par){
  RTMB::getAll(par,data)
  sigma = exp(log_sigma);
  kappa = exp(log_kappa);
  mu = exp(log_mu);
  c_mmpp = exp(log_c_mmpp)+1;
  bHazard = exp(logB); 
  
  nll = 0.0;
  "[<-" = RTMB::ADoverload("[<-") #Problem with R's byte compiler.
  
  #Precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  Q_intensity = kappa[1]^4 * spdeMatrices$M0 + 2 * kappa[1]^2 * spdeMatrices$M1 + spdeMatrices$M2
  
  if(matern_intensity != 0){
    nll = nll - RTMB::dgmrf(x_intensity,0,Q_intensity, TRUE);
    if(usePCpriors==1){
      d = 2; #Part of spatial pc-prior
      R = -log(pcPriorsRange_intensity[2])*pcPriorsRange_intensity[1]^(d/2)
      S = -log(pcPriorsSD_intensity[2])/pcPriorsSD_intensity[1];
      rhoP = sqrt(8)/kappa[1];
      nll = nll- log( d/2 * R *S *rhoP^(-1-d/2)* exp(-R* rhoP^(-d/2) -S* sigma[1] )); #pc-prior contribution
    }
  }
  
  scaleS = 1 /((4*base::pi)*kappa[1]*kappa[1]); #needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
  x_intensity = x_intensity/sqrt(scaleS)*sigma[1];
    
  #----- Classical line transect likelihood----------
  nObs = length(distObs)
  for(i in 1:nObs){
    z_i = X_g_obs[i,];
    y_i = distObs[i];
    # We do NOT include normalizing constant "eshw" (included later)
    if(g_function==1){
      nll = nll - Term(log(g_half_normal(y_i,z_i,beta_g)));
    }else if(g_function==2){
      nll = nll - Term(log(g_hazard(y_i,z_i,beta_g, bHazard)));
    }
  }
  
  #------ Log Gaussian Cox process ----------
  #Structure needed for effort constribution
  Z_transect =  X_z*beta_z+  AalongLines%*%x_intensity;
  Z_transect_endpoints =  X_z*beta_z+  AalongLinesEndpoints%*%x_intensity; 
  
  # Transect codes:
  # 0: start transect leg
  # 1: node used for numerical integration (including change g(y) or end transect)
  # 2: detection of animal
  #Likelihood contribution from effort-----------------
  P = matrix(0,nrow = 1,ncol = 2) 
  for(i in 1:length(Z_transect)){
    tmp_left = X_g_left[i,];
    tmp_right = X_g_right[i,];
    esw = 0;  # Effective strip width
    
    if(code[i]!=0){
      if(g_function==1){
        if(abortLeft[i]==0){
          esw = eshw_half_normal(tmp_left,beta_g, detectionTrunc);
        }
        if(abortRight[i]==0){
          esw = esw +  eshw_half_normal(tmp_right,beta_g, detectionTrunc);
        }
      }else if(g_function==2){
        if(abortLeft[i]==0){
          esw = eshw_hazard(tmp_left,beta_g, bHazard,detectionTrunc);
        }
        if(abortRight[i]==0){
          esw = esw +  eshw_hazard(tmp_right,beta_g, bHazard,detectionTrunc);
        }
      }
    }
    
    if(useMMPP==1){
      # MMPP quantities
      lambda = c(0,0)
      lambda[1] = exp(Z_transect[i])*esw;  # Low Poisson rate
      lambda[2] = c_mmpp*lambda[1];        # High Poisson rate
      
      if(code[i]==0){
        P[1,1] = mu[2]/sum(mu); # // Eq. (3) in (2006)
        P[1,2] = mu[1]/sum(mu); # // Eq. (3) in (2006)
      }else if(code[i]==1){# Not having detected anything since last node
        P = P%*%expm_mmpp2(lambda,mu,lineIntegralDelta[i]);
      }else if(code[i]==2){# Observation node
        P = P%*%expm_mmpp2(lambda,mu,lineIntegralDelta[i]);
        P[1,1] = P[1,1]* exp(Z_transect_endpoints[i]);
        P[1,2] = P[1,2]* c_mmpp*exp(Z_transect_endpoints[i]);
      }
      if(code[i]!=0){
        nll = nll- Term(log(sum(P)));
        P = P/sum(P); #  Renormalize probability after event.
      }
    }else{
      if(code[i]>0){
        lambda = exp(Z_transect[i])*esw;  # Low Poisson rate
        nll = nll + Term(lambda*lineIntegralDelta[i]);
      }
      if(code[i]==2){
        nll = nll - Term(Z_transect_endpoints[i]); 
      }
    }
  }
  

  
  pi_1 = mu[2]/sum(mu);
  pi_2 = mu[1]/sum(mu);
  k_psi = pi_1 + c_mmpp*pi_2;
  RTMB::ADREPORT(k_psi);
  #Penalize
  if(penalize[1]==1){
    nll = nll- RTMB::dexp(c_mmpp,penalize[2],TRUE);
  }
  
  
  #Abundance
  linPred =   exp(X_z_pred*beta_z + Apred%*%x_intensity) * k_psi;
  
  #Size part
  if(applyPodSize==1){
    if(matern_size != 0){
      Q_size = kappa[2]^4 * spdeMatrices$M0 + 2 * kappa[2]^2 * spdeMatrices$M1 + spdeMatrices$M2
      nll = nll - RTMB::dgmrf(x_size,0,Q_size,TRUE);
      scaleS_size = 1/(4*base::pi*kappa[2]*kappa[2]); #Needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
      x_size = x_size/sqrt(scaleS_size) *sigma[2];
    }
    
    if(independentPodSize==1){
      linPredSize =  exp(beta_size[1] + AObs%*%x_size);
      linPredSizePredicion =   exp(beta_size[1] + Apred%*%x_size);
    }else{
      #TODO: Dependence between spatial effect of group size and intensity
    }
    for(i in 1:nObs){
      if(podSizeDist==1){
        nll = nll - Term(dpois(size[i],linPredSize[i], TRUE) - log(1-dpois(0,linPredSize[i])));
      }else{
        sizeNB = exp(logSizeNB[1]);
        varNB = linPredSize[i] + linPredSize[i]*linPredSize[i]/sizeNB;
        nll = nll - Term((RTMB::dnbinom2(size[i], linPredSize[i],varNB ,log = TRUE)- log(1-RTMB::dnbinom2(0,linPredSize[i],varNB, log = FALSE))));
      }
    }
    for(i in 1:length(linPred)){
      if(podSizeDist==1){
        linPred[i] = linPred[i]*linPredSizePredicion[i]/(1-exp(-linPredSizePredicion[i]));
      }else{
        sizeNB = exp(logSizeNB[1]);
        varNB = linPredSizePredicion[i] + linPredSizePredicion[i]*linPredSizePredicion[i]/sizeNB;
        #Something not working with using dnbinom2, write the zero probability manually
        #linPred[i] = linPred[i]*(linPredSizePredicion[i]/(1-RTMB::dnbinom2(0, linPredSizePredicion[i], varNB,log = FALSE)));
        p = linPredSizePredicion[i]/ varNB
        r = linPredSizePredicion[i]^2 / (varNB-linPredSizePredicion[i])
        linPred[i] = linPred[i]*(linPredSizePredicion[i]/(1-p^(r))); 
      }
    }
  }
  
  if(spatialBiasCorFigure==1){ #Needed when producing spatial bias corrected plots in paper, removed by defauls because requires a couple of minutes computation time
     linPredFigure =   exp(beta_z[1] + x_intensity) * k_psi;
     linPredFigureSize =  exp(beta_size[1] + x_size);
    for(i in 1:length(linPredFigure)){
      sizeNB = exp(logSizeNB[1]);
      varNB = linPredFigureSize[i] + linPredFigureSize[i]*linPredFigureSize[i]/sizeNB;
      #Something not working with using dnbinom2, write the zero probability manually
      #linPredFigure[i] = linPredFigure[i]*(linPredFigureSize[i]/(1-RTMB::dnbinom2(0, linPredFigureSize[i], varNB,log = FALSE)));
      p = linPredFigureSize[i]/ varNB
      r = linPredFigureSize[i]^2 / (varNB-linPredFigureSize[i])
      linPredFigure[i] = linPredFigure[i]*(linPredFigureSize[i]/(1-p^(r))); 
    }
     RTMB::ADREPORT(linPredFigure); 
  }
  
  
  abundance = sum(areas*linPred);
  logAbundance = log(abundance);
  RTMB::ADREPORT(logAbundance);
  
  log_range_psi = log(-log(0.1)/sum(mu));
  RTMB::ADREPORT(log_range_psi);
  
  if(meanGroupFigure==1){
    if(applyPodSize==1){
      sizeNB = exp(logSizeNB[1]);
      meanGroup = exp(beta_size[1]+ x_size);
      for(i in 1: length(x_size)){
        varNB = meanGroup[i] + meanGroup[i]*meanGroup[i]/sizeNB;
        #Something not working with using dnbinom2, write the zero probability manually
        #meanGroup[i] = meanGroup[i]/(1-RTMB::dnbinom2(0, meanGroup[i], varNB));
        p = meanGroup[i]/ varNB
        r = meanGroup[i]^2 / (varNB-meanGroup[i])
        meanGroup[i] = meanGroup[i]/(1-p^(r)); 
      }
      RTMB::ADREPORT(meanGroup);
    }
  }
  
  #Adreport esw 
  logEswReport = rep(0,length(beta_g))
  for(i in 1:length(beta_g)){
    tmp = rep(0,length(beta_g))
    tmp[i]=1
    if(g_function==1){
      logEswReport[i] = log(eshw_half_normal(tmp,beta_g, detectionTrunc));
    }else if(g_function ==2){
      logEswReport[i] =   log(eshw_hazard(tmp,beta_g, bHazard, detectionTrunc));
    }
  }
  RTMB::ADREPORT(logEswReport);
  
  nll
}
