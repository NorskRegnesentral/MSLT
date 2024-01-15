#' MSLT model
#' @description MSLT model 
#' @param par parameters in model 
#' @return Liklihood of MSLT model given parameters
#' @export
#' @examples
mslt = function(par){
  getAll(par,data)
  sigma = exp(log_sigma);
  kappa = exp(log_kappa);
  mu = exp(log_mu);
  c_mmpp = exp(log_c_mmpp)+1;
  b = exp(logB); 
  
  nll = 0.0;
  Q_intensity = Q_spde(spdeMatrices,kappa[1]);
  REPORT(Q_intensity);
  if(matern_intensity != 0){
    nll = nll - dgmrf(x_intensity,0,Q_intensity, TRUE);
    if(usePCpriors==1){
      d = 2; #Part of spatial pc-prior
      R = -log(pcPriorsRange_intensity[2])*pcPriorsRange_intensity[1]^(d/2)
      S = -log(pcPriorsSD_intensity[2])/pcPriorsSD_intensity[1];
      rhoP = sqrt(8)/kappa[1];
      nll = nll- log( d/2 * R *S *rhoP^(-1-d/2)* exp(-R* rhoP^(-d/2) -S* sigma[1] )); #pc-prior contribution
    }
  }
  
  scaleS = 1 /((4*3.141592)*kappa[1]*kappa[1]); #needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
  x_intensity = x_intensity/sqrt(scaleS)*sigma[1];
  
  #----- Classical line transect likelihood----------
  nObs = length(distObs)
  for(i in 1:nObs){
    z_i = X_g_obs[i,];
    y_i = distObs[i];
    # We do NOT include normalizing constant "eshw" (included later)
    if(g_function==1){
      nll = nll - log(g_half_normal(y_i,z_i,beta_g));
    }else if(g_function==2){
      nll = nll - log(g_hazard(y_i,z_i,beta_g, b));
    }
  }
  
  #------ Log Gaussian Cox process ----------
  #Structure needed for effort constribution
  Z_transect =  X_z*beta_z+  AalongLines%*%x_intensity;
  
    
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
          esw = eshw_hazard(tmp_left,beta_g, b,detectionTrunc);
        }
        if(abortRight[i]==0){
          esw = esw +  eshw_hazard(tmp_right,beta_g, b,detectionTrunc);
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
      }else if(code[i]==1){
        # Not having detected anything since last node
        P = P%*%expm_mmpp2(lambda,mu,lineIntegralDelta[i]);
      }else if(code[i]==2){
        # Observation node
        P = P%*%expm_mmpp2(lambda,mu,lineIntegralDelta[i]);
        P[1,1] = P[1,1]*lambda[1];
        P[1,2] = P[1,2]*lambda[2];
        nll = nll+log(esw);  #This is the normalizing constant that was skipped above in the detection function
      }
      if(code[i]!=0){
        nll = nll- log(sum(P));
        P = P/sum(P); #  Renormalize probability after event.
      }
    }else{
      if(code[i]>0){
        lambda = exp(Z_transect[i])*esw;  # Low Poisson rate
        nll = nll + lambda*lineIntegralDelta[i];
      }
      if(code[i]==2){
        nll = nll - Z_transect[i];
      }
    }
  }

  pi_1 = mu[2]/sum(mu);
  pi_2 = mu[1]/sum(mu);
  k_psi = pi_1 + c_mmpp*pi_2;
  ADREPORT(k_psi);
  #Penalize
  if(penalize[1]==1){
    nll = nll- dexp(c_mmpp,penalize[2],TRUE);
  }
  
    
  #Abundance
  linPred =   exp(X_z_pred*beta_z + Apred%*%x_intensity) * k_psi;
  
  #Size part
  if(applyPodSize==1){
    if(matern_size != 0){
      Q_size = Q_spde(spdeMatrices,kappa[2]);
      nll = nll - dgmrf(x_size,0,Q_size,TRUE);
      scaleS_size = 1/(4*3.141592*kappa[2]*kappa[2]); #Needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
      x_size = x_size/sqrt(scaleS_size) *sigma[2];
    }
    
    if(independentPodSize==1){
      linPredSize =  exp(beta_size[1] + AObs%*%x_size);
      linPredSizePredicion =   exp(beta_size[1] + Apred%*%x_size);
    }else{
      #TODO
      #     linPredSize =  exp(beta_size(0) + beta_size(1)*(AObs*x_intensity));
      #    linPredSizePredicion =   exp(beta_size(0) + beta_size(1)*(Apred*x_intensity));
    }
    for(i in 1:nObs){
      if(podSizeDist==1){
        nll = nll - dpois(size[i],linPredSize[i], TRUE) - log(1-dpois(0,linPredSize[i]));
      }else{
        sizeNB = exp(logSizeNB[1]);
        varNB = linPredSize[i] + linPredSize[i]*linPredSize[i]/sizeNB;
        nll = nll - (dnbinom2(size[i], linPredSize[i],varNB ,log = TRUE)- log(1-dnbinom2(0,linPredSize[i],varNB, log = FALSE)));
      }
    }
    for(i in 1:length(linPred)){
      if(podSizeDist==1){
        linPred[i] = linPred[i]*linPredSizePredicion[i]/(1-exp(-linPredSizePredicion[i]));
      }else{
        sizeNB = exp(logSizeNB[1]);
        varNB = linPredSizePredicion[i] + linPredSizePredicion[i]*linPredSizePredicion[i]/sizeNB;
#       linPred[i] = linPred[i]*(linPredSizePredicion[i]/(1-dnbinom2(0, linPredSizePredicion[i], varNB, log = FALSE))); #Something not working with using dnbinom2 here, write the zero probability manually
        p = linPredSizePredicion[i]/ varNB
        r = linPredSizePredicion[i]^2 / (varNB-linPredSizePredicion[i])
        linPred[i] = linPred[i]*(linPredSizePredicion[i]/(1-p^(r))); 
      }
    }
  }
  
 
  
  abundance = area*sum(linPred)/length(linPred);
  logAbundance = log(abundance);
  ADREPORT(logAbundance);
  
  log_range_psi = log(-log(0.1)/sum(mu));
  ADREPORT(log_range_psi);
  
  if(applyPodSize==1){
    sizeNB = exp(logSizeNB[1]);
    linPredMeanSize = exp(beta_size[1]+ 0.5*sigma[2]*sigma[2]);
    varNB = linPredMeanSize + linPredMeanSize*linPredMeanSize/sizeNB;
    meanGroupSize = linPredMeanSize/(1-dnbinom2(0, linPredMeanSize, varNB, log = FALSE));
    ADREPORT(meanGroupSize);
  }

  nll
}
