#' MSLT model
#' @description MSLT model 
#' @param par parameters in model 
#' @return Liklihood of MSLT model given parameters
#' @export
mslt = function(par,data){
  
  if(exists("shared_env_simulations")){
    data = shared_env_simulations$data
  }
  
  if(sum(data$ridgeCorrectLine)==0){#Investigating solving correlation issues using ridge.correct, This is work in progress
    Term = function(f)f 
  }else{
    Term = RTMB:::Term #Using ridge.correct for some of the observations
  }
  
  RTMB::getAll(par,data)
  sigma_intensity = exp(log_sigma_intensity);
  sigma_size = exp(log_sigma_size);
  kappa_intensity = exp(log_kappa_intensity);
  kappa_size = exp(log_kappa_size);
  rho_intensity = -1 + 2 * (1 / (1 + exp(-logit_rho_intensity)))
  rho_size = -1 + 2 * (1 / (1 + exp(-logit_rho_size)))
  #rho = -1 +2*plogis(logit_rho)
  mu = exp(log_mu);
  c_mmpp = exp(log_c_mmpp)+1;
  bHazard = exp(logB); 
  
  nll = 0.0;
  "[<-" = RTMB::ADoverload("[<-") #Problem with R's byte compiler.
  RTMB::TapeConfig(atomic="disable") #Makes %*% more computationally efficient
  
  #Precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  Q_intensity_S = kappa_intensity[1]^4 * spdeMatrices$M0 + 2 * kappa_intensity[1]^2 * spdeMatrices$M1 + spdeMatrices$M2
  Q_intensity_ST = kappa_intensity[2]^4 * spdeMatrices$M0 + 2 * kappa_intensity[2]^2 * spdeMatrices$M1 + spdeMatrices$M2
  
  if(space_intensity_matern != 0){
    nll = nll - RTMB::dgmrf(x_intensity_S,0,Q_intensity_S, TRUE);
    if(usePCpriors==1){
      d = 2; #Part of spatial pc-prior
      R = -log(pcPriorsRange_intensity[2])*pcPriorsRange_intensity[1]^(d/2)
      S = -log(pcPriorsSD_intensity[2])/pcPriorsSD_intensity[1];
      rhoP = sqrt(8)/kappa_intensity[1];
      nll = nll- log( d/2 * R *S *rhoP^(-1-d/2)* exp(-R* rhoP^(-d/2) -S* sigma_intensity[1] )); #pc-prior contribution
    }
  }
  
  if(spaceTime_intensity_maternAR1 != 0){
    ## GMRF prior spatio-temporal effect
    f1 <- function(x) {
      RTMB::dgmrf(x,
            mu = 0,
            Q = Q_intensity_ST,
            log = TRUE)
    }
    f2 <- function(x) {
      RTMB::dautoreg(x,
               mu = 0,
               phi = rho_intensity,
               log = TRUE)
    }
    nll <- nll - RTMB::dseparable(f1, f2)(x_intensity_ST)
  
    if(usePCpriors==1){
      d = 2; #Part of spatial pc-prior
      R = -log(pcPriorsRange_intensity[2])*pcPriorsRange_intensity[1]^(d/2)
      S = -log(pcPriorsSD_intensity[2])/pcPriorsSD_intensity[1];
      rhoP = sqrt(8)/kappa_intensity[1];
      nll = nll- log( d/2 * R *S *rhoP^(-1-d/2)* exp(-R* rhoP^(-d/2) -S* sigma_intensity[1] )); #pc-prior contribution
    }
  }
  
  scaleS = 1 /((4*base::pi)*kappa_intensity[1]*kappa_intensity[1]); #needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
  scaleST = 1 /((4*base::pi)*kappa_intensity[2]*kappa_intensity[2]); #needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
  x_intensity_S = x_intensity_S/sqrt(scaleS)*sigma_intensity[1];
  x_intensity_ST = x_intensity_ST/sqrt(scaleST)*sigma_intensity[2];

  
  #----- Classical line transect likelihood----------
  nObs = length(distObs)
  for(i in 1:nObs){
    z_i = X_g_obs[i,];
    y_i = distObs[i];
    # We do NOT include normalizing constant "eshw" (included later)
    if(ridgeCorrectObservations[i]==1){
      if(g_function==1){
        nll = nll - Term(log(g_half_normal(y_i,z_i,beta_g)));
      }else if(g_function==2){
        nll = nll - Term(log(g_hazard(y_i,z_i,beta_g, bHazard)));
      }
    }else{
      if(g_function==1){
        nll = nll - log(g_half_normal(y_i,z_i,beta_g));
      }else if(g_function==2){
        nll = nll - log(g_hazard(y_i,z_i,beta_g, bHazard));
      } 
    }
  }
  

  #------ Log Gaussian Cox process ----------
  #Structure needed for effort constribution
  Z_transect =  X_z*beta_z+  AalongLines_S%*%x_intensity_S; 
  Z_transect_endpoints = X_z*beta_z+ AalongLinesEndpoints_S%*%x_intensity_S; 
  if(spaceTime_intensity_maternAR1 != 0){
    for(t in 1:nYear){
      indeks = which(year== min(year) + t-1)
      Z_transect[indeks] =  Z_transect[indeks] +  AalongLines_ST[[t]]%*%x_intensity_ST[,t]; 
      Z_transect_endpoints[indeks] =  Z_transect_endpoints[indeks] +  AalongLinesEndpoints_ST[[t]]%*%x_intensity_ST[,t]; 
    }
  }

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
      
      lambda_endpoints = c(0,0)
      lambda_endpoints[1] = exp(Z_transect_endpoints[i])*esw;  # Low Poisson rate
      lambda_endpoints[2] = c_mmpp*lambda_endpoints[1];        # High Poisson rate
      
      if(code[i]==0){
        P[1,1] = mu[2]/sum(mu); # // Eq. (3) in (2006)
        P[1,2] = mu[1]/sum(mu); # // Eq. (3) in (2006)
      }else if(code[i]==1){# Not having detected anything since last node
        P = P%*%expm_mmpp2(lambda,mu,lineIntegralDelta[i]);
      }else if(code[i]==2){# Observation node
        P = P%*%expm_mmpp2(lambda,mu,lineIntegralDelta[i]);
        P[1,1] = P[1,1]* lambda_endpoints[1];
        P[1,2] = P[1,2]* lambda_endpoints[2];
        if(ridgeCorrectLine[i]==1){
          nll = nll -Term( -log(esw));  # This is the normalizing constant that was skipped above in the detection function
        }else{
          nll = nll +log(esw);  # This is the normalizing constant that was skipped above in the detection function
        }
      }
      if(code[i]!=0){
        if(ridgeCorrectLine[i]==1){
          nll = nll- Term(log(sum(P)));
        }else{
          nll = nll- log(sum(P));
        }
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
  if(spaceTime_intensity_maternAR1==0){
    linPred =   exp(X_z_pred*beta_z + Apred%*%x_intensity_S) * k_psi;
  }else{
    linPred = RTMB::AD(matrix(NA,nrow = nYear,ncol = dim(Apred)[1]))
    for(t in 1:nYear){
      linPred[t,] =    RTMB::AD(as.vector(exp(X_z_pred*beta_z + Apred%*%x_intensity_S +  Apred%*%x_intensity_ST[,t]))) * k_psi;
    }
  }

  
  
  
  #Size part
  if(applyPodSize==1){
    if(matern_size != 0){
      Q_size_S = kappa_size[1]^4 * spdeMatrices$M0 + 2 * kappa_size[1]^2 * spdeMatrices$M1 + spdeMatrices$M2
      nll = nll - RTMB::dgmrf(x_size,0,Q_size_S,TRUE);
      scaleS_size = 1/(4*base::pi*kappa_size[1]*kappa_size[1]); #Needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
      x_size = x_size/sqrt(scaleS_size) *sigma_size[1];
    }
    
    linPredSize =  exp(beta_size[1] + AObs_S%*%x_size);
    linPredSizePredicion =   exp(beta_size[1] + Apred%*%x_size);

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
        p = linPredSizePredicion[i]/ varNB
        r = linPredSizePredicion[i]^2 / (varNB-linPredSizePredicion[i])
        linPred[i] = linPred[i]*(linPredSizePredicion[i]/(1-p^(r))); 
      }
    }
  }
  if(spaceTime_intensity_maternAR1==0){
    abundance = sum(areas*linPred);
    logAbundance = log(abundance);
    RTMB::ADREPORT(logAbundance);
    
    panel_names = unique(cell_integration$panel)
    abundance_panels = RTMB::AD(rep(0,length(panel_names)))
    for(i in 1:dim(cell_integration)[1]){
      abundance_panels[which(panel_names==cell_integration$panel[i])] = abundance_panels[which(panel_names==cell_integration$panel[i])] +
        as.numeric(cell_integration$panel_area)[i]*as.vector(linPred)[cell_integration$cell_id[i]]
    }
    logAbundance_panels = log(abundance_panels)
    RTMB::ADREPORT(logAbundance_panels);
    RTMB::REPORT(panel_names)
  }else{
    abundance_vector = RTMB::AD(rep(0,nYear))
    for(t in 1:nYear){
      abundance_vector[t] = sum(areas*linPred[t,]);
    }
    logAbundance_vector = log(abundance_vector);
    RTMB::ADREPORT(logAbundance_vector);
    
    panel_names = unique(cell_integration$panel)
    abundance_panels_matrix = RTMB::AD(matrix(0,nrow = nYear,ncol = length(panel_names)))
    for(t in 1:nYear){
      for(i in 1:dim(cell_integration)[1]){
        abundance_panels_matrix[t,which(panel_names==cell_integration$panel[i])] = abundance_panels_matrix[t,which(panel_names==cell_integration$panel[i])] +
          as.numeric(cell_integration$panel_area)[i]*as.vector(linPred[t,])[cell_integration$cell_id[i]]
      }
    }
    logAbundance_panels_matrix = log(abundance_panels_matrix)
    RTMB::ADREPORT(logAbundance_panels_matrix);
    RTMB::REPORT(panel_names)
  }
  log_range_psi = log(-log(0.1)/sum(mu));
  RTMB::ADREPORT(log_range_psi);
  
  if(spatialBiasCorFigure==1){ #Needed when producing spatial bias corrected plots in paper, removed by defauls because requires a couple of minutes computation time
    linPredFigure =   exp(beta_z[1] + x_intensity_S) * k_psi;
    linPredFigureSize =  exp(beta_size[1] + x_size);
    for(i in 1:length(linPredFigure)){
      sizeNB = exp(logSizeNB[1]);
      varNB = linPredFigureSize[i] + linPredFigureSize[i]*linPredFigureSize[i]/sizeNB;
      p = linPredFigureSize[i]/ varNB
      r = linPredFigureSize[i]^2 / (varNB-linPredFigureSize[i])
      linPredFigure[i] = linPredFigure[i]*(linPredFigureSize[i]/(1-p^(r))); 
    }
    RTMB::ADREPORT(linPredFigure); 
  }
  
  if(meanGroupFigure==1){
    if(applyPodSize==1){
      sizeNB = exp(logSizeNB[1]);
      meanGroup = exp(beta_size[1]+ x_size);
      for(i in 1: length(x_size)){
        varNB = meanGroup[i] + meanGroup[i]*meanGroup[i]/sizeNB;
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
