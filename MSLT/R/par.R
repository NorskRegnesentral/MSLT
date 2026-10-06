#' setPar
#' 
#' @param data data list returned by \code{\link{setData}}
#' @param conf configuration list returned by \code{\link{defConf}}
#' @return Initial model parameters
#' @export
setPar = function(data,conf){
  par = list(log_sigma =c(0.5,0.5),
             log_kappa = c(-5,-5),
             logit_rho = 0,
             beta_g = rep(0.2,dim(data$X_g_obs)[2]),
             beta_z = c(-6,rep(0,dim(data$X_z)[2]-1)),
             beta_size = 0,
             logB = 1,
             log_mu =c(-3,-3),
             log_c_mmpp = -1,
             logSizeNB = numeric(),
             x_intensity_S = rep(0,attributes(data)$mesh$n),
             x_intensity_ST =  matrix(0, nrow = attributes(data)$mesh$n, ncol = data$nYear),
             x_size = rep(0,attributes(data)$mesh$n))

  if(conf$applyPodSize==1){
    if(conf$podSizeDist == 2){#NB-distributed
      par$logSizeNB = 0
    }else{
      par$logSizeNB = numeric()
    }
    if(conf$independentPodSize==1){
      par$beta_size = 0
    }else{
      tt = stats::uniroot(function(x){
        mean(data$size) - x/(1-exp(-x))
      }, interval = c(0.01,10))$root
      par$beta_size = c(log(tt),0);
    }
  }else{
    par$beta_size = numeric()
    par$logSizeNB = numeric()
  }

  if(conf$mmpp!=1){
    par$log_c_mmpp = -10
  }


  return(par)
}
