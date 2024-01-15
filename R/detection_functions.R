#' Half normal detection function
#' @description Half normal detection function
#' @param y perpendicular distance
#' @param z row from (covariate) design matrix
#' @param beta parameter vector (all parameters needed for g)
#' @return Half normal detection density
#' @export
#' @examples
g_half_normal = function(y, z, beta){
  sigma = exp(sum(beta*z));
  g =  dnorm(y,0,sigma,FALSE)*sqrt(2*3.141592)*sigma;
  g
}

#' Half normal effective half width
#' @description Half normal effective half width
#' @param z row from (covariate) design matrix
#' @param beta parameter vector (all parameters needed for g)
#' @param w truncation distance
#' @return Half normal effective width
#' @export
#' @examples
eshw_half_normal = function(z, beta,  w){
  sigma = exp(sum(beta*z));
  if(w>0){
    eshw = sqrt(2*3.141592)*sigma*(pnorm(w/sigma,0,1) -0.5) ;
  }else{
    eshw = sqrt(2*3.141592)*sigma*0.5;
  }
  eshw
}


#' Hazard detection function
#' @description Hazard detection function
#' @param y perpendicular distance
#' @param z row from (covariate) design matrix
#' @param beta parameter vector (all parameters needed for g)
#' @param b parameter in hazard detection
#' @return Hazard detection density
#' @export
#' @examples
g_hazard = function(y, z, beta, b){
  sigma = exp(sum(beta*z));
  g =  1 - exp(-(y/sigma)^(-b));
  g;
}


#' Hazard effective half width
#' @description Hazard effective half width
#' @param z row from (covariate) design matrix
#' @param beta parameter vector (all parameters needed for g)
#' @param b parameter in hazard detection function
#' @param w truncation distance
#' @return Hazard effective width
#' @export
#' @examples
eshw_hazard = function(z, beta, b,w){
  sigma = exp(sum(beta*z));
  
  if(w>0){
    eshw = exp(lgamma((b-1)/b))*sigma;
 #   dTrunc = (0.5*w)/20
 #   truncHelper = seq(w + dTrunc/2, w*1.5 + dTrunc/2, by = dTrunc)
 #   truncSub = sum(g_hazard(truncHelper,z,beta,b)*dTrunc)
 #   eshw = eshw- truncSub
  }else{
    eshw = exp(lgamma((b-1)/b))*sigma;
  }
  
  eshw;
}


