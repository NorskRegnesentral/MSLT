#' Q_spde function
#' @description Set up precision matrix for SPDE-procedure
#' @param spdeMatricesS SPDE-matrices given by INLA needed to set up precision matrix 
#' @param kappa spatial scale parameter
#' @return Precision matrix for SPDE-procedure
#' @export
#' @examples
Q_spde <- function(spdeMatricesS, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spdeMatricesS$M0 + 2 * kappa_pow2 * spdeMatricesS$M1 + spdeMatricesS$M2    ## M0=G0, M1=G1, M2=G2
}


#' expm_mmpp2 function
#' @description  Direct implementation for the matrix exponential in dim 2 for MMPP
#' @param lamdba Vector with low and high intensity
#' @param mu parameters in MMPP
#' @param l distance moved allong transect
#' @return Matrix exponential for MMPP with two states
#' @export
#' @examples
expm_mmpp2 = function(lambda, mu, l){
  e = 1e-10;
  
  S = lambda+mu;
  K = prod(lambda) + lambda[1]*mu[2] + lambda[2]*mu[1];
  sumS = sum(S);
  U = sqrt(sumS*sumS-4*K + e);
  theta = c(0,0);
  theta[1] = (sumS+U)/2;
  theta[2] = (sumS-U)/2;
  
  M1 = matrix(0,nrow = 2,ncol=2)
  M2 = matrix(0,nrow = 2,ncol=2)
  M1[1,1] = S[2]-theta[2];
  M1[1,2] = mu[1];
  M1[2,1] = mu[2];
  M1[2,2] = S[1]-theta[2];
  
  M2 = M1;
  M2[1,1] = M2[1,1] + theta[2] - theta[1];
  M2[2,2] = M2[2,2] + theta[2] - theta[1];
  
  (exp(-theta[2]*l)*M1 - exp(-theta[1]*l)*M2)/U;
}
