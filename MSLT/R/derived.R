#' expm_mmpp2 function
#' @description  Direct implementation for the matrix exponential in dim 2 for MMPP
#' @param lamdba Vector with low and high intensity
#' @param mu parameters in MMPP
#' @param l distance moved allong transect
#' @return Matrix exponential for MMPP with two states
#' @export
expm_mmpp2 = function(lambda, mu, l){
  
  "[<-" = RTMB::ADoverload("[<-")#Problem with R's byte compiler.
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
