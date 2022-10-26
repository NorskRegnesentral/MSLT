#define TMB_LIB_INIT R_init_mslt

#include <TMB.hpp>
using namespace Eigen;//Needed for utilisation of sparse structures
using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
using namespace density;
#include "../inst/include/barrier.hpp"
#include "../inst/include/detection_functions.hpp"
#include "../inst/include/expm_mmpp.hpp"

template<class Type>
Type objective_function<Type>::operator() (){

  //Load data -------
  DATA_STRUCT(spdeMatrices,spde_t); //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(AalongLines);  //A-matrix for integration points along line
  DATA_SPARSE_MATRIX(Apred);  //A-matrix for prediction
  DATA_SPARSE_MATRIX(AObs);  //A-matrix as sightings
  DATA_SCALAR(area); //Area of area of interest

  DATA_MATRIX(X_g_obs); //Design matrix for observations (detection)
  DATA_MATRIX(X_g_left); //Design matrix for left side along transect lines
  DATA_MATRIX(X_g_right); //Design matrix for right side along transect lines
  DATA_MATRIX(X_z); //Design matrix for rate along transect
  DATA_MATRIX(X_z_pred); //Design matrix for rate in prediction points

  DATA_VECTOR(lineIntegralDelta); //Step length along lines for integration
  DATA_IVECTOR(code); // Codes 0,1,2 for node type. See details below
  DATA_VECTOR(distObs); //Perpendicular distances for observations
  DATA_INTEGER(g_function); //Which detection function used, 1: half normal 2:
  DATA_VECTOR(penalize); //First element is 1 if likelihood is penalized, second element is a penalize konstant
  DATA_VECTOR(pcPriorsRange_intensity); //pcPrior for spatial range
  DATA_VECTOR(pcPriorsSD_intensity); //pcPrior for spatial marginal variance
  DATA_VECTOR(pcPriorsRange_size); //pcPrior for spatial range
  DATA_VECTOR(pcPriorsSD_size); //pcPrior for spatial marginal variance
  DATA_INTEGER(usePCpriors); //1: apply pcpriors

  DATA_IVECTOR(abortLeft);
  DATA_IVECTOR(abortRight);
  DATA_SCALAR(detectionTrunc); //longest distance in truncated observation distribution, not used if negative

  DATA_INTEGER(matern_intensity); //0: Not include GMRF, 1: Include GMRF
  DATA_INTEGER(matern_size); //0: Not include GMRF, 1: Include GMRF

  DATA_VECTOR(size); //Size of groups
  DATA_INTEGER(applyPodSize); //1: Apply pod size
  DATA_INTEGER(podSizeDist); //1: poisson, 2: negative binomial
  DATA_INTEGER(independentPodSize); //1: pod size is independent of intensity

  PARAMETER_VECTOR(log_sigma); //marginal standard deviation
  PARAMETER_VECTOR(log_kappa); //spatial scaling parameter
  PARAMETER_VECTOR(beta_g); //covariates for detection
  PARAMETER_VECTOR(beta_z); //Spline regression parameters
  PARAMETER_VECTOR(beta_size); //covariates for group size
  PARAMETER_VECTOR(logB); //hazard
  PARAMETER_VECTOR(log_mu); // mu(1),mu(2) are swiching rates
  PARAMETER(log_c_mmpp);      // lambda2/lambda1
  PARAMETER_VECTOR(logSizeNB); //Size parameter in negative binomial pod size

  PARAMETER_VECTOR(x_intensity); //latent spatial effect intensity
  PARAMETER_VECTOR(x_size); //latent spatial effect group size

  //Transform parameters
  vector <Type> sigma = exp(log_sigma);
  vector <Type> kappa = exp(log_kappa);
  vector <Type> mu = exp(log_mu);
  Type c_mmpp = exp(log_c_mmpp)+1;
  Type bHazard = exp(logB(0));

  //Calculates nll
  Type nll = 0.0;

  Eigen::SparseMatrix<Type> Q_intensity;
  Q_intensity = Q_spde(spdeMatrices,kappa(0));

  if(matern_intensity != 0){
    nll += GMRF(Q_intensity)(x_intensity);
    SIMULATE{
      GMRF(Q_intensity).simulate(x_intensity);
    }

    Type d = 2; //Part of spatial pc-prior
    Type rhoP;
    Type R = -log(pcPriorsRange_intensity(1))*pow(pcPriorsRange_intensity(0),d/2);
    Type S = -log(pcPriorsSD_intensity(1))/pcPriorsSD_intensity(0);
    if(usePCpriors==1){
      rhoP = sqrt(8)/kappa(0);
      nll -= log( d/2 * R *S * pow(rhoP,(-1-d/2))* exp(-R* pow(rhoP,(-d/2)) -S* sigma(0) )); //pc-prior contribution
    }
  }
  SIMULATE{
    REPORT(x_intensity);
  }

  Type scaleS = Type(1)/((4*3.14159265)*kappa(0)*kappa(0)); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
  x_intensity = x_intensity/sqrt(scaleS) *sigma(0);



  //----- Classical line transect likelihood----------
  for(int i=0; i<distObs.size(); i++){
    vector<Type> z_i = X_g_obs.row(i);
    Type y_i = distObs(i);
    switch(g_function){  // We do NOT include normalizing constant "eshw" (included later)
    case 1:
      nll -= log(g_half_normal(y_i,z_i,beta_g));
      break;
    case 2:
      nll -= log(g_hazard(y_i,z_i,beta_g,bHazard));  // Here we put hazard rate
      break;
    }
  }

  //------ Log Gaussian Cox process ----------

  //Structure needed for effort constribution
  vector<Type> Z_transect =  X_z*beta_z+  AalongLines*x_intensity;



  //Penalize
  if(penalize(0)==1){
    nll -= dexp(c_mmpp,penalize(1),true);
  }

  //Abundance
  Type aveMMPP = (1/mu(1))/(1/mu(0) + 1/mu(1)); //TODO, double check Hans
  vector<Type> linPred =   exp(X_z_pred*beta_z + Apred*x_intensity) * (1 + exp(log_c_mmpp)*aveMMPP);


  //Size part
  Type sizeNB;
  Type varNB;


  if(applyPodSize==1){

    if(matern_size != 0){
      Eigen::SparseMatrix<Type> Q_size;
      Q_size = Q_spde(spdeMatrices,kappa(1));

      nll += GMRF(Q_size)(x_size);
      SIMULATE{
        GMRF(Q_size).simulate(x_size);
      }


      Type d = 2; //Part of spatial pc-prior
      Type rhoP;
      Type R = -log(pcPriorsRange_size(1))*pow(pcPriorsRange_size(0),d/2);
      Type S = -log(pcPriorsSD_size(1))/pcPriorsSD_size(0);
      if(usePCpriors==1){
        rhoP = sqrt(8)/kappa(1);
        nll -= log( d/2 * R *S * pow(rhoP,(-1-d/2))* exp(-R* pow(rhoP,(-d/2)) -S* sigma(1) )); //pc-prior contribution
      }
      SIMULATE{
        REPORT(x_size);
      }

      Type scaleS_size = Type(1)/((4*3.14159265)*kappa(1)*kappa(1)); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
      x_size = x_size/sqrt(scaleS_size) *sigma(1);

    }


    vector<Type> linPredSize;
    vector<Type> linPredSizePredicion;
    if(independentPodSize==1){
      linPredSize =  exp(beta_size(0) + AObs*x_size);
      linPredSizePredicion =   exp(beta_size(0) + Apred*x_size);
    }else{
      linPredSize =  exp(beta_size(0) + beta_size(1)*(AObs*x_intensity));
      linPredSizePredicion =   exp(beta_size(0) + beta_size(1)*(Apred*x_intensity));
    }
    for(int i = 0; i<size.size();++i){
      if(podSizeDist==1){
        nll -= dpois(size(i),linPredSize(i), true) - log(1-dpois(Type(0),linPredSize(i)));
      }else{
       sizeNB = exp(logSizeNB(0));
       varNB = linPredSize(i) + linPredSize(i)*linPredSize(i)/sizeNB;
       nll -= dnbinom2(size(i), linPredSize(i),varNB ,true)- log(1-dnbinom2(Type(0),linPredSize(i),varNB));
      }
    }
    for(int i = 0; i<linPred.size(); ++i){
      if(podSizeDist==1){
        linPred(i) = linPred(i)*linPredSizePredicion(i)/(1-exp(-linPredSizePredicion(i)));
      }else{
        sizeNB = exp(logSizeNB(0));
        varNB = linPredSizePredicion(i) + linPredSizePredicion(i)*linPredSizePredicion(i)/sizeNB;
        linPred(i) = linPred(i)*(linPredSizePredicion(i)/(1-dnbinom2(Type(0), linPredSizePredicion(i), varNB)));

      }
    }
  }


  Type abundance = area*linPred.sum()/linPred.size();
  Type logAbundance = log(abundance);
  ADREPORT(logAbundance);

  Type range_psi = -log(0.1)/mu.sum();
  Type k_psi =   1 + exp(log_c_mmpp)*aveMMPP;
  ADREPORT(range_psi);
  ADREPORT(k_psi);

  if(applyPodSize==1){
    sizeNB = exp(logSizeNB(0));
    Type linPredMeanSize = exp(beta_size(0)+ 0.5*sigma(1)*sigma(1));
    varNB = linPredMeanSize + linPredMeanSize*linPredMeanSize/sizeNB;
    Type meanGroupSize = linPredMeanSize/(1-dnbinom2(Type(0), linPredMeanSize, varNB));
    ADREPORT(meanGroupSize);
  }


  return nll;
}
