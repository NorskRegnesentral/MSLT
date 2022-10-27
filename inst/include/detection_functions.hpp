// Different detection functions g(y) and corresponding ESHW

// Notation:
// y = perpendicular distance
// z = row from (covariate) design matrix
// w_left = truncation distance left hand side
// w_right = truncation distance right hand side
// beta = parameter vector (all parameters needed for g)
// b is a shape parameter (only for some detection functions)


// ------------ Half normal class ------------------------------

template<class Type>
Type g_half_normal(Type y, vector<Type> z, vector<Type> beta){
  Type pi = 3.141592;
  Type sigma = exp((beta*z).sum());
  Type g =  dnorm(y,Type(0),sigma,false)*sqrt(2*pi)*sigma;

  return g;
}

template<class Type>
Type eshw_half_normal(vector<Type> z, vector<Type> beta, Type w){
  Type pi = 3.141592;
  Type sigma = exp((beta*z).sum());
  Type eshw;
  if(w>0){
    eshw = sqrt(2*pi)*sigma*(pnorm(w/sigma,Type(0),Type(1)) -0.5) ;
  }else{
    eshw = sqrt(2*pi)*sigma*0.5;
  }

  return eshw;
}

// ------------ Hazard rate class ------------------------------

template<class Type>
Type g_hazard(Type y, vector<Type> z, vector<Type> beta, Type b){
  Type sigma = exp((beta*z).sum());
  Type g =  Type(1.0) - exp(-pow(y/sigma,-b));

  return g;
}

template<class Type>
Type eshw_hazard(vector<Type> z, vector<Type> beta, Type b){
  Type sigma = exp((beta*z).sum());
  Type eshw = exp(lgamma((b-1)/b))*sigma;

  return eshw;
}
