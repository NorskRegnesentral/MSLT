// Direct implementation for the matrix exponential in dim 2 for MMPP

template<class Type>
matrix<Type> expm_mmpp2(vector<Type> lambda, vector<Type> mu, Type l){
  Type e = Type(1e-10);

  vector<Type> S = lambda+mu;
  Type K = lambda.prod() + lambda(0)*mu(1) + lambda(1)*mu(0);
  Type sumS = S.sum();
  Type U = sqrt(sumS*sumS-4*K + e);
  vector<Type> theta(2);
  theta(0) = (sumS+U)/2;
  theta(1) = (sumS-U)/2;

  matrix<Type> M1(2,2), M2(2,2);
  M1(0,0) = S(1)-theta(1);
  M1(0,1) = mu(0);
  M1(1,0) = mu(1);
  M1(1,1) = S(0)-theta(1);

  M2 = M1;
  M2(0,0) += theta(1) - theta(0);
  M2(1,1) += theta(1) - theta(0);

  return (exp(-theta(1)*l)*M1 - exp(-theta(0)*l)*M2)/U;
}

