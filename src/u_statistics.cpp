#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::NumericMatrix mmult( Rcpp::NumericMatrix m , Rcpp::NumericMatrix v)
{
  if( ! (m.ncol() == v.nrow()) ) stop("Non-conformable arrays") ;

  Rcpp::NumericMatrix out(m.nrow(),v.ncol()) ;

  for (int i = 0; i < m.nrow(); i++) {
    for (int j = 0; j < v.ncol(); j++) {
      for(int k = 0; k < m.ncol(); k++) {
        out(i,j) += m(i,k) * v(k,j) ;
      }
    }
  }
  return out ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix transpose( Rcpp::NumericMatrix m)
{
  Rcpp::NumericMatrix out(m.ncol(),m.nrow()) ;

  for (int i = 0; i < out.nrow(); i++) {
    for (int j = 0; j < out.ncol(); j++) {
      out(i,j) = m(j,i) ;
    }
  }
  return out ;
}

// [[Rcpp::export]]
float inner( Rcpp::NumericVector a, Rcpp::NumericVector b)
{
  if( ! (a.size() == b.size()) ) stop("Non-conformable arrays") ;

  float out = 0;

  for (int i = 0; i < a.size(); i++) {
    out += a[i]*b[i] ;
  }
  return out ;
}

// [[Rcpp::export]]
Rcpp::NumericVector calcUCpp(Rcpp::NumericMatrix dat, float n, Rcpp::NumericMatrix K, Rcpp::NumericVector colmeans) {

  Rcpp::NumericVector Q(2);
  Rcpp::NumericMatrix data = clone(dat);
  int j = 0;
  int dim = colmeans.size();
  Rcpp::NumericVector Z1(dim);
  Rcpp::NumericVector Z2(dim);
  Rcpp::NumericMatrix result(n,dim);
  float corr =1-1/n;

  for(int i = 0; i < n; i++) {
    data(i,_) = data(i,_) - colmeans;
  }

  result = mmult(data, K);

  for(int i = 0; i < n; i++){
    j = i + 1;
    while(j < n){

      Z1 = result(i,_);
      Z2 = result(j,_);

      Q(0) += inner(Z1,Z1)*inner(Z2,Z2);
      Q(1) += pow(inner(Z1,Z2), 2);

      j++;
    }
  }
  Q(1) = Q(1)*2/n*1/(n-1)*1/corr*1/corr;
  Q(0) = Q(0)*2/n*1/(n-1)*1/corr*1/corr;

  return Q;
}


// [[Rcpp::export]]
Rcpp::NumericVector calcUCppV(Rcpp::NumericVector dat, float n, float mean) {

  Rcpp::NumericVector Q(2);
  Rcpp::NumericVector data = clone(dat);
  int j = 0;
  int dim = 1;
  Rcpp::NumericVector Z1(1);
  Rcpp::NumericVector Z2(1);
  Rcpp::NumericVector result(n);
  float corr =1-1/n;

  for(int i = 0; i < n; i++) {
    data[i] = data[i] - mean;
  }

  result = data;


  for(int i = 0; i < n; i++){
    j = i + 1;
    while(j < n){

      float z1 = result(i);
      float z2 = result(j);
      Q(0) += pow(z1,2)*pow(z2,2);
      Q(1) += pow(z1*z2,2);

      j++;
    }
  }
  Q(1) = Q(1)*2/n*1/(n-1)*1/corr*1/corr;
  Q(0) = Q(0)*2/n*1/(n-1)*1/corr*1/corr;

  return Q;
}

