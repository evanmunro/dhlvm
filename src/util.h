#ifndef UTIL_H
#define UTIL_H 

#include <Rcpp.h>

using namespace Rcpp; 

namespace util { 

  NumericVector vectorRange(int K); 
  NumericMatrix initializeMultinomial(NumericMatrix prior); 
  int whichNonZero(NumericVector x); 
  IntegerVector callRMultinom(NumericVector prob); 
  NumericVector divide(NumericVector x, NumericVector y); 
  double sum_cpp(NumericVector x); 
  int sumEquals(NumericVector x,int k); 
  long double ldmultinom_cpp(NumericVector x,NumericVector prob); 
  NumericVector rdirichlet_cpp(NumericVector alpha_m); 
  NumericVector softmax(NumericVector x); 
  NumericMatrix softmax(NumericMatrix x); 

}

#endif 