#include "util.h"
#include "Array3d.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List lda_cpp(NumericMatrix data, List eta,
                    NumericMatrix alpha,int K, int steps) {
  
  int N = data.nrow(); 
  int J = data.ncol(); 
  Array3d g_samp = Array3d(N,K,steps);
  NumericMatrix z_samp = NumericMatrix(N,steps);  
  List beta_samp; 
  
  //initialize each of the parameter samples based on draws from prior
  List initBeta; 
  NumericMatrix initG(N,K); 
  NumericVector initZ(N); 
  
  for (int j=0; j< J ; j++) { 
    NumericMatrix eta_j = eta[j]; 
    initBeta.push_back(NumericMatrix(K,eta_j.ncol())); 
    for (int k=0; k < K; k++) {
      NumericMatrix initBeta_j = initBeta[j]; 
      initBeta_j(k,_) = util::rdirichlet_cpp(eta_j(k,_)); 
    } 
  }

  beta_samp.push_back(initBeta); 
  g_samp.setMatrix(0,initG); 
  NumericVector randomProb(K,(double) 1/K); 
  
  // TODO: make sampling a function or find existing function 
  for (int n=0; n<N; n++) {
    NumericVector mult_sample = as<NumericVector>(util::callRMultinom(randomProb)); 
    initZ[n] = util::whichNonZero(mult_sample); 
  }
  
  z_samp(0,_)= initZ ; 
  
  // run gibbs sampler
  for (int s=1; s< steps; s++) {
    NumericVector newZ = sampleZ_cpp(data,g_samp.getMatrix(s-1),
                                         z_samp(s-1,_),beta_samp[s-1]); 
    z_samp(s,_) = newZ; 
    P_samp.setMatrix(s,sampleTransition_cpp(S_samp.getMatrixRow(0,s),alpha)); //THIS IS FINE
    beta_samp.setMatrix(s,sampleBetaStatic(data,S_samp.getMatrixRow(0,s),eta)); 
  }
  List posteriors; 
  posteriors["S"] = S_samp.getVector(); 
  posteriors["P"] = P_samp.getVector(); 
  posteriors["beta"] = beta_samp.getVector(); 
  return posteriors ; 
}


List lda_cpp(NumericMatrix data, List eta,
             NumericMatrix alpha,int K, int steps) {
  
} 


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
