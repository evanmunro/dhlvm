#include "util.h"
#include "Array3d.h"
#include <Rcpp.h>

using namespace Rcpp;

List sampleBeta_cpp(NumericMatrix data,NumericVector Z,
                    List eta) { 
  int J = data.ncol(); 
  int N = data.nrow(); 
  List beta(J); 
  for(int j=0; j< J; j++) { 
    NumericMatrix eta_j = as<NumericMatrix>(eta[j]); 
    int Lj = eta_j.ncol(); 
    int K = eta_j.nrow(); 
    NumericMatrix profile_counts(K,Lj); 
    NumericMatrix beta_j(K,Lj); 
    NumericVector data_j = data(j,_); 
    for (int n=0; n< N; n++) { 
      int k = Z[n]; 
      int response = data_j[n]; 
      profile_counts(k,response)++; 
    }
    for (int k=0; k< K; k++) {
      beta_j(k,_) = util::rdirichlet_cpp(eta_j(k,_)+profile_counts(k,_)); 
    }
    beta[j] = beta_j; 
  }
  return(beta); 
}

NumericVector sampleZ_cpp(NumericMatrix data, NumericVector groups,NumericMatrix pi, List beta) {
  int N = data.nrow(); 
  int K = pi.ncol(); 
  NumericVector z(N); 
  int J = data.ncol(); 

  for (int i =0; i <N; i++) { 
    NumericVector X_i = data(i,_); 
    NumericVector prob(K); 
    int normalization = 0; 
    NumericVector pi_c = pi(groups[i],_); 
    for (int k =0; k < K; k++) { 
      double prob_k = log(pi_c[k]); 
      for (int j=0; j<J; j++) { 
        NumericMatrix beta_j = as<NumericMatrix>(beta[j]); 
        prob_k += log(beta_j(k,X_i[j])); 
      }
      prob_k = exp(prob_k); 
      normalization+=prob_k; 
    }
    prob = prob/normalization; 
    z[i] = sample(util::vectorRange(K),1,true,prob)[0];  
  }
  return z; 
} 

NumericMatrix samplepi_cpp(NumericVector groups,NumericVector z,NumericMatrix alpha) {
  int K = alpha.ncol(); 
  int C = alpha.nrow(); 
  int N = z.size(); 
  NumericMatrix counts(C,K); 
  NumericMatrix pi(C,K); 
  
  for (int i = 0; i < N; i++) { 
    counts(groups[i],z[i])++; 
  }
  for(int c= 0; c < C; c++) { 
    pi(c,_) = util::rdirichlet_cpp(alpha(c,_)+counts(c,_)); 
  }
  return pi; 
}

// [[Rcpp::export]]
List lda_cpp(NumericMatrix data,NumericVector groups, List eta,
                    NumericMatrix alpha, int steps) {
  
  int K = alpha.ncol(); 
  int C = alpha.nrow(); 
  int N = data.nrow(); 
  int J = data.ncol(); 
  Array3d pi_samp = Array3d(N,K,steps);
  NumericMatrix z_samp = NumericMatrix(N,steps);  
  List beta_samp(steps); 
  
  //initialize each of the parameter samples based on draws from prior
  List initBeta(J); 
  NumericMatrix initPi(C,K); 

  for (int j=0; j< J ; j++) { 
    NumericMatrix eta_j = as<NumericMatrix>(eta[j]); 
    initBeta[j] = util::initializeMultinomial(eta_j); 
  }
  beta_samp[0] = initBeta; 
  pi_samp.setMatrix(0,util::initializeMultinomial(alpha));
  
  // run gibbs sampler
  for (int s=1; s< steps; s++) {
    z_samp(s,_) = sampleZ_cpp(data,groups,pi_samp.getMatrix(s-1),beta_samp[s-1]); 
    pi_samp.setMatrix(s,samplepi_cpp(groups,z_samp(s,_),alpha)); 
    beta_samp[s] = sampleBeta_cpp(data,pi_samp.getMatrix(s-1),eta); 
  }
  List posteriors; 
  posteriors["Z"] = z_samp; 
  posteriors["pi"] = pi_samp.getVector(); 
  posteriors["beta"] = beta_samp; 
  return posteriors; 
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
