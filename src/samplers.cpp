#include "Array3d.h"
#include <cmath>

//returns 1-indexed index of non-zero element in a vector 
int whichNonZero(NumericVector x) {
  int index = -1; 
  for (int i =0; i < x.size(); i++) { 
    if (x[i] >0) { 
      index = i+1; 
    }
  }
  return index; 
}

IntegerVector callRMultinom(NumericVector prob) {
  int K = prob.size(); 
  IntegerVector d(K); 
  R::rmultinom(1,prob.begin(),K,d.begin()); 
  return d; 
}

NumericVector divide(NumericVector x, NumericVector y) {
  NumericVector result(x.size()); 
  for (int i =0; i < x.size(); ++i) {
    result[i] = x[i]/y[i]; 
  }
  return result; 
}

double sum_cpp(NumericVector x) {
  double total =0; 
  for (int i =0; i < x.size(); i++) {
    total += x[i]; 
  }
  return total; 
}

//calculates log multinomial density 
// [[Rcpp::export]]
long double ldmultinom_cpp(NumericVector x,NumericVector prob) {
  long double logp_sum = lgamma(sum_cpp(x)+1);  
  for (int i =0; i < x.size(); i++) { 
    logp_sum += (x[i]*log(prob[i]) - lgamma(x[i]+1)); 
  }
  return logp_sum; 
}

// This function is adapted from http://www.mjdenny.com/blog.html
// Takes a single sample from a dirichlet distribution 
// [[Rcpp::export]]
NumericVector rdirichlet_cpp(NumericVector alpha_m) {
  int distribution_size = alpha_m.size();
  // each row will be a draw from a Dirichlet
  NumericVector distribution(distribution_size);
  double sum_term = 0;
  // loop through the distribution and draw Gamma variables
  for (int j = 0; j < distribution_size; ++j) {
    double cur = R::rgamma(alpha_m[j],1.0);
    distribution[j] = cur;
    sum_term += cur;
    }
    // now normalize
    for (int j = 0; j < distribution_size; ++j) {
      distribution(j) = distribution(j)/sum_term;
    }
  return(distribution);
}

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
NumericMatrix sampleState_cpp(NumericMatrix data,NumericMatrix P,
                             NumericVector S_prev, NumericMatrix beta) {
  int K = P.nrow(); 
  int Time = data.nrow(); 
  NumericVector state_prob(K);
  NumericMatrix S_upd(K,Time);

  for (int t=0; t< Time; t++) {
    for (int k= 0; k< K; k++) {
      double last_per_prob;
      double next_per_prob;
      if (t==0) {
        last_per_prob = 0;
      } else {
        last_per_prob = log(P(k,S_upd(0,t-1)-1));
      }
      if (t==(Time-1)) {
        next_per_prob = 0;
      } else {
        next_per_prob = log(P(S_prev[t+1]-1,k));
      }
      
// should I start taking logs here to deal with the 0 issues..
      state_prob[k] = exp(last_per_prob+ldmultinom_cpp(data(t,_),beta(k,_))+next_per_prob);
    }
    if(sum_cpp(state_prob)==0) {
         NumericVector uniform_vec(K,(double) 1/K);
         state_prob = uniform_vec;
    }
    else {
      NumericVector denom(K,sum_cpp(state_prob));
      state_prob = divide(state_prob,denom);
    }
    NumericVector mult_sample = as<NumericVector>(callRMultinom(state_prob));
    int state = whichNonZero(mult_sample);
    S_upd(0,t) = state;
    for (int k= 1; k< K; k++) {
      S_upd(k,t) = state_prob[k-1]; 
    }
  }
  return S_upd;
}

// [[Rcpp::export]]
NumericMatrix sampleTransition_cpp(NumericVector S,NumericMatrix alpha) {
  int K = alpha.nrow(); 
  int Time = S.size(); 
  NumericMatrix N(K,K); 
  NumericMatrix P(K,K); 
  for (int t=1; t<Time; t++ ) {
    N(S[t]-1,S[t-1]-1)+=1; 
  }
  
  for (int k=0; k<K; k++) {
    P(_,k) = rdirichlet_cpp(alpha+N(_,k)); 
  }
  return P; 
}

// [[Rcpp::export]]
NumericMatrix sampleBeta_cpp(NumericMatrix data,NumericVector S,
                             NumericMatrix eta) {
  int Time = data.nrow(); 
  int K = eta.nrow(); 
  int V = eta.ncol(); 
  NumericMatrix beta(K,V); 
  NumericMatrix m(K,V); 
  for (int t=0; t<Time; t++) {
    m(S[t]-1,_) = m(S[t]-1,_) + data(t,_); 
  }
  for (int k=0; k< K; k++) {
    beta(k,_) = rdirichlet_cpp(eta(k,_)+m(k,_)); 
  }
  return beta; 
}

//Run gibbs sampler for markov switching model
//Eventually should port this to its own class 

List discreteMS(NumericMatrix data,NumericMatrix eta,
                NumericMatrix alpha,int K,int steps, int burn, int thin) {
  
  int T = data.nrow(); 
  int V = data.ncol(); 
  
  Array3d P_samp(K,K,steps) ;
  Array3d S_samp(K,T,steps); 
  Array3d beta_samp(K,V,steps); 
  
  //initialize each of the parameter samples based on draws from prior
  NumericMatrix initBeta(K,V); 
  NumericMatrix initP(K,K); 
  for (int k=0; k < K; k++) {
    initBeta(k,_) = rdirichlet_cpp(eta(k,_)); 
    initP(_,k) = rdirichlet_cpp(alpha(_,k)); 
  }
  beta_samp.setMatrix(1,initBeta); 
  P_samp.setMatrix(1,initP); 
  NumericMatrix initS(K,T); 
  NumericVector randomProb(K,0.5); 
  for (int t=0; t<T; t++) {
    NumericVector mult_sample = as<NumericVector>(callRMultinom(randomProb)); 
    initS[t] = whichNonZero(mult_sample); 
  }
  S_samp.setMatrix(1,initS); 
        
 // run gibbs sampler
  for (int s=1; s< steps; s++) {
    NumericMatrix newS = sampleState_cpp(data,P_samp.getMatrix(s-1),
                        S_samp.getMatrixRow(1,s-1),beta_samp.getMatrix(s-1)); 
    S_samp.setMatrix(s,newS); 
    P_samp.setMatrix(s,sampleTransition_cpp(S_samp.getMatrixRow(1,s-1),alpha)); 
    beta_samp.setMatrix(s,sampleBeta_cpp(data,S_samp.getMatrixRow(1,s-1),eta)); 
  }
  List posteriors; 
  posteriors["S"] = S_samp.getVector(); 
  posteriors["P"] = P_samp.getVector(); 
  posteriors["beta"] = beta_samp.getVector(); 
  return posteriors ; 
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
