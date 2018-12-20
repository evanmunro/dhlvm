#include "util.h"
#include "Array3d.h"
#include <Rcpp.h>

using namespace Rcpp; 

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
        last_per_prob = log(P(k,S_upd(0,t-1)));
      }
      if (t==(Time-1)) {
        next_per_prob = 0;
      } else {
        next_per_prob = log(P(S_prev[t+1],k));
      }
      
      // should I start taking logs here to deal with the 0 issues..
      state_prob[k] = exp(last_per_prob+util::ldmultinom_cpp(data(t,_),beta(k,_))+next_per_prob);
    }
    if(util::sum_cpp(state_prob)==0) {
      NumericVector uniform_vec(K,(double) 1/K);
      state_prob = uniform_vec;
    }
    else {
      NumericVector denom(K,util::sum_cpp(state_prob));
      state_prob = util::divide(state_prob,denom);
    }
// Rcout <<state_prob << std::endl; 
    NumericVector mult_sample = as<NumericVector>(util::callRMultinom(state_prob));
    int state = util::whichNonZero(mult_sample);
    S_upd(0,t) = state;
    for (int k= 1; k< K; k++) {
      S_upd(k,t) = state_prob[k-1]; 
    }
  }
  return S_upd;
}

NumericMatrix sampleTransition_cpp(NumericVector S,NumericMatrix alpha) {
  int K = alpha.nrow(); 
  int Time = S.size(); 
  NumericMatrix N(K,K); 
  NumericMatrix P(K,K); 
  for (int t=1; t<Time; t++ ) {
    N(S[t],S[t-1])+=1; 
  }
  
  for (int k=0; k<K; k++) {
    P(_,k) = util::rdirichlet_cpp(alpha(_,k)+N(_,k)); 
  }
  return P; 
}

NumericMatrix sampleBeta_cpp(NumericMatrix data,NumericVector S,
                             NumericMatrix eta) {
  
  int Time = data.nrow(); 
  int K = eta.nrow(); 
  int V = eta.ncol(); 
  NumericMatrix beta(K,V); 
  NumericMatrix m(K,V); 
  for (int t=0; t<Time; t++) {
    m(S[t],_) = m(S[t],_) + data(t,_); 
  }
  for (int k=0; k< K; k++) {
    beta(k,_) = util::rdirichlet_cpp(eta(k,_)+m(k,_)); 
  }
  return beta; 
}

/*
* @brief Run gibbs sampler for markov switching model 
* 
*/
// [[Rcpp::export]]
List discreteMS_cpp(NumericMatrix data,NumericMatrix eta,
                    NumericMatrix alpha,int K, int steps) {
 
  int T = data.nrow(); 
  int V = data.ncol(); 
  Array3d P_samp = Array3d(K,K,steps);
  Array3d S_samp = Array3d(K,T,steps); 
  Array3d beta_samp = Array3d(K,V,steps); 

  
  //initialize each of the parameter samples based on draws from prior
  NumericMatrix initBeta(K,V); 
  NumericMatrix initP(K,K); 
  for (int k=0; k < K; k++) {
    initBeta(k,_) = util::rdirichlet_cpp(eta(k,_)); 
    initP(_,k) = util::rdirichlet_cpp(alpha(_,k)); 
  }
  beta_samp.setMatrix(0,initBeta); 
  P_samp.setMatrix(0,initP); 
  NumericMatrix initS(K,T); 
  NumericVector randomProb(K,(double) 1/K); 
  
  // TODO: make sampling a function or find existing function 
  for (int t=0; t<T; t++) {
    NumericVector mult_sample = as<NumericVector>(util::callRMultinom(randomProb)); 
    initS[t] = util::whichNonZero(mult_sample); 
  }
  S_samp.setMatrix(0,initS); 
  
  // run gibbs sampler
  for (int s=1; s< steps; s++) {
    NumericMatrix newS = sampleState_cpp(data,P_samp.getMatrix(s-1),
                                         S_samp.getMatrixRow(0,s-1),beta_samp.getMatrix(s-1)); 
    S_samp.setMatrix(s,newS); 
    P_samp.setMatrix(s,sampleTransition_cpp(S_samp.getMatrixRow(0,s),alpha)); //THIS IS FINE
    beta_samp.setMatrix(s,sampleBeta_cpp(data,S_samp.getMatrixRow(0,s),eta)); 
  }
  List posteriors; 
  posteriors["S"] = S_samp.getVector(); 
  posteriors["P"] = P_samp.getVector(); 
  posteriors["beta"] = beta_samp.getVector(); 
  return posteriors ; 
}

