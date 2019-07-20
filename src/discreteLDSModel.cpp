#include "util.h"
#include "Array3d.h"
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List sampleZ_old(List data, NumericMatrix theta, NumericMatrix beta) { 
  int T = data.size(); 
  int K = theta.nrow(); 
  List z(T); 
  for (int t = 0; t < T; t++) {
    NumericVector data_t = as<NumericVector>(data[t]); 
    int Nt = data_t.size(); 
    NumericVector tStates(Nt); 
    for(int i = 0; i < Nt;i++ ) {
      int v = data_t[i]-1; 
      NumericVector prob(K); 
      double sum = 0; 
      for (int k =0; k < K; k++) {
        prob[k] = theta(k,t)*beta(k,v); 
        sum = sum + prob[k]; 
      }
      for (int k =0; k < K; k++) {
        prob[k] = prob[k]/sum; 
      }
      NumericVector mult_sample = as<NumericVector>(util::callRMultinom(prob));
      tStates[i] = util::whichNonZero(mult_sample);
    }
    z[t] = tStates; 
  }
  return z; 
}

NumericMatrix sampleSigmaOld(NumericMatrix gamma, double v0, double s0) { 
  int K = gamma.nrow(); 
  int T = gamma.ncol(); 
  NumericMatrix sigma(K,K); 
  for (int k = 0 ; k < K; k ++ ) {
    NumericVector gk = gamma(k,_); 
    
    double ssq = 0; 
    for (int t=1; t<T; t++) { 
      double diff = gk[t] - gk[t-1]; 
      ssq += pow(diff,2); 
    }
    double v1 = v0 + T; 
    double s1 = s0 + ssq; 
    // Rcout<<"s1: "<< s1 << std::endl; 
    //sigma(k,k) = 1/R::rgamma(v1/2,2/s1); 
    sigma(k,k) = 1/rgamma(1,v1/2,2/s1)[0]; 

  } 

  return sigma; 
}

// [[Rcpp::export]]
NumericMatrix sampleGamma_old(List z, NumericMatrix sigma,
                        NumericMatrix gammaprev,double shrink) { 
  int K = gammaprev.nrow(); 
  int T = gammaprev.ncol(); 
  NumericMatrix gammanew(K,T); 

  for (int k=0; k < K; k++){
    for (int t=0; t<T; t++) {
      double noise = R::rnorm(0,shrink);
      NumericVector zt = as<NumericVector>(z[t]);
      int Nt = zt.size();
      int ckt = util::sumEquals(zt,k);
      double term1 = 0;
      double term2 = 0;
      if(t!=0){
        term1 = -1/sigma(k,k)*(gammaprev(k,t) - gammaprev(k,t-1)); //should be gammanew on RHS
      }
      if(t!=(T-1)) {
        term2 = -1/sigma(k,k)*(gammaprev(k,t+1) -gammaprev(k,t));
      } 
      
      NumericVector thetat = util::softmax(gammaprev(_,t));
      double term3 = ckt - Nt*thetat(k);
      double gradlogprob = term1 + term2 + term3;
      double delta = shrink/2*gradlogprob + noise; //delta is growing expotentitally
   
      gammanew(k,t) = gammaprev(k,t) + delta;
    }
  }
  //Rcout << gammanew << std::endl; 
  return gammanew; 
}

// [[Rcpp::export]]
NumericMatrix sampleBetaLDA(List data, List z,NumericMatrix eta) { 
  int K = eta.nrow(); 
  int V = eta.ncol(); 
  int T = data.size(); 
  NumericMatrix beta(K,V); 
  NumericMatrix m(K,V); 
  for (int t= 0 ; t < T; t++) {
      NumericVector datat = as<NumericVector>(data[t]); 
      NumericVector zt = as<NumericVector>(z[t]); 
      int Nt = datat.size(); 
      for (int n = 0; n < Nt; n++ ) {
        int v = datat[n]-1;
        int k = zt[n];
        m(k,v) = m(k,v) +1;
      }
  }
  for (int k =0; k < K; k++) {
    beta(k,_) = util::rdirichlet_cpp(eta(k,_) + m(k,_)); 
  }
  return beta; 
  
}

// [[Rcpp::export]]
List discreteLDS_cpp(List data,NumericMatrix eta, double v0, double s0, 
                    double tune, int K,int V, int steps) {
  
  int T = data.size(); 
  Array3d gammaDraw = Array3d(K,T,steps);
  Array3d thetaDraw = Array3d(K,T,steps); 
  Array3d sigmaDraw = Array3d(K,K,steps); 
  Array3d betaDraw = Array3d(K,V,steps); 
  
  //initialize each of the parameter samples based on draws from prior
  NumericMatrix initBeta(K,V); 
  NumericMatrix initSigma(K,K); 
  NumericMatrix initGamma(K,T); 
  NumericMatrix initTheta(K,T); 
  for (int k=0; k < K; k++) {
    initBeta(k,_) = util::rdirichlet_cpp(eta(k,_));
    initSigma(k,k) = 1/R::rgamma(v0/2,2/s0);
    initGamma(k,0) = R::rnorm(0,initSigma(k,k)); 
  }
  initTheta(_,0) =util::softmax(initGamma(_,0)); 
  for (int t=1; t<T; t++) {
    for (int k = 0; k < K; k++) {
      initGamma(k,t) = R::rnorm(initGamma(k,t-1),initSigma(k,k)); 
    }
    initTheta(_,t) = util::softmax(initGamma(_,t)); 
  }
  betaDraw.setMatrix(0,initBeta); 
  sigmaDraw.setMatrix(0,initSigma); 
  gammaDraw.setMatrix(0,initGamma); 
  thetaDraw.setMatrix(0,initTheta); 

  //run gibbs sampler
  for (int s=1; s< steps; s++) {
    double shrink = tune*pow(1+s,-0.5);
    List z = sampleZ_old(data,thetaDraw.getMatrix(s-1),betaDraw.getMatrix(s-1));
    NumericVector z1  = as<NumericVector>(z[1]); 
    gammaDraw.setMatrix(s,sampleGamma_old(z,sigmaDraw.getMatrix(s-1),
                                    gammaDraw.getMatrix(s-1),shrink));
    thetaDraw.setMatrix(s,util::softmax_col(gammaDraw.getMatrix(s)));
    sigmaDraw.setMatrix(s,sampleSigmaOld(gammaDraw.getMatrix(s),v0,s0));
    betaDraw.setMatrix(s,sampleBetaLDA(data,z,eta));
  }
  
  List posteriors; 
  posteriors["theta"] = thetaDraw.getVector(); 
  posteriors["gamma"] = gammaDraw.getVector(); 
  posteriors["sigma"] = sigmaDraw.getVector(); 
  posteriors["beta"] = betaDraw.getVector(); 
  return posteriors ; 
}
