#include "util.h"
#include "Array3d.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double posteriorLikelihood(NumericMatrix data, NumericVector z, List beta) { 
  int J = data.ncol(); 
  int N = data.nrow(); 
  double likelihood = 0.0; 
  for (int i =0; i < N; i++) { 
    int k = z[i]; 
    for (int j =0; j <J; j++) { 
      NumericMatrix beta_j = as<NumericMatrix>(beta[j]);
      int response = data(i,j); 
      likelihood += log(beta_j(k,response)); 
    }
  }
  return(likelihood); 
}

List sampleBeta(NumericMatrix data,NumericVector Z,
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
    NumericVector data_j = data(_,j); 
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

NumericVector sampleZ(NumericMatrix data, NumericVector groups,NumericMatrix pi, List beta) {
  int N = data.nrow(); 
  int K = pi.ncol(); 
  NumericVector z(N); 
  int J = data.ncol(); 

  for (int i =0; i <N; i++) { 
    NumericVector X_i = data(i,_); 
    NumericVector prob(K); 
    NumericVector pi_c = pi(groups[i],_); 
    for (int k =0; k < K; k++) { 
      double prob_k = log(pi_c[k]); 
      for (int j=0; j<J; j++) { 
        NumericMatrix beta_j = as<NumericMatrix>(beta[j]);
        prob_k += log(beta_j(k,X_i[j])); 
      }
      prob[k] = exp(prob_k); 
    }
    prob = prob/sum(prob); 
    z[i] = sample(util::vectorRange(K),1,true,prob)[0];  
  }
  return z; 
} 

NumericMatrix samplePi(NumericVector groups,NumericVector z,NumericMatrix alpha) {
  int K = alpha.ncol(); 
  int G = alpha.nrow(); 
  int N = z.size(); 
  NumericMatrix counts(G,K); 
  NumericMatrix pi(G,K); 
  
  for (int i = 0; i < N; i++) { 
    counts(groups[i],z[i])++; 
  }
  for(int g= 0; g < G; g++) { 
    pi(g,_) = util::rdirichlet_cpp(alpha(g,_)+counts(g,_)); 
  }
  return pi; 
}



// [[Rcpp::export]]
List hlc_cpp(NumericMatrix data,NumericVector groups, List eta,
                    NumericMatrix alpha, int steps) {
  
  int K = alpha.ncol(); 
  int C = alpha.nrow(); 
  int N = data.nrow(); 
  int J = data.ncol(); 
  Array3d pi_samp = Array3d(C,K,steps); 
  NumericMatrix z_samp = NumericMatrix(N,steps);  
  List beta_samp(steps); 
  NumericVector likelihood = NumericVector(steps); 
  
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
    z_samp(_,s) = sampleZ(data,groups,pi_samp.getMatrix(s-1),beta_samp[s-1]); 
    pi_samp.setMatrix(s,samplePi(groups,z_samp(_,s),alpha));
    beta_samp[s] = sampleBeta(data,z_samp(_,s),eta); 
    likelihood[s] = posteriorLikelihood(data,z_samp(_,s),beta_samp[s]); 
  }
  List posteriors; 
  posteriors["Z"] = z_samp; 
  posteriors["pi"] = pi_samp.getVector(); 
  posteriors["beta"] = beta_samp; 
  posteriors["likelihood"] = likelihood; 
  return posteriors; 
}

NumericMatrix sampleSigma(NumericMatrix pi, double v0, double s0) { 
  int K = pi.ncol(); 
  int T = pi.nrow(); 
  NumericMatrix sigma(K,K); 
  for (int k = 0 ; k < K; k ++ ) {
    NumericVector pi_k = pi(_,k); 
    double ssq = 0; 
    for (int t=1; t<T; t++) { 
      double diff = pi_k[t] - pi_k[t-1]; 
      ssq += pow(diff,2); 
    }
    double v1 = v0 + T; 
    double s1 = s0 + ssq; 
    sigma(k,k) = sqrt((1/rgamma(1,v1/2,2/s1))[0]); 
  } 
  return sigma; 
}


NumericMatrix sampleGamma(NumericVector groups, NumericVector z, 
                          NumericMatrix gammaPrev,NumericMatrix sigma,double shrink) {
  int K = gammaPrev.ncol(); 
  int T = gammaPrev.nrow(); 
  int N = z.size(); 
  NumericMatrix counts(T,K); 
  NumericMatrix gammaNew(T,K); 
  NumericVector timeCounts(T); 
  NumericMatrix pi(T,K);
  for (int i=0; i<N; i++) { 
    counts(groups[i],z[i])++; 
    timeCounts(groups[i])++; 
  }

  for (int k=0; k < K; k++){
    for (int t=0; t<T; t++) {
      double noise = rnorm(1,0,shrink)[0];
      double term1 = 0;
      double term2 = 0;
      if(t!=0){
        term1 = -1/pow(sigma(k,k),2)*(gammaPrev(t,k) - gammaPrev(t-1,k)); //should be gammanew on RHS
      }
      if(t!=(T-1)) {
        term2 = -1/pow(sigma(k,k),2)*(gammaPrev(t+1,k) -gammaPrev(t,k));
      } 
      
      NumericVector thetat = util::softmax(gammaPrev(t,_));
      double term3 = counts(t,k) - timeCounts(t)*thetat(k);
      double gradlogprob = term1 + term2 + term3;
      double delta = shrink/2*gradlogprob + noise; //delta is growing expotentitally
      
      gammaNew(t,k) = gammaPrev(t,k) + delta;
    }
  }
  return gammaNew; 
}

// [[Rcpp::export]]
List dhlc_cpp(NumericMatrix data,NumericVector groups, List eta,
            int v0, int s0,double tune, int K, int T, int steps) {
  
  int N = data.nrow(); 
  int J = data.ncol(); 
  Array3d gammaDraw = Array3d(T,K,steps); 
  Array3d piDraw = Array3d(T,K,steps); 
  Array3d sigmaDraw = Array3d(K,K,steps); 
  NumericMatrix zDraw = NumericMatrix(N,steps);  
  List betaDraw(steps); 
  
  //initialize each of the parameter samples based on draws from prior
  List initBeta(J); 
  NumericMatrix initPi(T,K); 
  NumericMatrix initSigma(K,K); 
  NumericMatrix initGamma(T,K); 
  for (int k=0; k< K; k++ ) { 
    initSigma(k,k) = sqrt((1/rgamma(1,v0,1/s0))[0]);  
    initGamma(0,k) = rnorm(1,0,initSigma(k,k))[0]; 
    for (int t=1; t<T; t++) { 
      initGamma(t,k) = rnorm(1,initGamma(t-1,k),initSigma(k,k))[0]; 
    }
  }
  
  for (int t=0; t<T; t++) { 
    initPi(t,_) = util::softmax(initGamma(t,_)); 
  }

  for (int j=0; j< J ; j++) { 
    NumericMatrix eta_j = as<NumericMatrix>(eta[j]); 
    initBeta[j] = util::initializeMultinomial(eta_j); 
  }
  betaDraw[0] = initBeta; 
  gammaDraw.setMatrix(0,initGamma); 
  piDraw.setMatrix(0,initPi); 
  sigmaDraw.setMatrix(0,initSigma);  
  
  // run gibbs sampler
  for (int s=1; s< steps; s++) {
    double shrink = tune*pow(1+s,-0.5);
    zDraw(_,s) = sampleZ(data,groups,piDraw.getMatrix(s-1),betaDraw[s-1]); 
    gammaDraw.setMatrix(s,sampleGamma(groups,zDraw(_,s),gammaDraw.getMatrix(s-1), 
                                      sigmaDraw.getMatrix(s-1),shrink)); 
    piDraw.setMatrix(s,util::softmax(gammaDraw.getMatrix(s))); 
    sigmaDraw.setMatrix(s,sampleSigma(gammaDraw.getMatrix(s),v0,s0)); 
    betaDraw[s] = sampleBeta(data,zDraw(_,s),eta); 
  }
  
  List posteriors; 
  posteriors["Z"] = zDraw; 
  posteriors["sigma"] = sigmaDraw.getVector(); 
  posteriors["pi"] = piDraw.getVector(); 
  posteriors["beta"] = betaDraw; 
  return posteriors; 
}

/*** R
#mean(rgamma_cpp(1000,10,1)) 
*/
