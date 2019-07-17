#include "util.h"

#include <cmath>

namespace util { 

  /* 
   * @brief looks for non-zero element in a vector 
   * @param x vector of numbers   
   * 
   * @return 0-indexed location of last non-zero element in vector
   */
  NumericVector vectorRange(int K) { 
    NumericVector rangeVector(K); 
    for (int i = 0; i < K; i++) { 
      rangeVector[i] = i;  
    }
    return rangeVector; 
  }

  /* 
  * @brief looks for non-zero element in a vector 
  * @param x vector of numbers   
  * 
  * @return 0-indexed location of last non-zero element in vector
  */
  int whichNonZero(NumericVector x) {
    int index = -1; 
    for (int i =0; i < x.size(); i++) { 
      if (x[i] >0) { 
        index = i; 
      }
    }
    return index; 
  }

  /* 
   * @brief function for single multinomial draw 
   * @param prob vector of multinomial probabilities 
   * @return IntegerVector of 0s with single 1 at index drawn 
   * 
   */
  IntegerVector callRMultinom(NumericVector prob) {
    int K = prob.size(); 
    IntegerVector d(K); 
    R::rmultinom(1,prob.begin(),K,d.begin()); 
    return d; 
  }
  
  /* 
   * @brief Elementwise vector division 
   * @param x the numerator 
   * @param y the denominator 
   * 
   * @return NumericVector result of elementwise division x/y 
   */
  NumericVector divide(NumericVector x, NumericVector y) {
    NumericVector result(x.size()); 
    for (int i =0; i < x.size(); ++i) {
      result[i] = x[i]/y[i]; 
    }
    return result; 
  }
  
  /*
   * @brief vector sum 
   * @param x the input vector 
   * 
   * @return double the sum of the input vector 
   */
  double sum_cpp(NumericVector x) {
    double total =0; 
    for (int i =0; i < x.size(); i++) {
      total += x[i]; 
    }
    return total; 
  }

/*
 * @brief check how many elements of vector equal to integer 
 * @param x the input vector 
 * @param k the input integer 
 * 
 * @return int the number of elements of x equal to k 
 */
  
  int sumEquals(NumericVector x,int k) {
    IntegerVector y = as<IntegerVector>(x); 
    int sumeq = 0; 
    for (int i = 0 ; i < y.size(); i ++) {
      if(y[i]==k) {
        sumeq = sumeq+1; 
      }
    }
    
    return sumeq; 
  }
  
  
  /* 
   * @brief Calculates log multinomial density for a vector of counts
   * 
   * @param x vector of  counts
   * @param prob multinomial probability parameters 
   * 
   * @return log multinomial density of data given probability 
   */
  long double ldmultinom_cpp(NumericVector x,NumericVector prob) {
    long double logp_sum = lgamma(sum_cpp(x)+1);  
    for (int i =0; i < x.size(); i++) { 
      logp_sum += (x[i]*log(prob[i]) - lgamma(x[i]+1)); 
    }
    return logp_sum; 
  }
  
  
  
  /* 
   * credit: adapted from http://www.mjdenny.com/blog.html
   * @brief Takes a single sample from a dirichlet distribution 
   * @param alpha_m dirichlet parameters 
   * @ return vector of probabilities drawn from dirichlet distribution
  */ 
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


  /* 
   * @brief looks for non-zero element in a vector 
   * @param x NumericMatrix Dirichlet prior    
   * 
   * @return Multinomial matrix initialized from Dirichlet prior matrix 
   */
  NumericMatrix initializeMultinomial(NumericMatrix prior) { 
    int rows = prior.nrow(); 
    NumericMatrix initialized(prior.nrow(), prior.ncol()); 
    for (int r = 0; r < rows; r++) { 
      initialized(r,_) = rdirichlet_cpp(prior(r,_)); 
    }
    return initialized; 
  }

  NumericVector softmax(NumericVector x) {
    int n = x.size(); 
    NumericVector result(n); 
    double sumX = 0 ; 
    for (int i =0; i < n; i ++) {
      sumX = sumX + exp(x[i]); 
    }
    for (int i =0; i<n; i ++ ){
      result[i] = exp(x[i])/sumX; 
    }
    return result; 
  }

  NumericMatrix softmax(NumericMatrix x) {
    int nrow = x.nrow(); 
    int ncol = x.ncol(); 
    NumericMatrix result(nrow,ncol); 
    for (int i = 0; i < ncol; i ++) {
      result(_,i) = softmax(x(_,i)); 
    }
    return result; 
  }

} 