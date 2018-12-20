#ifndef DISCRETEMSMODEL_H
#define DISCRETEMSMODEL_H 

#include "Array3d.h"
#include <Rcpp.h>

using namespace Rcpp; 

class DiscreteMSModel { 
  Array3d P_samp; //samples of transition matrix
  Array3d S_samp; //samples of state
  Array3d beta_samp; //samples of beta
  int K; //< number of topics 
  int T; //< number of time periods 
  int V; //< number of categorical outcomes 
  NumericMatrix data; 
  NumericMatrix eta; //< prior for beta 
  NumericMatrix alpha; //< prior for P 
  int steps; 
  int burn; 
  int thin; 
  public: 
    DiscreteMSModel(NumericMatrix data,NumericMatrix eta,
                   NumericMatrix alpha,int K,int steps, int burn, int thin); 
}; 

#endif 