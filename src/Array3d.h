#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <Rcpp.h>
using namespace Rcpp;

class Array3d {
  NumericVector data; 
  int nrow; 
  int ncol; 
  int nslice; 
  int indexConvert(int, int, int); 
  public: 
    Array3d(void); 
    Array3d(int, int, int); 
    NumericMatrix getMatrix(int); 
    NumericVector getMatrixRow(int,int); 
    void setMatrix(int,NumericMatrix); 
    NumericVector getVector(void); 
};

#endif 

