#ifndef ARRAY3D_H
#define ARRAY3D_H

#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

class Array3d {
  NumericVector data; 
  int nrow; 
  int ncol; 
  int nslice; 
  int indexConvert(int, int, int); 
  public: 
    Array3d(int, int, int); 
    NumericMatrix getMatrix(int); 
    NumericVector getMatrixRow(int,int); 
    void setMatrix(int,NumericMatrix); 
    NumericVector getVector(void); 
};

#endif 

