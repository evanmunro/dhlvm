#include "Array3d.h"
  
//constructors 
Array3d::Array3d(int row, int col, int slice) {
  nrow = row; 
  ncol = col; 
  nslice = slice; 
  NumericVector vec(nrow*ncol*nslice); 
  data = vec; 
}

//convert 3D array to numeric vector
int Array3d::indexConvert(int row,int col,int slice) {
  int result = slice*nrow*ncol + col*nrow + row; 
  return result; 
}

void Array3d::setMatrix(int slice, NumericMatrix insertMat) { 
  for (int j = 0 ; j <ncol; j++ ) {
    for (int i = 0 ; i < nrow; i++) {
      data(indexConvert(i,j,slice)) = insertMat(i,j); 
    }
  }
}

NumericMatrix Array3d::getMatrix(int slice) {
  NumericMatrix matAtSlice(nrow,ncol); 
  for (int j = 0; j < ncol; j++) {
    for (int i = 0 ; i < nrow; i++) {
      matAtSlice(i,j) = data(indexConvert(i,j,slice)); 
    }
  }
  return matAtSlice; 
} 

NumericVector Array3d::getMatrixRow(int row, int slice) {
  NumericVector rowAtSlice(ncol); 
  for (int j = 0; j < ncol; j++) {
      rowAtSlice(j) = data(indexConvert(row,j,slice)); 
    }
  return rowAtSlice; 
} 

NumericVector Array3d::getVector() {
  data.attr("dim") = IntegerVector::create(nrow,ncol,nslice); 
  return data; 
} 





