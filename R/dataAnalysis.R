
#' Cleans up raw survey data matrix of numerical scores  by 1) converting all missing or not-applicable codes 
#' to a single response code L_j and 2) converting raw scores (that may not be contiguous) to 1:L_j, where the order of 
#' the original scores is maintained in the cleaned data 
#'
#' @param raw N x J matrix of raw numeric scores, where each score in a column corresponds to a different question response 
#' @param na.values vector of values in raw that should be considered missing or N/A responses.  
#' Any values in the matrix that are negative or NA are automatically included. 
#'
#' @return An N x J matrix where each column has values only in 1 to L_j 
#' @export
clean_data <- function(raw,na.values=NULL) { 
  na.code = 1000000
  raw[is.na(raw)] = na.code 
  raw[raw<0] = na.code 
  which.na <- t(apply(raw,MARGIN=1,FUN=function(x){return(x%in%na.values)})) 
  raw[which.na] = na.code 
  data_clean = apply(raw,MARGIN=2,FUN=function(x) {return(as.numeric(factor(x)))} )
  return(data_clean) 
}

#' Converts N x J response matrix to G x L response frequency matrix, when J>1 
#'
#' @param X the N x J response matrix 
#' @param groups N length vector of group membership indictors in 1,..,G 
#'
#' @return Y, a G x L response frequency matrix 
#' @export
#'
xtoAdjacency <- function(X,groups){ 
  N = nrow(X)
  J = ncol(X) 
  if (J <= 1) { 
    error("X must have more than 1 column")
  }
  G = length(unique(groups)) 
  L = apply(X,MARGIN=2,FUN=function(x){return(length(unique(x)))})
  V = sum(L)
  Y = matrix(0,nrow=G,ncol=V)
  for (i in 1:N) {
    for (j in 1:J) {
      Y[groups[i],sum(L[1:(j-1)])+X[i,j]] = Y[groups[i],sum(L[1:(j-1)])+X[i,j]] + 1 
    } 
  }
  for (j in 1:J){
    for (g in 1:G) {
      Y[g,(sum(L[1:(j-1)])+1):(sum(L[1:j]))] =Y[g,(sum(L[1:(j-1)])+1):(sum(L[1:j]))]/sum(Y[g,(sum(L[1:(j-1)])+1):(sum(L[1:j]))])
    }  
  }
  
  return(Y)
  
}