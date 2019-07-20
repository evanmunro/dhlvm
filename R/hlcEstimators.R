#Author: Evan Munro, Stanford University 
#This file contains the Gibbs Samplers for the various latent variable models for discrete data 
#that are being developed as part of a working paper with Serena Ng

# #convert responses count to list format
countToListInput <- function(X) {
  Time <- nrow(X)
  V <- ncol(X)
  result <- list()
  for (t in 1:Time) {
    result[[t]] = rep(1:V,times=X[t,])
  }
  return(result)
}

#' Indicate which elements of the profile specific multinomial distributions have the highest weight
#' 
#' @param beta KxV data frame of posterior estimates of 
#' @param identifier V-length vector of names for each v-index of beta
#' @param top number of elements of beta to return 
#'
#' @return table with top k responses for each profile-specific multinomial distribution 
#' @export
topBetas <- function(beta,identifier,top=5) {
  library(dplyr)
  beta = as.data.frame(beta)
  K <- ncol(beta)
  beta$Perm = identifier
  result = beta[order(-beta[,"V1"]),c("Perm","V1")] 
  colnames(result) = c("State 1", "Probabilities 1")
  for (k in 2:K) {
    stateName = paste("V",k,sep="")
    result = bind_cols(result,beta[order(-beta[,stateName]),c("Perm",stateName)])
  }
  return(result[1:top,])
}

#' Get posterior mean for set of posterior estimates from sampler 
#'
#' @param estimate 3d array from sampler 
#' @param out vector of iterations to take mean over for gibbs sampler
#' @return matrix of means of posterior over sampler slices
#' @export
posteriorMean <- function(estimate,out) { 
  return(apply(estimate[,,out],MARGIN=c(1,2),FUN=mean))
}

posteriorMeanForList <- function(estimate,out) { 
  
}

#' Title
#'
#' @param data XxV matrix of counts 
#' @param eta prior for dirichlet dist for beta
#' @param v0 prior for invgamma dist for sigma
#' @param s0 prior rate for invgamma
#' @param tune tuning parameter for SGLD step 
#' @param K number of states
#' @param steps number of steps for sampler 
#' @param burn burn-in iterations for sampler 
#' @param skip thinning iterations for sampler
#'
#' @return List with posterior samples from 1) theta 2) sigma 3) beta
#' @export
discreteLDSModel <- function(data,eta,v0,s0,tune,K,steps,burn,skip) {
  V <- ncol(data)
  data <- countToListInput(data)
  posterior <- discreteLDS_cpp(data,eta,v0,s0,tune,K,V,steps)
  out <- seq(from=burn,to=steps,by=skip)
  posterior$out = out 
  return(posterior)
}

#' Runs Gibbs Sampler for markov switching model for discrete time series data with c++ loop
#' 
#' @param data TimexV data frame of counts of data at time t with value v
#' @param eta V-length prior for state-specific multinomial probabilities beta
#' @param alpha K-length prior for transition matrix
#' @param K Number of states
#' @param steps Total number of gibbs sampling iterations 
#' @param burn Number of iterations to discard
#' @param skip Thinning parameter for gibbs sampler
#'
#' @return List with posterior mean of 1) state-specific multinomial probabilities beta 2) transition matrix 3) states
#'
#' @export
discreteMSM <- function(data,eta,alpha,K,steps,burn,skip) {
  posterior <- discreteMS_cpp(data,eta,alpha,K,steps)
  out <- seq(from=burn,to=steps,by=skip)
  posterior$out = out 
  return(posterior)
} 

#' Title
#'
#' @param X 
#' @param groups 
#' @param eta 
#' @param v0 
#' @param s0 
#' @param tune 
#' @param steps 
#' @param burn 
#' @param skip 
#'
#' @return
#' @export
#'
#' @examples
dhlcModel <- function(X,groups,eta,v0,s0,tune,K,steps,burn,skip) { 
  X = X - 1 
  groups = groups-1 
  G = length(unique(groups)) 
  posterior <- dhlc_cpp(as.matrix(X),as.vector(groups),eta,v0,s0,tune,K,G,steps)
  out <- seq(from=burn,to=steps,by=skip)
  posterior$out = out 
  return(posterior)
}

#' Title
#'
#' @param data 
#' @param groups 
#' @param eta 
#' @param alpha 
#' @param steps 
#' @param burn 
#' @param skip 
#'
#' @return
#' @export
#'
#' @examples
hlcModel <- function(X,groups,eta,alpha,steps,burn,skip) { 
  X = X - 1 
  groups = groups-1 
  posterior <- hlc_cpp(as.matrix(X),as.vector(groups),eta,as.matrix(alpha),steps)
  out <- seq(from=burn,to=steps,by=skip)
  posterior$out = out 
  return(posterior)
}





