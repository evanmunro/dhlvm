#Author: Evan Munro, Stanford University 
#This file contains the Gibbs Samplers for Bayesian hierarchical latent variable models for categorical survey data 
#that are described in a working paper with Serena Ng 

# helper utility: convert responses count to list format
countToListInput <- function(X) {
  Time <- nrow(X)
  V <- ncol(X)
  result <- list()
  for (t in 1:Time) {
    result[[t]] = rep(1:V,times=X[t,])
  }
  return(result)
}

#' Old HLC model with J=1
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
discreteMSM <- function(data,eta,alpha,K,steps,burn,skip) {
  posterior <- discreteMS_cpp(data,eta,alpha,K,steps)
  out <- seq(from=burn,to=steps,by=skip)
  posterior$out = out 
  return(posterior)
} 

#' Run MCMC sampler for LDA-DS Model 
#'
#' @param X N x J matrix of responses 
#' @param groups N x 1 vector of group membership value in 1 to G 
#' @param eta J length list of K x L_j matrices, Dirichlet prior for beta
#' @param v0 Gamma shape prior 
#' @param s0 Gamma rate prior 
#' @param tune Parameter a in Welling and Teh 2011 for shrinkage of gradient descent steps and variance of noise 
#' @param K Number of classes 
#' @param steps Number of MCMC steps 
#' @param burn Burn-in iterations for sampler
#' @param skip Thinning iterations for sampler 
#'
#' @return List of full posterior sample for each parameter in the model 
#' @export
dhlcModel <- function(X,groups,eta,v0,s0,tune,K,steps,burn,skip) { 
  X = X - 1 
  groups = groups-1 
  G = length(unique(groups)) 
  posterior <- dhlc_cpp(as.matrix(X),as.vector(groups),eta,v0,s0,tune,K,G,steps)
  out <- seq(from=burn,to=steps,by=skip)
  posterior$out = out 
  return(posterior)
}

#' Run MCMC sampler for LDA-S Model 
#'
#' @param X N x J matrix of responses
#' @param groups N x 1 vector of group membership value in 1 to G 
#' @param eta J length list of K x L_j matrices, Dirichlet prior for beta
#' @param alpha G x K matrix, Dirichlet prior for pi 
#' @param steps Number of MCMC steps 
#' @param burn Burn-in iterations for sampler
#' @param skip Thinning iterations for sampler
#'
#' @return List of full posterior sample for each parameter in the model 
#' @export
hlcModel <- function(X,groups,eta,alpha,steps,burn,skip) { 
  X = X - 1 
  groups = groups-1 
  posterior <- hlc_cpp(as.matrix(X),as.vector(groups),eta,as.matrix(alpha),steps)
  out <- seq(from=burn,to=steps,by=skip)
  posterior$out = out 
  return(posterior)
}
