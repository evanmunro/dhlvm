# Author: Evan Munro, Stanford University
# Date: July 25, 2018
# This file contains functions for each step of the gibbs sampler 
# to estimate two discrete time series models detailed in 
# an upcoming working paper: 1) the discrete markov switching model 
# and 2) the discrete state space model 

library(MCMCpack)

#' Updates state vector given transition matrix, data, beta, and current estimate of state 
#' Follows single-move sampler from Albert and Chib (1996)
#'
#' @param data Time x V counts (word counts, or survey respondent counts, for example)
#' @param P K x K Transition matrix
#' @param S_prev Time-length vector of state estimates in 1:K 
#' @param beta K V-length probability vectors for K multinomial distributions
#' @param Time number of time periods
#' @param K number of markov states
#'
#' @return Time-length vector of new state estimates 
#' @export

sampleState <- function (data,P,S,beta,Time,K) {
  state_prob <- rep(0,K)
  for (t in 1:Time) {
    for (k in 1:K) {
      if (t==1) {
        last_per_prob <- 1
      } else {
        last_per_prob <- P[k,S[t-1]]
      } 
      if (t==Time) {
        next_per_prob <- 1
      } else { 
        next_per_prob <- P[S[t+1],k]
      }
      state_prob[k] <- last_per_prob*dmultinom(data[t,],prob=beta[k,])*next_per_prob   
    }
    if(sum(state_prob)==0) {
      state_prob <-rep(1/K,K)
    }
    else {
      state_prob <- state_prob/sum(state_prob)
    }
    S[t] <- which(rmultinom(1,size=1,prob=state_prob)>0)
  }
  return(S)
}  


sampleStateSS <- function(z,a,sigma, gamma_prev,shrink,Time,K) { 
  gamma_new <- matrix(0,Time,K)
  n <- matrix(0,Time,K)
  print(a)
  for (k in 1:K){ 
    #Count the number of respondents assigned to sentiment k for each doc
    for (t in 1:Time) { 
      #second part of equation 1.3
      noise <- rnorm(n=1,mean=0,sd=shrink)
      Nt <- length(z[[t]])
      n[t,k] <- sum(z[[t]]==k)
      if(t==1) term1 = 0 else term1 = -1/sigma[k,k]*(gamma_prev[t,k] - a[k]*gamma_new[t-1,k]) 
      if(t==Time) term2 = 0 else term2 = -1/sigma[k,k]*(gamma_prev[t+1,k] - a[k]*gamma_prev[t,k])
      gradlogprob <- term1 + term2 + n[t,k] -Nt*exp(gamma_prev[t,k])/sum(exp(gamma_prev[t,]))
      delta <- shrink/2*gradlogprob + noise 
      gamma_new[t,k] = gamma_prev[t,k] + delta
    }
  }
  return(gamma_new)
}

#' Updates state assignments for each response based on current estimate for 
#' state-specific multinomial distributions and state mixtures for each date
#'
#' @param data T-length list of N_t length vectors of response indices for each t
#' @param theta TxK matrix of state mixture probabilities for each t
#' @param beta KxV matrix of state-specific multinomial probabilities 
#' @param Time Number of time periods in data 
#'
#' @return T-length list of Nt vectors containing state assignments for each response
#' @export
#'

sampleZ <- function(data, theta,beta,Time) { 
  z = list()
  for (t in 1:Time)
  {
    Nt <- length(data[[t]])
    z[[t]] <- rep(0,Nt)
    for (n in 1:Nt) {
      v <- data[[t]][n]
      prob <- theta[t,]*beta[,v]/sum(theta[t,]*beta[,v])
      if(sum(prob)==0) {
        prob <-rep(1/K,K)
      }
      z[[t]][n] <- sample(1:K,size=1,prob=prob)
    }
  }
  return(z)   
}

#' Updates transition matrix given current estimate of states 
#'
#' @param S Time length vector of states
#' @param K Number of states
#' @param alpha dirichlet prior for transition matrix columns
#' @param Time number of time periods in the data
#'
#' @return KxK state transition probability matrix
#' @export

sampleTransition <- function(S,alpha,K,Time) {
  N <- matrix(0,nrow=K,ncol=K)
  P <- matrix(0,nrow=K,ncol=K)
  for (t in 2:Time) {
    N[S[t],S[t-1]] <- N[S[t],S[t-1]]+1
  }
  
  for (k in 1:K) {
    P[,k] <- rdirichlet(1,alpha+N[,k])
  }
  return(P)
}

  
#' Updates state-specific parameters for multinomial distribution given data, prior, and 
#' state estimates for each time period 
#' 
#' @param data T x V matrix of counts 
#' @param S T length vector of states
#' @param eta Dirichlet prior on beta, parameters for multinomial distribution
#' @param K Number of states 
#' @param V number of discrete options that data can take
#' @param Time 
#'
#' @return updated estimate of parameter for state-specific multinomial distributions
#' @export 

sampleBeta <- function(data,S,eta,K,V,Time) {
  beta <- matrix(0,nrow=K,ncol=V)
  m <- matrix(0,nrow=K,ncol=V)
  for (t in 1:Time) {
    m[S[t],] = m[S[t],] + data[t,]
  }
  for (k in 1:K) {
    beta[k,] = rdirichlet(1,eta+m[k,])
  }
  return(beta)
}

#' Updates state-specific parameters for multinomial distribution given data, prior, and 
#' state estimates for each response in the data 
#'
#' @param data T-length list of Nt vectors of responses (integers in 1:V)
#' @param z T-length list of Nt vectors of state assignments for each response
#' @param eta V-length vector of prior on beta
#' @param K Number of states
#' @param V Number of unique outcomes
#' @param Time Number of time periods 
#'
#' @return a KxV matrix, updated estimate of parameter for state-specific multinomial dist.
#' @export
#'
sampleBetaSS <- function(data,z,eta,K,V,Time) {
  beta <- matrix(0,nrow=K,ncol=V)
  m = matrix(0,K,V)
  for (k in 1:K) {
    for (t in 1:Time) {
      for (v in 1:V) {
        #second part of equation 1.4
        m[k,v] = m[k,v] + sum((z[[t]]==k)*(data[[t]]==v))			
      }
    }
    #First part of equation 1.4 
    beta[k,] = rdirichlet(1,eta[k,] + m[k,])
  }
  return(beta)
}

sampleSigma <- function(g,a,v0,del0,K,Time) {
  sigma_new <- diag(rep(1,K))
  for (k in 1:K ) {
    y <- g[2:Time,k] 
    x <- g[1:(Time-1),k]
    v1 = v0 + Time -1 
    del1 = del0 + (y-x*a[k])%*%(y-x*a[k]) 
    sigma_new[k,k] = rinvgamma(1,v1/2,del1/2)
  } 
  print(sigma_new)
  return(sigma_new)
}

sampleCoef <- function(g,sigma,b0,s0,K,Time) {
  coef <- rep(1,K)
  for (k in 1:K) {
    print(k)
    y <- g[2:Time,k] 
    x <- g[1:(Time-1),k]
    s1 <- (s0^(-1) + sigma[k,k]^(-1)*(x%*%x))^-1 
    b1 <- s1*(s0^(-1)*b0 + sigma[k,k]^(-1)*x%*%y)
    coef[k] = unitCircleDraw(b1,s1)
  }
  return(coef)
}

unitCircleDraw <- function(mean,sd) { 
  result = rnorm(1,mean,sd)
    while (result > 1) { 
      result = rnorm(1,mean,sd)  
    }
  return(result)
}