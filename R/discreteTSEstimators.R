#Author: Evan Munro, Stanford University 
#This file contains the Gibbs Samplers for the various latent variable models for discrete time series
#that are being developed as part of a working paper with Serena Ng

#' Runs Gibbs Sampler for markov switching model for discrete time series data 
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
discreteMS <- function(data,eta,alpha,K,steps,burn,skip) {
  Time <- nrow(data)
  V <- ncol(data)
  
  #containers for gibbs sampled parameters
  P_samp <- array(0,dim=c(steps,K,K))
  S_samp <- array(0,dim=c(steps,K,Time)) 
  beta_samp <- array(0,dim=c(steps,K,V))
  
  #initialize each of the parameter samples based on draws from prior
  for (k in 1:K) {
    beta_samp[1,k,] <- rdirichlet(1,eta[k,])
  }
  P_samp[1,,] <- rdirichlet(K,alpha) 
  for (t in 1:Time) {
    S_samp[1,1,t] <- which(rmultinom(1,size=1,prob=c(rep(0.5,K)))>0)
  }
  
  #run gibbs sampler
  for (s in 2:steps) {
    S_samp[s,,] <- sampleState_cpp(data,P_samp[s-1,,],S_samp[s-1,1,],beta_samp[s-1,,])
    P_samp[s,,] <- sampleTransition_cpp(S_samp[s,1,],alpha)
    beta_samp[s,,] <- sampleBeta_cpp(data,S_samp[s,1,],eta)
  }
  
  out <- seq(from=burn,to=steps,by=skip)
  beta_mean <- apply(beta_samp[out,,],MARGIN=c(2,3),FUN=mean)
  P_mean <-apply(P_samp[out,,],MARGIN=c(2,3),FUN=mean) 
  S_mean <- apply(S_samp[out,,],MARGIN=c(2,3),FUN=mean)
#  for (k in 1:K) {
#    for (t in 1:Time) {
#      S_mean[k,t] <- sum(S_samp[out,t]==k)
#    }
#  }
#  S_mean = t(S_mean)/colSums(S_mean)
  return(list(beta_mean,P_mean,S_mean))
}


discreteStateSpace <- function(data,eta,v0,del0,b0,s0,tune,K,steps,burn,skip) {
  V <- ncol(data)
  data <- countToListInput(data)
  Time <- length(data)
  
  #containers for gibbs sampled parameters
  beta_samp <- array(0,dim=c(steps,K,V))
  gamma_samp <- array(0,dim=c(steps,Time,K))
  theta_samp <- array(0,dim=c(steps,Time,K))
  sigma_samp <- array(0,dim=c(steps,K,K))
  a_samp <- array(0,dim=c(steps,K))
  
  sigma_samp[1,,] <- diag(1/rgamma(K,v0,del0))
  a_samp[1,] <- c(1,1)
  for (k in 1:K) {
    #a_samp[1,k] = unitCircleDraw(b0,s0) 
    beta_samp[1,k,] <- rdirichlet(1,eta[k,])
  }
  gamma_samp[1,,] <- mvrnorm(Time,rep(0,K),diag(rep(1,K)))
  theta_samp[1,,] <- softmax(gamma_samp[1,,])
  plot(theta_samp[1,,1],type="l")
  for (s in 2:steps) {
    z = sampleZ(data,theta_samp[s-1,,],beta_samp[s-1,,],Time)
    shrink <- tune*(1+s)^-0.5
    gamma_samp[s,,] <- sampleStateSS(z,a_samp[s-1,],sigma_samp[s-1,,],gamma_samp[s-1,,],shrink,Time,K)
    #a_samp[s,] <- sampleCoef(gamma_samp[s,,],sigma_samp[s-1,,],b0,s0,K,Time)
    a_samp[s,] = c(1,1)
    sigma_samp[s,,] <- sampleSigma(gamma_samp[s,,],a_samp[s,],v0,del0,K,Time)
    print(sigma_samp[s,,])
    theta_samp[s,,] <- softmax(gamma_samp[s,,])
    beta_samp[s,,] <- sampleBetaSS(data,z,eta,K,V,Time)
    plot(theta_samp[s,,1],type="l")
  }

  out <- seq(from=burn,to=steps,by=skip)
  beta_mean <- apply(beta_samp[out,,],MARGIN=c(2,3),FUN=mean)
  theta_mean <- apply(theta_samp[out,,],MARGIN=c(2,3),FUN=mean)
  sigma_mean <- apply(sigma_samp[out,,],MARGIN=c(2,3),FUN=mean)
  a_mean <- apply(a_samp[out,],MARGIN=c(2),FUN=mean)
  return(list(beta_mean,theta_mean,sigma_mean,a_mean))
  
}

#convert responses count to list format 
countToListInput <- function(X) {
  Time <- nrow(X)
  V <- ncol(X)
  result <- list() 
  for (t in 1:Time) {
    result[[t]] = rep(1:V,times=X[t,]) 
  }
  return(result)
}





