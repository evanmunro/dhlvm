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

posteriorMeanForBeta <- function(betas,out) { 
  posterior_to_avg <- betas[out]
  J = length(betas[[1]])  
  expected_values <- list(J) 
  for (j in 1:J) { 
    get_j_matrices <- lapply(posterior_to_avg,FUN = function(x,j) {return(x[[j]])},j)
    expected_values[[j]] <- Reduce("+",get_j_matrices)/length(get_j_matrices) 
  }
  return(expected_values)
}

countInstance <- function(k,x) { 
  return(sum(x==k))
}

countMultipleInstances<- function(x,K) { 
  return(sapply(0:(K-1),countInstance,x))
}

posteriorProbZ <- function(Z,out) { 
  K=2
  postZ <- Z[1:10,out]
  Z_counts <- t(apply(Z,MARGIN=1,FUN=countMultipleInstances,K)) 
  return(Z_counts/apply(Z_counts,MARGIN=1,FUN=sum)) 
}

#' Title
#'
#' @param post 
#' @param plotPis 
#'
#' @return
#' @export
#'
#' @examples
posteriorMeans <- function(post) { 
  posterior_mean <- list() 
  if (!is.null(post$sigma)){ 
    print("Sigma Posterior Mean:") 
    print(posterior_mean$sigma <- posteriorMean(post$sigma,post$out)) 
  }
  posterior_mean$beta <- posteriorMeanForBeta(post$beta,post$out)
  posterior_mean$pi <- posteriorMean(post$pi,post$out)
  posterior_mean$z_prob <- posteriorProbZ(post$Z,post$out)
  posterior_mean$z_assign <- apply(posterior_mean$z_prob,MARGIN=1,FUN=which.max)
  return(posterior_mean) 
}

plotScatter <- function(x,y,z)  { 
  df <- data.frame(x=x,y=y,z=z) 
  pdf(file="scatter.pdf",width=400,height=300)
  ggplot2::ggplot(df,ggplot2::aes(x=x,y=y,color=z))+ggplot2::geom_point()
  dev.off() 
}

#' Title
#'
#' @param pi.ev 
#' @param dates 
#' @param gt 
#' @param withRecessions 
#'
#' @return
#' @export
#'
#' @examples
plotPis <- function(pi.ev,dates,withRecessions=F) { 
  df <- data.frame(x=dates,y=pi.ev) 
  df <- reshape2::melt(df,id.vars=c("x"))
  ggplot2::ggplot(df,ggplot2::aes(x=x,y=value,color=variable))+ggplot2::geom_line() 
  ggplot2::ggsave("pi_ev.pdf")
}

#' Title
#'
#' @param pi.ev 
#' @param dates 
#' @param comp 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
plotPisComp <- function(pi.ev,dates,comp,k) { 
  df <- data.frame(date=dates,y1=pi.ev,y2=comp) 
  df <- reshape2::melt(df,id.vars=c("date"))
  ggplot2::ggplot(df,ggplot2::aes(x=date,y=value,color=variable))+ggplot2::geom_line() 
  ggsave("pi_comp.pdf") 
}

#' Title
#'
#' @param betas 
#' @param rsponse_codes 
#' @param questions 
#'
#' @return
#' @export
#'
#' @examples
plotBetas <- function(betas,response_codes=NULL,questions=NULL) { 
  J = length(betas)
  K = nrow(betas[[1]])
  for (j in 1:J) { 
    beta_mat = t(betas[[j]]) 
    L_j = nrow(beta_mat)
    beta_df = data.frame(beta_mat)
    colnames(beta_df) = paste("K",1:K,sep="") 
    if(!is.null(response_codes)) { 
      x = response_codes[[j]]
    }
    else { 
      x= 1:L_j
    }
    beta_df$response = factor(x)
    data_long = reshape2::melt(beta_df,id.vars=c("response"))
    colnames(data_long) = c("response","class","prob")
    filename  = paste("beta",j,sep="") 
    filename = paste(filename,".pdf",sep="")
    ggplot2::ggplot(data_long,ggplot2::aes(x=response,y=prob,fill=class)) + 
          ggplot2::geom_bar(stat='identity', position='dodge')+ggplot2::ggtitle(paste("Q",j,sep=""))
    ggplot2::ggsave(filename) 
  }
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





