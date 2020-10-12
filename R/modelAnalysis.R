
#' Calculate likelihood for a single group of responses 
#' @param data N xJ response matrix 
#' @param pi 1x K vector of group-specific distribtion over classes
#' @param beta J-length list of K x L_j class specific response distributions 
#' 
#' @return Model ikelihood for provided group of data conditional on provided parameters
#' @export 
likelihood <- function(data,groups,pi,beta,select=NULL) { 
  if (!is.null(select)) {
    return(posteriorLikelihood(as.matrix(data[groups==select,])-1,as.vector(groups[groups==select])-1,pi,beta))
  } else {
    return(posteriorLikelihood(data-1,groups-1,pi,beta)) 
  }
}

#' Calculate approximate BIC for the model 
#'
#' @param data N x J response matrix  
#' @param groups N length vector of group membership indicators in 1,..,G  
#' @param pi  G x K matrix of group-specific distributions over classes, usually posterior mean
#' @param beta J-length list of K x L_j class-specific response distributions, usually posterior mean 
#' @param dynamic Adjusts number of parameters if this is a DHLC model 
#' 
#' @return BIC value for model with provided data and parameters 
#' @export
bic <- function(data,groups,pi,beta,dynamic=F) { 
  data <- as.matrix(data) -1 
  groups <- as.vector(groups)-1 
  lik <- posteriorLikelihood(data,groups,pi,beta)
  print(paste0("likelihood: ", lik))
  J <- length(beta)
  K <- dim(pi)[2]
  G <- dim(pi)[1]
  L <- unlist(lapply(beta, FUN=function(x) return(ncol(x))))
   
  nparams <- G*(K-1) + K*sum(L-1)
  nobs <- nrow(data)
  if(dynamic) {
    nparams = nparams+K 
   }
    
  bic <- -lik + 1/2*(nparams*(log(nobs)))
  return(bic) 
}


#' Check necessary condition for identification 
#'
#' @param G number of groups 
#' @param Lj J length vector with number of responses for each question 
#'
#' @return The maximum number of classes K that meet the identification condition given G, Lj
#' @export
checkIdCondition<- function(G,Lj) { 
  JL_1 = sum(Lj-1) 
  condition = G*(JL_1+1)/(G+JL_1) 
  print(paste("max K is ",floor(condition),sep="")) 
}

#' Check for which responses and questions the grouped response matrix is sparse 
#'
#' @param aux.data 
#' @param limit 
#'
#' @return prints which responses are spares 
#' @export
checkSparsity<- function(aux.data,limit=0.05) { 
  J=ncol(aux.data)
  N = nrow(aux.data) 
  L = apply(aux.data,MARGIN=2,FUN=function(x) return(length(unique(x))))
  for (j in 1:J ) { 
    question = colnames(aux.data)[j]
    for (i in 1:L[j]) { 
      frac = sum(aux.data[,j]==i)/N 
      if(frac < limit ){
        print(paste("Small fraction of",frac,"at",question,"response",i,"out of",L[j])) 
      }
    }
  }
}
#' Get posterior mean of parameters 
#'
#' @param post Posterior sample from running MCMC for LDA-S or LDA-DS 
#'
#' @return List with attribute containing posterior mean for each parameter in the model
#' @export
posteriorMeans <- function(post) { 
  posterior_mean <- list() 
  if (!is.null(post$sigma)){ 
   posterior_mean$sigma <- posteriorMean(post$sigma,post$out) 
  }
  if(!is.null(post$gamma)){
    posterior_mean$gamma <- posteriorMean(post$gamma,post$out)
  }
  posterior_mean$beta <- posteriorMeanForBeta(post$beta,post$out)
  posterior_mean$pi <- posteriorMean(post$pi,post$out)
  K <- ncol(posterior_mean$pi)
  posterior_mean$z_prob <- posteriorProbZ(post$Z,post$out,K)
  posterior_mean$z_assign <- apply(posterior_mean$z_prob,MARGIN=1,FUN=which.max)
  return(posterior_mean) 
}

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

posteriorProbZ <- function(Z,out,K) { 
  Z_counts <- t(apply(Z[,out],MARGIN=1,FUN=countMultipleInstances,K)) 
  return(Z_counts/apply(Z_counts,MARGIN=1,FUN=sum)) 
}
