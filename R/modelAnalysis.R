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
  data <- data -1 
  groups <- groups-1 
  lik <- posteriorLikelihood(data,groups,pi,beta)
  print(lik)
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
  posterior_mean$beta <- posteriorMeanForBeta(post$beta,post$out)
  posterior_mean$pi <- posteriorMean(post$pi,post$out)
  K <- ncol(posterior_mean$pi)
  posterior_mean$z_prob <- posteriorProbZ(post$Z,post$out,K)
  posterior_mean$z_assign <- apply(posterior_mean$z_prob,MARGIN=1,FUN=which.max)
  return(posterior_mean) 
}

#' Plot time specific class probabilies, possibly vs another time series for LDA-DS
#'
#' @param data TxM DataFrame to plot over time 
#' @param dates T-length vector of dates
#' @param path location to save plots 
#' @param save whether or not to save the plot 
#'
#' @return Displays plot 
#' @export
plotPis <- function(data,dates,path="") { 
  data$date <- dates 
  df <- reshape2::melt(df,id.vars=c("date"))
  ggplot2::ggplot(df,ggplot2::aes(x="date",y="value",color="variable"))+ggplot2::geom_line() 
  if(save) { 
    ggplot2::ggsave(paste(path,"pi_comp.pdf",sep="")) 
  } 
}

#' Plot class-specific distributions over responses in bar chart 
#'
#' @param betas J-length list of K x L_j matrices of class-specific distributions over responses 
#' @param response_codes Optional J length list, each containing a vector of response codes for each question to plot
#' @param questions Optional J length vector of question labels to title each graph 
#' @param path location to save plots 
#'
#' @return Saves plots of beta for each question to path 
#' @export 
plotBetas <- function(betas,response_codes=NULL,questions=NULL,path="") { 
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
      if(L_j==2) { 
        x = 0:1
      }
      else{ 
        x= 1:L_j 
      }
    }
    if(!is.null(questions)) { 
      title = questions[j] 
    }
    else{
      title=""
    }
    beta_df$response = factor(x)
    data_long = reshape2::melt(beta_df,id.vars=c("response"))
    colnames(data_long) = c("Response","Class","Probability")
    filename  = paste("beta",j,sep="") 
    filename = paste(title,".pdf",sep="")
    ggplot2::ggplot(data_long,ggplot2::aes(x="Response",y="Probability",fill="Class")) + 
      ggplot2::geom_bar(stat='identity', position='dodge')+ggplot2::ggtitle(title)
    ggplot2::ggsave(filename) 
  }
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
