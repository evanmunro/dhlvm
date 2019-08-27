#' Title
#'
#' @param posterior 
#'
#' @return
#' @export
#'
#' @examples
#' @export 
calculateLikelihood <- function(posterior,data) { 
  out <- posterior$out 
  zs <- posterior$Z[,out]
  zs <- lapply(1:ncol(zs), function(j) zs[,j]) 
  betas <- posterior$beta[out] 
  callCPPLik <- function(z,beta) { 
    return(posteriorLikelihood(data,z,beta))
  }
  lik <- mapply(FUN=callCPPLik,zs,betas)
  posterior$likelihood <- lik 
 
  print(2*(mean(lik)-var(lik)))
  return(lik)
}


#' Title
#'
#' @param posterior 
#'
#' @return
#' @export
#'
#' @examples
#' @export 
aicm <- function(posterior) {
  
  likelihood <- posterior$likelihood[posterior$out] 
  aic <- 2*(mean(likelihood)-var(likelihood)) 
  return(aic) 
}


#' Title
#'
#' @param G 
#' @param Lj 
#'
#' @return
#' @export
#'
#' @examples
#' @export 
checkIdCondition<- function(G,Lj) { 
  JL_1 = sum(Lj-1) 
  condition = G*(JL_1+1)/(G+JL_1) 
  print(paste("max K is ",floor(condition),sep="")) 
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

posteriorProbZ <- function(Z,out,K) { 
  Z_counts <- t(apply(Z[,out],MARGIN=1,FUN=countMultipleInstances,K)) 
  print(K)
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
  K <- ncol(posterior_mean$pi)
  posterior_mean$z_prob <- posteriorProbZ(post$Z,post$out,K)
  posterior_mean$z_assign <- apply(posterior_mean$z_prob,MARGIN=1,FUN=which.max)
  return(posterior_mean) 
}

#' Title
#'
#' @param x 
#' @param y 
#' @param z 
#'
#' @return
#' @export
#'
#' @examples
plotScatter <- function(x,y,z)  { 
  df <- data.frame(x=x,y=y,z=factor(z)) 
  ggplot2::ggplot(df,ggplot2::aes(x=x,y=y,color=z))+ggplot2::geom_point()
 # ggplot2::ggsave("scatter.pdf")
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
plotPis <- function(data,withRecessions=F) { 
  data <- reshape2::melt(data,id.vars=c("dates"))
  g <- ggplot2::ggplot(data)+ggplot2::geom_line(aes(x=dates,y=value,color=variable)) 
  if(withRecessions==T) { 
    recessions.df = read.table(textConnection(
        "Peak, Trough
        1857-06-01, 1858-12-01
        1860-10-01, 1861-06-01
        1865-04-01, 1867-12-01
        1869-06-01, 1870-12-01
        1873-10-01, 1879-03-01
        1882-03-01, 1885-05-01
        1887-03-01, 1888-04-01
        1890-07-01, 1891-05-01
        1893-01-01, 1894-06-01
        1895-12-01, 1897-06-01
        1899-06-01, 1900-12-01
        1902-09-01, 1904-08-01
        1907-05-01, 1908-06-01
        1910-01-01, 1912-01-01
        1913-01-01, 1914-12-01
        1918-08-01, 1919-03-01
        1920-01-01, 1921-07-01
        1923-05-01, 1924-07-01
        1926-10-01, 1927-11-01
        1929-08-01, 1933-03-01
        1937-05-01, 1938-06-01
        1945-02-01, 1945-10-01
        1948-11-01, 1949-10-01
        1953-07-01, 1954-05-01
        1957-08-01, 1958-04-01
        1960-04-01, 1961-02-01
        1969-12-01, 1970-11-01
        1973-11-01, 1975-03-01
        1980-01-01, 1980-07-01
        1981-07-01, 1982-11-01
        1990-07-01, 1991-03-01
        2001-03-01, 2001-11-01
        2007-12-01, 2009-06-01"), sep=',',
        colClasses=c('Date', 'Date'), header=TRUE)
        recessions.trim = subset(recessions.df, Peak >= min(data$dates) )
        g<- g+ geom_rect(data=recessions.trim, 
              aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.2)
  }
  g +theme_bw() 
  ggplot2::ggsave("pi_ev.pdf",width=7,height=3)
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
plotPisComp <- function(pi.ev,dates,comp) { 
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
    ggplot2::ggplot(data_long,ggplot2::aes(x=Response,y=Probability,fill=Class)) + 
      ggplot2::geom_bar(stat='identity', position='dodge')+ggplot2::ggtitle(title)
    ggplot2::ggsave(filename) 
  }
}