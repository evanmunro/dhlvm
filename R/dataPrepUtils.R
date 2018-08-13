#Author: Evan Munro, Stanford University
#This file contains utility functions useful for generating simulated data and preparing
#actual survey data for discrete time series models 

plotWithDatesandR <- function(data,lab) {
  library(quantmod)
  library(ggplot2)
  library(reshape2)
  umcsent = getSymbols('UMCSENT',src='FRED', auto.assign=F) 
  umcsent.df = data.frame(date=time(umcsent), umcsent = coredata(umcsent) )
  
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
  
  recessions.trim = subset(recessions.df, Peak >= min(data$date) )
  sentiment.trim  = subset(umcsent.df,date >=min(data$date))
  sentiment.trim  = subset(umcsent.df,date <=max(data$date))
  sentiment.trim$UMCSENT = (sentiment.trim$UMCSENT - min(sentiment.trim$UMCSENT))/(max(sentiment.trim$UMCSENT)- min(sentiment.trim$UMCSENT))
  g = ggplot(sentiment.trim) + geom_line(aes(x=date, y=UMCSENT,colour='red'),alpha=0.5) + theme_bw() +labs(UMCSENT="ICS")
  g = g + geom_rect(data=recessions.trim, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='gray', alpha=0.2) 
  g = g+ geom_line(data=data, aes(x=date,y=toplot,colour='blue'))
  g = g+labs(y="sentiment")
  g = g+  scale_color_discrete(name = "", labels = c(lab, "Michigan ICS"))
  return(g)
}

softmax <- function(x) { 
Time = nrow(x) 
cols = ncol(x)
result = matrix(0,Time,cols) 
for (t in 1:Time) {
  result[t,] = exp(x[t,])/sum(exp(x[t,]))
  
}
return(result)
}
topBetas <- function(beta,identifier,top=5) {
  library(dplyr)
  K <- ncol(beta)
  colnames(identifier) <- c("Permutation")
  beta<- bind_cols(identifier,beta)
  result = beta[order(-beta[,2]),c(1,2)] 
  colnames(result) = c("State 1", "Probabilities 1")
  for (k in 2:K) {
    result = bind_cols(result,beta[order(-beta[,k+1]),c(1,k+1)])
    colnames(result)[(k+1):(k+2)] = c(paste("State ",k,sep=""),
                                        paste("Probabilities ",k,sep=""))
  }
  return(result[1:top,])
}