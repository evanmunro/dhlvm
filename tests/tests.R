library(dhlvm) 
J=3 
L_j = 2 
N=1000
K=2 
ngroups = 2
groups = c(rep(1,N/2),rep(2,N/2))
beta=list() 
eta= list() 

eta[[1]] = matrix(1,nrow=K,ncol=L_j)
eta[[2]] = eta[[1]]
eta[[3]] = eta[[1]]
beta[[1]] = matrix(c(0.7,0.3,0.3,0.7),nrow=K,ncol=L_j) 
beta[[2]] = beta[[1]]
beta[[3]] = beta[[1]]

alpha = matrix(1,nrow=ngroups,ncol=K)
pi = matrix(c(0.2,0.8,0.8,0.2),nrow=ngroups,ncol=K)

X = matrix(0,nrow=N,ncol=J) 

for (i in 1:N) { 
  group = groups[i] 
  z_i = sample(1:K,1,prob=pi[group,])
  for (j in 1:J) { 
    beta_j = beta[[j]]
    prob = beta_j[z_i,] 
    X[i,j] = sample(1:2,1,prob=prob)
  }
}

steps = 1000 
burn = 100 
skip = 10 

posterior = hlcModel(X,groups,eta,alpha,steps,burn,skip)

posteriorMean(posterior$pi,posterior$out)

#test to see if Rcpp and R have same rgamma 

#as scale goes up, then the mean rgamma goes up
#as shape goes up, mean rgamma goes up 

#THE RGAMMA FUNCTION IN CPP IS TAKING Scale, not RATE, as the second parameter! CHECK VS THE ALGO FOR SIGMA  
#sampling from invgamma can sample from 1/gamma(shape,scale)

# I think it is scale in the gibbs sampling procedure 
mean(rgamma(100000,10,scale=1)) 


#Next test that the sigma sampling procedure works (DONE)

sigma = 5

g0 = 0 
Time = 100000 
g = rep(0,Time) 

for (t in 2:Time) { 
  g[t] = g[t-1] + rnorm(1,mean=0,sd=sigma)
}
sum((g[2:length(g)]- g[1:(length(g)-1)])^2)
sourceCpp("src/Models.cpp")
sampleSigma(as.matrix(g),1,1) 

#next test the full dynamic sampling procedure 

library(dhlvm)
library(Rcpp)
J=3 
L_j = 2 
Nt = 1000
N=Time*Nt
K=2 
Time = 100
sigma = 0.1
groups = c(rep(1:Time,each=Nt))
beta=list() 
eta= list() 

eta[[1]] = matrix(1,nrow=K,ncol=L_j)
eta[[2]] = eta[[1]]
eta[[3]] = eta[[1]]
beta[[1]] = matrix(c(0.7,0.3,0.3,0.7),nrow=K,ncol=L_j) 
beta[[2]] = beta[[1]]
beta[[3]] = beta[[1]]

z = rep(0,N)
pi = matrix(0,nrow=Time,ncol=K)
gammas = matrix(0,nrow=Time,ncol=K)

gammas[1,] = rnorm(K,0,sigma)
pi[1,]= exp(gammas[1,])/sum(exp(gammas[1,]))

X = matrix(0,nrow=N,ncol=J) 

for (t in 2:Time) { 
  for (k in 1:K) { 
    gammas[t,k] = rnorm(1,gammas[t-1,k],sigma)
  }
  pi[t,] = exp(gammas[t,])/sum(exp(gammas[t,]))

}

for (i in 1:N) { 
  time = groups[i] 
  z[i] = sample(1:K,1,prob=pi[time,])
  for (j in 1:J) { 
    beta_j = beta[[j]]
    prob = beta_j[z[i],] 
    X[i,j] = sample(1:2,1,prob=prob)
  }
}

v0=10
s0=1
steps = 1000
burn = 100
skip = 10
tune=0.05
posterior = dhlcModel(X,groups,eta,v0,s0,tune,K,steps,burn,skip)

post.pi = posteriorMean(posterior$pi,posterior$out)
plot(post.pi[,1],type="l")
post.sigma = posteriorMean(posterior$sigma,posterior$out) 
#post.gamma = posteriorMean(posterior$gamma,posterior$out) 


plot(pi[,1]) 

post.pi[1:10,]
pi[1:10,]

#test beta
sampleBeta(X-1,z-1,eta) 

#test Z 
sampleZ(X-1,groups-1,pi,beta)

#test gamma 
sourceCpp("src/Models.cpp")
sgam =sampleGamma(groups-1,z-1,gammas,diag(rep(sigma,K)),0.05)

#test sigma 
sampleSigma(gammas,v0,s0)

pi[1,]
exp(sgam[1,])/sum(exp(sgam[1,]))
pi[20,]
exp(sgam[20,])/sum(exp(sgam[20,]))
pi[50,]
exp(sgam[50,])/sum(exp(sgam[50,]))

pi[100,]
exp(sgam[100,])/sum(exp(sgam[100,]))


# test old function 
J=1
L_j = 3
Nt = 100
N=Time*Nt
K=2 
Time = 100
sigma = 0.2
groups = c(rep(1:Time,each=Nt))
beta = matrix(0,nrow=K,ncol=L_j) 
eta = matrix(1,nrow=K,ncol=L_j)
beta = matrix(c(0.7,0.2,0.1,0.1,0.2,0.7),nrow=K,ncol=L_j)


z = rep(0,N)
pi = matrix(0,nrow=Time,ncol=K)
gammas = matrix(0,nrow=Time,ncol=K)

gammas[1,] = rnorm(K,0,sigma)
pi[1,]= exp(gammas[1,])/sum(exp(gammas[1,]))

X = matrix(0,nrow=Time,ncol=L_j) 

for (t in 2:Time) { 
  for (k in 1:K) { 
    gammas[t,k] = rnorm(1,gammas[t-1,k],sigma)
  }
  pi[t,] = exp(gammas[t,])/sum(exp(gammas[t,]))
}

for (i in 1:N) { 
  time = groups[i] 
  z[i] = sample(1:K,1,prob=pi[time,])
  prob = beta[z[i],] 
  response = sample(1:L_j,1,prob=prob) 
  X[time,response]= X[time,response]+1 
}

v0=10
s0=1
steps = 1000
burn = 200
skip = 10
tune=0.1

posterior = discreteLDSModel(X,eta,v0,s0,tune,K,steps,burn,skip)

pi.out = t(posteriorMean(posterior$theta,posterior$out))
matplot(data.frame(truth=pi[,1],sampled=pi.out[,1]),type="l")
