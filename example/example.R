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

posterior = ldaModel(X,groups,eta,alpha,steps,burn,skip)


