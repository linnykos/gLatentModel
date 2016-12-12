set.seed(10)

d = 5
y = rnorm(d)
X = matrix(rnorm(d^2), d, d); X = X + t(X)

gamma.vec = seq(0.86, 0.87, length.out = 100)
z.mat = sapply(gamma.vec, function(x){
  res = ((2*t(y)%*%X - x)%*%solve(t(X)%*%X))/2
  idx = which(res > 0)

  tmp = (y%*%X[,idx] - x/2) %*% solve(t(X[,idx])%*%X[,idx])
  vec = rep(0, d)
  vec[idx] = tmp
  vec
})
z.vec = apply(z.mat, 2, function(x){sum(x[x>0])})
plot(z.vec)

idx = which.min(abs(z.vec - 1))
gamma = gamma.vec[idx]
z.vec[idx]
beta = z.mat[,idx]
#beta[beta < 0] = 0
#beta = beta/sum(beta)

#########################
sum((y - beta%*%X)^2)
beta2 <- beta/sum(beta)
sum((y - beta2%*%X)^2)
sum((y -  X%*%solve(X%*%X)%*%X%*%y)^2)

###########################
gamma = gamma.vec[idx]
lambda = -2*(y - beta%*%X)%*%X + gamma

abs(beta*lambda) < 1e-4
abs(gamma - 2*(y-beta%*%X)%*%X) >= -1e-4
#abs(beta - (2*t(y)%*%X - gamma)%*%solve(t(X)%*%X)/2) >= -1e-4
lambda > -1e-4

lambda.manual = sapply(1:d, function(x){
  -2*(y-beta%*%X)%*%X[,x] + gamma
})

#################################

## side test
bz = (y%*%X - gamma/2)%*%solve(t(X)%*%X)
-2*(y-bz%*%X)%*%X + gamma

idx = c(1,2)
bz =  (y%*%X[,idx] - gamma/2) %*% solve(t(X[,idx])%*%X[,idx])
bzz = rep(0, d); bzz[c(1,2)] = bz
-2*(y-bzz%*%X)%*%X + gamma

###################################

#given beta, solve for the optimal lambda and gamma
zz = -2*(y-beta%*%X)%*%X
idx = which(beta != 0)
gamma = mean(zz[idx])
lambda = rep(0, d)
lambda[-idx] = zz[-idx] + gamma
