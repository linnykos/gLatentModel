set.seed(10)

r1 <- 5; r2 <- 0.2; n <- 1000
ind <- runif(n, min = 0, max = 2*pi)
x <- c(cos(ind)*r1, cos(ind)*r2) + 0.1*rnorm(2*n)
y <- c(sin(ind)*r1, sin(ind)*r2) + 0.1*rnorm(2*n)

plot(x,y, asp = T)

dat <- cbind(x, y)

#do spectral clustering
eig <- eigen(cov(dat))
clus <- kmeans(eig$vectors)

#hahha failed experiment. right. spectral clustering usually takes in a weighted laplacian matrix based on the distances between points