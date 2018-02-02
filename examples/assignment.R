
library(auction.Rcpp)

n <- 3
m <- 2

cost <- matrix(c(3,6,7,2,8,2), nrow = n)

p  <-  rep(1,n)
q  <-  rep(1,m)

q[1] <- 2

auction(cost,q,p,1,FALSE)
