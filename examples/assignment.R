
library(auction.Rcpp)

n <- 3
m <- 2

cost <- matrix(c(3,6,7,2,8,2), nrow = n)

p  <-  rep(1,n)
q  <-  rep(1,m)

q[1] <- 2

auction(cost,q,p,1,FALSE)

auction(cost,q,p,4,FALSE, eps_init = 0.5, eps_min = 0.01)
auction(cost,q,p,3,FALSE, eps_init = 0.5, eps_min = 0.01)
auction(cost,q,p,1,FALSE)
