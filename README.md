# Auction.Rcpp

An R interface to the Auction C++ library, written by Joseph D. Walsh III, which can be obtained [here](http://jdwalsh03.com/coding.html).

This version appends an additional step to compute lattice bounds using a residual network approach.

## Installation

The quickest way to install Auction.Rcpp is via the devtools package:
``` R
install.packages("devtools")
devtools:::install_github("TraME-Project/auction")
```

Note that Auction.Rcpp requires compilation, so an appropriate development environment is necessary to install the package.
* Windows users should get [Rtools](https://cran.r-project.org/bin/windows/Rtools/). Please ensure that R and Rtools are installed to `C:\` (and not `C:\Program Files`), and that the PATH variables are set correctly (described during the Rtools installation process).
* Mac users should install Xcode and then check [here](https://cran.r-project.org/bin/macosx/tools/) for additional tools (Clang6 and gfortran).

## Example

Solving a simple assignment problem using the General Auction algorithm:
``` R
library(auction.Rcpp)

n <- 3
m <- 2

cost <- matrix(c(3,6,7,2,8,2), nrow = n)

p  <-  rep(1,n)
q  <-  rep(1,m)

q[1] <- 2

auction(cost,q,p,1,FALSE)
```

A full description of the `auction` function can be found on the help page `?auction`.

### R Package Authors

Keith O'Hara

### License

GPL (>= 2) 
