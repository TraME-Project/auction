################################################################################
##
##   Copyright (C) 2017 Keith O'Hara
##
##   This file is part of the R package auction.
##
##   The R package auction is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 2 of the License, or
##   (at your option) any later version.
##
##   The R package auction is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with auction. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

auction <- function(Phi,demand,supply,algorithm=1,max_prob=TRUE)
{
    #
    res <- .Call("auction_R", Phi,demand,supply,algorithm,max_prob, PACKAGE = "auction.Rcpp")
    #
    return(res)
}
