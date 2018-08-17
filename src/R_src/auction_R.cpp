/*################################################################################
  ##
  ##   Copyright (C) 2017-2018 Keith O'Hara
  ##
  ##   This file is part of auction.
  ##
  ##   auction is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   auction is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with auction. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

#include "auction_R.hpp"
using namespace Rcpp;

SEXP auction_R(SEXP Phi_R, SEXP demand_R, SEXP supply_R, SEXP algorithm_R, SEXP max_prob_R, SEXP eps_init_R, SEXP eps_min_R,
               SEXP run_bellman_ford_R, SEXP bf_upper_prices_R, SEXP verbose_R)
{
    try
    {
        double primal_cost, dual_cost;
        arma::mat solution_mat, dual_mat;
        
        auction(as<arma::mat>(Phi_R),as<arma::vec>(demand_R),as<arma::vec>(supply_R),
                solution_mat,primal_cost,dual_mat,dual_cost,
                as<int>(algorithm_R),as<bool>(max_prob_R),
                as<mfloat>(eps_init_R),as<mfloat>(eps_min_R),
                as<bool>(run_bellman_ford_R),as<bool>(bf_upper_prices_R),as<bool>(verbose_R));
        
        return Rcpp::List::create(Rcpp::Named("solution") = solution_mat, Rcpp::Named("dual") = dual_mat, 
                                  Rcpp::Named("primal_cost") = primal_cost, Rcpp::Named("dual_cost") = dual_cost);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "auction: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
