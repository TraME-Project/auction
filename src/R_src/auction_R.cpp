/*################################################################################
  ##
  ##   Copyright (C) 2017 Keith O'Hara
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

SEXP auction_R(SEXP Phi_R, SEXP demand_R, SEXP supply_R, SEXP algorithm_R, SEXP max_prob_R)
{
    try {
        arma::mat Phi = as<arma::mat>(Phi_R);
        arma::vec demand = as<arma::vec>(demand_R);
        arma::vec supply = as<arma::vec>(supply_R);

        int algo_choice = as<int>(algorithm_R);
        bool max_prob = as<bool>(max_prob_R);

        arma::mat solution_mat, dual_mat;
        
        auction(Phi,demand,supply,solution_mat,dual_mat,algo_choice,max_prob);

        double primal_cost = arma::accu(solution_mat%Phi);
        double dual_cost = arma::as_scalar(dual_mat.t()*arma::join_cols(supply,demand));
        //
        return Rcpp::List::create(Rcpp::Named("solution") = solution_mat, Rcpp::Named("dual") = dual_mat, Rcpp::Named("primal_cost") = primal_cost, Rcpp::Named("dual_cost") = dual_cost);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
