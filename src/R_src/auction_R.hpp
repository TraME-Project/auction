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

#ifndef _auction_R_HPP
#define _auction_R_HPP

#include "auction.hpp"

RcppExport SEXP auction_R(SEXP Phi_R, SEXP demand_R, SEXP supply_R, 
                          SEXP algorithm_R, SEXP max_prob_R, 
                          SEXP eps_init_R, SEXP eps_min_R);

#endif