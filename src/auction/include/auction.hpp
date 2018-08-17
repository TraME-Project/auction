/** ----------------------------------------------------------------------------
Copyright (C) 2016 Joseph D Walsh III <math@jdwalsh03.com>

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.

This file is part of the "AUCTION ALGORITHMS IN C++" software project. See the
document <filelist.txt> for a full list of the project files. If
<filelist.txt> or any other file is missing, go to <http://www.jdwalsh03.com/>
to download the complete project.

This material is based upon work supported by the National Science Foundation
Graduate Research Fellowship under Grant No. DGE-1148903. Any opinion,
findings, and conclusions or recommendations expressed in this material are
those of the author and do not necessarily reflect the views of the National
Science Foundation.
---------------------------------------------------------------------------- **/

//
// Modifications of functions in 'main.cpp'
//
// Keith O'Hara
//

#include <algorithm>  // std::find
#include <cmath>      // std::abs, std::sqrt
#include <cstdlib>    // exit, EXIT_FAILURE
#include <cstdio>     // std::sprintf
#include <fstream>    // std::ifstream, std::ofstream
#include <iostream>   // std::cout, std::endl, std::ios::out, std::ios::trunc
#include <limits>     // std::numeric_limits
#include <random>     // std::mt19937_64, std::uniform_int_distribution
#include <sstream>    // std::stringstream
#include <string>     // std::string, std::getline
#include <vector>     // std::vector
#include "glob.hpp"   // mfloat, mfvec, mint, uint_t
#include "object.hpp" // objlist, voblist
#include "apmap.hpp"  // APmap
#include "gamap.hpp"  // GAmap
#include "sopmap.hpp" // SOPmap
#include "somap.hpp"  // SOmap

#ifdef USE_RCPP_ARMADILLO
    #include <RcppArmadillo.h>
#else
    #ifndef ARMA_DONT_USE_WRAPPER
        #define ARMA_DONT_USE_WRAPPER
    #endif
    #include "armadillo"
#endif

int auction(const arma::mat& Phi, const arma::vec& demand_vec, const arma::vec& supply_vec, 
            arma::mat& solution_mat, double& primal_cost, arma::mat& dual_mat, double& dual_cost, 
            const int algo_choice, const bool max_prob, mfloat eps = -1.0, mfloat eps_min = -1.0,
            const bool run_bellman_ford = true, const bool bf_upper_prices = true, const bool verbose = true);

//

inline
void
APrun(const mfvec &DWT, const mfvec &SWT, const voblist &A, const mfloat MX, const mfloat MN, const mfloat ST, objlist &T, mfvec &PR)
{
    APmap apslv(DWT, SWT, A, MX, MN, ST);
    apslv.Solve(T, PR);
}

inline
void
GArun(const mfvec &DWT, const mfvec &SWT, const voblist &A, const mfloat MX, const mfloat MN, const mfloat ST, objlist &T, mfvec &PR)
{
    GAmap gslv(DWT, SWT, A, MX, MN, ST);
    gslv.Solve(T, PR);
}

inline
void
SOPrun(const mfvec &DWT, const mfvec &SWT, const voblist &A, const mfloat MX, const mfloat MN, const mfloat ST, objlist &T, mfvec &PR)
{
    SOPmap sopslv(DWT, SWT, A, MX, MN, ST);
    sopslv.Solve(T, PR);
}

inline
void
SOrun(const mfvec &DWT, const mfvec &SWT, const voblist &A, const mfloat MX, const mfloat MN, const mfloat ST, objlist &T, mfvec &PR)
{
    SOmap soslv(DWT, SWT, A, MX, MN, ST);
    soslv.Solve(T, PR);
    return;
}

inline
mfloat
Dual(const mfvec &DWT, const mfvec &SWT, const voblist &ARX, const mfvec &PR, arma::mat& dual_mat, const bool max_prob)
{
    mfloat dcst = 0.0;
    mfvec M;

    for (size_t i = 0; i < DWT.size(); i++)
    {
        M.clear();

        for (size_t j = 0; j < ARX[i].size(); j++)
        {
            M.push_back(ARX[i][j].c - PR[uint_t(ARX[i][j].j)]);
        }

        dual_mat(i,0) = *std::max_element(M.begin(), M.end());

        dcst += DWT[i] * dual_mat(i,0);
    }

    for (size_t it = 0; it < SWT.size(); it++)
    {
        dual_mat(it + DWT.size(),0) = PR[it];

        dcst += SWT[it] * PR[it];
    }

    if (max_prob == false)
    {
        dcst = - dcst;
        dual_mat = - dual_mat;
    }

    return dcst;
}

inline
mfloat
Primal(const voblist& ARX, const objlist& T, arma::mat& solution_mat, const bool max_prob)
{
    mfloat pcst = 0.0;
    uint_t tst;

    for (uint_t it = 0; it < T.size(); it++)
    {
        tst = 0;
        while ((tst < ARX[uint_t(T[it].i)].size()) && (ARX[uint_t(T[it].i)][tst].j != T[it].j))
        {
            tst++;
        }

        if (tst < ARX[uint_t(T[it].i)].size())
        {
            pcst += ARX[uint_t(T[it].i)][tst].c * T[it].c;

            solution_mat(T[it].i,T[it].j) = T[it].c;
        }
        else
        {
            std::cout << "Overflow in cost calculation" << std::endl;
        }
    }

    if (max_prob == false) {
        pcst = - pcst;
    }

    return pcst;
}
