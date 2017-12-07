
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
#include "glob.hpp"   // mfloat, mfvec, mint, muint
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

int auction(const arma::mat& Phi, const arma::vec& demand_vec, const arma::vec& supply_vec, arma::mat& solution_mat, arma::mat& dual_mat, const int algo_choice, const bool max_prob);

//

inline
void APrun(const mfvec &DWT, const mfvec &SWT, const voblist &A, const mfloat MX, const mfloat MN, const mfloat ST, objlist &T, mfvec &PR)
{
    APmap apslv(DWT, SWT, A, MX, MN, ST);
    apslv.Solve(T, PR);
    return;
}

inline
void GArun(const mfvec &DWT, const mfvec &SWT, const voblist &A, const mfloat MX, const mfloat MN, const mfloat ST, objlist &T, mfvec &PR)
{
    GAmap gslv(DWT, SWT, A, MX, MN, ST);
    gslv.Solve(T, PR);
    return;
}

inline
void SOPrun(const mfvec &DWT, const mfvec &SWT, const voblist &A, const mfloat MX, const mfloat MN, const mfloat ST, objlist &T, mfvec &PR)
{
    SOPmap sopslv(DWT, SWT, A, MX, MN, ST);
    sopslv.Solve(T, PR);
    return;
}

inline
void SOrun(const mfvec &DWT, const mfvec &SWT, const voblist &A, const mfloat MX, const mfloat MN, const mfloat ST, objlist &T, mfvec &PR)
{
    SOmap soslv(DWT, SWT, A, MX, MN, ST);
    soslv.Solve(T, PR);
    return;
}

inline
mfloat Dual(const mfvec &DWT, const mfvec &SWT, const voblist &ARX, const mfvec &PR, arma::mat& dual_mat, const bool max_prob)
{
    mfloat dcst = 0.0;
    mfvec M;

    for (muint i = 0; i < DWT.size(); i++)
    {
        M.clear();

        for (muint j = 0; j < ARX[i].size(); j++)
        {
            M.push_back(ARX[i][j].c - PR[muint(ARX[i][j].j)]);
        }

        dual_mat(i,0) = *std::max_element(M.begin(), M.end());

        dcst += DWT[i] * dual_mat(i,0);
    }

    for (muint it = 0; it < SWT.size(); it++)
    {
        dual_mat(it + DWT.size(),0) = PR[it];

        dcst += SWT[it] * PR[it];
    }

    if (max_prob == false) {
        dcst = - dcst;
        dual_mat = - dual_mat;
    }

    return dcst;
}

inline
mfloat Primal(const voblist& ARX, const objlist& T, arma::mat& solution_mat, const bool max_prob)
{
    mfloat pcst = 0.0;
    muint tst;

    for (muint it = 0; it < T.size(); it++)
    {
        tst = 0;
        while ((tst < ARX[muint(T[it].i)].size()) && (ARX[muint(T[it].i)][tst].j != T[it].j))
        {
            tst++;
        }

        if (tst < ARX[muint(T[it].i)].size())
        {
            pcst += ARX[muint(T[it].i)][tst].c * T[it].c;

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

