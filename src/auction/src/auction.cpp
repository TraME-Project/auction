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
// Modification of 'main.cpp'
//

#include "auction.hpp"
#include "sp.hpp"

//

uint_t gVBS = 0; // verbosity (0, 1, or 2)

int auction(const arma::mat& Phi, const arma::vec& demand_vec, const arma::vec& supply_vec, 
            arma::mat& solution_mat, double& primal_cost, arma::mat& dual_mat, double& dual_cost, 
            const int algo_choice, const bool max_prob, mfloat eps, mfloat eps_min,
            const bool run_bellman_ford, const bool bf_upper_prices, const bool verbose)
{
    const long int NSINK = demand_vec.n_elem;
    const long int NSORC = supply_vec.n_elem;

    bool rev = false; // if NSORC > NSINK

    char str[255];

    int err = 0;

    if (std::abs(arma::accu(demand_vec) - arma::accu(supply_vec)) > gEPS)
    {
        std::cout << "auction error: sum(demand) != sum(supply)" << std::endl;
        return 1;
    }
    
    // mfloat eps = -1.0;
    // mfloat eps_min = -1.0;

    mfloat stp = 0.25;

    mfvec DWT, SWT;
    voblist ARX, RARX;

    dual_mat.zeros(NSORC+NSINK,1);

    // uint_t gVBS = 0; // verbosity (0, 1, or 2)

    //

    if (NSORC > NSINK)
    {
        rev = true;
        DWT.reserve(NSORC);
        SWT.reserve(NSINK);
        ARX.resize(NSORC);
        RARX.resize(NSINK);

        solution_mat.zeros(NSORC,NSINK);
    }
    else
    {
        DWT.reserve(NSINK);
        SWT.reserve(NSORC);
        ARX.resize(NSINK);
        RARX.resize(NSORC);

        solution_mat.zeros(NSINK,NSORC);
    }

    for (uint_t i=0; i < NSORC; i++)
    {
        mfloat AMT = supply_vec(i);
        if (rev) {
            DWT.push_back(AMT);
        } else {
            SWT.push_back(AMT);
        }
    }

    for (uint_t j=0; j < NSINK; j++)
    {
        mfloat AMT = demand_vec(j);
        if (rev) {
            SWT.push_back(AMT);
        } else {
            DWT.push_back(AMT);
        }
    }

    for (uint_t i=0; i < NSORC; i++)
    {
        for (uint_t j=0; j < NSINK; j++)
        {
            if (Phi(i,j) != 0.0)
            {
                mfloat CST = (max_prob) ? Phi(i,j) : - Phi(i,j);

                if (rev) {
                    ARX[i].emplace_back(CST, i, j);
                    RARX[j].emplace_back(CST, j, i);
                } else {
                    ARX[j].emplace_back(CST, j, i);
                    RARX[i].emplace_back(CST, i, j);
                }
            }
        }
    }

    //

    arma::uvec n_nzcost = arma::find( Phi );

    uint_t ar = n_nzcost.n_elem;
    mfloat C = arma::abs(Phi.elem( n_nzcost )).max();

    if (eps < gEPS)
    {
        eps = C / 5.0;
    }
    if (eps_min < gEPS)
    {
        eps_min = 1.0 / std::min( mfloat(SWT.size()), mfloat(DWT.size()) );
    }

    objlist T;   // transport plan 
    mfvec PR;    // price vector

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<mfloat> dur;

    if (verbose)
    {
        std::cout << "  ----------------------------------------------------------" << std::endl;
        std::sprintf(str, "  GRAPH: %lu sinks, %lu sources, %lu arcs", NSINK, NSORC, ar);
        std::cout << str << std::endl;
        std::sprintf(str, "  EPS  : %f starting, %e minimum", eps, eps_min);
        std::cout << str << std::endl;
    }

    //
    // run solvers

    t1 = std::chrono::high_resolution_clock::now();

    if (algo_choice == 1) // general auction
    {
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "  General Auction:" << std::endl;
        }

        GArun(DWT, SWT, ARX, eps, eps_min, stp, T, PR);
    }
    else if (algo_choice == 2) // assignment auction
    {
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "  Assignment Auction:" << std::endl;
        }

        APrun(DWT, SWT, ARX, eps, eps_min, stp, T, PR);
    }
    else if (algo_choice == 3) // auction - similar objects
    {
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "  Auction-SO:" << std::endl;
        }

        SOrun(DWT, SWT, ARX, eps, eps_min, stp, T, PR);
    }
    else if (algo_choice == 4) // auction - similar objects and persons
    {
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "  Auction-SOP:" << std::endl;
        }

        objlist S;
        SOPrun(SWT, DWT, RARX, eps, eps_min, stp, S, PR);

        if (S.empty())
        {
            std::cout << "error in Auction-SOP solution" << std::endl;
        }
        else
        {
            T.clear();
            for (uint_t it = 0; it < S.size(); it++)
            {
                T.emplace_back(S[it].c, S[it].j, S[it].i);
            }
        }
    }
    else
    {
        std::cout << "auction error: unrecognized algorithm choice" << std::endl;
        return 1;
    }

    t2 = std::chrono::high_resolution_clock::now();
    dur = std::chrono::duration_cast<std::chrono::duration<mfloat>>(t2 - t1);

    //
    // collection solution

    if (T.empty())
    {
        std::cout << "auction error: empty solution" << std::endl;
        return 1;
    }
    else
    {
        primal_cost = Primal(ARX, T, solution_mat, max_prob);
        dual_cost   = Dual(DWT, SWT, ARX, PR, dual_mat, max_prob);

        if (verbose)
        {
            std::cout << "     Primal cost = " << primal_cost << ". Dual cost = " << dual_cost << "." << std::endl;
        }
    }

    if (verbose)
    {
        std::cout << std::endl;
        std::sprintf(str, "  Runtime: %13.3f sec", dur.count());
        std::cout << str << std::endl;
    }

    if (!rev)
    {
        solution_mat = arma::trans(solution_mat);
    }

    //
    // additional step to compute dual prices using the residual network

    double dual_gap = std::abs(primal_cost - dual_cost);

    if (run_bellman_ford && dual_gap > gEPS)
    {
        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "  Duality gap detected: |primal cost - dual cost| = " << dual_gap << std::endl;
            std::cout << "     Running Bellman-Ford on the residual network." << std::endl;
        }

        int nbX = NSINK;
        int nbY = NSORC;

        arma::vec mu_x0 = supply_vec - arma::sum(solution_mat,1);
        arma::vec mu_0y = demand_vec - arma::trans(arma::sum(solution_mat,0));

        //

        arma::mat arcs_xy(nbX*nbY,3);

        for (int i=0; i < nbY; i++)
        {
            for (int j=0; j < nbX; j++)
            {
                arcs_xy(i*nbY + j,0) = j;
                arcs_xy(i*nbY + j,1) = i + nbX;
            }
        }

        arcs_xy.col(2) = - arma::vectorise(Phi);

        //

        arma::mat arcs_x0 = arma::zeros(nbX,3);
        for (int i=0; i < nbX; i++)
        {
            arcs_x0(i,0) = i;
            arcs_x0(i,1) = nbX + nbY;
        }

        arma::mat arcs_0y = arma::zeros(nbY,3);
        for (int i=0; i < nbY; i++)
        {
            arcs_0y(i,0) = nbX + nbY;
            arcs_0y(i,1) = i + nbX;
        }

        arma::mat reduced_graph = arma::join_cols(arma::join_cols(arcs_xy,arcs_x0),arcs_0y);

        //

        arma::uvec select_xy = arma::find( arma::vectorise(solution_mat) );
        arma::uvec select_x0 = arma::find( mu_x0 );
        arma::uvec select_0y = arma::find( mu_0y );

        //

        arma::mat add_arcs_xy = arcs_xy.rows(select_xy);
        add_arcs_xy.swap_cols(0,1);
        add_arcs_xy.col(2) = - add_arcs_xy.col(2);

        reduced_graph = arma::join_cols(reduced_graph,add_arcs_xy);

        if (select_x0.n_elem > 0U)
        {
            arma::mat add_arcs_tmp = arcs_x0.rows(select_x0);
            add_arcs_tmp.swap_cols(0,1);
            add_arcs_tmp.col(2) = - add_arcs_tmp.col(2);

            reduced_graph = arma::join_cols(reduced_graph,add_arcs_tmp);
        }

        if (select_0y.n_elem > 0U)
        {
            arma::mat add_arcs_tmp = arcs_0y.rows(select_0y);
            add_arcs_tmp.swap_cols(0,1);
            add_arcs_tmp.col(2) = - add_arcs_tmp.col(2);

            reduced_graph = arma::join_cols(reduced_graph,add_arcs_tmp);
        }

        if (!bf_upper_prices)
        {
            reduced_graph.swap_cols(0,1);
        }

        //

        int n_nodes = reduced_graph.col(0).max() + 1;
        sp::graph_t arc_list;

        sp::arma_to_graph(n_nodes,reduced_graph,arc_list);

        //

        int source_ind = nbX + nbY;

        std::vector<double> min_distance;
        std::vector<int> path_list;

        bool bf_success = sp::bellman_ford::compute_paths(source_ind, arc_list, min_distance, path_list, nullptr);

        //

        if (bf_success)
        {
            for (int i=0; i < nbX; i++)
            {
                dual_mat.row(i) = min_distance[i];
            }
            for (int i=0; i < nbY; i++)
            {
                dual_mat.row(i+nbX) = -min_distance[nbX+i];
            }

            if (!bf_upper_prices)
            {
                dual_mat = -dual_mat;
            }

            dual_cost = arma::dot(arma::join_cols(demand_vec,supply_vec),dual_mat);

            if (verbose)
            {
                std::cout << "     Bellman-Ford indicated successful completion." << std::endl;
                std::cout << "     New dual cost = " << dual_cost << std::endl;
            }
        }
        else
        {
            if (verbose)
            {
                std::cout << "     Bellman-Ford did not complete successfully." << std::endl;
            }
        }
    }

    if (verbose)
    {
        std::cout << "  ----------------------------------------------------------" << std::endl;
    }

    //

    return err;
}
