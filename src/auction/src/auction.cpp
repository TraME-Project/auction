
//
// Modification of 'main.cpp'
//

#include "auction.hpp"

//

muint gVBS = 0; // verbosity (0, 1, or 2)

int auction(const arma::mat& Phi, const arma::vec& demand_vec, const arma::vec& supply_vec, arma::mat& solution_mat, double& primal_cost, arma::mat& dual_mat, double& dual_cost, const int algo_choice, const bool max_prob)
{
    const long int NSINK = demand_vec.n_elem;
    const long int NSORC = supply_vec.n_elem;

    bool rev = false; // if NSORC > NSINK

    char str[255];

    int err = 0;
    
    mfloat eps = -1.0;
    mfloat min = -1.0;

    mfloat stp = 0.25;

    mfvec DWT, SWT;
    voblist ARX, RARX;

    dual_mat.zeros(NSORC+NSINK,1);

    // muint gVBS = 0; // verbosity (0, 1, or 2)

    //

    if (NSORC > NSINK) {
        rev = true;
        DWT.reserve(NSORC);
        SWT.reserve(NSINK);
        ARX.resize(NSORC);
        RARX.resize(NSINK);

        solution_mat.zeros(NSORC,NSINK);
    } else {
        DWT.reserve(NSINK);
        SWT.reserve(NSORC);
        ARX.resize(NSINK);
        RARX.resize(NSORC);

        solution_mat.zeros(NSINK,NSORC);
    }

    for (int i=0; i < NSORC; i++) {
        mfloat AMT = supply_vec(i);
        if (rev) {
            DWT.push_back(AMT);
        } else {
            SWT.push_back(AMT);
        }
    }

    for (int j=0; j < NSINK; j++) {
        mfloat AMT = demand_vec(j);
        if (rev) {
            SWT.push_back(AMT);
        } else {
            DWT.push_back(AMT);
        }
    }

    for (int i=0; i < NSORC; i++) {
        for (int j=0; j < NSINK; j++) {
            if (Phi(i,j) != 0.0) {
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

    if (!rev) {

    }

    //

    arma::uvec n_nzcost = arma::find( Phi );

    muint ar = n_nzcost.n_elem;
    mfloat C = arma::abs(Phi.elem( n_nzcost )).max();

    if (eps < gEPS)
    {
        eps = C / 5.0;
    }
    if (min < gEPS)
    {
        min = 1.0 / mfloat(SWT.size());
    }

    objlist T;   // transport plan 
    mfvec PR;    // price vector

    std::chrono::high_resolution_clock::time_point t1, t2;
    std::chrono::duration<mfloat> dur;

    std::cout << "  ----------------------------------------------------------"
              << std::endl;
    std::sprintf(str, "  GRAPH: %lu sinks, %lu sources, %lu arcs", NSINK,
                 NSORC, ar);
    std::cout << str << std::endl;
    std::sprintf(str, "  EPS  : %f starting, %e minimum", eps, min);
    std::cout << str << std::endl;

    //
    // run solvers

    t1 = std::chrono::high_resolution_clock::now();

    if (algo_choice == 1) // general auction
    {
        std::cout << std::endl;
        std::cout << "  General Auction:" << std::endl;

        GArun(DWT, SWT, ARX, eps, min, stp, T, PR);
    }
    else if (algo_choice == 2) // assignment auction
    {
        std::cout << std::endl;
        std::cout << "  Assignment Auction:" << std::endl;

        APrun(DWT, SWT, ARX, eps, min, stp, T, PR);
    }
    else if (algo_choice == 3) // auction - similar objects
    {
        std::cout << std::endl;
        std::cout << "  Auction-SO:" << std::endl;

        SOrun(DWT, SWT, ARX, eps, min, stp, T, PR);
    }
    else if (algo_choice == 4) // auction - similar objects and persons
    {
        std::cout << std::endl;
        std::cout << "  Auction-SOP:" << std::endl;

        objlist S;
        SOPrun(SWT, DWT, RARX, eps, min, stp, S, PR);

        if (S.empty())
        {
            std::cout << "error in Auction-SOP solution" << std::endl;
        }
        else
        {
            T.clear();
            for (muint it = 0; it < S.size(); it++)
            {
                T.emplace_back(S[it].c, S[it].j, S[it].i);
            }
        }
    }
    else
    {
        std::cout << "unrecognized algorithm choice" << std::endl;
        return 1;
    }

    t2 = std::chrono::high_resolution_clock::now();
    dur = std::chrono::duration_cast<std::chrono::duration<mfloat>>(t2 - t1);

    //

    if (T.empty())
    {
        std::cout << "error in solution" << std::endl;
        return 1;
    }
    else
    {
        primal_cost = Primal(ARX, T, solution_mat, max_prob);

        dual_cost = Dual(DWT, SWT, ARX, PR, dual_mat, max_prob);

        std::cout << "     Primal cost = " << primal_cost << ". Dual cost = " << dual_cost << "." << std::endl;
    }

    if (!rev) {
        solution_mat = arma::trans(solution_mat);
    }

    std::cout << std::endl;
    std::sprintf(str, "  Runtime: %13.3f sec", dur.count());
    std::cout << str << std::endl;
    std::cout << "  ----------------------------------------------------------"
              << std::endl;

    // 

    return err;
}
