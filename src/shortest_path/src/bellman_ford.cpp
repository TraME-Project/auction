/*################################################################################
  ##
  ##   Copyright (C) 2017-2018 Keith O'Hara
  ##
  ##   This file is part of the Shortest Path C++ library (SPLib).
  ##
  ##   SPLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   SPLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

#include "sp.hpp"

//

bool
sp::bellman_ford::compute_paths(const int source_ind, const graph_t& node_list, 
                                std::vector<double>& min_distance, std::vector<int>& path_list,
                                comptime_t* time_spent)
{
    const int n = node_list.size(); // number of nodes

    std::vector<std::vector<double>> dist_mat(n + 2, std::vector<double>(n, max_weight));
    std::vector<int> path_mat(n, -1);

    dist_mat[0][source_ind] = 0;

    //

    bool should_cont = true;
    int end_ind = n;

    clocktime_t start_time;
    if (time_spent) {
        start_time = tic();
    }

    for (int i=1; i < n + 2; i++)
    {   // relax edges repeatedly
        should_cont = false;

        for (int j=0; j < n; j++)
        {   // iterate over nodes
            if (dist_mat[i-1][j] != max_weight) 
            {   // avoid unvisited nodes
                if (dist_mat[i-1][j] < dist_mat[i][j]) 
                {   // update with results from previous relaxation
                    dist_mat[i][j] = dist_mat[i-1][j];
                }

                for (auto &node_iter : node_list[j])
                {   // iterate over connected nodes
                    if (dist_mat[i][node_iter.index] == max_weight && dist_mat[i-1][node_iter.index] != max_weight)
                    {
                        dist_mat[i][node_iter.index] = dist_mat[i-1][node_iter.index];
                    }

                    if (dist_mat[i-1][j] + node_iter.weight < dist_mat[i][node_iter.index]) 
                    {
                        dist_mat[i][node_iter.index] = dist_mat[i-1][j] + node_iter.weight;
                        path_mat[node_iter.index] = j;

                        should_cont = true;
                    }
                }
                //
            }
        }

        if (should_cont == false)
        {
            end_ind = i-1;
            break;
        }
    }

    // check for negative weight cycles

    for (int j=0; j < n; j++)
    {
        if (dist_mat[end_ind][j] != dist_mat[end_ind+1][j]) 
        {
            std::cout << "error: graph contains negative weight cycle!" << std::endl;
            return false;
        }
    }

    if (time_spent) {
        *time_spent = tic() - start_time;
    }

    //

    min_distance = std::move(dist_mat[end_ind]);
    path_list = std::move(path_mat);

    return true;
}
