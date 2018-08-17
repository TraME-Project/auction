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

// structures

struct node_data {
    int index;
    double weight;
    
    node_data(int arg_index, double arg_weight) : index(arg_index), weight(arg_weight) { }
};

// typedef convenience

const double max_weight = std::numeric_limits<double>::infinity();

typedef std::vector<std::vector<node_data>> graph_t;
