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

#ifndef __GLOB_HPP_INCLUDED
#define __GLOB_HPP_INCLUDED

#include <cmath>     // std::abs, std::max, std::sqrt
#include <vector>    // std::vector
#include <iostream>  // std::cout, std::endl, std::sprintf
#include <limits>     // std::numeric_limits

typedef long int          mint;
typedef unsigned long int uint_t;
typedef double            mfloat;

typedef std::vector<mfloat> mfvec;
typedef std::vector<mint>   mivec;
typedef std::vector<uint_t>  muvec;

static const mfloat gEPS = std::sqrt(std::numeric_limits<mfloat>::epsilon());
static const mfloat gINF = std::numeric_limits<mfloat>::infinity();

// extern mfloat gEPS;
// extern mfloat gINF;
extern uint_t gVBS;

inline
bool
equal(const mfloat x, const mfloat y)
{
  return (std::abs (x-y) <
          std::max (gEPS, gEPS * std::max (std::abs (x), std::abs (y))));
}

// inline mfloat dist (const mfloat y1, const mfloat x1,
//                     const mfloat y2, const mfloat x2) {
//   return mfloat (std::sqrt ((y2 - y1)*(y2 - y1) + (x2 - x1)*(x2 - x1)));
// }

#endif // __GLOB_HPP_INCLUDED
