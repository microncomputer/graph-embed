// Copyright 2019 Lawrence Livermore National Security. 
// Produced at the Lawrence Livermore National Laboratory.
// LLNL-CODE-781679. All Rights reserved. See file LICENSE for details.
//
// This file is part of graph-embed. For more information and source code
// availability, see github.com/LLNL/graph-embed
//
// SPDX-License-Identifier: LGPL-2.1


#ifndef MICRON_EMBED_HPP
#define MICRON_EMBED_HPP

#include "partitioner.hpp"
#include <set>

namespace partition {

  /*
   * NEEDS COMMENTS
   */
    std::pair<double, double> functional_with_derivative(const SparseMatrix& A,
                                                        const int i, 
                                                        const std::vector<std::vector<double>>& coords, 
                                                        const std::vector<double>& direction, 
                                                        const std::set<int>& random_nonneighbors,
                                                        double epsilon,
                                                        double t );
 
    std::vector<double> random_direction (const int D);

    std::vector<std::vector<double>>embedViaMinimization_newprototype_micron(const SparseMatrix& A, 
                                                        const int d, 
                                                        std::vector<std::vector<double>>& new_coords, 
                                                        const int num_dirs_per_coord, const int ITER, 
                                                        const int size_S);

}
#endif // MICRON_EMBED_HPP
