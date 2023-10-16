
// Copyright 2019 Lawrence Livermore National Security. 
// Produced at the Lawrence Livermore National Laboratory.
// LLNL-CODE-781679. All Rights reserved. See file LICENSE for details.
//
// This file is part of graph-embed. For more information and source code
// availability, see github.com/LLNL/graph-embed
//
// SPDX-License-Identifier: LGPL-2.1


#ifndef MATRIXUTILS_HPP
#define MATRIXUTILS_HPP

#include "sparsematrix.hpp"

using SparseMatrix = linalgcpp::SparseMatrix<double>;
using coord = std::vector<double>;
using coordinates = std::vector<coord>;

namespace partition {

    /* 
       @param n the desired dimension
       Returns the n*n identity matrix
       */
    SparseMatrix identity (int n);  

    /* 
       @param A the adjacency matrix  
       Given an adjacency matrix A, returns the corresponding graph Laplacian matrix L
       */
    SparseMatrix toLaplacian (const SparseMatrix& A);

    /* 
       @param L the graph Laplacian matrix
       Given a graph Laplacian matrix L, returns the corresponding adjacency matrix A
       */
    SparseMatrix fromLaplacian (const SparseMatrix& L); 
    double abs (double val);

    double distance (const std::vector<double>& v1, const std::vector<double>& v2);
      
    double magnitude (const std::vector<double>& v);
    void normalize(std::vector<double> & v);
     
    std::vector<double> addVectors(const std::vector<double>& vec1, const std::vector<double>& vec2);
    std::vector<double> multScalarVector(double scalar, const std::vector<double>& vec);
    

    double innerproduct_euclidean(const std::vector<double>& vec1, const std::vector<double>& vec2);
}

#endif // MATRIXUTILS_HPP
