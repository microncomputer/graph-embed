
// Copyright 2019 Lawrence Livermore National Security. 
// Produced at the Lawrence Livermore National Laboratory.
// LLNL-CODE-781679. All Rights reserved. See file LICENSE for details.
//
// This file is part of graph-embed. For more information and source code
// availability, see github.com/LLNL/graph-embed
//
// SPDX-License-Identifier: LGPL-2.1


#include "matrixutils.hpp"

namespace partition {

    SparseMatrix identity(int n){
        /*std::vector<int> nI (n+1);
          nI[0] = 0;
          std::vector<int> nJ (n);
          std::vector<double> nD (n, 1.0);

          for (int i=0; i<n; i++) {
          nI[i+1] = i+1;
          nJ[i] = i;
          }

          return SparseMatrix(nI, nJ, nD, n, n);*/
        return SparseMatrix(std::vector<double>(n, 1.0));
    }

    SparseMatrix toLaplacian(const SparseMatrix& A) {
        int numRows = A.Rows();
        const std::vector<int>& I = A.GetIndptr();
        const std::vector<int>& J = A.GetIndices();
        const std::vector<double>& D = A.GetData();

        std::vector<int> nI (numRows+1);
        for (int i=1; i<numRows+1; i++) 
        {nI[i] = I[i]+i;}
        std::vector<int> nJ(A.nnz()+numRows);
        std::vector<double> nD(A.nnz()+numRows);

        int count=0;
        for (int i=0; i<numRows; i++) {
            bool passed=false;
            for (int k=I[i]; k<I[i+1]; k++) {
                if (!passed && J[k]>i){
                    double sum=0; 
                    for (int k0=I[i]; k0<I[i+1]; k0++) 
                    {sum+=D[k0];}
                    nJ[k+count] = i; 
                    nD[k+count] = sum;
                    count++;
                    passed=true;
                }
                nJ[k+count] = J[k]; 
                nD[k+count] = -D[k];
            }
            if (!passed) {
                double sum=0; 
                for (int k0=I[i]; k0<I[i+1]; k0++) 
                {sum+=D[k0];}

                nJ[I[i+1]+count] = i; 
                nD[I[i+1]+count] = sum;
                count++;
            }
        }
        return SparseMatrix(nI, nJ, nD, numRows, A.Cols());
    }

    SparseMatrix fromLaplacian(const SparseMatrix& L) {
        int numRows = L.Rows();
        const std::vector<int>& I = L.GetIndptr();
        const std::vector<int>& J = L.GetIndices();
        const std::vector<double>& D = L.GetData();

        std::vector<int> nI (numRows+1);
        std::vector<int> nJ (L.nnz() - numRows);
        std::vector<double> nD (L.nnz() - numRows);
        for (int i=1; i<numRows+1; i++) 
        {nI[i] = I[i]-i;}

        int count=0;
        for (int i=0;i<numRows; i++) {
            bool passed=false;
            for(int k=nI[i]; k<nI[i+1]; k++) {
                if (!passed && J[k+count]>=i) {
                    count++;
                    passed=true;
                }
                nJ[k] = J[k+count]; nD[k] = -D[k+count];
            }
            if (!passed)
            {count++;}
        }
        return SparseMatrix(nI, nJ, nD, numRows, L.Cols());
    }
    double abs (double val) {
        return  (val < 0) ? -val : val;
    }

    double distance (const std::vector<double>& v1, const std::vector<double>& v2) {
        assert(v1.size() == v2.size());
        double sum = 0.0;
        for (int i=0; i<v1.size(); i++) {
            double d = v2[i] - v1[i];
            sum += d * d;
        }
        return sqrt(sum);
    }

    double magnitude (const std::vector<double>& v) {
        double sum = 0.0;
        for (int i=0; i<v.size(); i++) {
            double d = v[i];
            sum += d * d;
        }
        return sqrt(sum);
    }

    void normalize(std::vector<double> & v){
        double m = magnitude(v);
        for(int i =0; i< v.size(); i++){
            v[i] /= m;
        }
        return;
    }

    std::vector<double> addVectors(const std::vector<double>& vec1, const std::vector<double>& vec2)  
    {  
        // Check if the vectors have the same size  
        if (vec1.size() != vec2.size())  
        {  
            std::cout << "Error: Vectors must have the same size to add them." << std::endl;  
            return std::vector<double>();  
        }  

        // Create a vector to store the result  
        std::vector<double> result(vec1.size());  

        // Add the elements of vec1 and vec2 and store the result in result  
        for (int i = 0; i < vec1.size(); i++)  
        {  
            result[i] = vec1[i] + vec2[i];  
        }  

        return result;  
    }  

    std::vector<double> multScalarVector(double scalar, const std::vector<double>& vec)  
    {                                         
        // Create a vector to store the result  
        std::vector<double> result(vec.size());  

        for (int i = 0; i < vec.size(); i++)  
        {  
            result[i] = scalar * vec[i];  
        }  

        return result;  
    }

    double innerproduct_euclidean(const std::vector<double>& vec1, const std::vector<double>& vec2){
        // Check if the vectors have the same size  
        if (vec1.size() != vec2.size())  
        {  
            std::cout << "Error: Vectors must have the same size to do inner product." << std::endl;  
            return 0;
        }  

        double result = 0;  

        // multiply the elements of vec1 and vec2 and store the result in result  
        for (int i = 0; i < vec1.size(); i++)  
        {  
            result += vec1[i] * vec2[i];  
        }  

        return result;  
    }

}
