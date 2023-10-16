#include <cstdlib>
#include <iostream>
#include <ctime>
#include <algorithm>

#include "parser.hpp"
#include "export.hpp"
#include "micron_embed.hpp"

#include <iomanip>
#include <fstream>

#include "partitioner.hpp"

using SparseMatrix = linalgcpp::SparseMatrix<double>;

/* TO DO:
 *      add functionality for inputting other formats (mtx, etc)
 * 
 */
int main (int argc, char* argv[]) {
    /* inputs:
     * 1. -f followed by filename of CooList graph to use.
     * 2. (optional) -d followed by integer dimension to embed to. if none given, it is 3.
     * 3. (optional) -c followed by coarsening factor. if none given, it is 0.1 .
     *
     * note: assumes symmetric. if needed, this can be easily changed.
     */


// Processing Inputs from Program Execution
    std::string inputpath;
    std::string outputpath;
    int dimension = 3;
    int numdirs = 1;
    int size_S = 0; //size of random subset of neighbors for the minimization scheme. if none, does standard process (all neighbors)
    double coarseningFactor = 1.0 / 10.0; 


    for (int i=0; i<argc-1; i++) {
        if (std::string(argv[i]) == "-f") {
            inputpath = std::string(argv[i+1]);
        } else if (std::string(argv[i]) == "-o") {
            outputpath = std::string(argv[i+1]);
        }else if (std::string(argv[i]) == "-d") {
            dimension = std::stoi(std::string(argv[i+1]));
        } else if (std::string(argv[i]) == "-c") {
            coarseningFactor = std::stoi(std::string(argv[i+1]));
        } else if (std::string(argv[i]) == "-nd") {
            numdirs = std::stoi(std::string(argv[i+1]));
        } else if (std::string(argv[i]) == "-ss") {
            size_S = std::stoi(std::string(argv[i+1]));
        }
    }

    if (inputpath == "" or outputpath == "") {
        std::cerr << "the required inputs are as follows:" << std::endl;
        std::cerr << "-f graph_input_filename -o outputpath_for_write_coords" << std::endl;
        std::cerr << "follow that with optional inputs:" << std::endl;
        std::cerr << "-d dimension_to_embed_into -c coarsening_factor -nd number_of_directions -ss size_random_subset_neighbors" << std::endl;
        return 1;
    }

    SparseMatrix A = linalgcpp::ReadCooList(inputpath, true);
    int n = A.Rows();


    // Partitioning 
    std::vector<SparseMatrix> hierarchy = partition::partition(A, coarseningFactor, false, true, 1.0, 1, false);
    
    //Embedding via minimization of new functional
    std::vector<std::vector<double>> coords (A.Rows(), std::vector<double>(dimension, 0.0));
    coords = partition::embedViaMinimization_newprototype_micron (A, dimension, coords, numdirs, 10, 10);


    //std::vector<std::vector<double>> coords = partition::embed(As, hierarchy, dimension);

    partition::writeCoords(coords, outputpath);

}
    
