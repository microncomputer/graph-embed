#include <vector>
#include <numeric>
#include <set>
#include <algorithm>
#include "micron_embed.hpp"
//#include "forceatlas.hpp"


/* prototyping, notes and questions:
 * at each level of the hierarchy during coarsening? or do we wait till we get the coarsened graph and then go through and embed each level?
 * would it be useful to get METIS to do the partitioning?
 * 
 * different approaches (do i need to implement them all?)
 *  -one vertex at a time(with gauss seidel to minimize)  or random small batch of vertices (with stochastic gradient to minimize)
 *  -
 *  
 *  Q: if no random neighbors set, is the base case that we process each neighbor in that second summation(is this what is meant by sum over j in phi_of_t on p.2 section 2.1? or all nodes not in the neighbors as seen on p.2 1.1 in PSV?
 * ok wow, so i was doing random set of neighbors but its actually random subset of non neighbors
 */


namespace partition {

    //evaluate the functional at t for the i'th node, return result and derivative in given direction
    std::pair<double, double> functional_with_derivative(const SparseMatrix& A,
            const int i, 
            const std::vector<std::vector<double>>& coords, 
            const std::vector<double>& direction, 
            const std::set<int>& random_nonneighbors,
            double epsilon,
            double t ){ //ben thinks t should come in but I think the whole range of t's should be processed here on each iteration
        /* 
         * @param A, graph in CSR(compressed sparse row) format
         * @param i, index of node we are processing
         *  
         *
         */

        int n = A.Rows();
        const std::vector<int>& INDPTR = A.GetIndptr(); //used to be I
        const std::vector<int>& INDICES = A.GetIndices(); //J
        const std::vector<double>& DATA = A.GetData(); //D


        // here's an outline of what's to happen:
        // f_summation_1: Sum over the neighbors of i: norm of (sum of coord x_i, t*g, and negative coord x_j)
        // f_summation_2: epsilon *(1/norm(coord x_i + t*direction) + 1/(1-norm(coord x_i + t*direction)) 
        //              + sum over the nodes in the random subset of nonneighbors (1/norm(coord x_i + t*direction - coord x_j)))
        // entire functional(t): f_summation_1 + epsilon(f_oneterm + f_twoterm + f_summation_2)
        // these pieces used in function and differential
        std::vector<double> xi = coords[i];
        std::vector<double> tg = multScalarVector(t, direction);
        std::vector<double> add_xitg = addVectors(xi, tg);
        double norm_add_xitg = magnitude(add_xitg);

        double f_summation_1 = 0.0;
        double f_summation_2 = 0.0; 
        double f_oneterm = 1/ norm_add_xitg;
        double f_twoterm = 1/(1 - norm_add_xitg);

        //these pieces used in differential
        auto xi_g = innerproduct_euclidean(xi, direction);
        auto norm_g_squared = innerproduct_euclidean(direction, direction);

        double df_summation_1 = 0.0;
        double df_oneterm = -(xi_g + t * norm_g_squared) / pow(norm_add_xitg,3);
        double df_twoterm = (xi_g + t * norm_g_squared)/(norm_add_xitg * pow((1-norm_add_xitg),2));
        double df_summation_2 = 0.0; 

        // REMINDER: indices[indptr[i]:indptr[i+1]] has column indices of row i and data[indptr[i]:indptr[i+1]] has the data values for those elements
        for (int col_index=INDPTR[i]; col_index<INDPTR[i+1]; col_index++) {   
            int neighbor = INDICES[col_index];
            auto xj = coords[neighbor];
            if (i != neighbor) {
                double normsumofterms = magnitude(addVectors(addVectors(xi , multScalarVector(-1, xj)), tg));
                f_summation_1 += normsumofterms;
                df_summation_1 += (innerproduct_euclidean(addVectors(xi, multScalarVector(-1, xj)), direction) + t * norm_g_squared) / magnitude(addVectors(addVectors(xi , multScalarVector(-1, xj)), tg));
            }
        }

        //processing nonneighbors for second summations in both f and df
        //if nonneighbor in the random subset, it is considered in summation_2 
        for (const auto& nonneighbor : random_nonneighbors) {
            auto xj = coords[nonneighbor];
            double normsumofterms = magnitude(addVectors(addVectors(xi , multScalarVector(-1, xj)), tg));
            f_summation_2 += 1/normsumofterms;
            df_summation_2 += (innerproduct_euclidean(addVectors(xi, multScalarVector(-1, xj)), direction) 
                    + t * norm_g_squared) 
                / pow(magnitude(addVectors(addVectors(xi , multScalarVector(-1, xj)), tg)),3);
        }


        //add terms together:
        double phi_of_t = f_summation_1 + epsilon * (f_oneterm + f_twoterm + f_summation_2);
        double phi_of_t_prime = df_summation_1 + epsilon * (df_oneterm + df_twoterm - df_summation_2);
        auto f_df = std::make_pair(phi_of_t, phi_of_t_prime);
        return f_df;
    }

    std::vector<double> random_direction (const int D){
        /* D is the dimension to embed into
         * return g, random direction
         *
         * add input of random gen instead so that the seed is reproducible
         */

        std::random_device rd;
        std::mt19937 gen (rd());
        std::uniform_real_distribution<double> random (-1.0, 1.0);
        auto g = std::vector<double>(D);
        for(int i=0; i<D; i++){
            g[i] = (random(gen));
        }
        normalize(g); //normalize g

        return g;
    }

#if 0
    std::vector<std::vector<double>> embedViaMinimization_newprototype_micron (const SparseMatrix& A, const std::vector<SparseMatrix>& P_Ts,  const int d, const int size_S, const int subintervals_t) {
        /* wrapper function for the embedViaMinimization function that follows this one.
         * Here is where we initialize an empty vector of coordinates to hold the right amount of nodes,
         * each node coordinate being its own vector of size d, the dimension to embed into.
         *  
         * Also, we process the hierarchy and pass each level to the inner function.
         * 
         *
         * Inputs:
         *  SparseMatrix A: graph
         *  vector of 'SparseMatri'ces P_Ts: the Hierarchy of partition matrices in aggregate_vertex form 
         *                                  (NOTE: this is the transpose of what you may expect, hence P_Ts)
         *  int d: dimension to embed into
         *  int size_S: size of the random set of indices to be evaluated over on each level of the hierarchy
         *  int subintervals_t: for t in [a,b], partition [a,b] into subintervals_t sub-intervals. eg, 10,20,50
         *   
         *   new direction and S every iteration
         *  try with neighbors and without
         *  TODO: if size_s is size of all the neihbors of xi or if it is 0, then don't make random subset and run on all neighbors
         *
         *
         * 
         * email from ben:
         So yeah, the goal is to find the t that minimizes the functional. The derivative tells you what direction to move (positive -> decrease t, negative -> increase t). You can probably find such a t via binary search: evaluate the function at t=1/2 and take steps of +/- 1/4, 1/8, etc each iteration in the direction the derivative instructs. But as you noticed, you need to evaluate the main term to to keep track of the minimum value. Evaluating both at the same time is definitely fine.


*/


        /* not using yet
           assert(As.size() == P_Ts.size() + 1);
           for (int i=0; i<P_Ts.size(); i++) {
           assert(As[i].Rows() == P_Ts[i].Cols());
           }
           for (int i=0; i<P_Ts.size(); i++) {
           assert(As[i+1].Rows() == P_Ts[i].Rows());
           }
           std::vector<double> none (0);:
           std::vector<std::vector<double>> none2 (0);
           return embedMultilevel(As, P_Ts, d, 0, none, none2);

           std::vector<std::vector<double>> coords (A.Rows(), std::vector<double>(d, 0.0));
           embedViaMinimization(A, d, coords, 1000);
           return coords;
           */
        return;
    }
#endif

    std::vector<std::vector<double>> embedViaMinimization_newprototype_micron (const SparseMatrix& A, const int d, std::vector<std::vector<double>>& new_coords, const int num_dirs_per_coord, const int ITER, const int size_S ) {
        /* function in prototyping stage. 
         * TODO: turn this into return void, alter coords implicitly and pass the coords from wrapper
         *
         * assume A is one level of hierarchy (for now)
         * 
         * subintervals_t not being used. instead t is updated with a step size that gets smaller. could be way more efficient..
         *
         */

        double inf = std::numeric_limits<double>::infinity();
        double epsilon = 10e-3;

        std::random_device rd;
        std::mt19937 gen (rd());
        std::uniform_real_distribution<double> random (-1.0, 1.0);

        int n = A.Rows();

        const std::vector<int>& I = A.GetIndptr();
        const std::vector<int>& J = A.GetIndices();
        const std::vector<double>& D = A.GetData();

        if (new_coords.size() == 0) {
            new_coords = std::vector<std::vector<double>> (n, std::vector<double> (d));
        }

        // initialize x_i to random values
        for (int i=0; i<n; i++) {
                new_coords[i] = random_direction(d);
        }


        //iterate for a number of directions for each coord 
        //on each direction take xi some amount that way.
        // iterate to find approximate solution for all nodes
        for (int iter=0; iter<ITER; iter++) 
        {
            for (int i=0; i<n; i++) 
            {
                const std::vector<double>& x_i = new_coords[i];
                int neighborcount = A.RowSize(i);
                //only do if has neighbors TODO: confirm that this doesn't break something
                if (neighborcount == 0)
                    break;

                std::uniform_int_distribution<int> uni(0,n-1);
                int nonneighborcount = n-neighborcount; 

                int size_S_confirmed;
                //initialize the random set of nonneighbors
                if (nonneighborcount < size_S)
                    size_S_confirmed = nonneighborcount;
                else
                    size_S_confirmed = size_S;
                auto S = std::set<int>(); 
                const std::vector<int> neighbors = A.GetIndices(i);
                while (S.size() < size_S_confirmed){
                    auto s = uni(gen);
                    // add to set if not i or any neighbors of i. could this be optimized? 
                    // yes, i think since the indices are in order, a binary search would work. maybe later
                    if (s != i){
                        if(std::find(neighbors.begin(), neighbors.end(), s) == neighbors.end()) {
                            /*s is not a neighbor, can be added to set */
                            S.insert(s);
                        }
                    }
                }

                for(int dirs = 0; dirs < num_dirs_per_coord; dirs++){
                    //initialize random normalized direction
                    const std::vector<double> g = random_direction(d);

                    // the constraint on the functional is that norm(x_i + t*g) < 1
                    // This implies the closed interval bounds for t are -<x_i,g> +- sqrt(<x_i,g>2 + 1 - norm(x_i)^2)
                    // ie. t is in [a,b] taken at subintervals_t equally spaces values within [a,b]
                    double xdotg = innerproduct_euclidean(x_i, g);
                    double a = -xdotg - sqrt(pow(xdotg,2) + 1 - pow(magnitude(x_i),2)); 
                    double b = -xdotg + sqrt(pow(xdotg,2) + 1 - pow(magnitude(x_i),2)); 
                    //double oneportion_t = (b-a)/subintervals_t;
                    //starting t in middle of interval
                    double where_in_interval = (b-a)/2;
                    double curr_t ; 
                    double t_for_min_f;
                    double min_f = inf;



                    /*
                     * find minimum t for the direction and then take xi to be xi+tg
                     */
                    do{
                        //for(double t = a; t < b; t+= oneportion_t)
                        //linear interpolation to get current t value
                        curr_t = std::lerp(a, b, where_in_interval);

                        //as defined above, this returns a pair with the functional at t and the differential at t
                        auto f_df_at_t = functional_with_derivative(A, i, new_coords, g, S, epsilon, curr_t);

                        double f = std::get<0>(f_df_at_t);
                        double df = std::get<1>(f_df_at_t);

                        if(f < min_f)
                        {
                            min_f = f;
                            t_for_min_f = curr_t;

                            //check if df +,-,0
                            //Q: if approaching a saddlepoint, do I take something into account?
                            if (df < 0.0) //take t up a step
                                where_in_interval += where_in_interval/2;
                            else if (df > 0.0) //take t down a step
                                where_in_interval -= where_in_interval/2;
                            else //df == 0 or really close
                                continue; //this could be wrong

                        }
                        //else f is larger than the last f, safe to say that was min
                        else break;

                    }while(where_in_interval > 1e-4);

                    //found f minimum, update coords with respect to the t value for f_min
                    //if(curr_t)
                    std::vector<double> temp = addVectors(new_coords[i], multScalarVector(curr_t, g));
                    //after removing bugs in iteration process, assert the following norm
                    if (magnitude(temp) > 1)
                        normalize(temp);
                    new_coords[i] = temp;
                }
            }
        }
        return new_coords;

        // center the coords at origin
        if (n > 1) {
            std::vector<double> avg (d);
            for (int i=1; i<n; i++) {
                for (int k=0; k<d; k++) {
                    avg[k] = avg[k] + new_coords[i][k];
                }
            }
            for (int k=0; k<d; k++) {
                avg[k] = avg[k] / (n - 1);
            }
            for (int i=0; i<n; i++) {
                for (int k=0; k<d; k++) {
                    new_coords[i][k] -= avg[k];
                }
            }

            double max_length = 0.0;
            for (int i=1; i<n; i++) {
                double length = magnitude(new_coords[i]);
                if (max_length < length) {
                    max_length = length;
                }
            }
            for (int i=0; i<n; i++) {
                for (int k=0; k<d; k++) {
                    new_coords[i][k] = new_coords[i][k] / max_length;
                }
            }
        }
    return new_coords;
    }
}
/* following is Ben's code copied from embedViaMinimization in embed.cpp
 * using as needed
 if (count != 0) {
 double min_J_loc = inf;
 double min_t = 0.0f;
 double min_s = -1;
 double w = 1000000.0;
 for (int s=0; s<dirs.size(); s++) {
// find the minimum for the interval (x_i, x_s)
const std::vector<double>& x_s = dirs[s];

double t = 0.5;
double jump = 0.25;

do {
double dJ_loc_dt = 0.0;
for (int r=0; r<n; r++) {
if (i != r) {
const auto& x_r = new_coords[r];
double term1 = 0.0;
double term2 = 0.0;
for (int k=0; k<d; k++) {
double u_k = x_s[k] - x_i[k];
double v_k = x_i[k] - x_r[k];
double z_k = (u_k * t + v_k);
term1 = term1 + z_k * z_k;
term2 = term2 + z_k * u_k;
}
if (term1 < epsilon) {
term1 = epsilon;
}
dJ_loc_dt += - ((1.0 / sqrt(term1 * term1 * term1)) * term2);
}
}

for (int kk=I[i]; kk<I[i+1]; kk++) {
int r = J[kk];
if (i != r) {
const auto& x_r = new_coords[r];
double term = 0.0;
for (int k=0; k<d; k++) {
double a = (1 - t) * x_i[k] + t * x_s[k] - x_r[k];
term += w * 2.0 * a * (x_s[k] - x_i[k]);
}
dJ_loc_dt += term;
}
}

if (dJ_loc_dt < 0.0) {
t = t + jump;
} else {
t = t - jump;
}
jump = jump / 2.0;

} while (jump > 1.e-4); // now the step size is small enough to move on

double J_loc = 0.0;
for (int r=0; r<n; r++) {
if (i != r) {
const auto& x_r = new_coords[r];
double term1 = 0.0;
for (int k=0; k<d; k++) {
double u_k = x_s[k] - x_i[k];
double v_k = x_i[k] - x_r[k];
double z_k = (u_k * t + v_k);
term1 = term1 + z_k * z_k;
}
if (term1 < epsilon) {
term1 = epsilon;
}
J_loc = J_loc + 1.0 / sqrt(term1);
}
}
for (int kk=I[i]; kk<I[i+1]; kk++) {
    int r = J[kk];
    if (i != r) {
        const auto& x_r = new_coords[r];
        double term = 0.0;
        for (int k=0; k<d; k++) {
            double a = (1 - t) * x_i[k] + t * x_s[k] - x_r[k];
            term += a * a;
        }
        J_loc += w * term;
    }
}
if (J_loc < min_J_loc) {
    min_J_loc = J_loc;
    min_t = t;
    min_s = s;
}

//int num = 100;
//double dt = 1.0 / num;
//for (int a=1; a<num; a++) {
//	double t = a * dt;
//	double J_loc = 0.0;
//	for (int r=0; r<s_num; r++) {
//	  auto& x_r = x_j[r];
//	  double bottom = 0.0;
//	  for (int k=0; k<d; k++) {
//	double z_k = (1 - t) * x_i[k] + t * x_s[k] - x_r[k];
//		bottom = bottom + z_k * z_k;
//	  }
//	  J_loc = J_loc + 1.0 / sqrt(bottom);
//	}
//	if (J_loc > min_J_loc) {
//	  min_J_loc = J_loc;
//	  min_t = t;
//	  min_s = s;
//	}
//}

}
//std::cout << "i: " << i << " | t: " << min_t << " | s: " << min_s << std::endl;
if (min_s >= 0) {
    for (int k=0; k<d; k++) {
        new_coords[i][k] = new_coords[i][k] * (1 - min_t) + dirs[min_s][k] * min_t;
    }
}
}
}
}

if (n > 1) {
    // center coords at 0
    std::vector<double> avg (d);
    for (int i=1; i<n; i++) {
        for (int k=0; k<d; k++) {
            avg[k] = avg[k] + new_coords[i][k];
        }
    }
    for (int k=0; k<d; k++) {
        avg[k] = avg[k] / (n - 1);
    }
    for (int i=0; i<n; i++) {
        for (int k=0; k<d; k++) {
            new_coords[i][k] -= avg[k];
        }
    }

    double max_length = 0.0;
    for (int i=1; i<n; i++) {
        double length = magnitude(new_coords[i]);
        if (max_length < length) {
            max_length = length;
        }
    }
    for (int i=0; i<n; i++) {
        for (int k=0; k<d; k++) {
            new_coords[i][k] = new_coords[i][k] / max_length;
        }
    }
    */

