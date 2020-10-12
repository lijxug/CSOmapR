#include <Rcpp.h>
#include "mytsne.h"
#include <vector>
using namespace Rcpp;

// Function that runs the exact t-SNE
// [[Rcpp::export]]
Rcpp::List runExactTSNE_wrapper(
    NumericMatrix X, int no_dims,   
    bool verbose, int max_iter,  
    NumericMatrix Y_in, bool init,  
    int rand_seed, bool skip_random_init, double max_step_norm,  
    int mom_switch_iter, double momentum, double final_momentum, double df ) {

    if (verbose) Rprintf("Wrapper started\n");
    int N = X.ncol(), D = X.nrow();
    double * P =X.begin();
    
    // NumericVector costs_vec(max_iter);
    // double* costs = costs_vec.begin();
    
    if (verbose) Rprintf("Read the %i x %i data matrix successfully!\n", N, D);
    // std::vector<double> Y(N * no_dims), costs(N), itercosts(static_cast<int>(std::ceil(max_iter/50.0)));
    std::vector<double> Y(N * no_dims), costs(max_iter+1);
    
    // NumericMatrix Y_mt(N, no_dims);
    // double *Y = Y_mt.begin();
  
    // Providing user-supplied solution.
    if (init) {
        // for (int i = 0; i < Y_in.size(); i++) Y[i] = Y_in[i];
        double * Y = Y_in.begin();
        if (verbose) Rprintf("Using user supplied starting positions\n");
    }
    
    
    // Run tsne
    int exit_code = runExactTSNE(P, N, D, Y.data(), no_dims, rand_seed, skip_random_init, max_iter, 
                mom_switch_iter, momentum, final_momentum, costs.data(), df, max_step_norm, verbose);
    
    if(verbose){
      if(exit_code < 0){
        Rprintf("Error occured, exit_code = %i\n", exit_code);
      } else {
        Rprintf("Run successful! Returning values now.\n");
      }
    }

    return Rcpp::List::create(Rcpp::_["Y"]=Rcpp::NumericMatrix(no_dims, N, Y.data()), 
                              Rcpp::_["costs"]=Rcpp::NumericVector(costs.begin(), costs.end()),
                              Rcpp::_["N"]=N, 
                              Rcpp::_["D"]=D);
}

