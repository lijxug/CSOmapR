#include <Rcpp.h>
#include <thread>
#include "tsne.h"
using namespace Rcpp;


// Function that runs as a wrapper for connection of R and tsne.cpp. Need to be modified
// [[Rcpp::export]]
Rcpp::List Rtsne_cpp(
  double X, int N, int D, double Y, int no_dims, double perplexity, double theta, int rand_seed,
  bool skip_random_init, int max_iter, int stop_lying_iter, int mom_switch_iter, 
  double momentum, double final_momentum, double learning_rate, int K, double sigma,
  int knn_algo, double early_exag_coeff, 
  bool no_momentum_during_exag, int start_late_exag_iter, double late_exag_coeff, int n_trees, int search_k,
  int nterms, double intervals_per_integer, int min_num_intervals, unsigned int nthreads, 
  int load_affinities, int perplexity_list_length, double perplexity_list, double df,
  double max_step_norm
) {
  
  // parse args
  bool no_momentum_during_exag_bool = true;
  if (no_momentum_during_exag == 0) no_momentum_during_exag_bool = false;
  
  
  // Run tsne
  TSNE* tsne = new TSNE();
  // Now fire up the SNE implementation
  double* costs = (double*) calloc(max_iter, sizeof(double));
  if(costs == NULL) { printf("Memory allocation failed!\n"); exit(1); }
  int error_code = 0;
  error_code = tsne->run(&X, N, D, &Y, no_dims, perplexity, theta, rand_seed, skip_random_init, max_iter, 
                         stop_lying_iter, mom_switch_iter, momentum, final_momentum, learning_rate, K, sigma, knn_algo, 
                         early_exag_coeff, costs, no_momentum_during_exag_bool, start_late_exag_iter, late_exag_coeff, n_trees,search_k, 
                         nterms, intervals_per_integer, min_num_intervals, nthreads, load_affinities,
                         perplexity_list_length, &perplexity_list, df, max_step_norm);
  
  if (error_code < 0) {
    exit(error_code);
  }
  
  // arrange outputs
  
  return Rcpp::List::create(Rcpp::_["Y"]=Y,
                            Rcpp::_["costs"]= *costs,
                            Rcpp::_["max_iter"]=max_iter,
                            Rcpp::_["N"]=N,
                            Rcpp::_["D"]=D);
}


