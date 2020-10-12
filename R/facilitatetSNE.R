
load_FItSNE = function(path2fast_tsneR = NULL){
  source(path2fast_tsneR)
}

#' Wrapper function for FItSNE: fast_tsne.R 
#' 
#' @param path2fast_tsne a string specify the fast_tsne.R from FIt-SNE
#' @param data_path a string specify the data_path passed to FIt-SNE
#' @param fast_tsne_path a string specify the path of executable binary fast_tsne
#' @param load_affinities
#'            If 'precomputed', input data X will be regarded as precomputed similarities and passed to fast_tsne
#'            If 'load', input similarities will be loaded from the file.
#'            If 'save', input similarities are saved into a file.
#'            If 0 or NULL(default), affinities are neither saved nor loaded
#' @param ... include all the following fields that will be passed to fast_tsne
#' @param X data matrix 
#' @param dims dimensionality of the embedding. Default 2. 
#' @param perplexity perplexity is used to determine the bandwidth of the Gaussian kernel in the input space. Default 30.                                                   
#' @param theta Set to 0 for exact.  If non-zero, then will use either Barnes Hut or FIt-SNE based on nbody_algo. If Barnes Hut, then                                       
#'           this determins the accuracy of BH approximation. Default 0.5.                                                                                                  
#' @param max_iter Number of iterations of t-SNE to run. Default 1000.
#' @param ann_not_vptree use vp-trees (as in bhtsne) or approximate nearest neighbors (default).                                                                            
#'            set to be True for approximate nearest neighbors                                                                                                              
#' @param exaggeration_factor coefficient for early exaggeration (>1). Default = 1, not used.
#' @param no_momentum_during_exag Set to 0 to use momentum and other optimization tricks. 1 to do plain,vanilla                                                             
#'           gradient descent (useful for testing large exaggeration coefficients)                                                                                          
#' @param stop_early_exag_iter When to switch off early exaggeration. Default 250.                                                                                          
#' @param start_late_exag_iter When to start late exaggeration. 'auto' means that late exaggeration is not used, unless late_exag_coeff>0. In that                          
#'        case, start_late_exag_iter is set to stop_early_exag_iter. Otherwise, set to equal the iteration at which late exaggeration should begin. Default 'auto'.         
#' @param late_exag_coeff Late exaggeration coefficient. Set to -1 to not use late exaggeration. Default -1          
#' @param learning_rate Set to desired learning rate or 'auto', which sets learning rate to N/exaggeration_factor where N is the sample size, or to 200 if                  
#'       N/exaggeration_factor < 200. Default 'auto'                                                                                                                        
#' @param max_step_norm Maximum distance that a point is allowed to move on                                                                                                 
#'       one iteration. Larger steps are clipped to this value. This prevents                                                                                               
#'       possible instabilities during gradient descent.  Set to -1 to switch it                                                                                            
#'       off. Default: 5                                                                                                                                                    
#' @param mom_switch_iter Default 250
#' @param momentum Default 0.5 
#' @param final_momentum Default 0.8
#' @param nterms If using FIt-SNE, this is the number of interpolation points per sub-interval                                                                              
#' @param intervals_per_integer See min_num_intervals                                                                                                                       
#' @param min_num_intervals Let maxloc = ceil(max(max(X))) and minloc = floor(min(min(X))). i.e. the points are in a [minloc]^no_dims by [maxloc]^no_dims interval/square.  
#'           The number of intervals in each dimension is either min_num_intervals or ceil((maxloc - minloc)/intervals_per_integer), whichever is                           
#'           larger. min_num_intervals must be an integer >0, and intervals_per_integer must be >0. Default: min_num_intervals=50, intervals_per_integer = 1                
#' @param sigma Fixed sigma value to use when perplexity==-1 Default -1 (None)                                                                                              
#' @param K Number of nearest neighbours to get when using fixed sigma Default -30 (None)                                                                                   
#' @param initialization 'random', or N x no_dims array to intialize the solution. Default: 'random'.                                                                       
#' @param perplexity_list if perplexity==0 then perplexity combination will                                                                                                 
#'            be used with values taken from perplexity_list. Default: NULL                                                                                                 
#' @param df Degree of freedom of t-distribution, must be greater than 0.                                                                                                   
#'       Values smaller than 1 correspond to heavier tails, which can often                                                                                                 
#'       resolve substructure in the embedding. See Kobak et al. (2019) for                                                                                                 
#'       details. Default is 1.0                                                                                                                                            
#' @param fft_not_bh logical, default FALSE, we only use tSNE for interpreting 3D coordinates
#' @param verbose Print running infos for debugging.                                                                                                                        
#' @export
#' 
run_tSNE = function(path2fast_tsneR = NULL,
                    fast_tsne_path = NULL,
                    verbose = T, ...) {
  load_FItSNE(path2fast_tsneR)
  args = list(...)
  
  args$dims = ifelse(is.null(args$dims), 2, args$dims)
  args$perplexity = ifelse(is.null(args$perplexity), 30, args$perplexity)
  args$theta = ifelse(is.null(args$theta), 0.5, args$theta)
  args$max_iter = ifelse(is.null(args$max_iter), 1000, args$max_iter)
  args$fft_not_bh = ifelse(is.null(args$fft_not_bh), FALSE, args$fft_not_bh)
  args$ann_not_vptree = ifelse(is.null(args$ann_not_vptree), TRUE, args$ann_not_vptree)
  args$stop_early_exag_iter = ifelse(is.null(args$stop_early_exag_iter),
                                     250,
                                     args$stop_early_exag_iter)
  args$exaggeration_factor = ifelse(is.null(args$exaggeration_factor),
                                    1,
                                    args$exaggeration_factor)
  args$no_momentum_during_exag = ifelse(is.null(args$no_momentum_during_exag),
                                        FALSE,
                                        args$no_momentum_during_exag)
  args$start_late_exag_iter = ifelse(is.null(args$start_late_exag_iter),
                                     -1 ,
                                     args$start_late_exag_iter)
  args$late_exag_coeff = ifelse(is.null(args$late_exag_coeff), -1, args$late_exag_coeff)
  args$mom_switch_iter = ifelse(is.null(args$mom_switch_iter), 250 , args$mom_switch_iter)
  args$momentum = ifelse(is.null(args$momentum), 0.5 , args$momentum)
  args$final_momentum = ifelse(is.null(args$final_momentum), 0.8 , args$final_momentum)
  args$learning_rate = ifelse(is.null(args$learning_rate), 'auto', args$learning_rate)
  args$n_trees = ifelse(is.null(args$n_trees), 50 , args$n_trees)
  args$search_k = ifelse(is.null(args$search_k),-1 , args$search_k)
  args$rand_seed = ifelse(is.null(args$rand_seed),-1, args$rand_seed)
  args$nterms = ifelse(is.null(args$nterms), 3 , args$nterms)
  args$intervals_per_integer = ifelse(is.null(args$intervals_per_integer),
                                      1 ,
                                      args$intervals_per_integer)
  args$min_num_intervals = ifelse(is.null(args$min_num_intervals), 50 , args$min_num_intervals)
  args$K = ifelse(is.null(args$K),-1 , args$K)
  args$sigma = ifelse(is.null(args$sigma),-30 , args$sigma)
  args$initialization = ifelse(is.null(args$initialization), 'random', args$initialization)
  args$max_step_norm = ifelse(is.null(args$max_step_norm), 5, args$max_step_norm)
  # args$data_path = ifelse(is.null(args$data_path), NULL , args$data_path)
  # args$result_path = ifelse(is.null(args$result_path), NULL, args$result_path)
  # args$load_affinities = ifelse(is.null(args$load_affinities), NULL, args$load_affinities)
  # args$fast_tsne_path = ifelse(is.null(args$fast_tsne_path), NULL , args$fast_tsne_path)
  args$nthreads = ifelse(is.null(args$nthreads), 0 , args$nthreads)
  # args$perplexity_list = ifelse(is.null(args$perplexity_list), NULL , args$perplexity_list)
  args$get_costs = ifelse(is.null(args$get_costs), FALSE , args$get_costs)
  args$df = ifelse(is.null(args$df), 1.0, args$df)
  
  if (!is.null(args$load_affinities)) {
    # regarded input X as precomputed affinity matrix and output as the fast_tsne requested
    if (args$load_affinities == "precomputed") {
      if (is.null(args$X)) {
        stop("Empty data slot.")
      }
      if (args$theta == 0) {
        # only support exact tSNE
        if (verbose) {
          message(paste0("Saving affinity matrix to:", "P.dat"))
        }
        f <- file("P.dat", "wb")
        writeBin(as.numeric(args$X), f)
        close(f)
        args$load_affinities = "load" # rewrite for passing to fast_tsne
      } else {
       
        # row_P = ones([size(P,1)+1,1],'uint32');
        # col_P = ones([1,1],'uint32');
        # val_P = zeros([1,1], 'double');
        # k = 0;
        # for i = 1:size(P,1)
        #   row_P(i, 1) = k;
        #   for j = 1:size(P,1)
        #     if P(i,j) ~= 0 
        #       col_P(k+1, 1) = j;
        #       val_P(k+1, 1) = P(i, j);
        #       k = k+1;
        #     end
        #   end
        # end
        # row_P(size(P,1)+1,1) = k;
        
        if (verbose) {
          message(paste0("Saving affinity matrix to:", "val/row/col_P.dat"))
        }
        P = args$X
        
        tic()
        min_P = min(P) # min_P can be specified in the future, to be modified
        num_notZero = sum(P > min_P) # the minimum values were ignored and recognized as zero
        row_P = rep(0L, nrow(P) + 1)
        col_P = rep(0L, num_notZero)
        val_P = rep(0, num_notZero)
        
        k = 0
        for(i in 1:nrow(P)){
          row_P[i] = k
          for(j in 1:ncol(P)){
            if(P[i,j] != min_P){
              col_P[k] = j
              val_P[k] = P[i, j]
              k = k+1
            }
          }
        }
        row_P[i+1] = k
        toc()
        
        row_P = as.integer(row_P)
        col_P = as.integer(col_P)
        
        # write .dat
        f <- file("P_row.dat", "wb")
        writeBin(row_P, f)
        close(f)
        
        f <- file("P_col.dat", "wb")
        writeBin(col_P, f)
        close(f)
        
        f <- file("P_val.dat", "wb")
        writeBin(val_P, f)
        close(f)
        
        args$load_affinities = "load"
      }
    }
  }
  args$fast_tsne_path = fast_tsne_path
  
  out =  do.call(fftRtsne, args)
}


#' Run exact tsne, wrapper for integrated Exact TSNE calculation cpp
#' 
#' @param X affinity matrix to input
#' @param no_dims integer; Output dimensionality; Default = 3
#' @param verbose logical; Whether to print out debug information; Default = TRUE
#' @param max_iter integer; maximum iteration; Default = 1000
#' @param Y_in user-defined intiate coordinates; Default = NULL
#' @param init TRUE if Y_in were specified
#' @param rand_seed integer; random seed default = -1, set by time.
#' @param max_step_norm Maximum distance that a point is allowed to move on one iteration. Larger steps are clipped to this value. 
#' This prevents possible instabilities during gradient descent. Set to -1 to switch it off. (Default: 5) #' 
#' @param mom_switch_iter Numeric; (Default: 250)
#' @param momentum numeric; (Default 0.5)
#' @param final_momentum numeric; (Default 0.8)
#' @param df Degree of freedom of t-distribution, must be greater than 0. Values smaller than 1 correspond to heavier tails, 
#' which can often resolve substructure in the embedding. See Kobak et al. (2019) for details. Default is 1.0
#' @useDynLib CSOmapR
#' @import Rcpp
#' @export
#' 
runExactTSNE_R = function(X, no_dims = 3, ...){
  # NumericMatrix X, int no_dims, 
  # bool verbose, int max_iter, 
  # NumericMatrix Y_in, bool init, 
  # int rand_seed, bool skip_random_init, double max_step_norm, 
  # int mom_switch_iter, double momentum, double final_momentum, double df
  args = list(...)
  args$verbose = ifelse(is.null(args$verbose), T, args$verbose)
  args$max_iter = ifelse(is.null(args$max_iter), 1000, args$max_iter)
  args$rand_seed = ifelse(is.null(args$rand_seed), -1, args$rand_seed)
  args$init = ifelse(is.null(args$Y_in), F, T)
  if(!args$init){args$Y_in = matrix(0, 1, 1)} # default for input
  args$skip_random_init = args$init
  args$max_step_norm = ifelse(is.null(args$max_step_norm), 5, args$max_step_norm)
  args$mom_switch_iter = ifelse(is.null(args$mom_switch_iter), 250 , args$mom_switch_iter)
  args$momentum = ifelse(is.null(args$momentum), 0.5 , args$momentum)
  args$final_momentum = ifelse(is.null(args$final_momentum), 0.8 , args$final_momentum)
  args$df = ifelse(is.null(args$df), 1, args$df)
  if(args$verbose) message("Passing values to cpp...")
  out = do.call(runExactTSNE_wrapper, c(list(X = X, no_dims = no_dims), args))
  out$Y = t(out$Y)
  return(out)
}
