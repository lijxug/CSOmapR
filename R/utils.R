# suppressPackageStartupMessages(library("Rcpp"))
# suppressPackageStartupMessages(library("RcppEigen"))
# sourceCpp("utils/functions.cpp") # hard coded, should be modified in the future

# Main interface ----

#' Test for significance level
#' 
#' @param coordinate a 3D matrix
#' @param labels a vector of celltype labels, correspond to the coordinates matrix
#' @param k a integer for top k connections
#' @export
#' @return A list contains coordinates, counts, p values and q values
#' 
getSignificance = function(coordinates, labels, k = 3, adjusted.method = "fdr", verbose = F) {
  stopifnot(is.matrix(coordinates))
  if(ncol(coordinates) < 2 | ncol(coordinates) > 3)
    warning("Abnormal number of dimensions of the coordinates")
  # preprocess
  labels = setNames(labels, nm = rownames(coordinates))
  standards <- unique(labels)
  # labelIx <- match(labels, standards)
  cellCounts <- table(labels)
  
  # Calc dist_mtance
  dist_mt <- as.matrix(dist(coordinates))
  
  # identify topK as a cutoff
  if(verbose) loginfo("identify topK")
  topKs <- c()
  diag(dist_mt) <- Inf
  topKs = apply(dist_mt, 1, function(dist_mt_row_i){
    dist_mtSorted <- sort(dist_mt_row_i)
    dist_mtSorted[k]
  })
  topK <- median(topKs)
  
  # initiate counts
  counts <-
    matrix(0, nrow = length(standards), ncol = length(standards))
  colnames(counts) <- standards
  rownames(counts) <- standards
  
  # Cells within topK range are recognized as connected
  
  # if(verbose) loginfo("calculate connection")
  # for (i in 1:nrow(dist_mt)) { # explore cells one by one
  #   connects <- which(dist_mt[i,] <= topK)
  #   for (j in connects) {
  #     counts[labelIx[i], labelIx[j]] = counts[labelIx[i], labelIx[j]] + 1
  #   }
  # }
  
  # an alternative way
  if(verbose) loginfo("calculate detailed connections")
  connects_mt = dist_mt <= topK
  counts =
    apply(apply(connects_mt, 1, function(x) {
      tapply(x, labels, sum)
    }), 1, function(x) {
      tapply(x, labels, sum)
    })
  
  diag(counts) <- diag(counts) / 2 # diag elements were counted twice
  
  # Store detailed connected cell pairs
  detailed_connections = list()
  clusterPairs2run = rbind(t(combn(standards, 2)), matrix(rep(standards, 2), ncol = 2))
  for(i_row in 1:nrow(clusterPairs2run)) {
    cluster1 = clusterPairs2run[i_row, 1]
    cluster2 = clusterPairs2run[i_row, 2]
    
    cellsfrom1 = names(labels)[labels == cluster1]
    cellsfrom2 = names(labels)[labels == cluster2]
    
    sub_connects_mt = connects_mt[cellsfrom1, cellsfrom2, drop = F]
    
    if(cluster1 == cluster2){ # only count once if self2self
      sub_connects_mt[lower.tri(sub_connects_mt)] = F
    }
    sub_connects_df = data.frame(
      cell1 = cellsfrom1[row(sub_connects_mt)][sub_connects_mt], 
      cell2 = cellsfrom2[col(sub_connects_mt)][sub_connects_mt], stringsAsFactors = F
    )
    detailed_connections[[paste0(cluster1, "---", cluster2)]] = sub_connects_df
  }
  
  # calculate pvalue using hypergeometric distribution
  if(verbose) loginfo("calculate pvalues")
  K <- (sum(counts) + sum(diag(counts))) / 2
  p_value <- counts
  
  assertthat::assert_that(all(rownames(counts) == names(cellCounts)))
  assertthat::assert_that(all(colnames(counts) == names(cellCounts)))
  for (i in 1:nrow(counts)) {
    for (j in 1:ncol(counts)) {
      if (i == j) {
        M <- as.numeric(cellCounts[rownames(counts)[i]]) * (as.numeric(cellCounts[colnames(counts)[j]]) - 1) / 2
      } else {
        M <- as.numeric(cellCounts[rownames(counts)[i]]) * (as.numeric(cellCounts[colnames(counts)[j]]))
      }
      N <- sum(cellCounts) * (sum(cellCounts) - 1) / 2 - M
      p_value[i, j] <-
        phyper(counts[i, j], M, N, K, lower.tail = FALSE)
    }
  }
  
  # p adjust
  clusters = colnames(p_value)
  cluster_pair = paste(clusters[row(p_value)], clusters[col(p_value)], sep="---")
  low_idx = lower.tri(p_value, diag = T)
  p_value_df = 
    data.frame(cluster_pair = cluster_pair[low_idx], 
               p.value = p_value[low_idx], 
               q.value = p.adjust(p_value[low_idx], method = adjusted.method), 
               stringsAsFactors = F)
  
  q_value = p_value
  q_value[low_idx] = p_value_df$q.value
  q_value = t(q_value)
  q_value[low_idx] = p_value_df$q.value
  
  result = list()
  result$connections = counts
  result$pvalue = p_value
  result$qvalue = q_value
  result$pvalue_tbl = p_value_df
  result$detailed_connections = detailed_connections
  result$topK = topK
  return(result)
}


#' Optimize the 3D coordinates(cpp)
#' This function is inspired from tsne algorithm. We use similar gradient descent 
#' method to optimize our target function specified in our paper. 
#' condition can be loose or tight, we suggest using "loose" condition 
#' for dataset with over 10000 cells 
#' @param affinityMat affinity matrix
#' @param initial_config initial configuration
#' @param k k cells 
#' @param max_iter Maximum iteration time
#' @param min_cost Minimum cost
#' @param condition A string, either 'loss' or 'tight'
#' @param momentum initial momentum, default = 0.5
#' @param final_momentum final momentum, default = 0.8
#' @param mom_switch_iter value to which momentum is changed, default = 250
#' @param epsilon initial learning rate, default = 1000
#' @param min_gain minimum gain for delta-bar-delta, default = 0.01
#' @param eps Minimum distances between cells
#' @param epoch numeric, print out lost funciton cost after every *epoch* iterations
#' @param verbose logical. If TRUE, print out the progress information
#' @return a matrix of optimized 3D coordinates
#' @export
#' 
optimization <-
  function (affinityMat,
            initial_config = NULL,
            k = 3,
            max_iter = 1000,
            min_cost = 0,
            condition = "tight",
            momentum = 0.5, 
            final_momentum = 0.8, 
            mom_switch_iter = 250, 
            epsilon = 1000, 
            min_gain = 0.01,
            eps = 2.2251e-308,
            epoch = 100,
            verbose = F) {
    n = nrow(affinityMat)
    
    if (!is.null(initial_config) && is.matrix(initial_config)) {
      if (nrow(initial_config) != n | ncol(initial_config) !=
          k) {
        stop("initial_config argument does not match necessary configuration for X")
      }
      ydata = initial_config
    }
    else {
      # ydata = matrix(rnorm(k * n), n)
      ydata = (matrix(runif(k * n), n) - 0.5) * 50
    }
    P = affinityMat
    # P = 0.5 * (affinityMat + t(affinityMat))
    # P[P < eps] <- eps
    # P = P / sum(P)
    grads = matrix(0, nrow(ydata), ncol(ydata))
    incs = matrix(0, nrow(ydata), ncol(ydata))
    gains = matrix(1, nrow(ydata), ncol(ydata))
    if(verbose) loginfo("Iteration started")
    iter01 = ".I1_fun"
    initiatePB(iter01)
    
    for (iter in 1:max_iter) {
      d = calc_d_rcpp(ydata)
      num = 1/(1+d)
      diag(num) = 0
      Q = num/sum(num)
      Q[Q < eps] = eps
      P_Q = P - Q
      P_Q[P_Q > 0 & d <= 0.01] = -0.01
      
      # stiffnesses = 4 * P_Q * num
      stiffnesses = 4 * (P-Q) * num
      
      grads = update_grads_rcpp(grads, ydata, stiffnesses)
      
      gains = ((gains + 0.2) * abs(sign(grads) != sign(incs)) +
                 gains * 0.8 * abs(sign(grads) == sign(incs)))
      gains[gains < min_gain] = min_gain
      incs = momentum * incs - epsilon * (gains * grads)
      ydata = ydata + incs
      ydata = sweep(ydata, 2, apply(ydata, 2, mean))
      if (iter == mom_switch_iter)
        momentum = final_momentum
      if (iter %% epoch == 0) {
        cost = sum(apply(P * log((P + eps) / (Q + eps)), 1,
                         sum))
        message("Iteration #", iter, " loss function cost is: ",
                cost)
        if (cost < min_cost)
          break
      }
      range = max(abs(ydata))
      if (condition == "tight") {
        if (range > 50 && iter %% 10 == 0) {
          ydata = ydata * 50 / range
        }
      } else {
        if (range > 50 && iter %% max_iter == 0) {
          ydata = ydata * 50 / range
        }
      }
    }
    ydata
  }

#' Calculate affinity matrix
#' @param TPM a TPM matrix with gene names as rownames and cell names as colnames
#' @param LR a dataframe/tibble record the information of ligand receptor pairs, 
#' have to have colnames "ligand", "receptor" and an optional third column with weights
#' 
#' @param denoise numeric value, 
#' @param eps Minimum distances between cells
#' @param verbose logical. If TRUE, print out the progress information
#' @export
#' 
getAffinityMat = function(TPM,
                          LR,
                          denoise = 50,
                          eps = 2.2251e-308,
                          verbose = F,
                          ...) {
  
  genenames = rownames(TPM)
  cellnames = colnames(TPM)
  
  # get the TPM of ligands and receptors
  if(verbose) loginfo("Extracting affinity matrix")
  # ligandsIndex <- match(LR[, 1, drop = T], genenames)
  # receptorIndex <- match(LR[, 2, drop = T], genenames)
  
  flt_LR = LR[(LR[, 1] %in% rownames(TPM)) & (LR[, 2] %in% rownames(TPM)), ]
  
  reverse_flag = flt_LR[, 1] != flt_LR[, 2]
  reverse_LR = flt_LR
  reverse_LR[, 1] = flt_LR[,2]
  reverse_LR[, 2] = flt_LR[,1]
  reverse_LR = reverse_LR[reverse_flag, ]
  
  combn_LR = rbind(
    flt_LR, 
    reverse_LR
  )
  
  ligandsTPM <-
    as.matrix(TPM[combn_LR[, 1], ])
  receptorTPM <-
    as.matrix(TPM[combn_LR[, 2], ])
  
  # determine weight scores
  if(ncol(combn_LR) > 3){
    LRscores = combn_LR[, 3]
  } else {
    LRscores <- rep(1, nrow(combn_LR))
  }
  
  if(verbose) loginfo("Extracting coordinates affinity matrix")
  affinityMat <- t(ligandsTPM) %*% diag(LRscores) %*% receptorTPM
  
  if(verbose) loginfo("Denoising ...")
  # get coordinates through affinity matrix
  for (i in 1:nrow(affinityMat)) {
    affinityArray <- affinityMat[i, ]
    affinityArraySorted <- sort(affinityArray, decreasing = TRUE)
    affinityArray[affinityArray <= affinityArraySorted[denoise]] = 0
    affinityMat[i, ] = affinityArray
  }
  
  # symmetrize P-values
  P = 0.5 * (affinityMat + t(affinityMat))
  P[P < eps] <- eps
  # normalize 
  P = P/sum(P)
  return(P)
}

#' Calculate 3D coordinates from expression
#' A wrapper function to get 3D coordinates directly from expression
#' @param TPM TPM matrix, with gene names as rownames and cell names as colnames
#' @param LR dataframe/tibble; record the information of ligand receptor pairs, have to have colnames "ligand" and "receptor"
#' @param method string; sepcify the optimization method to use. Can be one of 'Rcpp', 'tSNE' or 'BHtSNE'.
#' @param verbose logical. If TRUE, print out the progress information
#' @param ... arguments passsed to different optimization method
#' @export
#' 
#'
getCoordinates = function(TPM, LR, method = 'tSNE', verbose = F, ...) {
  # Get affinity
  affinityMat = getAffinityMat(TPM, LR)
  # optimization
  if(verbose) loginfo("Optimizing coordinates")
  if(method == 'Rcpp'){
    coords <- optimization(affinityMat, verbose = verbose, ...)
  } else if(method == 'tSNE'){
    coords_res = runExactTSNE_R(
      X = affinityMat,
      no_dims = 3,
      max_iter = 1000,
      verbose = verbose, ...
    )
    coords = coords_res$Y
  }
  rownames(coords) <- colnames(TPM)
  colnames(coords) <- c('x', 'y', 'z')
  coords
}


#' get LR contribution for all the listed cluster pairs
#' 
#' @param TPM a TPM matrix with gene names as rownames and cell names as colnames
#' @param LR a dataframe/tibble record the information of ligand receptor pairs, 
#' have to have colnames "ligand", "receptor" and an optional third column with weights
#' @param detailed_connections a list generated by `getSignificance`, 
#' which stored the connected cell pairs for each clutser pair
#' @param verbose logical; whether to print progress
#' @return a named list with sorted LR contributions
#' @export
#' 
getContribution = function(TPM, LR, detailed_connections, verbose = T){
  LR[, 1] = as.character(LR[, 1])
  LR[, 2] = as.character(LR[, 2])
  if(verbose) {loginfo("Extracting data matrix")}
  ligands_existed = intersect(rownames(TPM), LR[, 1])
  receptors_existed = intersect(rownames(TPM), LR[, 2])
  
  flt_LR = LR[(LR[, 1] %in% ligands_existed) & (LR[, 2]%in% receptors_existed), ]
  
  if(ncol(LR) > 2){
    LRscores = flt_LR[, 3]
  } else {
    LRscores = rep(1, nrow(LR))
  }
  
  if(verbose) {loginfo("Calculate contribution of ", length(detailed_connections), " cluster pairs.")}
  LR_contri_lst = list()
  # calculate contribution cell-pair by cell-pair
  for(target_clusterPair in names(detailed_connections)){
    if(nrow(detailed_connections[[target_clusterPair]]) < 3){
      warning("Number of connected cells in ", target_clusterPair, " is lower than 3.\n")
    }
    L1S = TPM[flt_LR[, 1], detailed_connections[[target_clusterPair]][, 1], drop = F]
    R1S = TPM[flt_LR[, 2], detailed_connections[[target_clusterPair]][, 1], drop = F]
    L1R = TPM[flt_LR[, 1], detailed_connections[[target_clusterPair]][, 2], drop = F]
    R1R = TPM[flt_LR[, 2], detailed_connections[[target_clusterPair]][, 2], drop = F]
    
    all_intensity = L1S * LRscores * R1R + R1S * LRscores * L1R 
    rownames(all_intensity) = paste0(flt_LR[, 1], "---", flt_LR[, 2])
    contribution_mt = t(t(all_intensity) / (colSums(all_intensity)))
    
    contribution_forCluster = sort(rowSums(contribution_mt) / ncol(contribution_mt), decreasing = T)
    # head(contribution_forCluster)
    LR_contri_lst[[target_clusterPair]] = contribution_forCluster
    if(which(names(detailed_connections) == target_clusterPair) %% 100 ==0 && verbose){
      loginfo(sprintf("%d/%d cluster pairs calculated\n", which(names(detailed_connections) == target_clusterPair), length(names(detailed_connections))))
    }
  }
  # clear memory every loop 
  rm(L1S, L1R, R1S, R1R, all_intensity, contribution_mt, contribution_forCluster)
  invisible(gc())
  
  return(LR_contri_lst)
}

#' Calculate normalized connection based on connection matrix and cell counts
#' 
#' @param connection_mt Named matrix
#' @param cell_count_table Cell count table generated by `table()` or a named vector recording the cell counts
#' @return Named matrix of normalized connection
#' @export
calcNormalizedConnection = function(connection_mt, cell_count_table){
  stopifnot(all(rownames(connection_mt) == colnames(connection_mt)))
  cell_count_table = cell_count_table[rownames(connection_mt)]
  count_mt = lapply(rownames(connection_mt), function(x){
    cell_count_table[x] * cell_count_table
  }) %>% do.call(what = rbind, .)
  rownames(count_mt) = rownames(connection_mt)
  count_mt = count_mt[, colnames(connection_mt)]
  diag(count_mt) = cell_count_table * (cell_count_table-1)
  normalized_connection = connection_mt/count_mt
  return(normalized_connection)
}

# Density ----
#' get 3D density
#' Estimate 3D density around each data point based on coordinates. 
#' @param x,y,z coordinates 
#' @param n numbers of grid points to use for each dimension; recycled if length is less than 3. Default 100.
#' @param ... other parameters passed to kde3d
#' @export
getDensity3D = function(x, y, z, n = 100, ...) {
  tryCatch({
    dens <- kde3d(x = x, y = y, z =z, n = n, ...)
  }, error = function(e) {
    print(e)
    warning("Swith bandwidth to h = 1")
    dens <<- kde2d(x = x,
                         y = y,
                         n = n,
                         h = 1)
  })
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  iz <- findInterval(z, dens$z)
  ii <- cbind(ix, iy, iz)
  return(dens$d[ii])
}


# function adapted from R package misc3d
# https://cran.r-project.org/web/packages/misc3d/index.html
kde3d = function (x, y, z, h, n = 20, lims = c(range(x), range(y), range(z))) 
{
  nx <- length(x)
  if (length(y) != nx || length(z) != nx) 
    stop("data vectors must be the same length")
  if (missing(h)) 
    h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y), 
           MASS::bandwidth.nrd(z))/6
  else if (length(h) != 3) 
    h <- rep(h, length = 3)
  if (length(n) != 3) 
    n <- rep(n, length = 3)
  if (length(lims) == 2) 
    lims <- rep(lims, length = 6)
  gx <- seq(lims[1], lims[2], length = n[1])
  gy <- seq(lims[3], lims[4], length = n[2])
  gz <- seq(lims[5], lims[6], length = n[3])
  mx <- matrix(outer(gx, x, dnorm, h[1]), n[1], nx)
  my <- matrix(outer(gy, y, dnorm, h[2]), n[2], nx)
  mz <- matrix(outer(gz, z, dnorm, h[3]), n[3], nx)
  v <- array(0, n)
  tmy.nx <- t(my)/nx
  for (k in 1:n[3]) {
    tmy.nz.zk <- tmy.nx * mz[k, ]
    v[, , k] <- mx %*% tmy.nz.zk
  }
  return(list(x = gx, y = gy, z = gz, d = v))
}

# function adapted from R package MASS
kde2d = function (x, y, h, n = 25, lims = c(range(x), range(y))) 
# https://cran.r-project.org/web/packages/MASS/index.html
{
  nx <- length(x)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  if (any(!is.finite(x)) || any(!is.finite(y))) 
    stop("missing or infinite values in the data are not allowed")
  if (any(!is.finite(lims))) 
    stop("only finite values are allowed in 'lims'")
  n <- rep(n, length.out = 2L)
  gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
  gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
  h <- if (missing(h)) 
    c(bandwidth.nrd(x), bandwidth.nrd(y))
  else rep(h, length.out = 2L)
  if (any(h <= 0)) 
    stop("bandwidths must be strictly positive")
  h <- h/4
  ax <- outer(gx, x, "-")/h[1L]
  ay <- outer(gy, y, "-")/h[2L]
  z <- tcrossprod(matrix(dnorm(ax), , nx), matrix(dnorm(ay), , nx))/(nx * h[1L] * h[2L])
  list(x = gx, y = gy, z = z)
}


# Inform functions ----
loginfo <- function(..., printnow = T) {
  msg = paste0(list(...), collapse = "")
  msg <- paste0("[",format(Sys.time()), "] ", msg,"\n")
  if(printnow)
    cat(msg)
  invisible(msg)
}

# Others ----
paste2columns = function(x, y, delim = "---") {
  stopifnot(length(x) == length(y))
  x = as.character(x)
  y = as.character(y)
  z = c()
  for (i in 1:length(x)) {
    z =  c(z, paste0(sort(c(x[i], y[i])), collapse = delim))
  }
  z
}

# 3DPlot ----
#' Plot 3D figure using plotly
#' A wrapper function to plot 3D plots using plotly. Not necessary for CSOmap's core functions.
#' 
#' @param plt_tbl data.frame/tibble; Should provide coordinates x,y,z.
#' @param color_by string; Specify that by which columns should the data points will be colored 
#' @param title string; Title
#' @param alpha numeirc; 0-1 specify the alpha of dots
#' @param save_path string; Speicfy the saving path of the output 3D plot. a `/lib/` will also be generated with the output html. Default = NULL.
#' @param ... Other arguments that will be passed to htmlwidgets::saveWidget
#' @return a plotly object
#' @export
plot3D = function(plt_tbl,
                  color_by = "density",
                  title = "3D density",
                  alpha = 0.8,
                  save_path = NULL,
                  ...) {
  
  if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("htmlwidgets")) {
    stop("Package \"plotly\" and \"htmlwidgets\" are needed for this function to work. Please install it.")
  }
  fig_density = plotly::plot_ly(
    plt_tbl,
    x = ~ x,
    y = ~ y,
    z = ~ z,
    alpha = alpha
  ) 
  fig_density = plotly::add_markers(fig_density, color = eval(parse(text = sprintf("~%s", color_by))))
  fig_density = plotly::layout(fig_density, title = title)
  
  if(!is.null(save_path)){
    htmlwidgets::saveWidget(
      fig_density,
      file = save_path,
      selfcontained = F,
      libdir = paste0(dirname(save_path), "/lib/")
    )
  }
  invisible(fig_density)
}