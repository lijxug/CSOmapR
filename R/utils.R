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
    
    
    sub_connects_mt = connects_mt[cellsfrom1, cellsfrom2]
    
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
  i = lower.tri(p_value, diag = T)
  p_value_df = 
    data.frame(cluster_pair = cluster_pair[i], 
               p.value = p_value[i], 
               q.value = p.adjust(p_value[i], method = adjusted.method))
  
  q_value = p_value
  q_value[i] = p_value_df$q.value
  q_value = t(q_value)
  q_value[i] = p_value_df$q.value
  
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
      if (verbose) processBar(iter01, iter, max_iter, tail = "ETA")
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
  if(method = 'Rcpp'){
    coords <- optimization(affinityMat, verbose = verbose, ...)
  } else if(method == 'tSNE'){
    
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
#' @return a named list with sorted LR contributions
#' @export
#' 
getContribution = function(TPM, LR, detailed_connections){
  LR[, 1] = as.character(LR[, 1])
  LR[, 2] = as.character(LR[, 2])
  
  ligands_existed = intersect(rownames(TPM), LR[, 1])
  receptors_existed = intersect(rownames(TPM), LR[, 2])
  
  flt_LR = LR[(LR[, 1]%in%ligands_existed) & (LR[, 2]%in% receptors_existed), ]
  
  if(ncol(LR) > 2){
    LRscores = flt_LR[, 3]
  } else {
    LRscores = rep(1, nrow(LR))
  }
  
  LR_contri_lst = list()
  # calculate contribution cell-pair by cell-pair
  for(target_clusterPair in names(detailed_connections)){
    L1S = TPM[flt_LR[, 1], detailed_connections[[target_clusterPair]][, 1]]
    R1S = TPM[flt_LR[, 2], detailed_connections[[target_clusterPair]][, 1]]
    L1R = TPM[flt_LR[, 1], detailed_connections[[target_clusterPair]][, 2]]
    R1R = TPM[flt_LR[, 2], detailed_connections[[target_clusterPair]][, 2]]
    
    all_intensity = L1S * LRscores * R1R + R1S * LRscores * L1R 
    rownames(all_intensity) = paste0(flt_LR[, 1], "---", flt_LR[, 2])
    contribution_mt = t(t(all_intensity) / (colSums(all_intensity)))
    
    contribution_forCluster = sort(rowSums(contribution_mt) / ncol(contribution_mt), decreasing = T)
    # head(contribution_forCluster)
    LR_contri_lst[[target_clusterPair]] = contribution_forCluster
  }
  
  return(LR_contri_lst)
}

# Density ----
#' get 3D density
#' 
#' Wrapper for misc3d::kde3d to get density estimation
#' 
#' @param x,y,z coordinates 
#' @param n numbers of grid points to use for each dimension; recycled if length is less than 3.
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
initiatePB = function(iterOBJ){
  .tic = sprintf("%s", paste0(".TIC_", iterOBJ))
  rm_list = c(iterOBJ, .tic)
  if(any(exists(rm_list, inherits = T)))
    rm(list = rm_list, inherits = T)
}

processBar = function(objName,
                      i,
                      cycles,
                      title = "Process",
                      scale = 40,
                      sign = "#", 
                      tail = "",
                      terminal = "R", # terminal could be R/log, others default to shell
                      final = "Work done!") {
  suppressPackageStartupMessages(require("iterators"))
  if (!exists(objName)) {
    if (terminal != "R")
      words_list = unlist(lapply(1:cycles, function(x) {
        sprintf(
          paste0("\033[?25l\r%s %5.1f%% | %-", scale, "s | "),
          title,
          x * 100 / cycles ,
          paste0(rep(sign, ceiling(x * scale / cycles)), collapse = "")
        )
      }))#\033[?25l hide the cursor - linux control code
    else
      words_list = unlist(lapply(1:cycles, function(x) {
        sprintf(
          paste0("\r%s %5.1f%% | %-", scale, "s | "),
          title,
          x * 100 / cycles ,
          paste0(rep(sign, ceiling(x * scale / cycles)), collapse = "")
        )
      }))
    eval(parse(text = sprintf("%s <<- iter(words_list)", objName)))
    eval(parse(text = sprintf("%s <<- Sys.time()", paste0(".TIC_",objName))))
    # if i didn't start at 1
    times = i
    while (times > 1) {
      msg = eval(parse(text = sprintf("nextElem(%s)", objName)))
      times = times - 1
    }
  }
  
  msg = eval(parse(text = sprintf("nextElem(%s)", objName)))
  if(tail == "ETA"){
    .tic = eval(parse(text = sprintf("%s", paste0(".TIC_",objName))))
    if(terminal != "R"){
      tail = paste0("ETA: ", format(round((Sys.time() - .tic) / i * (cycles - i), digits = 2)), "\033[K")
    }
    else {
      tail = paste0("ETA: ", format(round((Sys.time() - .tic) / i * (cycles - i), digits = 2)), "   ")
    }
  }
  if(terminal == "log")
    tail = paste0(tail, "\n")
  cat(paste0(msg, tail))
  if(i == cycles){
    if(nchar(final)) final = loginfo(final, printnow = F)
    if(terminal != "R") cat("\033[?25h")
    cat(paste0("\n", final))
    rm(list = objName, inherits = T)
  }
}


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
