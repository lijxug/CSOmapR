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
getSignificance = function(coordinates, labels, k = 3, verbose = T) {
  require(reshape2)
  require(dplyr)
  # preprocess
  labels = setNames(labels, nm = rownames(coordinates))
  standards <- unique(labels)
  labelIx <- match(labels, standards)
  cellCounts <- table(labelIx)
  
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
  for (i in 1:nrow(counts)) {
    for (j in 1:ncol(counts)) {
      if (i == j) {
        M <- as.numeric(cellCounts[i]) * (as.numeric(cellCounts[i]) - 1) / 2
      } else {
        M <- as.numeric(cellCounts[i]) * (as.numeric(cellCounts[j]))
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
               q.value = p.adjust(p_value[i]))
  
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


#' Optimize the 3D coordinates(origin)
#' 
#' this function is inspired from tsne algorithm. We use similar gradient descent 
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
#' 
optimization_origin <-
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
            verbose = F) 
{ # this function is inspired from tsne algorithm. We use similar gradient descent
  # method to optimize our target function specified in our paper.
  # condition can be loose or tight, we suggest using "loose" condition 
  # for dataset with over 10000 cells
  n = dim(affinityMat)[1]
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
  P = 0.5 * (affinityMat + t(affinityMat))
  P[P < eps] <- eps
  P = P/sum(P)
  grads = matrix(0, nrow(ydata), ncol(ydata))
  incs = matrix(0, nrow(ydata), ncol(ydata))
  gains = matrix(1, nrow(ydata), ncol(ydata))
  for (iter in 1:max_iter) {
    
    sum_ydata = apply(ydata^2, 1, sum)
    d = sum_ydata + sweep(-2 * ydata %*% t(ydata), 2, -t(sum_ydata))
    num = 1/(1 + d)
    diag(num) = 0
    Q = num/sum(num)
    if (any(is.nan(num))) 
      message("NaN in grad. descent")
    Q[Q < eps] = eps
    P_Q = P - Q
    P_Q[P_Q > 0 & d <= 0.01] = -0.01;
    stiffnesses = 4 * (P - Q) * num
    # stiffnesses = 4 * P_Q * num
    for (i in 1:n) {
      grads[i, ] = apply(sweep(-ydata, 2, -ydata[i, ]) * 
                           stiffnesses[, i], 2, sum)
    }
    gains = ((gains + 0.2) * abs(sign(grads) != sign(incs)) + 
               gains * 0.8 * abs(sign(grads) == sign(incs)))
    gains[gains < min_gain] = min_gain
    incs = momentum * incs - epsilon * (gains * grads)
    ydata = ydata + incs
    ydata = sweep(ydata, 2, apply(ydata, 2, mean))
    if (iter == mom_switch_iter) 
      momentum = final_momentum
    if (iter%%epoch == 0) {
      cost = sum(apply(P * log((P + eps)/(Q + eps)), 1, 
                       sum))
      message("Iteration #", iter, " loss function cost is: ", 
              cost)
      if (cost < min_cost) 
        break
    }
    range = max(abs(ydata))
    if (condition == "tight") {
      if (range > 50 && iter%%10 == 0) {
        ydata = ydata * 50/range
      }
    } else {
      if (range > 50 && iter%%max_iter == 0) {
        ydata = ydata * 50/range
      }
    }
  }
  ydata
}

#' Optimize the 3D coordinates(cpp)
#' 
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
    P = 0.5 * (affinityMat + t(affinityMat))
    P[P < eps] <- eps
    P = P / sum(P)
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

#' Get affinityMatrix
#' @param TPM a TPM matrix with gene names as rownames and cell names as colnames
#' @param LR a dataframe/tibble record the information of ligand receptor pairs, 
#' have to have colnames "ligand", "receptor" and an optional third column with weights
#' @param verbose logical. If TRUE, print out the progress information
#' 
getAffinityMat = function(TPM, LR, verbose = F, ...) {
  genenames = rownames(TPM)
  cellnames = colnames(TPM)
  
  # get the TPM of ligands and receptors
  if(verbose) loginfo("Extracting affinity matrix")
  ligandsIndex <- match(LR[, 1, drop = T], genenames)
  receptorIndex <- match(LR[, 2, drop = T], genenames)
  ligandsTPM <-
    as.matrix(TPM[ligandsIndex[!is.na(ligandsIndex) &
                                 !is.na(receptorIndex)], ])
  receptorTPM <-
    as.matrix(TPM[receptorIndex[!is.na(ligandsIndex) &
                                  !is.na(receptorIndex)], ])
  
  # determine weight scores
  if(ncol(LR) > 3){
    LRscores = LR[!is.na(ligandsIndex) & !is.na(receptorIndex), 3]
  } else {
    LRscores <- rep(1, nrow(LR))[!is.na(ligandsIndex) & !is.na(receptorIndex)]
  }
  
  affinityMat <- t(ligandsTPM) %*% diag(LRscores) %*% receptorTPM
  
  # get coordinates through affinity matrix
  if(verbose) loginfo("Extracting coordinates affinity matrix")
  for (i in 1:nrow(affinityMat)) {
    affinityArray <- affinityMat[i, ]
    affinityArraySorted <- sort(affinityArray, decreasing = TRUE)
    affinityArray[affinityArray <= affinityArraySorted[50]] = 0
    affinityMat[i, ] = affinityArray
  }
  affinityMat
}

#' Calculate 3D coordinates from expression
#'
#' @param TPM a TPM matrix with gene names as rownames and cell names as colnames
#' @param LR a dataframe/tibble record the information of ligand receptor pairs, have to have colnames "ligand" and "receptor"
#' @param optimiztion a string sepcify the optimization method to use. Can be one of 'origin.R','Rcpp', 'tSNE' or 'umap'.
#' @param verbose logical. If TRUE, print out the progress information
#' @param ... arguments passsed to optimization
#' @export
#' 
#'
getCoordinates = function(TPM, LR, optimization = 'Rcpp', verbose = F, ...) {
  # Get affinity
  affinityMat = getAffinityMat(TPM, LR)
  # optimization
  if(verbose) loginfo("Optimizing coordinates")
  if(version == 'origin'){
    coords <- optimization_origin(affinityMat, verbose = verbose, ...)
  } else if(version == "tSNE"){
    coords = Rtsne::Rtsne(affinityMat, ...)
  } else if(version == "umap"){
    coords = umap::umap(affinityMat, ...)
  } else {
    coords <- optimization(affinityMat, verbose = verbose, ...)
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

#' #' calcualte LR contribution
#' #' 
#' calculateContribution = function(TPM, )

# Calculations: R version ----
calc_d = function(ydata) {
  sum_ydata = apply(ydata ^ 2, 1, sum)
  d = sum_ydata + sweep(-2 * ydata %*% t(ydata), 2,-t(sum_ydata))
  # d = sum_ydata + sweep(-2*tcrossprod(ydata), 2,-sum_ydata)
  return(d)
}

update_grads = function(grads, ydata, stiffnesses) {
  for (i in 1:nrow(grads)) {
    grads[i, ] = apply(sweep(-ydata, 2, -ydata[i, ]) *
                         stiffnesses[, i], 2, sum)
  }
  return(grads)
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
  require(misc3d)
  tryCatch({
    dens <- misc3d::kde3d(x = x, y = y, z =z, n = n, ...)
  }, error = function(e) {
    print(e)
    warning("Swith bandwidth to h = 1")
    dens <<- MASS::kde2d(x = x,
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