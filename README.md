# CSOmapR
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/lijxug/CSOmapR/blob/master/LICENSE)

---

R package for CSOmap(developing)

Right now, CSOmapR is only available on the linux-based machines. Installation on windows may encounter errors.

# Installation

``` r
# install.packages("devtools")
devtools::install_github("lijxug/CSOmapR")

# install CSOmapR.demo package for easy-loading demo data
# devtools::install_github("lijxug/CSOmapR.demo")
```

# Usage 

## Load demo dataset
``` r
library(CSOmapR)
library(CSOmapR.demo)
invisible(TPM)
invisible(LR)
invisible(labelData)
```

## Calculate optimized 3D coordinates
``` r
affinityMat = getAffinityMat(TPM, LR, verbose = T)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)
coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')

```

## Visualization(by 3D density)
``` r
require(dplyr)
# arrange data
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")

```

## Get significance
``` r
signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)
```

# When dealing with large dataset
We provide two options: optimize coordinates through tSNE(BH algorithm), or downsample the original huge dataset first.

``` r
# under development
```

# Citation 
Ren, X., Zhong, G., Zhang, Q., Zhang, L., Sun, Y., and Zhang, Z. (2020). Reconstruction of cell spatial organization from single-cell RNA sequencing data based on ligand-receptor mediated self-assembly. Cell Res.


