# CSOmapR
[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](https://github.com/lijxug/CSOmapR/blob/master/LICENSE)

---

R package of CSOmap(developing

# Installation

``` r
# install.packages("devtools")
devtools::install_github("lijxug/CSOmapR")
```

# Usage 

## Load demo dataset
``` r
invisible(TPM)
invisible(LR)
invisible(labelData)
```

## Calculate optimized 3D coordinates
``` r
affinityMat = getAffinityMat(TPM, LR)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = F
)
coords = coords_res$Y

```

## Visualization(by 3D density)
``` r
require(dplyr)
# arrange data
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

p_3Ddensity = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")
```

# When dealing with large dataset
We provide two options: tSNE(BH algorithm) or downsampling
``` r

```