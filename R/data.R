#' TPM matrix of 1000 genes with the expression of 23,686 genes
#'
#' A dataset containing the expression data of 1000 cells.
#'
#' @format A matrix with 23686 rows and 1000 variables:
"TPM"

#' LR data frame
#' 
#' A data frame of 1961 ligand-receptor pairs.
#'
#' @format A data frame with 1961 rows and 3 variables:
#' \describe{
#'   \item{V1}{Gene names of ligands}
#'   \item{V2}{Gene names of receptors}
#'   \item{V3}{Weigths of the ligand-receptor pairs}
#' }
"LR"

#' label data frame
#' 
#' A data frame of the labels of cells, corresponding to cells in TPM.
#'
#' @format A data frame with 4645 rows and 2 variables:
#' \describe{
#'   \item{cells}{Cell names}
#'   \item{labels}{Cell labels}
#' }
"LR"