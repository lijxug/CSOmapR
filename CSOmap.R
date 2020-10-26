#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("methods")) # Rscript for CMD doesn't load this automatically

# ---- inputs ----
suppressWarnings(library("optparse"))

parser <- OptionParser()
option_list <- list( 
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false", 
              dest="verbose", help="Print little output"),
  make_option(c("-n", "--nCore"), type="integer", default=1,
              help="Number of cores to use. Under development. [default %default]",
              metavar="number"),
  make_option(c("-e", "--version"), type="character", default="cpp",
              help="version of optimization to use. origin or cpp [default %default]",
              metavar="number")
)

inputs = commandArgs(trailingOnly=FALSE)
script_dir = dirname(sub("--file=", "", inputs[grep("--file=",inputs)]))
if(!length(script_dir)){
  script_dir = getwd()
}
# Loading dependencies & functions ----


# get real args
parser = OptionParser(
  usage =
    "%prog [options] <TPMpath> <LRpath> <Labeltblpath> <OUTPUT_DIR>
Description: 
  Labeltblpath: a path pointed to a label table, with the first column being the corresponding cellID of TPM matrix. 
    Every following columns will be parsed as a label vector and passed to `getSignificance` one by one.
  This script is designed for running CSOmapR in console.",
  option_list = option_list
) # donot change the format of this doc.

itfs = parse_args(parser, positional_arguments = TRUE) 
opts = itfs$options
args = list(
  TPMpath = itfs$args[1],
  LRpath = itfs$args[2],
  Labeltblpath = itfs$args[3],
  output_dir = itfs$args[4]
)

if(!dir.exists(args$output_dir)) {
    dir.create(args$output_dir)
}

# if none arguments input
if(anyNA(args)){
  warning("Missing arguments!")
	parse_args(parser, args = c("--help"))
	q("no", status = 2)
}


# ---- loading dependencies ----
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("tictoc"))
suppressPackageStartupMessages(library("Rcpp"))
source("/lustre1/zeminz_pkuhpc/lijiesheng/00.data/05.COVID19/01.codes/CSOmapR/utils/utils.R") # hardcoded, should be modified in the future

# ---- main ----
loginfo("Loading input file ...")
tic()
TPM_tbl <- read_tsv(args$TPMpath)
LR <- read.table(args$LRpath, header = FALSE, sep = "\t")
labelData <- read_tsv(args$Labeltblpath)
toc()

# set variables & validation ----
TPM = as.matrix(TPM_tbl[,-1])
rownames(TPM) = TPM_tbl[, 1, drop = T] 
TPM[is.na(TPM)] = 0

# genenames <- TPM$X
# cellnames <- colnames(TPM)

# check labels
# stopifnot(all(colnames(TPM) %in% labelData[, 1, drop = T]))
  
loginfo("CSOmap started ...")
# get coordinates ----
coords = getCoordinates(TPM, LR, version = opts$version)
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))
coords_outdir = paste0(args$output_dir, "/coordinates.txt")
loginfo("Writing coordintates to ", coords_outdir)
write_tsv(coords_tbl, path = coords_outdir)

join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)
cellinfo_outdir = paste0(args$output_dir, "/cell_infos.tsv")
write_tsv(cellinfo_tbl, path = cellinfo_outdir)

# get significance ----
# for(iCol in 2:ncol(labelData)) {
iCol = 2 # ignore the extra columns
labels <-
  labelData[, iCol, drop = T][match(colnames(TPM), labelData[, 1, drop = T])]
labels[is.na(labels)] = "unlabeled"
standards <- unique(labels)
labelIx <- match(labels, standards)
cellCounts <- table(labelIx)

# Get significance
signif_res = getSignificance(coords, labels = labels)
write_csv(
  signif_res$pvalue_tbl,
  path = paste0(
    args$output_dir,
    "/signif_interacting_clusters_",
    colnames(labelData)[iCol],
    ".csv"
  )
)
# }

loginfo("Work Done!")
