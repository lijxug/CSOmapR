CSOmap <- function(DataSetName) {
  library(plotly)
  
  # load data ----
  TPMpath <- paste0("./data/", DataSetName, "/TPM.txt")
  LRpath <- paste0("./data/", DataSetName, "/LR_pairs.txt")
  Labelpath <- paste0("./data/", DataSetName, "/label.txt")
  TPM_tbl <- read.table(TPMpath, header = TRUE, sep = "\t")
  LR <- read.table(LRpath, header = FALSE, sep = "\t")
  labelData <- read.table(Labelpath, header = TRUE, sep = "\t")
  # create output path
  dir.create(paste0("./results/", DataSetName))
  
  
  # set variables properly ----
  TPM = as.matrix(TPM_tbl[,-1])
  TPM[is.na(TPM)] = 0
  
  rownames(TPM) = TPM_tbl$X
  # genenames <- TPM$X
  # cellnames <- colnames(TPM)
  labels <- labelData$labels[match(colnames(TPM), labelData$cells)]
  labels[is.na(labels)] = "unlabeled"
  standards <- unique(labels)
  labelIx <- match(labels, standards)
  cellCounts <- table(labelIx)
  
  coords = getCoordinates(TPM, LR)
  
  # coordsPath <-
  #   paste0("./results/", DataSetName, "/coordinates.txt")
  # write coords
  # write.table(coords, coordsPath, quote = FALSE, sep = "\t")
  
  # test for significance level----
  significance_result = getSignificance(coords, labels, 3)
  
  return(significance_result)
}

