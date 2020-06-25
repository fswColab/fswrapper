library(tools)
library(FlowSOM)
library(ggplot2)
library(pheatmap)

GetFlowFrame <- function(fcsFilePath, colsToUse, asinh_scale=150, seed=1) {
  #' Return flowFrame of specified .fcs file.
  #' 
  #' @param fcsFilePath file path of .fsc file.
  #' @param colsToUse column indices to prepare flowFrame matrix.
  #' @param asinh_scale scale of asinh transformation. Default to "150"
  #' @param seed seed int used to standardize result of different files. Default to "1"
  #'
  
  # set seed int
  set.seed(1)
  
  # build matrix
  mat=flowCore::exprs(read.FCS(fcsFilePath, transformation=FALSE,truncate_max_range=FALSE))
  
  # apply transformation to selected columns
  mat[, colsToUse]=asinh(mat[, colsToUse] / asinh_scale)
  
  return(flowFrame(mat))
}


SanitizeMfisData <- function(fSOM, seed=1) {
  #' Sanitize and return MFIs data
  #' 
  #' @param fSOM FlowSOM object.
  #' @param seed seed int used to standardize result of different files. Default to "1"
  #' 
  
  # set seed int
  set.seed(1)
  
  mfis <- MetaclusterMFIs(fSOM)
  
  data_mfis=as.data.frame(mfis)
  data_mfis$`FSC-A`=NULL
  data_mfis$`FSC-W`=NULL
  data_mfis$`FSC-H`=NULL
  data_mfis$`SSC-A`=NULL
  data_mfis$`SSC-W`=NULL
  data_mfis$`SSC-H`=NULL
  data_mfis$Time=NULL
  
  return(data_mfis)
}


GetFilenames <- function(targetPath) {
  #' Return list of filenames of a given path
  #' 
  #' @param targetPath directory path used for lookup
  #' 

  l=list.files(path = targetPath, pattern = NULL, all.files = FALSE,
                             full.names = FALSE, recursive = FALSE,
                             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  return(as.list(l))
}


GenPHM <- function(fSOM, fileName, phmPath) {
  #' Use pheatmap tp generate and store heat of fcs files
  #' 
  #' @param fSOM FlowSOM object.
  #' @param fileName file name of generated heatmap
  #' @param phmPath pheatmap files directory path
 
  pheatmap(SanitizeMfisData(fSOM), scale = "none", fontsize=5, file=sprintf("%s/%s.pdf", phmPath, fileName), bg="transparent", useDingbats=FALSE, cluster_rows = FALSE,
           cluster_cols = FALSE)
}

GetOutLiersPercent <- function(fSOM_new, fSOM) {
  #' Return percentage of outliers cell of new fsom
  #' 
  #' @param fSOM_new fsom added to the refference one
  #' @param fSOM reference FlowSOM object
  #' 
  
  test_outliers <- TestOutliers(fSOM_new, mad_allowed = 4, fsom_reference = fSOM)
  max_outliers <- max(test_outliers$Number_of_outliers)
  n_outliers <- sum(test_outliers$Number_of_outliers)
  
  return(round(n_outliers/nrow(fSOM_new[[1]]$data) *  100, 2))
}



# ==================
# Begining of script
# ==================

# some vars
fcsDir <- c("fcs")
phmDir <- c("phm")
marker_cols <- (7:19)
k <- 21

# get list of all fcs files
fcsList <- GetFilenames(fcsDir)

# create first flowSOM object
fFrame <- GetFlowFrame(sprintf("%s/%s", fcsDir, fcsList[[1]]), marker_cols)

set.seed(1)
fSOM_frame <- FlowSOM(fFrame,scale=FALSE,colsToUse=marker_cols,nClus=k)

# generate pheatmap
GenPHM(fSOM_frame, tools::file_path_sans_ext(fcsList[[1]]), phmDir)

# create metacluster dataframe
metaClusterTable = table(GetMetaclusters(fSOM_frame))
metaClusterDataFrame = as.data.frame(metaClusterTable)[,2]

# instantiate new objets
fSOM <- fSOM_frame
dataFrame <- metaClusterDataFrame

# instantiate outliers df
outliersPercentReport <- data.frame(matrix(ncol = length(fcsList), nrow = 1))

for(f in 2:length(fcsList)){
  
  # get base name of fcs file
  fBaseName <- tools::file_path_sans_ext(fcsList[[f]])
  
  # get new flowSOM
  fSOM <- NewData(fSOM_frame, GetFlowFrame(sprintf("%s/%s", fcsDir, fcsList[[f]])))
  
  # fill outliers list with percentage of outliers of current FlowSOM object
  outliersPercentReport[f] <- GetOutLiersPercent(fSOM, fSOM_frame)
  
  # Gen pheatmap
  GenPHM(fSOM, tools::file_path_sans_ext(fcsList[[f]]), phmDir)
  
  
  dataFrame <- cbind(dataFrame, as.data.frame(table(GetMetaclusters(fSOM)))[,2])
  
}

# rename some cols
colnames(dataFrame) <- fcsList
colnames(outliersPercentReport) <- fcsList



