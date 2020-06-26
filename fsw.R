library(flowCore)
library(tools)
library(FlowSOM)
library(gridExtra)
library(pheatmap)


DfToPdf <- function(df, filename, dirPath, title=NULL) {
  #' Create pdf based on given dataframe
  #' 
  #' @param df dataframe to use
  #' @param filename filename of pdf
  #' @param dirPath path of pdf directory
  #' @param title pdf title. Default to "NULL". If not specified, filename wil be used.
  #' 
  
  # set title with filename if NULL
  if (is.null(title)) {
    title <- filename
  }
  
  # generate PDF
  pdf(sprintf("%s/%s.pdf", dirPath, filename, title=title), height=11, width=10)
  grid.table(df)
  dev.off()
}


GenPheatmap <- function(fSOM, fcsFilename, pheatmapDirPath) {
  #' Use pheatmap to generate and store heat of fcs files
  #' 
  #' @param fSOM FlowSOM object.
  #' @param fcsFilename .fcs file name used to generate pheatmap
  #' @param phmPath pheatmap files directory fullpath
  
  fname = tools::file_path_sans_ext(fcsFilename)
  pheatmap(SanitizeMfisData(fSOM), scale = "none", fontsize=5, file=sprintf("%s/%s.pdf", phmDirPath, fname), bg="transparent", useDingbats=FALSE, cluster_rows = FALSE,
           cluster_cols = FALSE)
}


GetFcsFilenames <- function(targetPath) {
  #' Return list of filenames of a given path
  #' 
  #' @param targetPath directory path used for lookup
  #' 
  
  l=list.files(path = targetPath, pattern = "\\.fcs$", all.files = FALSE,
               full.names = FALSE, recursive = FALSE,
               ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  return(as.list(l))
}


GetFlowFrame <- function(fcsFilePath, colsToUse, asinh_scale=150, seed=1) {
  #' Return flowFrame of specified .fcs file.
  #' 
  #' @param fcsFilePath fullpath of .fsc file.
  #' @param colsToUse column indices to prepare flowFrame matrix.
  #' @param asinh_scale scale of asinh transformation. Default to "150"
  #' @param seed seed int used for reproducibility. Default to "1"
  #'
  
  # set seed int
  set.seed(seed)
  
  # build matrix
  mat=flowCore::exprs(read.FCS(fcsFilePath, transformation=FALSE,truncate_max_range=FALSE))
  
  # apply transformation to selected columns
  mat[, colsToUse]=asinh(mat[, colsToUse] / asinh_scale)
  
  return(flowFrame(mat))
}



GetOutLiersPercent <- function(fSOM, refFSOM) {
  #' Return percentage of outliers cell of new fsom
  #' 
  #' @param fSOM fsom added to the refference one
  #' @param refFSOM reference FlowSOM object
  #' 
  
  test_outliers <- TestOutliers(fSOM, mad_allowed = 4, fsom_reference = refFSOM)
  max_outliers <- max(test_outliers$Number_of_outliers)
  n_outliers <- sum(test_outliers$Number_of_outliers)
  
  return(round(n_outliers/nrow(fSOM[[1]]$data) *  100, 2))
}


GetRefFSOM <- function(refFcsPath, colsToUse, nClus, seed=1) {
  #' Return reference FlowSOM object
  #' 
  #' @param refFcsPath reference fcs file path
  #' @param colsToUse arg pass to FlowSOM
  #' @param nClus arg pass to FlowSOM
  #' @param seed seed int used for reproducibility. Default to "1"
  
  ff <- GetFlowFrame(refFcsPath, marker_cols)
  set.seed(seed)
  fs <- FlowSOM(ff,scale=FALSE,colsToUse=colsToUse,nClus=nClus)
  return(fs)
}


SanitizeMfisData <- function(fSOM, seed=1) {
  #' Sanitize and return MFIs data
  #' 
  #' @param fSOM FlowSOM object.
  #' @param seed seed int used for reproducibility. Default to "1"
  #' 
  
  # Todo: Could we use marker_cols to trim columns instead of hardcode?
  
  # set seed int
  set.seed(seed)
  
  mfis <- MetaclusterMFIs(fSOM)
  
  data_mfis=as.data.frame(mfis)
  data_mfis$`FSC-A` <- NULL
  data_mfis$`FSC-W` <- NULL
  data_mfis$`FSC-H` <- NULL
  data_mfis$`SSC-A` <- NULL
  data_mfis$`SSC-W` <- NULL
  data_mfis$`SSC-H` <- NULL
  data_mfis$Time <- NULL
  
  return(data_mfis)
}



# ==================
# Begining of script
# ==================

# import params
source('params.R')


# computed var
phmDirPath <- sprintf("%s/%s", fcsDirPath, phmDirName)
extraDirPath <- sprintf("%s/%s", fcsDirPath, extraDirName)

# create phm and extra directories
dir.create(phmDirPath, showWarnings = FALSE)
dir.create(extraDirPath, showWarnings = FALSE)

# get list of all fcs filenames
fcsList <- GetFcsFilenames(fcsDirPath)

# set reference FlowSOM object
refFSOM <- GetRefFSOM(sprintf("%s/%s", fcsDirPath, referenceFcsFileName),marker_cols, nClus)
fSOM <- refFSOM

# gen pheatmap
GenPheatmap(refFSOM, referenceFcsFileName, phmDirPath)

# create metacluster dataframe
metaClusterTable <- table(GetMetaclusters(refFSOM))
metaClusterDF <- as.data.frame(metaClusterTable)[,2]



# Main loop
orderedSampleList <- list(referenceFcsFileName)
orderedIndex <- 1
outliersPercentReport <- data.frame(matrix(ncol = length(fcsList), nrow = 1))
for(f in 1:length(fcsList)){
   
  if(fcsList[[f]]==referenceFcsFileName) {
    # pass iteration if hit the reference fcs file
    next
  }else {
    # increment ordered index and names
    orderedIndex = orderedIndex + 1
    orderedSampleList[[orderedIndex]] <- fcsList[[f]]
  
    # add new data to FlowSOM object
    fSOM <- NewData(refFSOM, GetFlowFrame(sprintf("%s/%s", fcsDirPath, fcsList[[f]]), marker_cols))
    
    # append outliers and metacluster data frame
    outliersPercentReport[f] <- GetOutLiersPercent(fSOM, refFSOM)
    metaClusterDF <- cbind(metaClusterDF, as.data.frame(table(GetMetaclusters(fSOM)))[,2])
    
    # gen pheatmap
    GenPheatmap(fSOM, fcsList[[f]], phmDirPath)
  }
  
 
}

# clean and export some dataframes
colnames(metaClusterDF) <- orderedSampleList
DfToPdf(metaClusterDF, 'metaClusterDF', extraDirPath)

colnames(outliersPercentReport) <- orderedSampleList
DfToPdf(outliersPercentReport, 'outilersReport', extraDirPath)


