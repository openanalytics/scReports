#' Read single cell data contained in 10X H5 or H5AD files
#'
#' @inheritParams commonParams
#'
#' @importFrom checkmate assertCharacter assertTRUE assertLogical 
#'  assertFileExists
#' @return [\code{list}] \cr
#'  A named list of sparse matrices of class dgCMatrix containing the 
#'  gene expression data and/or antibody data of the samples 
#'  supplied to \code{samplePaths}. 
#' @export
#' @author Trishanta Padayachee
#' @examples
#' \dontrun{
#' names <- c("pbmc_1k_totalseqB_SI_h5ad", "pbmc_1k_totalseqB_SI_h5")
#' paths <- c(
#'  "./testdata/pbmc_1k_totalseqB_SI/pbmc_1k_totalseqB_SI_raw_feature_bc_matrix.h5ad",
#'  "./testdata/pbmc_1k_totalseqB_SI/pbmc_1k_totalseqB_SI_raw_feature_bc_matrix.h5"
#' )
#' data <- readSingleCellData(
#'     sampleNames = names,
#'     samplePaths = paths, 
#'     verbose = TRUE
#' )
#' }
readSingleCellData <- function(sampleNames, samplePaths, 
                               verbose = FALSE,
                               makeFeaturesUnique = TRUE){
  
  # Check parameters
  if (verbose) message("Checking parameters...")
  checkmate::assertCharacter(x = sampleNames,
                            any.missing = FALSE,
                            unique = TRUE)
  checkmate::assertCharacter(x = samplePaths,
                            any.missing = FALSE,
                            unique = TRUE)
  checkmate::assertTRUE(length(x = sampleNames) == length(samplePaths))
  checkmate::assertLogical(x = verbose,
                            any.missing = FALSE,
                            len = 1)
  
  if (verbose) {
    
    missingFiles <- samplePaths[!sapply(samplePaths, file.exists)]
    
    if (length(missingFiles) != 0) {

      stop(
        "One of more files not found. Please check the following paths: \n",
        paste(missingFiles, collapse = ",\n")
      )
      
    }
    
  }
  checkmate::assertFileExists(x = samplePaths,
                              extension = c("h5", "h5ad"))
  
  ## Load data
  if (verbose) message("Reading in data...")
  scdata <- lapply(
    X = samplePaths, 
    FUN = readData, 
    verbose = verbose, 
    makeFeaturesUnique = makeFeaturesUnique
  )
  names(scdata) <- sampleNames
  
  return(scdata)
}



#' Read single cell data contained in 10X H5 or H5AD files
#'
#' @inheritParams commonParams
#'
#' @return [\code{list}] \cr
#'  A list of one or more dgCMatrix (sparse, column-oriented,
#'  numeric matrix) objects containing the gene expression data and/or antibody 
#'  data of a particular sample.
#' @author Trishanta Padayachee
readData <-  function(path, 
                      verbose = FALSE, 
                      makeFeaturesUnique = TRUE) {
  
  splitPath <- strsplit(basename(path), split="[.]")[[1]]
  ext <- splitPath[length(splitPath)]
  
  if (ext == "h5") {
    
    data <- read10xH5(path = path, 
                       verbose = verbose, 
                       makeFeaturesUnique = makeFeaturesUnique)
    
  } else if (ext == "h5ad") {
    
    data <- readH5AD(path = path, 
                      verbose = verbose, 
                      makeFeaturesUnique = makeFeaturesUnique)
    
  }
  
  return(data)
}



#' Read 10X h5 file 
#' 
#' Reads 10x Genomics h5 data files.
#'
#' @inheritParams commonParams
#' @return [\code{list}] \cr
#'  A list of one or more dgCMatrix (sparse, column-oriented,
#'  numeric matrix) objects containing gene expression data and/or antibody 
#'  data.
#' @author Trishanta Padayachee
read10xH5 <- function(path, 
                       verbose = FALSE, 
                       makeFeaturesUnique = TRUE){
  
  message("Reading ", path)
  h5 <- hdf5r::h5file(filename = path, mode = 'r')   
  
  on.exit(expr = {
    
    if (verbose) message("Closing H5 file...")
    if ("h5" %in% ls()) h5$close_all()
    
  })
  
  groupName <- hdf5r::list.groups(h5)[1]
  
  if (hdf5r::existsGroup(h5, 'matrix')) {
    
    features <- h5[[paste0(groupName, "/features/name" )]][]
    featureTypes <- h5[[paste0(groupName, "/features/feature_type")]][]
    
  } else {
    
    features <- h5[[paste0(groupName, "/gene_names")]][]
    featureTypes <- rep(x = "Gene Expression", 
                        times = h5[[paste0(groupName, "/shape")]][][1])
  }
  
  indices <- h5[[paste0(groupName, "/indices")]][]
  indptr <- h5[[paste0(groupName, "/indptr")]][]
  data <- h5[[paste0(groupName, "/data")]][]
  shape <- h5[[paste0(groupName, "/shape")]][]
  barcodes <- h5[[paste0(groupName, "/barcodes")]][]

  X <- Matrix::sparseMatrix(
    i = indices,
    p = indptr,
    x = data,
    dims = shape,
    dimnames = list(features, barcodes),
    index1 = FALSE
  )
  
  data <- list()
  for (type in unique(featureTypes)) {
    
    typeShort <- switch(EXPR = type,
                        `Gene Expression` = "RNA",
                        `Antibody Capture` = "ADT",
                        `CRISPR Guide Capture` = "CRISPR",
                        `Custom` = "Custom")
    
    if (verbose) message("Extracting ", type, " data...")
    data[[typeShort]] <- X[featureTypes == type, ]
    
    if (makeFeaturesUnique) {
      
      rownames(data[[typeShort]]) <- 
        make.unique(rownames(data[[typeShort]]))
      
    }
    
  }
  
  return(data)
}



#' Read h5ad file
#' 
#' \code{readH5AD} assumes that Gene Expression data is stored in the \code{X}
#'  slot of the anndata object and that if present, Antibody Capture data is 
#'  stored in the \code{obsm} slot of the anndata object under the name 
#'  \code{counts_antibody}.
#' 
#' @inheritParams commonParams
#' 
#' @return [\code{list}] \cr
#'  A list of one or more sparse matrices of class dgCMatrix containing 
#'  gene expression data and/or antibody data.
#' @importFrom reticulate import
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as
#' @author Trishanta Padayachee
readH5AD <- function(path, 
                      verbose = FALSE, 
                      makeFeaturesUnique = TRUE){
  
  ## Load python anndata module using reticulate
  anndata <- reticulate::import("anndata")
  
  ## Read data
  message("Reading ", path)
  
  aData <- anndata$read_h5ad(filename = path)
  dimNames <- list(aData$var_names$to_list(), 
                   aData$obs_names$to_list()) # since aData$X@Dimnames is NULL
  
  ## Transpose aData$X for GEX and convert from dgRMatrix to dgCMatrix
  data <- list()
  
  if (verbose) message("Extracting Gene Expression data...")
  
  if (inherits(aData$X, "dgRMatrix")) {
    
    data[["RNA"]] <- Matrix::sparseMatrix(
      i = aData$X@j,
      p = aData$X@p,
      x = aData$X@x,
      dims = rev(aData$X@Dim),
      dimnames = dimNames,
      index1 = FALSE
    )
    
  } else if (inherits(aData$X, "matrix")) {
    
    data[["RNA"]] <- Matrix::Matrix(
      data = t(aData$X),
      dimnames = dimNames,
      sparse = TRUE
    )
    
    if (!inherits(data[["RNA"]], "dgCMatrix")) {
      
      data[["RNA"]] <- methods::as(
        object = data[["RNA"]], 
        Class = "dgCMatrix"
      )
      
    }
    
  } else if (is.null(aData$X)) {
    
    data[["RNA"]] <- NULL
    
  } else {
    
    stop("Object is not of class `dgRMatrix` or `matrix`. ",
         "Please contact the package maintainer.")
    
  }
  
  if (makeFeaturesUnique & !is.null(data[["RNA"]])) {
    
    rownames(data[["RNA"]]) <- make.unique(
      names = as.character(rownames(data[["RNA"]]))
    )
    
  }
  
  ## Extract antibody data if available
  if (!is.null(aData$obsm) & "counts_antibody" %in%  aData$obsm_keys()) {
    
    if (inherits(aData$obsm[["counts_antibody"]], "data.frame") &
        (ncol(aData$obsm[["counts_antibody"]]) != 0)) {
      
      if (verbose) message("Extracting Antibody Capture data...")
      
      data[["ADT"]] <- Matrix::Matrix(
        data = t(as.matrix(aData$obsm[["counts_antibody"]])),
        sparse = TRUE
      )
      
      if (!inherits(data[["ADT"]], "dgCMatrix")) {
        
        data[["ADT"]] <- methods::as(
          object = data[["ADT"]], 
          Class = "dgCMatrix"
        )
        
      }
      
      if (makeFeaturesUnique) {
        
        rownames(data[["ADT"]]) <- make.unique(
          names = as.character(rownames(data[["ADT"]]))
        )
        
      }
      
    } else if (!inherits(aData$obsm[["counts_antibody"]], "data.frame")) {
      
      stop("Object is not of class `data.frame`. ", 
           "Please contact the package maintainer.")
      
    }
    
  }
  
  return(data)
  
}
