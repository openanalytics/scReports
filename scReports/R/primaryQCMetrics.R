#' Extract metrics from an h5ad file
#' 
#' @param nFeatureRNACol [\code{character(1)}] \cr
#'  Name of the column in the obs slot of the AnnData object that corresponds 
#'  to the number of genes expressed per cell.
#' @param nFeatureADTCol [\code{character(1)}] \cr 
#'  Name of the column in the obs slot of the AnnData object that corresponds 
#'  to the number of antibodies per cell.
#' @param nCountRNACol [\code{character(1)}] \cr 
#'  Name of the column in the obs slot of the AnnData object that corresponds 
#'  to the total number of UMI counts per cell.
#' @param nCountADTCol [\code{character(1)}] \cr 
#'  Name of the column in the obs slot of the AnnData object that corresponds 
#'  to the total number of ADT counts per cell.
#' @param percentMitoCol [\code{character(1)}] \cr 
#'  Name of the column in the obs slot of the AnnData object that corresponds 
#'  to the percentage of mitochondrial UMI counts per cell.
#' @param nCellRNACol [\code{character(1)}] \cr 
#'  Name of the column in the var slot of the AnnData object that corresponds 
#'  to the number of cells expressing each gene.
#' @param nCellADTCol [\code{character(1)}] \cr 
#'  Name of the column in the var slot of the AnnData object that corresponds 
#'  to the number of cells with each antibody.
#' @inheritParams commonParams
#' 
#' @return [\code{list}] \cr
#'  A list of one or more sparse matrices of class dgCMatrix containing 
#'  gene expression data and/or antibody data.
#' @importFrom reticulate import
#' @importFrom Matrix sparseMatrix
#' @importFrom dplyr rename
#' @author Trishanta Padayachee
extractMetrics <- function(path,
                           sampleName, 
                           verbose = FALSE, 
                           makeFeaturesUnique = TRUE,
                           nFeatureRNACol = "n_genes",
                           nFeatureADTCol = NULL,
                           nCountRNACol = "n_counts",
                           nCountADTCol = NULL,
                           percentMitoCol = "percent_mito",
                           nCellRNACol = "n_cells",
                           nCellADTCol = NULL){
  
  # load python anndata module using reticulate
  anndata <- reticulate::import("anndata")
  
  # read data
  if (verbose) message("Reading ", path)
  
  aData <- anndata$read_h5ad(filename = path)
  
  # extract metrics
  if (verbose) message("Extracting primary QC metrics...")
  
  metricObsNames <- c(nFeatureRNACol, 
                      nFeatureADTCol, 
                      nCountRNACol, 
                      nCountADTCol, 
                      percentMitoCol)
  metricVarRNANames <- c(nCellRNACol)
  metricVarADTNames <- c(nCellADTCol)
  
  ## check if metrics are present in the anndata object
  if (!all(metricObsNames %in% colnames(aData$obs))) {
    
    missing <- metricObsNames[!(metricObsNames %in% colnames(aData$obs))]
    if (length(missing) == length(metricObsNames) & !is.null(aData$X)) {
      stop("The following metrics are not present in the AnnData object: ",
           paste0(missing, sep = " "))
    } else {
      warning("The following metrics are not present in the AnnData object: ",
              paste0(missing, sep = " "))
      metricObsNames <- metricObsNames[metricObsNames %in% colnames(aData$obs)]
    }
    
  }
  
  missing <- metricVarRNANames[!(metricVarRNANames %in% colnames(aData$var))]
  if (!all(metricVarRNANames %in% colnames(aData$var)) & !is.null(aData$X)) {
    
    stop("The following metrics are not present in the AnnData object: ",
         paste0(missing, sep = ", "))
    
  } else {
    
    warning("The following metrics are not present in the AnnData object: ",
         paste0(missing, sep = ", "))
    metricVarRNANames <- metricVarRNANames[metricVarRNANames %in% colnames(aData$var)]
    
  }
  
  missing <- metricVarADTNames[!(metricVarADTNames %in% colnames(aData$var))]
  if (!all(metricVarADTNames %in% colnames(aData$var))) {
    
    stop("The following metrics are not present in the AnnData object: ",
         paste0(missing, sep = ", "))
    
  } else {
    
    warning("The following metrics are not present in the AnnData object: ",
            paste0(missing, sep = ", "))
    metricVarADTNames <- metricVarADTNames[metricVarADTNames %in% colnames(aData$var)]
    
  }
  
  ## extract and rename cell-level metrics
  nameKey <- c("nFeature_RNA" = nFeatureRNACol, 
               "nFeature_ADT" = nFeatureADTCol, 
               "nCount_RNA" = nCountRNACol, 
               "nCount_ADT" = nCountADTCol, 
               "percentMito" = percentMitoCol)
  
  qcMetrics <- NULL
  qcMetrics[["cells"]] <- aData$obs[, metricObsNames, drop = FALSE]
  qcMetrics[["cells"]] <- dplyr::rename(
    .data = qcMetrics[["cells"]],
    dplyr::any_of(nameKey)
  )
  if (nrow(qcMetrics[["cells"]]) != 0) {
    
    qcMetrics[["cells"]][,"sampleName"] <-  sampleName
    
  } else {
    
    ### include row of zeros if there are no cells
    qcMetrics[["cells"]][1, ] <- 0
    qcMetrics[["cells"]][, "sampleName"] <- sampleName
    
  }
  
  ## extract and rename feature-level RNA metrics
  nameKey <- c("nCells_RNA" = nCellRNACol)
  qcMetrics[["features"]][["RNA"]] <- aData$var[, metricVarRNANames, drop = FALSE]
  qcMetrics[["features"]][["RNA"]] <- dplyr::rename(
    .data = qcMetrics[["features"]][["RNA"]],
    dplyr::any_of(nameKey)
  )
  if (nrow(qcMetrics[["features"]][["RNA"]]) == 0) {
    
    qcMetrics[["features"]][["RNA"]][1, ] <- 0
    
  }
  
  if (makeFeaturesUnique) {
    
    rownames(qcMetrics[["features"]][["RNA"]]) <- 
      make.unique(rownames(qcMetrics[["features"]][["RNA"]]))
    
  }
  
  ## extract and rename feature-level ADT metrics
  if (!is.null(metricVarADTNames)) {
    
    nameKey <- c("nCells_ADT" = nCellADTCol)
    qcMetrics[["features"]][["ADT"]] <- aData$var[, metricVarADTNames, drop = FALSE]
    qcMetrics[["features"]][["ADT"]] <- dplyr::rename(
      qcMetrics[["features"]][["ADT"]],
      dplyr::any_of(nameKey)
    )
    if (nrow(qcMetrics[["features"]][["ADT"]]) == 0) {
      
      ### include row of NAs if there are no cells
      qcMetrics[["features"]][["ADT"]][1, ] <- 0
      
    }
    
    if (makeFeaturesUnique) {
      
      rownames(qcMetrics[["features"]][["ADT"]]) <- 
        make.unique(rownames(qcMetrics[["features"]][["ADT"]]))
      
    }
    
  }
  
  
  return(qcMetrics)
  
}


#' Compute the primary quality control metrics
#' 
#' Computes the primary quality control (QC) metrics for a dgCMatrix object
#'  or list of dgCMatrix objects for a particular sample. The following QC 
#'  metrics are computed:  
#'  
#'  - the total number of counts per cell (nCounts_RNA and/or nCounts_ADT),  
#'  
#'  - the number of features expressed per cell (nFeature_RNA and/or 
#'  nFeature_ADT), and
#'  
#'  - the percentage of counts that can be attributed to mitochondrial genes
#'  (percent.mito) 
#'
#' @param X [\code{list}] \cr
#'  A named list of dgCMatrix objects (features by cells) of raw counts.
#' @param replaceNaN [\code{logical(1)}, default: TRUE] \cr
#'  Should the \code{NaN} values in the percentMito be replaced with \code{0}?
#' @inheritParams commonParams 
#'
#' @return [\code{data.frame}] \cr
#'  A data frame of quality control metrics. 
#' @export
#' @author Trishanta Padayachee
computeMetrics <- function(X, 
                           sampleName,
                           verbose = FALSE,
                           replaceNaN = TRUE) {
  
  # check parameters
  if (verbose) message("Checking parameters...")
  
  checkmate::assertList(
    x = X,
    any.missing = FALSE, 
    names = "unique"
  )
  X <- lapply(
    X = X, 
    FUN = checkmate::assertClass, 
    classes = "dgCMatrix"
  )
  checkmate::assertCharacter(
    x = sampleName,
    len = 1,
    any.missing = FALSE
  )
  
  modes <- names(X)
  nTotalCells <- sapply(
    X = modes, 
    FUN = function(mode) ncol(X[[mode]]), 
    simplify = TRUE
  )
  
  if (any(nTotalCells > 0)) {
    
    # compute cell-level metrics
    if (verbose) message("Computing cell-level metrics...")
    
    # compute counts per cell
    nCounts <- sapply(X = X, FUN = Matrix::colSums, simplify = TRUE)
    colnames(nCounts) <- paste0("nCount_", colnames(nCounts))
    
    # features per cell
    nFeatures <- sapply(
      X = X,
      FUN = function(x) {
        Matrix::colSums(x > 0)
      }
    )
    colnames(nFeatures) <- paste0("nFeature_", colnames(nFeatures))
    
    cellMetrics <- data.frame(sampleName, nCounts, nFeatures)
    
    # percentage of mitochondrial counts per cell
    if ("RNA" %in% modes) {
      
      mitoGenes <- grep("^(MT-|mt-)", rownames(X[["RNA"]]))
      mitoSubset <- X[["RNA"]][mitoGenes, , drop = FALSE]
      percMito <- Matrix::colSums(mitoSubset)/nCounts[, "nCount_RNA"] * 100
      
      if (replaceNaN) percMito[is.nan(percMito)] <- 0 
      
      percMito <- data.frame("percentMito" = percMito)
      cellMetrics <- data.frame(cellMetrics, percMito)
      
    }
    
    ## cells per gene
    if (verbose) message("Computing feature-level metrics...")
    
    nCells <- sapply(
      X = X,
      FUN = function(x) {
        Matrix::rowSums(x > 0)
      },
      simplify = FALSE
    )
    
    featureMetrics <- NULL
    
    for (mode in modes) {
      
      featureMetrics[[mode]] <- data.frame(nCells[[mode]])
      colnames(featureMetrics[[mode]]) <- paste0("nCells_", mode)
      
    }
    
  } else {
    
    if (verbose) message("Matrix contains zero cells. ",
                         "Including a fake cell for plotting...")
    
    cellMetrics <- data.frame(
      "sampleName" = sampleName, 
      "nCount_RNA" = 0, 
      "nFeature_RNA" = 0
    )
    
    featureMetrics <- list(
      "RNA" = data.frame("nCells_RNA" = 0)
    )
      
    
  }
  
  return(list("cells" = cellMetrics, "features" = featureMetrics))
  
}


#' Combine a list of primary QC metrics into a single data frame
#'
#' @param X [\code{list}] \cr
#'  Data frames of primary QC metrics.
#' @param fillMissing [\code{logical(1)}, default: \code{TRUE}] \cr
#'  For data frames with missing metric column(s) (e.g., due to the absence of 
#'  an assay), should a single zero observation be included in the (missing) 
#'  metric columns, so that the sample's name still appears in figures?
#'
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' 
#' @return [\code{data.frame}] \cr
#'  A data.frame of primary QC metrics.
#' @export
#' @author Trishanta Padayachee
combineMetrics <- function(X, 
                           fillMissing = TRUE){
  
  metricNames <- unique(
    unlist(
      sapply(X = X, FUN = colnames, simplify = FALSE)
    )
  )
  
  if (fillMissing) {
    
    X <- lapply(
      X = seq_along(X), 
      FUN = function(i) {
        extendMetrics(X = X[[i]],
                      sampleName = names(X)[i],
                      metricNames = metricNames)
      }
    )
    
  }
  
  metrics <- X %>%
    dplyr::bind_rows()
  
  
  return(metrics)
  
}


#' Fill in missing columns 
#'
#' @param X [\code{data.frame}] \cr
#'  A data.frame of primary QC metrics.
#' @param sampleName [\code{character(1)}] \cr
#'  Name of sample associated with \code{X}.
#' @param metricNames [\code{character}] \cr
#'  A vector of the expected metrics. 
#' @return [\code{data.frame}] \cr
#'  A data.frame of primary QC metrics containing the columns listed in 
#'  \code{metricNames}.
#' @export
#' @author Trishanta Padayachee
extendMetrics <- function(X, sampleName, metricNames) {
  
  missingCols <- metricNames[!(metricNames %in% colnames(X))]
  
  if (!is.null(missingCols) & nrow(X) != 0) {
    
    # add missing columns (i.e., metrics)
    for (col in missingCols) {
      
      X[, col] <- c(0, rep(NA, nrow(X) - 1))
      
    }
    
  } else if (nrow(X) == 0) {
    
    # in case of zero cells after filtering
    countCols <- metricNames[grep("(^n|^percent)", metricNames)]
    logicalCols <- metricNames[grep("(^isSelected)", metricNames)]
    X[1, ] <- 0
    X[1, "sampleName"] <- sampleName
    X[1, countCols] <- 0
    X[1, logicalCols] <- NA
    
  }
  
  return(X)
  
}