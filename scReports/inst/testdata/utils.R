createExample10xH5 <- function(
    outputFile = "raw_feature_bc_matrix.h5",
    nRows = 500,
    nCols = 100,
    nNonZero = 1000,
    nMito = 100, 
    nCRISPR = 0,
    nADT = 0,
    nBinomMu = 10,
    nBinomSize = 0.5,
    groupName = "matrix",
    version = NA,
    seed = 555
)
{
  
  if (nRows < nMito + nADT + nCRISPR) 
    stop("nRows should be greater than the sum of nADT and nCRISPR")
  
  set.seed(seed = seed)
  
  message(">> Creating random sparse matrix ...")
  # zero-based dgCMatrix
  rawCounts <- Matrix::rsparsematrix(
    nrow = nRows-nADT,
    ncol = nCols, 
    nnz =  nNonZero,
    rand.x = function(n) rnbinom(n, mu = nBinomMu, size = nBinomSize) + 1,
    repr = "C"
  )
  
  if (nADT > 0) {
    
    # zero-based dgCMatrix
    rawCountsADT <- Matrix::rsparsematrix(
      nrow = nADT,
      ncol = nCols, 
      nnz =  2*nCols,
      rand.x = function(n) rnbinom(n, mu = 10, size = 0.5),
      repr = "C"
    )
    
    rawCounts <- rbind(rawCounts, rawCountsADT)
    
  }
  
  
  
  message(">> Creating test data h5 file ...")
  on.exit(
    expr = {
      message(">> Closing h5 files ...")
      if ("testdata_h5" %in% ls()) {
        testdata_h5$close_all()
      }
    }
  )
  
  if (!dir.exists(dirname(outputFile))){
    dir.create(dirname(outputFile), recursive = TRUE)
  }
  
  testdata_h5 <- H5File$new(outputFile, mode = "w")
  testdata_h5$create_group(groupName)
  
  if (is.na(version)) {
    
    testdata_h5$create_group(paste0(groupName, "/features"))
    
    times <- c(nRows-nMito-nADT-nCRISPR, nMito, nADT, nCRISPR)
    featureTypes <- rep(c("Gene Expression", "Gene Expression", 
                          "Antibody Capture", "CRISPR Guide Capture"),
                        times = times)
    featureIDs <- rep(c("GeneID", "MTGeneID", 
                        "ADTID", "CRISPRID"),
                      times = times)
    featureNames <- rep(c("Gene", "MTGene", 
                          "ADT", "CRISPR"),
                        times = times)
    
    datasets <- list(
      "indices" = rawCounts@i,
      "indptr" = rawCounts@p,
      "data" = rawCounts@x,
      "shape" = rawCounts@Dim,
      "features/_all_tag_keys" = c("genome", "read", "pattern", "sequence"),
      "features/feature_type" = featureTypes,
      "features/genome" = rep("GRCh38", times = nRows),
      "features/id" = paste0(featureIDs, 1:nRows),
      "features/name" = paste0(featureNames, 1:nRows),
      "features/read" = rep("", nRows),
      "features/pattern" = rep("", nRows),
      "features/sequence" = rep("", nRows),
      "barcodes" = paste0("CellBarcode", 1:nCols)
    )
    names(datasets) <- paste0(groupName, "/", names(datasets))
    
    attributes <- list(
      "filetype" = "matrix",
      "version" = "version 1",
      "software_version" = "version 1",
      "chemistry_description" = "version 1",
      "library_ids" = "test data A",
      "original_gem_groups" = 1
    )
    
    for (dname in names(datasets)){
      
      if (length(grep(pattern = "features", x = dname)) != 0) {
        
        dType <- H5T_STRING$new(size = 256)
        dType$set_strpad(strpad = 1)  
        
      } else  if (length(grep(pattern = "barcodes", x = dname)) != 0) {
        
        dType <- H5T_STRING$new(size = 24)
        dType$set_strpad(strpad = 1) 
        
      } else {
        
        # indices|indptr|data|shape should be H5T_STD_I64LE or H5T_STD_I32LE
        dType = NULL
        
      }
      
      testdata_h5$create_dataset(
        name = dname,
        robj = datasets[[dname]],
        dtype = dType,
        space=H5S$new("simple", dims=length(datasets[[dname]]), maxdims=length(datasets[[dname]])),
        chunk_dims = "auto"
      )
      
    }
    
    
    for (aname in names(attributes)) {
      
      testdata_h5$create_attr(
        attr_name = aname, 
        robj = attributes[[aname]]
      )
      
    }
    
  } else if (version == 2) {
    
    testdata_h5$create_group(paste0(groupName, "/features"))
    
    times <- c(nRows-nMito, nMito)
    genes <- rep(c("GeneID", "MTGeneID"), times = times)
    geneNames <- rep(c("Gene", "MTGene"), times = times)
    
    datasets <- list(
      "indices" = rawCounts@i,
      "indptr" = rawCounts@p,
      "data" = rawCounts@x,
      "shape" = rawCounts@Dim,
      "gene_names" = paste0(geneNames, 1:nRows),
      "genes" = paste0(genes, 1:nRows),
      "barcodes" = paste0("CellBarcode", 1:nCols)
    )
    names(datasets) <- paste0(groupName, "/", names(datasets))
    
    attributes <- list(
      "filetype" = "matrix",
      "VERSION" = "version 1",
      "chemistry_description" = "version 2",
      "library_ids" = "test data A",
      "original_gem_groups" = 1
    )
    
    for (dname in names(datasets)){
      
      if (length(grep(pattern = "gene", x = dname)) != 0) {
        
        dType <- H5T_STRING$new(size = 18)
        dType$set_strpad(strpad = 1)  
        
      } else  if (length(grep(pattern = "barcodes", x = dname)) != 0) {
        
        dType <- H5T_STRING$new(size = 18)
        dType$set_strpad(strpad = 1) 
        
      } else {
        
        dType = NULL
      }
      
      testdata_h5$create_dataset(
        name = dname,
        robj = datasets[[dname]],
        dtype = dType,
        space=H5S$new("simple", 
                      dims=length(datasets[[dname]]), 
                      maxdims=length(datasets[[dname]])),
        chunk_dims = "auto"
      )
      
    }
    
    
    for (aname in names(attributes)) {
      
      testdata_h5$create_attr(
        attr_name = aname, 
        robj = attributes[[aname]]
      )
      
    }
    
  }
  
  
  
  message(">> Dimension of raw h5 dataset (rows x cols): ", 
          paste0(testdata_h5[[paste0(groupName, "/shape")]][], collape = " "))
  
}


convertH5toH5AD <- function(
    filePath,
    outputFile = "raw_feature_bc_matrix.h5ad"
)
{
  
  scanpy <- import("scanpy")
  
  message(">> Reading H5 file ...")
  data <- scanpy$read_10x_h5(
    filename = filePath,
    gex_only = FALSE
  )
  
  if ("feature_types" %in% data$var_keys()) {
    
    X <- data$X[, data$var[,"feature_types"] == "Gene Expression"]
    X <- as(
      object = as.matrix(X), 
      Class = "dgRMatrix"
    )
    obs <- data$obs
    var <- data$var[data$var[, "feature_types"] == "Gene Expression", ]
    obsm <- py_none()
    
    if ("Antibody Capture" %in% unique(data$var[, "feature_types"])) {
      
      obsm <- NULL
      obsm[["counts_antibody"]] <- 
        data$to_df()[, data$var[, "feature_types"] == "Antibody Capture"]
      
      message(">> Dimension of raw ADT dataset: Obs:",
              nrow(obsm[["counts_antibody"]]),
              " Vars:", ncol(obsm[["counts_antibody"]]))
      
    }
    
    data <- scanpy$AnnData(
      X = X,
      obs = obs, 
      var = var,
      obsm = obsm
    )
    
  }

  
  message(">> Saving AnnData object ...")
  data$write_h5ad(filename = outputFile, compression = "gzip") 
  
  message(">> Dimension of raw GEX dataset: Obs:", 
          data$n_obs, " Vars:", data$n_vars)
  
}


filterH5AD <- function(
    filePath,
    outputFile = "filtered_feature_bc_matrix.h5ad",
    minCellsPerGene = 3,
    minGenesPerCell = 5,
    maxGenesPerCell = 10,
    minCountsPerCell = 30,
    maxCountsPerCell = 150,
    minFractionMito = 0,
    maxFractionMito = 1, 
    removePercentMito = FALSE
) 
{
  
  scanpy <- import("scanpy")
  numpy <- import("numpy")
  
  data <- scanpy$read_h5ad(filename = filePath)
  
  message(">> Creating the columns n_genes, n_cells and n_counts ...")
  scanpy$pp$filter_cells(data, min_genes=0)
  scanpy$pp$filter_cells(data, min_counts=0)
  scanpy$pp$filter_genes(data, min_cells=0)

  message(">> Creating the percent_mito column ...")
  mito_genes <- data$var_names$str$startswith("MT-")
  data$obs['percent_mito'] <-
    numpy$ravel(numpy$sum(data$to_df()[, mito_genes], axis=1)) / 
    numpy$ravel(numpy$sum(data$X, axis=1))
  
  print(summary(data$obs[,"n_genes"]))
  print(summary(data$obs[,"n_counts"]))
  print(summary(data$obs[,"percent_mito"]))
  
  message(">> Filtering genes based on (1) n_cells ... ")
  geneSelection <- data$var[, 'n_cells'] >= minCellsPerGene
  
  message(">> Filtering cells based on (1) n_genes, ",
          "(2) n_counts, and (3) percent_mito ...")
  cellSelection <- 
    (data$obs[, 'n_genes'] >= minGenesPerCell) &
    (data$obs[, 'n_genes'] <= maxGenesPerCell) &
    (data$obs[, 'n_counts'] >= minCountsPerCell) &
    (data$obs[, 'n_counts'] <= maxCountsPerCell) &
    (data$obs[, 'percent_mito'] <= maxFractionMito) &
    (data$obs[, 'percent_mito'] >= minFractionMito) 
  
  message(">> Filtering cells with zero counts after gene filtering ...")
  cellSelection <- cellSelection &
    (numpy$ravel(numpy$sum(data$to_df()[, geneSelection], axis=1)) > 0)
  
  keepCols <- TRUE
  if (removePercentMito) {
    keepCols <- data$obs_keys() != "percent_mito"
  }
  
  message(">> Creating filtered object ...")
  X <- data$to_df()[cellSelection, geneSelection]
  obs <- data$obs[cellSelection, keepCols]
  var <- data$var[geneSelection, ]
  obsm <- py_none()
  if ("counts_antibody" %in% data$obsm_keys()) {
    obsm <- NULL
    obsm[["counts_antibody"]] <- 
      data$obsm[["counts_antibody"]][cellSelection, ]
  }

  if (nrow(X) != 0 & ncol(X) != 0) {
    
    data <- scanpy$AnnData(
      X = X, 
      obs = obs, 
      var = var,
      obsm = obsm
    )
    data$raw <- data$copy()
    
  } else {
    
    data <- scanpy$AnnData()
    
  }
  
  message(">> Saving filtered object ...")
  data$write_h5ad(filename = outputFile, compression = "gzip") 
  
  message(">> Dimension of filtered dataset: Obs:", 
          data$n_obs, " Vars:", data$n_vars)
  
}

createExampleMetricsSummary <- function(
    names,
    values,
    outputFile = "metrics_summary.csv"
)
{
  
  if (length(names) != length(values))
    stop("`names` [",length(names), "] and `values` [", 
         length(values), "] should be of equal lengths")
  
  message(">> Generating summary metrics ...")
  data <- data.frame(
    matrix(values, nrow = 1, ncol = length(names))
  )
  colnames(data) <- names
  
  message(">> Saving file ...")
  if (!dir.exists(dirname(outputFile))){
    dir.create(dirname(outputFile), recursive = TRUE)
  }
  
  write.table(
    x = data, 
    file = outputFile, 
    sep = ",", 
    row.names = FALSE, 
    col.names = TRUE
  )
  
}


gexNames <- function() {
  
  x <- c(
    "Estimated Number of Cells",
    "Mean Reads per Cell",
    "Median Genes per Cell",
    "Number of Reads",
    "Valid Barcodes",
    "Sequencing Saturation",
    "Q30 Bases in Barcode",
    "Q30 Bases in RNA Read",
    "Q30 Bases in Sample Index",
    "Q30 Bases in UMI",
    "Reads Mapped to Genome",
    "Reads Mapped Confidently to Genome",
    "Reads Mapped Confidently to Intergenic Regions",
    "Reads Mapped Confidently to Intronic Regions",
    "Reads Mapped Confidently to Exonic Regions",
    "Reads Mapped Confidently to Transcriptome",
    "Reads Mapped Antisense to Gene",
    "Fraction Reads in Cells",
    "Total Genes Detected",
    "Median UMI Counts per Cell"
  )
  
  return(x)
  
}

adtNames <- function() {
  
  x <- c(
    "Antibody: Number of Reads",
    "Antibody: Mean Reads per Cell",
    "Antibody: Valid Barcodes",
    "Antibody: Sequencing Saturation",
    "Antibody: Q30 Bases in Barcode",
    "Antibody: Q30 Bases in Antibody Read",
    "Antibody: Q30 Bases in UMI",
    "Antibody: Q30 Bases in Sample Index",
    "Antibody: Fraction Antibody Reads",
    "Antibody: Fraction Antibody Reads Usable",
    "Antibody: Antibody Reads Usable per Cell",
    "Antibody: Fraction Reads in Barcodes with High UMI Counts",
    "Antibody: Fraction Unrecognized Antibody",
    "Antibody: Antibody Reads in Cells",
    "Antibody: Median UMIs per Cell (summed over all recognized antibody barcodes)"
  )
  
  return(x)
  
}

crisprNames <- function() {
  
  x <- c(
    "CRISPR: Number of Reads",
    "CRISPR: Mean Reads per Cell",
    "CRISPR: Valid Barcodes",
    "CRISPR: Sequencing Saturation",
    "CRISPR: Q30 Bases in Barcode",
    "CRISPR: Q30 Bases in RNA Read",
    "CRISPR: Q30 Bases in Sample Index",
    "CRISPR: Q30 Bases in UMI",
    "CRISPR: Fraction Reads with Putative Protospacer Sequence",
    "CRISPR: Fraction Guide Reads",
    "CRISPR: Fraction Guide Reads Usable",
    "CRISPR: Guide Reads Usable per Cell",
    "CRISPR: Fraction Protospacer Not Recognized",
    "CRISPR: Guide Reads in Cells",
    "CRISPR: Cells with 1 or more protospacers detected",
    "CRISPR: Cells with 2 or more protospacers detected",
    "CRISPR: Median UMIs per Cell"
  )
  
  return(x)
  
}
customNames <- function() {
  
  x <- c(
    "Custom: Number of Reads",
    "Custom: Mean Reads per Cell",
    "Custom: Valid Barcodes",
    "Custom: Sequencing Saturation",
    "Custom: Q30 Bases in Barcode",
    "Custom: Q30 Bases in Feature Read",
    "Custom: Q30 Bases in Sample Index",
    "Custom: Q30 Bases in UMI",
    "Custom: Fraction Feature Reads",
    "Custom: Fraction Feature Reads Usable",
    "Custom: Feature Reads Usable per Cell",
    "Custom: Fraction Unrecognized Feature",
    "Custom: Feature Reads in Cells",
    "Custom: Median UMIs per Cell"
  )
  
  return(x)
} 


createMetaFile <- function(nSamples, 
                           outputFile = "meta.tsv"){
  
  data <- data.frame(
    matrix(NA, nrow = nSamples, ncol = 6)
  )
  colnames(data) <- c(
    "LibraryID", "SampleID", "SubjectID", 
    "RunID", "Technology", "Species"
  )
  data[["LibraryID"]] <- paste0("LI000", 1:nSamples)
  data[["SampleID"]] <- paste0("SA000", 1:nSamples)
  data[["SubjectID"]] <- paste0("SU000", 1:nSamples)
  data[["RunID"]] <- paste0("RUN000", 1:nSamples)
  data[["Technology"]] <- "10x Chromium"
  data[["Species"]] <- "human"
  
  message(">> Saving file ...")
  if (!dir.exists(dirname(outputFile))){
    dir.create(dirname(outputFile), recursive = TRUE)
  }
  
  write.table(
    x = data, 
    file = outputFile, 
    sep = "\t", 
    row.names = FALSE, 
    col.names = TRUE
  )
  
}
