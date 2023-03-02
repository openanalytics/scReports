#' Filter features and/or cells
#' 
#' Filters cells based on user specified primary quality control metric cutoffs.
#'
#' @param X [\code{list}] \cr
#'  A named list of dgCMatrix objects (features by cells) of raw counts 
#'  for a single sample.
#' @param qcMetrics [\code{list}] \cr
#'  A named list of primary quality control metrics. 
#' @param applyMultiModalCellFiltering 
#'  [\code{logical(1)}, default: \code{FALSE}] \cr
#'  Whether to apply multi-modal filtering. If \code{TRUE}, cells with 
#'  nCounts_RNA and nCounts_ADT less than \code{cellFilteringMinCounts} and 
#'  greater than \code{cellFilteringMaxCounts} are filtered out. 
#' @inheritParams commonParams 
#' 
#' @return [\code{list}] \cr
#'  A named list containing the filtered dgCMatrix objects, an updated 
#'  data frame of the raw QC metrics with filtering additional columns and a 
#'  filtered QC metrics data frame. 
#' @export
#' @author Trishanta Padayachee
filterSCData <- function(X, qcMetrics, sampleName,
                         applyMultiModalCellFiltering = FALSE,
                         cellFilteringMinGenes = 0,
                         cellFilteringMaxGenes = Inf,
                         cellFilteringMinCounts = 10,
                         cellFilteringMaxCounts = Inf,
                         cellFilteringMaxPercentMito = 100,
                         geneFilteringMinCells = 3,
                         antibodyFilteringMinCells = 0, 
                         cellFilteringMinCountsRNA = 10,
                         cellFilteringMaxCountsRNA = Inf, 
                         cellFilteringMinCountsADT = 10,
                         cellFilteringMaxCountsADT = Inf, 
                         verbose = FALSE) {
  
  message("> Filtering ", sampleName, " data...")
  
  # check parameters
  if (verbose) message("Checking parameters...")
  
  checkmate::assertList(X, 
                        names = "unique")
  checkmate::assertNames(names(X), 
                         type = "unique",
                         subset.of = c("RNA", "ADT", "CRISPR", "CUSTOM"))
  checkmate::assertList(qcMetrics, 
                        names = "unique")
  checkmate::assertNames(names(qcMetrics), 
                         type = "unique",
                         subset.of = c("cells", "features"))
  checkmate::assertLogical(applyMultiModalCellFiltering,
                           any.missing = FALSE,
                           len = 1)
  checkmate::assertNumber(cellFilteringMinGenes, 
                          lower = 0)
  checkmate::assertNumber(cellFilteringMaxGenes, 
                          lower = 0)
  checkmate::assertNumber(cellFilteringMinCounts, 
                          lower = 0)
  checkmate::assertNumber(cellFilteringMaxCounts, 
                          lower = 0)
  checkmate::assertNumber(cellFilteringMaxPercentMito, 
                          lower = 0,
                          upper = 100)
  checkmate::assertNumber(geneFilteringMinCells, 
                          lower = 0)
  checkmate::assertNumber(antibodyFilteringMinCells, 
                          lower = 0)
  checkmate::assertNumber(cellFilteringMinCountsRNA, 
                          lower = 0)
  checkmate::assertNumber(cellFilteringMaxCountsRNA, 
                          lower = 0)
  checkmate::assertNumber(cellFilteringMinCountsADT, 
                          lower = 0)
  checkmate::assertNumber(cellFilteringMaxCountsADT, 
                          lower = 0)
  checkmate::assertLogical(verbose,
                           any.missing = FALSE,
                           len = 1)
  
  # initialize
  modes <- names(X)
  
  nTotalCells <- sapply(
    X = modes, 
    FUN = function(mode) ncol(X[[mode]]), 
    simplify = TRUE
  )
  
  cellMetrics <- qcMetrics[["cells"]]
  featureMetrics <- qcMetrics[["features"]]
  
  requiredMetricsRNA <- c("nCount_RNA", "nFeature_RNA")
  requiredMetricsADT <- c("nCount_ADT", "nFeature_ADT")
  
  # set filtering logicals
  
  filterRNA <- FALSE
  filterADT <- FALSE
  includesCountsRNA <- FALSE
  includesCountsADT <- FALSE
  includesMetricsRNA <- FALSE
  includesMetricsADT <- FALSE
  
  if ("RNA" %in% modes) {
    
    if (nTotalCells["RNA"] == 0) {
      
      message("RNA matrix contains zero cells. No RNA filtering applied.")
      
    }
    
    includesCountsRNA <- ("RNA" %in% modes) & (nTotalCells["RNA"] != 0)
    includesMetricsRNA <- 
      all(requiredMetricsRNA %in% colnames(cellMetrics))
      "nCells_RNA" %in% colnames(featureMetrics[["RNA"]])
    
    filterRNA <- includesCountsRNA & includesMetricsRNA
    
  }
  
  if ("ADT" %in% modes) {
    
    if (nTotalCells["ADT"] == 0) {
      
      message("RNA matrix contains zero cells. No ADT filtering applied.")
      
    }
    
    includesCountsADT <- "ADT" %in% names(X) & (nTotalCells["ADT"] != 0)
    includesMetricsADT <- 
      all(requiredMetricsADT %in% colnames(cellMetrics)) &
      "nCells_ADT" %in% colnames(featureMetrics[["ADT"]])
    
    filterADT <- includesCountsADT & includesMetricsADT
    
  }
  
  if (filterRNA) {
    
    nGenesPreFiltering <- nrow(X[["RNA"]])
    nCellsPreFiltering <- ncol(X[["RNA"]])
    
  }
  if (filterADT) {
    
    nADTPreFiltering <- nrow(X[["ADT"]])
    nCellsPreFiltering <- ncol(X[["ADT"]])
    
  }
  
  
  # warnings/messages based on filtering logicals
  if (filterRNA | filterADT) {
    
    if (applyMultiModalCellFiltering) {
      
      if (!(includesCountsRNA & includesCountsRNA))  {
        
        if (verbose) {
          
          warning(
            strwrap(
              "X should contain 'RNA' and 'ADT' data when the argument
          applyMultiModalCellFiltering is set to TRUE. Setting
          applyMultiModalCellFiltering to FALSE.",
              prefix = " ", initial = ""
            )
          )
          
        }
        applyMultiModalCellFiltering <- FALSE
        
      } else if (!(includesMetricsRNA & includesMetricsADT)) {
        
        if (verbose) {
          
          warning(
            strwrap(
              "qcMetrics should have the columns 
          'nCount_RNA', 'nCells_RNA', 'nCount_ADT' and 'nCells_ADT' when  
          applyMultiModalCellFiltering is TRUE.
          Setting applyMultiModalCellFiltering to FALSE.",
              prefix = " ", initial = ""
            )
          )
          
        }
        applyMultiModalCellFiltering <- FALSE
        
      }
      
    }
    
  }
  
  if (includesCountsRNA & !includesMetricsRNA) {
    
    warning(
      strwrap(
        "RNA counts are available. 
        However, the required RNA metrics are missing.
        No filtering is performed based on RNA primary QC metrics.",
        prefix = " ", initial = ""
      )
    )
    
  }
  
  if (includesCountsADT & !includesMetricsADT) {
    
    warning(
      strwrap(
        "ADT counts are available. 
        However, the required ADT metrics are missing.
        No filtering is performed based on ADT primary QC metrics.",
        prefix = " ", initial = ""
      )
    )
    
  }
  
  XFiltered <- X
  qcMetricsFiltered <- NULL
  qcMetricsFiltered[["cells"]] <- cellMetrics
  qcMetricsFiltered[["features"]] <- featureMetrics
  
  if (filterRNA | filterADT) {
    
    ## logicals for filtering
    if (verbose) message("Computing filtering variables...")
    
    applyNGeneCellFiltering <- FALSE
    applyNCountsMultiModalCellFiltering <- FALSE
    applyNCountsRNACellFiltering <- FALSE
    applyNCountsADTCellFiltering <- FALSE
    applyPercentMitoCellFiltering <- FALSE
    applyNCellsAntibodyFiltering <- FALSE
    
    if (filterRNA) {
      
      applyNGeneCellFiltering <- 
        cellFilteringMinGenes != 0 | cellFilteringMaxGenes != Inf
      
      applyNCountsRNACellFiltering <- 
        cellFilteringMinCountsRNA != 0 | cellFilteringMaxCountsRNA != Inf
      
      applyPercentMitoCellFiltering <-
        "percentMito" %in% colnames(cellMetrics) & 
        cellFilteringMaxPercentMito != 100 
      
      applyNCellsRNAGeneFiltering <- geneFilteringMinCells != 0
      
      featureMetrics[["RNA"]][, "isSelectedNCellsRNA"] <- TRUE
      
    }
    
    if (filterADT) {
      
      applyNCountsADTCellFiltering <- 
        cellFilteringMinCountsADT != 0 | cellFilteringMaxCountsADT != Inf
      
      applyNCellsAntibodyFiltering <- antibodyFilteringMinCells != 0
      
      featureMetrics[["ADT"]][, "isSelectedNCellsADT"] <- TRUE
      
    }
    
    if (filterRNA & filterADT) {
      
      if (applyMultiModalCellFiltering) {
        
        applyNCountsRNACellFiltering <- FALSE
        applyNCountsADTCellFiltering <- FALSE
        applyNCountsMultiModalCellFiltering <- 
          cellFilteringMinCounts != 0 | cellFilteringMaxCounts != Inf
        
      }
      
    }
    
    cellMetrics[, "isSelected"] <- TRUE  
    cellMetrics[, "isSelectedNGenes"] <- TRUE
    cellMetrics[, "isSelectedNCountsRNA"] <- TRUE
    cellMetrics[, "isSelectedNCountsADT"] <- TRUE
    cellMetrics[, "isSelectedPercMito"] <- TRUE
    
    if (applyNGeneCellFiltering) {
      
      cellMetrics[, "isSelectedNGenes"] <- 
        cellMetrics[, "nFeature_RNA"] >= cellFilteringMinGenes &
        cellMetrics[, "nFeature_RNA"] <= cellFilteringMaxGenes
      
    }
    
    if (applyNCountsRNACellFiltering) {
      
      cellMetrics[, "isSelectedNCountsRNA"] <- 
        cellMetrics[, "nCount_RNA"] >= cellFilteringMinCountsRNA &
        cellMetrics[, "nCount_RNA"] <= cellFilteringMaxCountsRNA
      
    }
    
    if (applyPercentMitoCellFiltering) {
      
      cellMetrics[, "isSelectedPercMito"] <- 
        cellMetrics[, "percentMito"] <= cellFilteringMaxPercentMito
      
    }
    
    if (applyNCountsADTCellFiltering) {
      
      cellMetrics[, "isSelectedNCountsADT"] <- 
        cellMetrics[, "nCount_ADT"] >= cellFilteringMinCountsADT &
        cellMetrics[, "nCount_ADT"] <= cellFilteringMaxCountsADT
      
    }
    
    if (applyNCountsMultiModalCellFiltering) {
      
      cellMetrics[, "isSelectedNCountsRNA"] <-
        cellMetrics[, "isSelectedNCountsADT"] <-
        (cellMetrics[, "nCount_RNA"] >= cellFilteringMinCounts &
           cellMetrics[, "nCount_RNA"] <= cellFilteringMaxCounts) |
        (cellMetrics[, "nCount_ADT"] >= cellFilteringMinCounts &
           cellMetrics[, "nCount_ADT"] <= cellFilteringMaxCounts)
      
    }
    
    cellMetrics[, "isSelected"] <- 
      cellMetrics[, "isSelectedNGenes"] & 
      cellMetrics[, "isSelectedNCountsRNA"] & 
      cellMetrics[, "isSelectedNCountsADT"] & 
      cellMetrics[, "isSelectedPercMito"]
    
    if (applyNCellsRNAGeneFiltering) {
      
      featureMetrics[["RNA"]][, "isSelectedNCellsRNA"] <- 
        featureMetrics[["RNA"]][, "nCells_RNA"] >= geneFilteringMinCells 
      
    }
    
    if (applyNCellsAntibodyFiltering) {
      
      featureMetrics[["ADT"]][, "isSelectedNCellsADT"] <- 
        featureMetrics[["ADT"]][, "nCells_ADT"] >= antibodyFilteringMinCells 
      
    }
    
    # filter RNA and ADT count data
    if (verbose) message("Filtering count matrices...")
    
    if (filterRNA) {
      
      XFiltered[["RNA"]] <- X[["RNA"]][
        featureMetrics[["RNA"]][, "isSelectedNCellsRNA"], 
        cellMetrics[, "isSelected"]
      ]
      
    }
    
    if (filterADT) {
      
      XFiltered[["ADT"]] <- X[["ADT"]][
        featureMetrics[["ADT"]][, "isSelectedNCellsADT"],
        cellMetrics[, "isSelected"]
      ]
      
    }
    
    # filter cell and feature metrics
    if (verbose) message("Filtering metrics...")
    
    if (filterRNA | filterADT) {
      
      qcMetricsFiltered[["cells"]] <- cellMetrics[cellMetrics[, "isSelected"], ]
      
    }
    if (filterRNA) {
      
      qcMetricsFiltered[["features"]][["RNA"]] <- 
        featureMetrics[["RNA"]][
          featureMetrics[["RNA"]][, "isSelectedNCellsRNA"], 
        ]
      
    }
    if (filterADT) {
      
      qcMetricsFiltered[["features"]][["ADT"]] <- 
        featureMetrics[["ADT"]][
          featureMetrics[["ADT"]][, "isSelectedNCellsADT"], 
        ]
      
    }
    
    # return filtering stats
    if (filterRNA) {
      
      nCellsPostFiltering <- sum(cellMetrics[, "isSelected"])
      
      nGenesPostFiltering <- sum(featureMetrics[["RNA"]][, "isSelectedNCellsRNA"])
      
      nCellsOutsideNGeneInterval <- sum(!cellMetrics[, "isSelectedNGenes"])
      nCellsOutsideNCountRNAInterval <- 
        sum(!cellMetrics[, "isSelectedNCountsRNA"])
      nCellsOutsideNCountADTInterval <- 
        sum(!cellMetrics[, "isSelectedNCountsADT"])
      nCellsAboveMaxPercentMito <- 
        sum(!cellMetrics[, "isSelectedPercMito"])
      
    } 
    
    if (filterADT) {
      
      nADTPostFiltering <- sum(featureMetrics[["ADT"]][, "isSelectedNCellsADT"])
      
    }
    
    if (verbose) {
      
      if (filterRNA) {
        
        message(
          "\n", "Gene filtering: ",
          nGenesPreFiltering-nGenesPostFiltering,
          " gene(s) were removed due to expression in fewer than ",
          geneFilteringMinCells, 
          " cell(s)"
        )
        message("Pre filtering: ", nGenesPreFiltering, " gene(s)" )
        message("Post filtering: ", nGenesPostFiltering, " gene(s)" )
        
      }
      
      if (filterADT) {
        
        message(
          "\n", "Antibody filtering: ",
          nADTPreFiltering-nADTPostFiltering,
          " antibodie(s) were removed due to observation in fewer than ", 
          antibodyFilteringMinCells, 
          " cell(s)"
        )
        message("Pre filtering: ", nADTPreFiltering, " antibodie(s)")
        message("Post filtering: ", nADTPostFiltering, " antibodie(s)")
        
      }
      
      message(
        "\n", "Cell filtering: ",
        nCellsPreFiltering-nCellsPostFiltering,
        " cell(s) were removed"
      )
      message("Pre filtering: ", nCellsPreFiltering, " cell(s)")
      message("Post filtering: ", nCellsPostFiltering, " cell(s)")
      
      if (applyMultiModalCellFiltering) {
        
        message(
          "\n", nCellsOutsideNCountRNAInterval,
          " cell(s) have an nCount_RNA and nCount_ADT value less than ", 
          cellFilteringMinCounts, 
          " or greater than ", cellFilteringMaxCounts, "."
        )
        
      } else {
        
        if (filterRNA) {
          
          message(
            "\n", nCellsOutsideNGeneInterval,
            " cell(s) have an nFeature_RNA value less than ", cellFilteringMinGenes, 
            " or greater than ", cellFilteringMaxGenes, "."
          )
          
          message(
            "\n", nCellsOutsideNCountRNAInterval,
            " cell(s) have an nCount_RNA value less than ", cellFilteringMinCountsRNA, 
            " or greater than ", cellFilteringMaxCountsRNA, "."
          )
          
        }
        
        if (filterADT) {
          
          message(
            "\n", nCellsOutsideNCountADTInterval,
            " cell(s) have an nCount_ADT value less than ", cellFilteringMinCountsADT, 
            " or greater than ", cellFilteringMaxCountsADT
          )
          
        }
        
      }
      
      if (filterRNA & "percentMito" %in% colnames(cellMetrics)) {
        
        message(
          "\n", nCellsAboveMaxPercentMito,
          " cell(s) have a percentMito value greater than ", 
          cellFilteringMaxPercentMito, "."
        )
        
      }
      
    }
    
  }
  
  return(
    list("XFiltered" = XFiltered,
         "qcMetricsFiltered" = qcMetricsFiltered)
  )
  
}
