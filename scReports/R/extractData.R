#' Assess gene expression assay metrics based on Cell Ranger alerts 
#'
#' Assess the Cell Ranger gene expression metrics for various assays/samples 
#' to establish whether they meet expectations. Expectations are based on the 
#' 10xGenomics Cell Ranger alerts. 
#' 
#' @param data [\code{data.frame}] \cr
#'  A data frame of Cell Ranger gene expression metrics. Must contain a 
#'  'Sample ID' column and metric columns with names corresponding to those used 
#'  by the 10xGenomics Cell Ranger pipeline 
#'  (see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/gex-metrics).
#'
#' @importFrom utils read.delim
#' @import checkmate
#' 
#' @return [\code{data.frame}] \cr
#'  A data frame of sample status and metric alerts. 
#' @export
assessMetrics <- function(data){
  
  checkmate::assertDataFrame(
    data, 
    types = c("numeric", "character")
  )
  checkmate::assertNames(
    colnames(data), 
    type = "unique", 
    must.include = "Sample ID" 
  )
  
  # Read in file of metric warning and error conditions
  conditions <- read.delim(
    system.file("template", "rmarkdown", "qc", "metricAlerts.csv", 
                package = "scReports"), 
    comment.char = "#", 
    sep = ","
  )
  # Select applicable subset of warnings and errors for which to test
  testConditions <- conditions[which(conditions$Metric %in% colnames(data)), ]
  
  # Initialize data frame of alerts
  alerts <- data.frame("Sample" = data[, "Sample ID"], "Status" = "",
                       "Alerts" = "", "Details" = "")
  
  # Populate data frame of alerts
  for (i in 1:nrow(testConditions)) {
    
    nameOfMetric <- testConditions[i, "Metric"]
    nameOfTest <- testConditions[i, "Name"]
    
    result <- 
      (data[, nameOfMetric] >= testConditions[i, "LowerLimit"]) & 
      (data[, nameOfMetric] < testConditions[i, "UpperLimit"])
    
    alerts[, nameOfTest] <- ifelse(result, testConditions[i, "Log"], "Pass") 
    alerts[which(result), "Alerts"] <- sapply(
      which(result), 
      function(x) {
        paste0(alerts[x, "Alerts"],
               testConditions[i, "Log"], ": ",
               testConditions[i, "Message"], "<br>")
      }
    )
    alerts[which(result), "Details"] <- sapply(
      which(result), 
      function(x) {
        paste0(alerts[x, "Details"], 
               testConditions[i, "Log"], ": ",
               testConditions[i, "Details"], "<br>")
      }
    )
  }

  alerts[, "Status"] <- apply(alerts[, -(1:4)], 1, 
                              function(x) {
                                if (all(x == "Pass", na.rm = TRUE)) {
                                  "OK"
                                } else if (any(x == "Error", na.rm = TRUE)) {
                                  "Error"
                                } else if (any(x == "Warning", na.rm = TRUE)) {
                                  "Warning"
                                }
                              }
  )
  alerts <- alerts[, 1:3]
  
  
  return(alerts)
}



#' Hierarchical clustering and ordering of observations
#' Uses the \code{stats::dist()} and \code{stats::hclust()} functions for 
#'  hierarchical clustering. Clustering is applied to the rows of the data \code{x}.
#'
#' @param x [\code{matrix} or \code{data.frame}] \cr
#'  A numeric matrix or data frame. Passed to \code{stats::dist}.
#' @param distMethod [\code{character(1)}] \cr
#'  The distance measure to be used, e.g., "euclidean".
#' @param clustMethod [\code{character(1)}] \cr
#'  The clustering (agglomeration) method, e.g., "ward.D", "ward.D2", or 
#'  "single". 
#'
#' @importFrom stats dist hclust
#'
#' @return [\code{list}] \cr
#'  A list of the \code{stats::hclust} output and \code{hclust} ordering of 
#'  observations.
#' @author Trishanta Padayachee  
hierarchicalOrdering <- function(x, 
                                 distMethod = "euclidean",
                                 clustMethod = "ward.D"){
  
  distData <- stats::dist(x, method = distMethod)
  distData[is.na(distData)] <- max(distData, na.rm = TRUE)
  hc <- stats::hclust(distData, method = clustMethod)
  ord <- attr(distData, "Labels")[hc$order]
  
  return(list("hclust" = hc, "order" = ord))
}



#' Compute the ADT score of an antibody for a sample
#' The antibody scores are computed by counting the number of unique UMI Counts 
#'  (excluding zero) per sample that are observed in intervals of 100 and summing 
#'  this value over the intervals.
#'  
#' @param x [\code{numeric}] \cr
#'  A vector of ADT UMI counts for a single sample.
#' 
#' @return [\code{numeric}] \cr
#'  ADT score for a particular antibody.
#' @author Trishanta Padayachee
adtScore <- function(x){
  
  uniqueCounts <- unique(x)
  intervalAllocation <- findInterval(
    uniqueCounts[uniqueCounts > 0],
    vec = seq(0, max(uniqueCounts) + (100 - max(uniqueCounts) %% 100), 100), 
    left.open = TRUE # intervals open on left and closed on the right
  )
  uniqueCountsPerInterval <- tabulate(x)
  score <- sum(uniqueCountsPerInterval)
  
  return(score)
}



#' ADT percentages of total scores per sample and antibody
#'
#' @param obj [\code{list}] \cr
#'  A named list of dgCMatrix objects. Exactly one of the names should be 'ADT'.
#'
#' @return [\code{data.frame}] \cr
#'  The percentage of total scores by sample and antibody.
#' @author Trishanta Padayachee
adtPctScores <- function(obj){
  
  counts <- obj[["ADT"]]
  abNames <- gsub("(-|_|totalseqa|totalseqb|totalseqc)", "", 
                   rownames(counts), ignore.case = TRUE)
  abScores <- unlist(sapply(seq_along(abNames),function(x) adtScore(counts[x,])))
  abPct <- data.frame("Antibody" = abNames, 
                      "AntibodyMetric" = abScores/sum(abScores)*100)
  return(abPct)
}



#' Arc-SinH cofactor
#'
#' @param x [\code{numeric}] A vector of counts for a single sample
#'
#' @return [\code{numeric}] \cr
#'  The arcsinh cofactor.
#' @author Trishanta Padayachee
arcSinHCofactor <- function(x){
  
  gm <- sinh(mean(asinh(x)))
  
  return(gm)
}



#' ADT arc-sinh cofactors per sample
#'
#' @param obj [\code{list}] \cr
#'  A named list of dgCMatrix objects. Exactly one of the names should be 'ADT'.
#'
#' @return [\code{data.frame}] \cr
#'  A data frame of arcsinh cofactors for each sample-antibody combination.
adtArcSinHCofactor <- function(obj){
  
  counts <- obj[["ADT"]]
  abNames <- gsub("(-|_|totalseqa|totalseqb|totalseqc)", "", 
                   rownames(counts), ignore.case = TRUE)
  abCofactors <- sapply(seq_along(abNames), function(x) arcSinHCofactor(counts[x,]))
  abCofactors <- data.frame("Antibody" = abNames, 
                            "AntibodyMetric" = abCofactors)
  
  return(abCofactors)
}



#' ADT counts per sample
#'
#' @param obj [\code{list}] \cr
#'  A named list of dgCMatrix objects. Exactly one of the names should be 'ADT'.
#' 
#' @importFrom Matrix rowSums
#' 
#' @return [\code{data.frame}] \cr
#'  A data frame of the counts by sample and antibody.
#' @author Trishanta Padayachee
adtCounts <- function(obj){
  
  rSums <- Matrix::rowSums(obj[["ADT"]])
  abNames <- gsub("(-|_|totalseqa|totalseqb|totalseqc)", "", 
                  names(rSums), ignore.case = TRUE)
  abCounts <- data.frame("Antibody" = abNames, 
                         "AntibodyMetric" = rSums)
  
  return(abCounts)
}


#' ADT percentages of total counts per sample
#'
#' @param obj [\code{list}] \cr
#'  A named list of dgCMatrix objects. Exactly one of the names should be 'ADT'.
#'
#' @importFrom Matrix rowSums
#' 
#' @return [\code{data.frame}] \cr
#'  A data frame of the percentage of total counts by sample and antibody.
#' @author Trishanta Padayachee
adtPctCounts <- function(obj){
  
  rSums <- Matrix::rowSums(obj[["ADT"]])
  percentages <- rSums/sum(rSums)*100
  abNames <- gsub("(-|_|totalseqa|totalseqb|totalseqc)", "", 
                  names(percentages), ignore.case = TRUE)
  abPctCounts <- data.frame("Antibody" = abNames, 
                           "AntibodyMetric" = percentages)
  
  return(abPctCounts)
}



#' Convert ADT metrics for heatmap data 
#'
#' @param rawData [\code{list}] \cr
#'  A named list of dgCMatrix objects (features by cells) of raw counts.
#' @param type [\code{character(1)}] \cr
#'  Should be either "score", "gm", "counts" or "pctCounts".
#'
#' @importFrom rlang .data :=
#' @importFrom purrr map_dfr
#' @importFrom stats median
#' @import tidyr
#' @import dplyr
#' 
#' @return [\code{list}] \cr
#'  A list of the heatmap data and hierarchical clustering results of 
#'  the samples and antibodies.
#' @author Trishanta Padayachee
#' @export
computeADTHeatmapData <- function(rawData, type){
  
  n <- length(rawData)
  func <- switch(type, 
                 score = "adtPctScores",
                 gm = "adtArcSinHCofactor",
                 counts = "adtCounts",
                 pctCounts = "adtPctCounts")
  metric <- switch(type, 
                   score = "AntibodyScores",
                   gm = "AntibodyGMs",
                   counts = "AntibodyCount",
                   pctCounts = "AntibodyPercentages")

  # Compute metric
  metricPerSample <- rawData %>%
    purrr::map(get(func)) %>%
    dplyr::bind_rows(.id = "Sample") %>%
    tidyr::complete(.data$Antibody, .data$Sample, fill = list("AntibodyMetric" = 0))
  
  # Obtain sample and antibody ordering
  metricPerSampleWide <- metricPerSample %>%
    tidyr::pivot_wider(id_cols = "Antibody", 
                       names_from = "Sample", 
                       values_from = "AntibodyMetric") %>%
    as.data.frame()
  rownames(metricPerSampleWide) <- metricPerSampleWide$Antibody
  metricPerSampleWide <- metricPerSampleWide[, -1, drop = FALSE]
  
  sampleClust <- NA
  if (n == 1) sampleOrder <- c(unique(metricPerSample$Sample), "Median")
  if (n >= 2) {
    sampleClust <- hierarchicalOrdering(t(metricPerSampleWide))
    sampleOrder <- c(sampleClust$order, "Median")
  }
  antibodyClust <- NA
  nAbs <- length(unique(metricPerSample$Antibody))
  if (nAbs == 1) antibodyOrder <- unique(metricPerSample$Antibody)
  if (nAbs >= 2) {
    antibodyClust <- hierarchicalOrdering(metricPerSampleWide)
    antibodyOrder <- antibodyClust$order
  }
  
  
  # Compute median count across samples and apply ordering
  heatmapData <- metricPerSample %>%
    dplyr::group_by(.data$Antibody) %>%
    dplyr::summarise("Median" = median(.data$AntibodyMetric, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Sample = "Median") %>%
    dplyr::rename(AntibodyMetric = "Median") %>%
    dplyr::bind_rows(metricPerSample) %>%
    dplyr::mutate("Sample" := factor(.data$Sample, ordered = TRUE, levels = sampleOrder)) %>%
    dplyr::mutate("Antibody" := factor(.data$Antibody, ordered = TRUE, levels = antibodyOrder)) %>%
    dplyr::rename(!!metric := "AntibodyMetric")
  
  return(list(
    "heatmapData" = heatmapData,
    "antibodyClust" = antibodyClust,
    "sampleClust" = sampleClust
  ))
}



#' Exclude the top N points based on two variables
#' 
#' @param topNCutoff [\code{logical(1)}, default: \code{TRUE}] \cr
#'  Number of top cells to exclude (a temporary solution
#'  to remove the effect of potential outliers). The top N based on 
#'  \code{featureX} and the top N based on \code{featureY} are removed.
#' @inheritParams commonParams
#'  
#' @author Trishanta Padayachee  
.scatterExcludeTopN <- function(df, xObs, yObs, 
                         topNCutoff = 10){
  
  dropLocations <- unique(
      c(order(df[, xObs], decreasing = TRUE)[1:topNCutoff],
        order(df[, yObs], decreasing = TRUE)[1:topNCutoff])
  )

  dfReduced <- df[-dropLocations, ]
  
  return(dfReduced)
}





