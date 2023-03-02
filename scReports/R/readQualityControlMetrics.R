#' Read in the Cell Ranger gene expression metrics
#'
#' \code{readCRQC} Reads in the Cell Ranger 3 gene expression metrics for
#'  a single sample or multiple samples.
#'  
#' @param metrics [\code{data.frame}] \cr
#'  A single data frame of the QC metrics for all samples in \code{sampleNames}. 
#'  Order of rows in the data frame should correspond to the order of names in 
#'  \code{sampleNames}.  
#' @inheritParams commonParams
#'  
#' @importFrom data.table fread rbindlist
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr mutate_all
#' @importFrom stringr str_replace_all 
#'
#' @return [\code{data.frame}] \cr
#'  A data frame containing the loaded gene expression metrics (columns) 
#'  for each sample (rows) in \code{sampleNames}.
#' @author Trishanta Padayachee
readCRQC <- function(sampleNames, 
                     sampleCRQC, 
                     metrics = NULL) {

  ## Read in data (samples X metrics)
  if (is.null(metrics)) {
    
    metrics <- lapply(sampleCRQC, data.table::fread)
    metrics <- data.table::rbindlist(metrics, fill = TRUE)
    
  }

  ## Remove commas and percentage signs and convert columns to numeric
  metrics <- tibble::as_tibble(metrics) %>%
    dplyr::mutate_all(stringr::str_replace_all, pattern = ",", replacement = "") %>%
    dplyr::mutate_all(stringr::str_replace_all, pattern = "%", replacement = "") %>%
    dplyr::mutate_all(as.numeric) %>%
    tibble::add_column("Sample ID" = sampleNames, .before = 1) %>% 
    as.data.frame

  return(metrics)
}



#' Convert Cell Ranger multi metric_summary.csv format to Cell Ranger count format
#'
#' @param path [\code{character}] \cr 
#'  Path to a metrics_summary.csv file obtained through the Cell Ranger multi 
#'  pipeline.
#' 
#' @importFrom purrr map_vec
#' @importFrom utils read.csv
#' 
#' @return A data frame of QC metrics in the layout of the Cell Ranger count
#'  metrics_summary.csv files. 
#' @export
convertCellRangerMultiToCellRangerCount <- function(path) {
  
  multi <- read.csv(path)
  mappings <- read.csv(
    system.file("template", "rmarkdown", "qc", "metricMappings.csv",
                package = "scReports")
  )
  multi["libraryMetric"] <- paste(
    multi[,"Library.Type"], multi[,"Metric.Name"], 
    sep = ": "
  )
  
  lookUp <- setNames(mappings[,"cellRangerCount"], mappings[,"cellRangerMulti"])
  multi["libraryMetricMapping"] <- purrr::map_vec(multi[,"libraryMetric"], ~ lookUp[.x])
  multi["libraryMetricRenamed"] <- multi["libraryMetricMapping"]
  multi[is.na(multi["libraryMetricMapping"]), "libraryMetricRenamed"] <- 
    multi[is.na(multi["libraryMetricMapping"]), "libraryMetric"]
  
  count <- pivot_wider(
    data = distinct(multi[c("libraryMetricRenamed", "Metric.Value")]), 
    names_from = "libraryMetricRenamed", 
    values_from = "Metric.Value"
  )
  
  return(count)
  
}

