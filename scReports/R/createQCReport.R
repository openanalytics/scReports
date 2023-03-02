#' Create a QC report for single cell data based on user specified parameters
#'
#' @inheritParams commonParams
#' @param TA [\code{character(1)}] \cr
#'  Therapeutic area.
#' @param project [\code{character(1)}] \cr
#'  Project code and project name.
#' @param experiment [\code{character(1)}] \cr
#'  Project code, experiment ID and experiment name.
#' @param analyst [\code{character(1)}] \cr
#'  Name of the individual generating the report.
#' @param autoRotate [\code{logical(1)}, default: \code{TRUE}] \cr
#'  If \code{TRUE}, figures will auto rotate (flip) based on the number of 
#'  samples analysed. Note, auto-rotation can only be turned off for less than 
#'  30 samples. 
#' @param cellRangerCount [\code{logical(1)}, default: \code{TRUE}] \cr
#'  Logical indicating whether metrics_summary.csv files specified in 
#'  \code{sampleCRQC} where obtained using the Cell Ranger count pipeline. 
#'  \code{FALSE} indicates that the Cell Ranger multi pipeline was used. 
#' @param saveRObjects [\code{logical(1)}, default: \code{FALSE}] \cr
#'  Logical value indicating whether to save the R objects. 
#'  Note, this can be extremely time consuming. 
#' @param workingDirectory [\code{character(1)}] \cr
#'  Working directory.
#' @param outputDirectory [\code{character(1)}] \cr
#'  Output directory.
#' @param quietRendering [\code{logical(1)}, default: \code{FALSE}] \cr
#'  An option to suppress printing of the pandoc command line.
#' @param isFilteredData [\code{logical(1)}, default: \code{FALSE}] \cr
#'  Whether \code{samplePaths} point to filtered h5ad files.
#'  If \code{TRUE}, the QC metrics stored in the object are used (i.e., the metrics 
#'  aren't recomputed). This is only possible for h5ad files with n_counts, 
#'  n_genes and/or percent_mito columns in the obs slot. 
#' 
#' @import Cairo 
#' @importFrom stringi stri_split_fixed
#' @importFrom rmarkdown render
#'
#' @author Trishanta Padayachee
#' @export
createQCReport <- function(TA, project, experiment,
                           sampleIDs, sampleNames, 
                           samplePaths, 
                           sampleCRQC = NULL, 
                           metaFile = NULL,
                           analyst = NULL,
                           autoRotate = TRUE,
                           horizontalBarCharts = FALSE,
                           horizontalViolinPlots = FALSE,
                           isFilteredData = FALSE,
                           cellRangerCount = TRUE,
                           saveRObjects = FALSE,
                           workingDirectory = "./wdir/",
                           outputDirectory = "./wdir/",
                           quietRendering = FALSE) {
  
  # Cairo isn't explicitly used in this function
  # import Cairo was included as a temp fix for devtools::check() which is 
  # producing Cairo errors
  inputParams <- as.list(environment())
  ## VIASH adaptation
  if (!is.null(sampleCRQC) & length(sampleCRQC) == 1) {
    
    if (sampleCRQC == "") inputParams[["sampleCRQC"]] <- NULL
    
  }
  if (!is.null(metaFile)) { 
    
    if (metaFile == "") inputParams[["metaFile"]] <- NULL
    
  }
  if (!is.null(analyst)) { 
    
    if (analyst == "") inputParams[["analyst"]] <- NULL
    
  }
  
  for (param in c("sampleIDs", "sampleNames", "samplePaths", "sampleCRQC")) {
    inputParams[[param]] = gsub(" ", "", inputParams[[param]])
  }
  remove <- which(names(inputParams) == "quietRendering")
  inputParams <- lapply(inputParams[-remove],
                        function(x) switch(
                          is.character(x) + 1, x, 
                          unlist(stringi::stri_split_fixed(x, ","))))
  
  ## Check that all the necessary input parameters are defined
  necessaryParams <- c("TA", "project", "experiment", "sampleIDs", 
                       "sampleNames",  "samplePaths")
  if (!all(necessaryParams %in% names(inputParams))) {
    
    missingParams <-  necessaryParams[!(necessaryParams %in% names(inputParams))]
    stop("Parameters missing. The following parameter(s) should be specified: ", 
         paste(missingParams, collapse = ", "))
    
  }
  
  ## Create working directory if it does not exist
  if (!dir.exists(workingDirectory))
    dir.create(workingDirectory, recursive = TRUE)
  ## Create output directory if it does not exist
  if (!dir.exists(outputDirectory))
    dir.create(outputDirectory, recursive = TRUE)
  
  ## Normalize paths
  inputParams[["workingDirectory"]] <- normalizePath(workingDirectory)
  inputParams[["outputDirectory"]] <- normalizePath(outputDirectory)
  inputParams[["samplePaths"]] <- sapply(inputParams[["samplePaths"]], normalizePath)
  if (!is.null(inputParams[["sampleCRQC"]])) {
    inputParams[["sampleCRQC"]] <- sapply(inputParams[["sampleCRQC"]], normalizePath)
  }
  if (!is.null(inputParams[["metaFile"]])) {
    inputParams[["metaFile"]] <- normalizePath(inputParams[["metaFile"]])
  }
  
  ## Copy rmarkdown file to the working directory
  file.copy(from = system.file("template", "rmarkdown", "qc", "qualityControl.Rmd",
                               package = "scReports"),
            to = workingDirectory,
            overwrite = TRUE)
  ## Copy style.css file to the workingDirectory
  file.copy(from = system.file("template", "rmarkdown", "style.css",
                               package = "scReports"),
            to = workingDirectory, 
            overwrite = TRUE)
  ## Copy 10xGenomics_QCDefinitions to the workingDirectory
  file.copy(from = system.file("template", "rmarkdown", "10xGenomics_QCDefinitions.yml", 
                               package = "scReports"),
            to = workingDirectory, 
            overwrite = TRUE)
  
  ## Render the report using user specified parameters
  rmarkdown::render(paste0(workingDirectory, "/qualityControl.Rmd"),
                    params = inputParams,
                    quiet = quietRendering)
  
  ## Copy output from working directory to output directory, if different.
  if (workingDirectory != outputDirectory) {
    
    file.copy(from = paste0(workingDirectory, "/qualityControl.html"), 
              to = outputDirectory, 
              overwrite = TRUE)
    
    if (saveRObjects) {
      file.copy(from = paste0(workingDirectory, "/qc.rds"), 
                to = outputDirectory, 
                overwrite = TRUE)
      file.copy(from = paste0(workingDirectory, "/rawData.rds"), 
                to = outputDirectory, 
                overwrite = TRUE)
    }
  }
}
