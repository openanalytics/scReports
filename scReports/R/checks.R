#' Initial parameter checks
#'
#' Checks that:
#' \itemize{
#'  \item{the length of sampleNames, samplePaths and/or sampleCRQC are the same}
#'  \item{the files in samplePaths and sampleCRQC exist}
#' }
#' 
#' @inheritParams commonParams
#'
#' @return None
#' @author Trishanta Padayachee
#' @export
initialParameterChecks <- function(sampleNames,
                                   samplePaths, 
                                   sampleCRQC,
                                   metaFile){
  
  ## Check that the lengths of sampleNames, samplePaths and/or sampleCRQC 
  ## are the same
  nNames <- length(sampleNames)
  nPaths <- length(samplePaths)
  nCRQC <- length(sampleCRQC)
  
  if (!is.null(sampleCRQC) & (nNames != nCRQC)) {
    
    stop("Length of sampleNames and sampleCRQC should be the same.")
    
  }
  
  if (!is.null(samplePaths) & (nNames != nPaths)) {
    
    stop("Length of sampleNames and samplePaths should be the same.")
    
  }
  
  ## Check that files exist
  paths <- c(metaFile, samplePaths, sampleCRQC)
  fileExists <- file.exists(paths)
  if (any(!fileExists)) {

    stop(
      "One or more files don't exist.",
      "Check that the following paths have been correctly specified:\n",
      paste0(paths[!fileExists], sep = "\n")
    )
    
  }
  
  
}





