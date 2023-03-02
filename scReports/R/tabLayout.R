#' Create an R markdown tab layout in which each tab can contains a paragraph of 
#' text and/or a figure or table.
#'
#' @param tabNames [\code{character}] \cr
#'  A vector of tab names with length equal to that of \code{contentNames}. 
#' @param tabLevel [\code{numeric(1)}] \cr
#'  An integer representing the level of the tab-heading. 
#' @param addContentNameToTextList [\code{logical(1)}, default: \code{TRUE}] \cr
#'  Whether to prepend the associated tab name to the \code{textList} entry.  
#' @param plotString [\code{character}] \cr
#'  Optionally, a string that contains the \code{sprintf} text for creating 
#'  the plots in \code{plotList}. 
#' @param grobList [\code{list}] \cr
#'  Optionally, a named list of grobs.
#' @inheritParams commonParams
#'
#' @importFrom knitr include_graphics
#' 
#' @return [\code{character}] \cr
#'  A string of code for creating the tab header and tab content.
#' @author Trishanta Padayachee
tabLayout <- function(tabNames, tabLevel, 
                      contentNames,
                      textList = NULL, 
                      addContentNameToTextList = TRUE,
                      plotList = NULL, plotPaths = NULL, plotString = NULL,
                      tableList = NULL,
                      grobList = NULL,
                      includeFigDimensions = FALSE, 
                      figWidth = NA, figHeight = NA,
                      includeOutDimensions = FALSE, 
                      outWidth = NA, outHeight = NA) {
  
  ## Check that tabNames and contentNames are of the same length
  if (length(tabNames) != length(contentNames)) { 
    
    stop("Length of tabNames should be equal to the length of contentNames.\n")
  
  }
  
  ## Check that contentNames correspond with names in textList, plotList, 
  ## plotPaths, grobList and tableList
  if (!is.null(textList)) {

    contentPresent <- all(contentNames %in% names(textList))

    if (contentPresent == FALSE)
      warning("Text for one or more contentNames is missing from textList.\n")

  }
  if (!is.null(plotList)) {

    contentPresent <- all(contentNames %in% names(plotList))

    if (contentPresent == FALSE)
      warning("Figure for one or more contentNames is missing from plotList.\n")

  }
  if (!is.null(plotPaths)) {

    contentPresent <- all(contentNames %in% names(plotPaths))

    if (contentPresent == FALSE)
      warning("Figure for one or more contentNames is missing from plotPaths.\n")

  }
  if (!is.null(tableList)) {

    contentPresent <- all(contentNames %in% names(tableList))

    if (contentPresent == FALSE)
      warning("Table for one or more contentNames is missing from tableList.\n")

  }
  if (!is.null(grobList)) {
    
    contentPresent <- all(contentNames %in% names(grobList))
    
    if (contentPresent == FALSE)
      warning("Grob for one or more contentNames is missing from grobList.\n")
    
  }

  
  
  k <- 0
  codeString <- NULL
  for (tabName in tabNames) {
    
    k <- k + 1
    hashString <- paste(rep("#", tabLevel), collapse = "")
    name <- contentNames[k]
    chunkName <- gsub("[[:punct:]]|[[:space:]]","", name)
    
    ## Tab name
    line1 <- sprintf(paste0(hashString, " %s \n"), 
                     paste0(tabName, "{-}"))
    
    ## Tab text
    line2 <- ""
    if (!is.null(textList)) {
      
      if (addContentNameToTextList) {
        line2 <- sprintf("**%s**: %s \n", 
                         name, textList[[name]])
      } else {
        
        line2 <- sprintf("%s \n", 
                         textList[[name]])
      }
      
    }
    
    ## Start r chunk
    if (includeFigDimensions == TRUE & includeOutDimensions == FALSE) {
      
      line3 <- sprintf("```{r %s, fig.width = %f, fig.height = %f} \n", 
                       chunkName, figWidth[k], figHeight[k]) 
      
    } else if (includeFigDimensions == FALSE & includeOutDimensions == TRUE) {
      
      line3 <- sprintf("```{r %s,  out.width = '%spx', out.height = '%spx'} \n", 
                       chunkName, outWidth[k], outHeight[k])
      
    } else if (includeFigDimensions == TRUE & includeOutDimensions == TRUE) {
      
      line3 <- sprintf("```{r %s,  out.width = '%spx', out.height = '%spx', fig.width = %f, fig.height = %f} \n", 
                       chunkName, outWidth[k], outHeight[k], figWidth[k], figHeight[k])
      
    } else {
      
      line3 <- sprintf("```{r %s} \n", 
                       chunkName)
      
    }
    
    ## Chunk content
    if (!is.null(plotPaths)) {
      
      line4 <- sprintf("knitr::include_graphics(%s[['%s']]) \n", 
                       "plotPaths", name)
      
    } else  if (!is.null(plotList)) {
      
      if (is.null(plotString)) {
        
        if (!is.null(plotList[[name]])) {
          
          line4 <- sprintf("%s[['%s']] \n", 
                           "plotList", name)
          
        } else {
          
          line4 <- ""
          
        }

      } else {
        
        line4 <- sprintf(paste0(plotString, " \n"),
                         "plotList", name)
        
      }
      
    } else if (!is.null(tableList)) {
      
      line4 <- sprintf("%s[['%s']] \n", 
                       "tableList", name)

    } else if (!is.null(grobList)) {
      
      line4 <- sprintf("grid::grid.draw(%s[['%s']]) \n", 
                       "grobList", name)
    } else {
      
      line4 = ""
      
    }
    
    ## End r chunk
    line5 <- "``` \n" 
    
    
    codeString[[k]] <- paste0(line1, line2, line3, line4, line5)
    
  }
  
  return(codeString)
}

