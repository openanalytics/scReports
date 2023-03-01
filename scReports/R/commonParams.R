#' Arguments that are frequently used in the scReports package
#'
#' @param sampleCRQC [\code{character}] \cr 
#'  File paths to Cell Ranger outputted csv files containing sample-specific 
#'  summary metrics. The order of the file paths in \code{sampleCRQC} should
#'  correspond to the order of the sample names in \code{sampleNames}.
#' @param sampleIDs [\code{character}] \cr 
#'  Sample IDs. The order of the IDs should correspond to the order of the file 
#'  paths provided to the argument \code{samplePaths}. These IDs are typically 
#'  less descriptive than those provided to the \code{sampleNames} argument.   
#' @param sampleNames [\code{character}] \cr 
#'  Sample names. The order of the names should correspond to the file paths 
#'  provided to the argument \code{samplePaths} and/or \code{sampleCRQC}. 
#'  These names are used to label the samples and need not be identical to the 
#'  actual file names.
#' @param sampleName [\code{character(1)}] \cr 
#'  Sample name.  
#' @param samplePaths [\code{character}] \cr 
#'  File paths to either 10xGenomics h5 files or to AnnData h5ad files. If 
#'  \code{sampleNames} are also specified, the order of the file paths should 
#'  correspond to the order of the sample names in \code{sampleNames}. 
#' @param path [\code{character(1)}] \cr 
#'  Absolute or relative file path.
#' @param metaFile [\code{character(1)}] \cr 
#'  Path to a meta.tsv file containing meta data on all the samples analysed.
#' @param verbose [\code{logical(1)}, default: \code{FALSE}] \cr 
#'  Should progress messages be displayed? The default value, \code{FALSE}, 
#'  indicates that logging messages should be kept to a minimum.   
#' @param makeFeaturesUnique [\code{logical(1)}, default: \code{TRUE}] \cr 
#'  Whether to make feature names unique. The default value, \code{TRUE}, 
#'  implies that if feature names are not unique, they should be made unique. 
#'  This is done using the \code{base::make.unique} function. If \code{FALSE}, 
#'  the original feature names are kept.
#' @param horizontalBarCharts [\code{logical(1)}, default: \code{FALSE}] \cr 
#'  If \code{TRUE}, samples are plotted on the y-axis resulting in horizontal 
#'  bars. This layout allows for figure heights to adjust as a function of the 
#'  number of samples. 
#' @param horizontalViolinPlots [\code{logical(1)}, default: \code{TRUE}] \cr 
#'  If \code{TRUE}, samples are plotted on the y-axis. This layout allows for 
#'  figure heights to be adjusted as a function of the number of samples.
#' @param geneFilteringMinCells [\code{numeric(1)}, default: \code{3}] \cr 
#'  Genes with fewer than \code{geneFilteringMinCells} cells will be filtered 
#'  out.
#' @param adtFilteringMinCells [\code{numeric(1)}, default: \code{0}] \cr 
#'  Antibodies with fewer than \code{adtFilteringMinCells} cells will be 
#'  filtered out. Since the absence of certain antibodies is also informative, 
#'  the default value is \code{0}. 
#' @param antibodyFilteringMinCells [\code{numeric(1)}, default: \code{0}] \cr 
#'  Antibodies with fewer than \code{antibodyFilteringMinCells} cells will be 
#'  filtered out. Since the absence of certain antibodies is also informative, 
#'  the default value is \code{0}. 
#' @param cellFilteringMinGenes [\code{numeric(1)}, default: \code{1}] \cr 
#'  Cells with counts for fewer than \code{cellFilteringMinGenes} genes will be 
#'  filtered out. This filtering is applied when only Gene Expression data is 
#'  available for a sample.
#' @param cellFilteringMaxGenes [\code{numeric(1)}, default: \code{1}] \cr 
#'  Cells with counts in greater than \code{cellFilteringMaxGenes} genes will 
#'  be filtered out. This filtering is applied when only Gene Expression data 
#'  is available for a sample.
#' @param cellFilteringMinCounts [\code{numeric(1)}, default: \code{10}] \cr 
#'  Cells with fewer than \code{cellFilteringMinCounts} total counts will be 
#'  filtered out. If the object contains both gene expression and antibody 
#'  capture data, cells with fewer than \code{cellFilteringMinCounts} RNA or 
#'  ADT counts will be filtered out.
#' @param cellFilteringMaxCounts [\code{numeric(1)}, default: \code{10}] \cr 
#'  Cells with fewer than \code{cellFilteringMaxCounts} total counts will be 
#'  filtered out. If the object contains both gene expression and antibody 
#'  capture data, cells with fewer than \code{cellFilteringMinCounts} RNA or 
#'  ADT counts will be filtered out. 
#' @param cellFilteringMinCountsRNA [\code{numeric(1)}, default: \code{10}] \cr 
#'  Cells with fewer than \code{cellFilteringMinCountsRNA} total gene counts 
#'  will be filtered out. 
#' @param cellFilteringMaxCountsRNA [\code{numeric(1)}, default: \code{10}] \cr 
#'  Cells with greater than \code{cellFilteringMaxCountsRNA} total gene counts 
#'  will be filtered out.  
#' @param cellFilteringMinCountsADT [\code{numeric(1)}, default: \code{10}] \cr 
#'  Cells with fewer than \code{cellFilteringMinCountsADT} total antibody 
#'  counts will be filtered out. 
#' @param cellFilteringMaxCountsADT [\code{numeric(1)}), default: \code{10}] \cr 
#'  Cells with greater than \code{cellFilteringMaxCountsADT} total antibody 
#'  counts will be filtered out. 
#' @param cellFilteringMaxPercentMito [\code{numeric(1)}, default: \code{100}] \cr 
#'  Cells with a percentage of mitochondrial gene counts greater than the 
#'  specified value will be filtered out.
#' @param colours [\code{character}] \cr
#'  A character vector of colours to be used for plotting.
#' @param palette [\code{character(1)}] \cr
#'  Name of an \code{RColourBrewer} palette. 
#' @param title [\code{character(1)}] \cr
#'  Plot title.
#' @param xLabel [\code{character(1)}] \cr
#'  X-axis label.
#' @param yLabel [\code{character(1)}] \cr
#'  Y-axis label.
#' @param lLabel \cr
#'  Legend name.
#' @param xAngle [\code{numeric(1)}, default: 0] \cr
#'  Angle of text on the x-axis. A value between 0 and 360.
#' @param xAxisTextAngle [\code{numeric(1)}, default: 0] \cr
#'  Angle of text on the x-axis. A value between 0 and 360.
#' @param textSize [\code{numeric(1)}] \cr
#'  Value representing the font size.  
#' @param addText [\code{logical(1)}, default: \code{TRUE}] \cr
#'  Whether to add text to the plot.
#' @param makeInteractive [\code{logical(1)}, default: \code{TRUE}] \cr
#'  Should an interactive plot be returned?
#' @param df [\code{list}] \cr
#'  A \code{data.frame}. 
#' @param obs [\code{character(1)}] \cr
#'  Name of the observation to plot. Should correspond to a column in \code{df}. 
#' @param obsX Name of the observation to plot on the X-axis. Should correspond 
#'  to a column in \code{df}. 
#' @param obsY Name of the observation to plot on the Y-axis. Should correspond 
#'  to a column in \code{df}. 
#' @param xObs Name of the observation to plot on the X-axis. Should correspond 
#'  to a column in \code{df}. 
#' @param yObs Name of the observation to plot on the Y-axis. Should correspond 
#'  to a column in \code{df}.
#' @param xLimits [\code{numeric(2)}] \cr
#'  X-axis limits.
#' @param yLimits [\code{numeric(2)}] \cr
#'  Y-axis limits.   
#' @param fillBy Name of a column in \code{df} that colours should be based on.
#' @param includeFigDimensions [\code{logical(1)}, default: \code{FALSE}] \cr
#'  Whether to include figure dimensions in the R chunk options.
#' @param figWidth [\code{numeric}] \cr
#'  A vector of figure widths with length equal to that of \code{plotList}.
#' @param figHeight [\code{numeric}] \cr
#'  A vector of figure heights with length equal to that of \code{plotList}.
#' @param includeOutDimensions [\code{logical(1)}] \cr
#'  Logical indicating whether to include output dimensions in the R chunk 
#'  options.   
#' @param outWidth [\code{numeric}] \cr
#'  A vector of output widths with length equal to that of \code{plotList}.
#' @param outHeight [\code{numeric}] \cr
#'  A vector of output heights with length equal to that of \code{plotList}.
#' @param tableList [\code{list}] \cr
#'  Optionally, a named list of tables.
#' @param plotList [\code{list}] \cr
#'  Optionally, a named list of plots. 
#' @param plotPaths [\code{list}] \cr
#'  Optionally, a named list of paths to saved figures. \code{plotPaths} can be 
#'  specified instead of a \code{plotList}.
#' @param contentNames [\code{character}] \cr
#'  A vector of names corresponding to the names in \code{textList}, 
#'  \code{plotList}, \code{plotPaths} and \code{tableList}, if specified.
#' @param textList [\code{list}] \cr
#'  Optionally, a named list of paragraphs. If an element of the list is 
#'  \code{NULL} for a particular tab, no text will be displayed in that tab.
#' @param performChecks \code{logical(1)}, default: \code{TRUE}] \cr
#'  Whether to perform routine parameter checks before proceeding.
#' @name commonParams
NULL