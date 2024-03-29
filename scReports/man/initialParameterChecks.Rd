% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checks.R
\name{initialParameterChecks}
\alias{initialParameterChecks}
\title{Initial parameter checks}
\usage{
initialParameterChecks(sampleNames, samplePaths, sampleCRQC, metaFile)
}
\arguments{
\item{sampleNames}{[\code{character}] \cr 
Sample names. The order of the names should correspond to the file paths 
provided to the argument \code{samplePaths} and/or \code{sampleCRQC}. 
These names are used to label the samples and need not be identical to the 
actual file names.}

\item{samplePaths}{[\code{character}] \cr 
File paths to either 10xGenomics h5 files or to AnnData h5ad files. If 
\code{sampleNames} are also specified, the order of the file paths should 
correspond to the order of the sample names in \code{sampleNames}.}

\item{sampleCRQC}{[\code{character}] \cr 
File paths to Cell Ranger outputted csv files containing sample-specific 
summary metrics. The order of the file paths in \code{sampleCRQC} should
correspond to the order of the sample names in \code{sampleNames}.}

\item{metaFile}{[\code{character(1)}] \cr 
Path to a meta.tsv file containing meta data on all the samples analysed.}
}
\value{
None
}
\description{
Checks that:
\itemize{
 \item{the length of sampleNames, samplePaths and/or sampleCRQC are the same}
 \item{the files in samplePaths and sampleCRQC exist}
}
}
\author{
Trishanta Padayachee
}
