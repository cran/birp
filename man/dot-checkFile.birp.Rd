% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/birp.R
\name{.checkFile.birp}
\alias{.checkFile.birp}
\title{Function to check if a file exists and generate error message if it was not found}
\usage{
.checkFile.birp(
  path,
  files,
  patterns,
  allowMultiMatch = FALSE,
  sep = "\\t",
  mustExist = TRUE
)
}
\arguments{
\item{path}{A file path}

\item{files}{A vector of character strings corresponding to file names found in the path}

\item{patterns}{A vector of patterns to search for within 'files'}

\item{allowMultiMatch}{Logical. If \code{TRUE}, multiple matches are allowed}

\item{sep}{The field separator character of the input file}

\item{mustExist}{Logical. If \code{TRUE}, an exception is raised if the file does not exist}
}
\value{
A string denoting the filename(s)
}
\description{
Function to check if a file exists and generate error message if it was not found
}
\keyword{internal}
