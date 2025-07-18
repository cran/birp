% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/birp_data.R
\name{birp_data_from_data_frame}
\alias{birp_data_from_data_frame}
\title{Create a 'birp_data' Object from a Data Frame or List of Data Frames}
\usage{
birp_data_from_data_frame(data)
}
\arguments{
\item{data}{A data frame or a list of data frames, each representing one method. Each data frame must include the following columns: 'timepoint', 'location', 'counts', 'effort', and 'CI_group'. Each row represents a survey conducted at a specific timepoint and location, and for a specific control-intervention (CI) group.}
}
\value{
An object of type \link{birp_data}
}
\description{
Constructs a 'birp_data' object from a single data frame or a list of data frames, with one per method.
}
\examples{
df <- data.frame(
  timepoint = 1:10,
  location = rep(1, 10),
  counts = runif(10, 0, 100),
  effort = rexp(10),
  CI_group = "intervention"
)
data <- birp_data_from_data_frame(df)
}
