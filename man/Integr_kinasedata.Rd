% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinases.R
\name{Integr_kinasedata}
\alias{Integr_kinasedata}
\title{Integrate the kinase library, infer the SiKs and calculate the EI}
\usage{
Integr_kinasedata(
  object,
  kldata,
  mean_method = NULL,
  cell_group,
  fold_threshold = 10
)
}
\arguments{
\item{object}{Incytr object}

\item{kldata}{data frame of the kinase library with variables 'gene', 'site_pos', and 'motif.geneName'}

\item{mean_method}{the method name used to calculate the average expressed value. NULL is the default value, and the arithmetic mean is used if it is "mean".}

\item{cell_group}{the group/cluster to calculate the EI}

\item{fold_threshold}{if the difference between the selected value and the second highest value passes the threshold, then we think it is "highly exclusive", EI = 1. The default setting is fold_threshold = 10.}
}
\value{
an Incytr object
}
\description{
Integrate the kinase library, infer the SiKs and calculate the EI
}
