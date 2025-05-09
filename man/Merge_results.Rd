% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{Merge_results}
\alias{Merge_results}
\title{Merge the tables from the pairwise analysis (i.e. multiple rounds of analysis)}
\usage{
Merge_results(
  results.list,
  cutoff_PDS = NULL,
  group = "Both",
  select.top = NULL
)
}
\arguments{
\item{results.list}{a list of tables for multiple rounds of analysis}

\item{cutoff_PDS}{selecting pathways whose PDS is higher than the cutoff value, the default value is NULL}

\item{group}{same pathway may appear in different Sender-Receiver communications, grouping the pathways by "Sender.group", "Receiver.group", or "Both". The default value is "Both".}

\item{select.top}{select the top pathways based on the PDS value, in default select.top = NULL}
}
\value{
a data frame
}
\description{
Merge the tables from the pairwise analysis (i.e. multiple rounds of analysis)
}
