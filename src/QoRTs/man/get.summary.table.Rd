\name{get.summary.table}
\docType{methods}
\alias{get.summary.table}
\title{
   Get summary data table
}
\description{
   Retrieves and compiles a summary data table.
}
\usage{
  get.summary.table(res, outfile, debugMode)
}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  \item{outfile}{
    Optional. A file name where the table should be written.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
}
\details{
  Returns summary data in table form.
}

\seealso{
  \code{\link{read.qc.results.data}} 
}