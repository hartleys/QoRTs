\name{get.summary.table}
\docType{methods}
\alias{get.summary.table}
\alias{get.size.factors}

\title{
   Get summary data tables
}
\description{
   Retrieves and compiles a summary data table.
}
\usage{
  get.summary.table(res, outfile, debugMode);
  
  get.size.factors(res, 
       sf.method = c("DESeq2","DESeq2_GEO","TC",
                     "edgeR","edgeR_TMM","edgeR_UQ","edgeR_RLE"), 
       outfile, debugMode)
}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  \item{outfile}{
    Optional. A file name where the table should be written.
  }
  \item{sf.method}{
    The size factor method to use. Note that most of these methods (except "TC") are reliant
    on external packages not included with QoRTs. You will need to install DESeq2 and edgeR
    to use these methods.
  }
  
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
}
\details{
  Returns summary data in table form.
}
\examples{
  data(res,package="QoRTsExampleData");
  get.summary.table(res);
}

\seealso{
  \code{\link{read.qc.results.data}} 
}