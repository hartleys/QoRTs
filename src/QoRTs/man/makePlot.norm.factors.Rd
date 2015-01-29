\name{makePlot.norm.factors}
\docType{methods}
\alias{makePlot.norm.factors}
\alias{makePlot.norm.factors.vs.TC}
\title{
   Plot Alignment Clipping
}
\description{
   Plots the rate at which the aligner soft-clips off portions of
   aligned reads.
}
\usage{
  makePlot.norm.factors(plotter, by.sample = TRUE, 
                return.table = FALSE,
                debugMode, singleEndMode, ...)
  makePlot.norm.factors.vs.TC(plotter, 
                by.sample = TRUE, 
                return.table = FALSE,
                debugMode, singleEndMode, ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{by.sample}{
    Logical. Whether to combine all lanebams for each sample before 
    calculating normalization factors. By default, normalization factors
    for each lanebam will be calculated seperately.
  }
  \item{return.table}{
    Logical. Whether to return a data table containing the data that is plotted.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  makePlot.norm.factors plots the normalization factors for each sample/lanebam. 
  Note that unless DESeq2 and edgeR are installed, it will only plot total-count
  normalization by default. Also note that unless calc.DESeq2 = TRUE and/or
  calc.edgeR = TRUE in the execution of the read.qc.results.data that produced the
  QC results, only the total counts normalization factors will be calculated.
  
  makePlot.norm.factors.vs.TC plots the ratio of alternative normalization factors against
  the total count normalization.
}
\value{
  Usually returns nothing, unless return.table is TRUE, in which case it returns a data.frame containing the plotted data for each sample. 
}
\seealso{
  \code{\link{build.plotter}}

}