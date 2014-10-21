\name{makePlot.dropped.rates}
\docType{methods}
\alias{makePlot.dropped.rates}
\title{
   Read Drop Plot
}
\description{
   Plots the rates at which reads are dropped from
   analysis for various causes.
}
\usage{
makePlot.dropped.rates(plotter, debugMode, \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   For each bam file, this function plots the rates and reasons for reads being dropped from QC analysis. 
   
   Note that in the example dataset reads were never dropped. This is a consequence of the pre-processing steps in the example pipeline.
}

\seealso{
  \code{\link{build.plotter}}
}