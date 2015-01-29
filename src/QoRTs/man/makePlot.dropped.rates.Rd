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
makePlot.dropped.rates(plotter, dropAlwaysZeroRows = FALSE, 
                       debugMode, singleEndMode, \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{dropAlwaysZeroRows}{
    Logical. If TRUE, drop-reasons that never occur in the dataset will not be plotted. This often cleans up the plot somewhat, since in many production pipelines reads that fail many of the filtering steps may have already been filtered out.
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
   For each bam file, this function plots the rates and reasons for reads being dropped from QC analysis. 
   
   Note that in the example dataset reads were never dropped. This is a consequence of the pre-processing steps in the example pipeline.
}

\seealso{
  \code{\link{build.plotter}}
}