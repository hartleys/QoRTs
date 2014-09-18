\name{makePlot.mapping.rates}
\docType{methods}
\alias{makePlot.mapping.rates}
\title{
   Plot mapping rates
}
\description{
   Plots the rates at which reads map to the genome.
}
\usage{
makePlot.mapping.rates(plotter, plot.mm = TRUE, 
                   y.counts.in.millions = TRUE, 
                   debugMode,
                   \dots)
}
\arguments{
  \item{plotter}{
    A QoRT_Plotter object.
  }
  \item{plot.mm}{
    Plot multi-mapping rates.
  }
  \item{y.counts.in.millions}{
    Label/scale the y-axis in millions of reads.
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
   TODO!
}

\seealso{

}