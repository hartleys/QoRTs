\name{makePlot.missingness.rate}
\docType{methods}
\alias{makePlot.missingness.rate}
\title{
   Plot N-Rate by read cycle
}
\description{
   Plots rate by which the sequencer cannot determine the base, 
   by read cycle.
}
\usage{
  makePlot.missingness.rate(plotter, r2.buffer = NULL, 
                        debugMode, ...)
}
\arguments{
  \item{plotter}{
    A QoRT_Plotter object.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1 and read 2.
    By default this will choose a reasonable value.
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
  x-axis: Read Cycle
  
  y-axis: Rate at which the sequencer assigns an 'N' base.
}

\seealso{

}