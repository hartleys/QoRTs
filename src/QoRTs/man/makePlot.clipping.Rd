\name{makePlot.clipping}
\docType{methods}
\alias{makePlot.clipping}
\title{
   Plot Alignment Clipping
}
\description{
   Plots the rate at which the aligner soft-clips off portions of
   aligned reads.
}
\usage{
  makePlot.clipping(plotter, rate.per.million = FALSE, 
                r2.buffer = NULL, 
                debugMode,
                ...)
}
\arguments{
  \item{plotter}{
    A QoRTs_Plotter object.
  }
  \item{rate.per.million}{
    A logical value indicating whether or not to scale the y 
    axis to rate-per-million-reads, rather than rate-per-read.
    Some people may find the results more readable this way, even
    though the plots themselves will appear the same.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
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
  x-axis: The read cycle (ie. the base-pair position in the read).
  
  y-axis: The rate at which the bases at the given read-cycle is
  clipped off.
}

\seealso{

}