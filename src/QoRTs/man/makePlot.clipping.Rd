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
                debugMode, singleEndMode,
                ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
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
  x-axis: The read cycle (ie. the base-pair position in the read).
  
  y-axis: The rate at which the bases at the given read-cycle is
  clipped off.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.clipping(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}