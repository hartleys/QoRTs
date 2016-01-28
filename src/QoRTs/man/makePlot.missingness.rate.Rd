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
                        debugMode, singleEndMode,...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1 and read 2.
    By default this will choose a reasonable value.
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
  x-axis: Read Cycle
  
  y-axis: Rate at which the sequencer assigns an 'N' base.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.missingness.rate(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}