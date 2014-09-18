\name{makePlot.genebody.coverage}
\docType{methods}
\alias{makePlot.genebody.coverage}
\alias{makePlot.genebody.coverage.UMQuartile}
\alias{makePlot.genebody.coverage.lowExpress}
\title{
   Plot Gene-Body coverage distribution
}
\description{
   Plots Gene-Body coverage distribution
}
\usage{
  makePlot.genebody.coverage(plotter, plot.medians,
                       plot.means = TRUE,  
                       debugMode, ...)
  makePlot.genebody.coverage.UMQuartile(plotter, plot.medians, 
                       plot.means = TRUE, 
                       debugMode, ...)
  makePlot.genebody.coverage.lowExpress(plotter, plot.medians, 
                       plot.means = TRUE, 
                       debugMode, ...)
}
\arguments{
  \item{plotter}{
    A QoRTs_Plotter object.
  }
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average for each bamfile at the bottom of the plot. 
    Overrides \code{plot.means}.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average for each bamfile at the bottom of the plot.
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
  x-axis: Gene-body quantile. By default this is broken up into
  40 quantiles containing 2.5% each.
  
  y-axis: Rate at which reads falls into the given quantile of the
  genes' bodies.
  
  INCOMPLETE!!!!!
}

\seealso{

}