\name{makePlot.gene.cdf}
\docType{methods}
\alias{makePlot.gene.cdf}
\title{
   Plot the cumulative gene diversity curve
}
\description{
   Plots the cumulative gene diversity curve
}
\usage{
  makePlot.gene.cdf(plotter, 
                sampleWise = FALSE,
                plot.intercepts = TRUE,
                label.intercepts = FALSE,
                debugMode, ...)
}
\arguments{
  \item{plotter}{
    A QoRT_Plotter object.
  }
  \item{sampleWise}{
    A logical value indicating whether to compile each sample
    together across all runs.
  }
  \item{plot.intercepts}{
    A logical value indicating whether or not to plot the
    intercepts with the round numbers on the x-axis, in 
    dotted lines.
  }
  \item{label.intercepts}{
    A logical value indicating whether or not to label the
    intercepts. No effect if plot.intercepts is FALSE.
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