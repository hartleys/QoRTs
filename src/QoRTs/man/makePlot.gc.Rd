\name{makePlot.gc}
\docType{methods}
\alias{makePlot.gc}
\title{
   Plot GC content
}
\description{
   Plots GC content.
}
\usage{
  makePlot.gc(plotter, plot.medians = NULL, plot.means = TRUE, 
              debugMode, paired = FALSE, ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average GC content for each bam file at the bottom 
    of the plot. Overrides \code{plot.means}.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average GC content for each bam file at the bottom 
    of the plot.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{paired}{
    Logical. If FALSE, GC content rates will be calculated for
    all reads individually. If TRUE, GC content rates will be
    calculated for all read-pairs. Note that when the insert
    size is small, the paired plot can sometimes appear jagged,
    as completely-overlapped read-pairs will always have an even
    number of G/C nucleotides.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: Percent of the read-pairs that is made up of G's or C's.
  
  y-axis: Rate at which the given percent occurs.
}

\seealso{
  \code{\link{build.plotter}}
}