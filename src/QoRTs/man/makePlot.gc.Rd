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
              plotRate = FALSE, byPair = FALSE,
              debugMode, singleEndMode, ...)
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
  \item{plotRate}{
    A logical value indicating whether or not the X-axis should be 
    the raw number of nucleotides that are G/C, vs the rate G/C.
  }
  \item{byPair}{
    Logical. If FALSE, GC content rates will be calculated for
    all reads individually. If TRUE, GC content rates will be
    calculated for all read-pairs. Note that when the insert
    size is small, the paired plot can sometimes appear jagged,
    as completely-overlapped read-pairs will always have an even
    number of G/C nucleotides.
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
  x-axis: Percent of the read-pairs that is made up of G's or C's.
  
  y-axis: Rate at which the given percent occurs.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.gc(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}