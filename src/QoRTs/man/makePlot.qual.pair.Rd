\name{makePlot.qual.pair}
\docType{methods}
\alias{makePlot.qual.pair}
\title{
   Plot quality score by read cycle
}
\description{
   Plots the Phred quality score as a function of the read cycle
   for both reads.
}
\usage{
  makePlot.qual.pair(plotter, y.name, r2.buffer = NULL, 
                 debugMode, singleEndMode, ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{y.name}{
    The name of the quality score metric to plot. Must be one of:
      \itemize{
        \item \code{"min"}
        \item \code{"lowerQuartile"} 
        \item \code{"median"} 
        \item \code{"upperQuartile"} 
        \item \code{"max"}
      }
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
These plots display information about the phred quality score 
(y-axis) as a function of the position in the read (x-axis). Five 
statistics can be plotted: minimum, maximum, upper and lower 
quartiles, and median. These statistics are calculated individually
for each bam file and each read position (ie, each plotted line 
corresponds to a bam file).

Note that the Phred score is always an integer, and as such these 
plots would normally be very difficult to read because lines 
would be plotted directly on top of one another. To reduce this 
problem, the lines are vertically offset from one another. Most
plotters offset each line by lane.ID.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.qual.pair(plotter,y.name="min");
  makePlot.qual.pair(plotter,y.name="median");
  makePlot.qual.pair(plotter,y.name="max");
}
\seealso{
  \code{\link{build.plotter}}
}