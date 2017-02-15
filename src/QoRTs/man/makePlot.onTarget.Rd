\name{makePlot.onTarget}
\docType{methods}
\alias{makePlot.targetDistribution}
\alias{makePlot.onTarget.counts}
\alias{makePlot.onTarget.rates}

\title{
   On-Target Rates
}
\description{
   On-Target Rates.
}
\usage{
makePlot.targetDistribution(plotter,
       plot.rates = TRUE, 
       plot.hist = TRUE, 
       log.y = TRUE, 
       singleEndMode,
       byPair = !singleEndMode, 
       rasterize.plotting.area = FALSE, 
       raster.height = 1000, 
       raster.width = 1000,
       debugMode,
       plot = TRUE,
       ...)

makePlot.onTarget.counts(plotter, 
      y.counts.in.millions = TRUE, 
      debugMode, 
      singleEndMode,
      plot = TRUE, ...)

makePlot.onTarget.rates(plotter,  
      debugMode, 
      singleEndMode,
      plot = TRUE, ...)



}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{y.counts.in.millions}{
    Logical. If TRUE, y axis labels are labelled in millions.
  }
  \item{plot.rates}{
    Logical. If TRUE, plot rates, otherwise plot raw counts.
  }
  \item{plot.hist}{
    Plot a histogram.
  }
  \item{log.y}{
    Logical. If TRUE, the y-axis should be log-scaled.
  }
  \item{byPair}{
    Logical. If TRUE, count read-pairs rather than reads. In other words, if two reads overlap, the overlapping region is counted only once.
  }

  \item{rasterize.plotting.area}{
    Logical. If \code{TRUE}, the plotting area will be written to a temp
    png file then drawn to the current device as a raster image.
    This option is generally preferred for vector devices, because
    NVC plots can be very large when drawn in vector format. 
    \emph{Note: this requires the \code{png} package!}
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels.
  }
  
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
      Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }

}
\details{
   For each sample run, indicates the amount of time spent running the QoRTs QC data processing tool
}
\value{
   By default, this function returns nothing. If the return.table parameter is TRUE, then it returns a data.frame with the data that was plotted.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.runTimePerformance(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}
