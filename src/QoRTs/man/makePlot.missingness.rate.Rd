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
                        debugMode, singleEndMode,
                        rasterize.plotting.area = FALSE, 
                        raster.height = 1000, 
                        raster.width = 1000,
                        log.y = FALSE, 
                        ylim = NULL,
                        plot = TRUE,
                        ...)
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
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{log.y}{
      Logical. If TRUE, the y-axis will be log-scale.
  }
  \item{ylim}{
      Numeric of length 2. Sets the y-axis limits.
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