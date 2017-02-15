\name{makePlot.readLengthDist}
\docType{methods}
\alias{makePlot.readLengthDist}
\title{
   Plot the distribution of read lengths.
}
\description{
   Plots the distribution of read lengths. Only useful for data with variable trimming (which is generally not recommended at least for RNA-Seq data).
}
\usage{
makePlot.readLengthDist(plotter,
                        plot.rates = TRUE, 
                        plot.means = TRUE, 
                        plot.medians = NULL,
                        include.full.length = FALSE, 
                        cumulative = TRUE,
                        singleEndMode,
                        rasterize.plotting.area = FALSE, 
                        raster.height = 1000, 
                        raster.width = 1000,
                        debugMode,
                        r2.buffer,
                        plot = TRUE,
                        ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.rates}{
    A logical value indicating whether or not the X-axis should be 
    the raw number of nucleotides that are G/C, vs the rate G/C.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average GC content for each bam file at the bottom 
    of the plot.
  }
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average GC content for each bam file at the bottom 
    of the plot. Overrides \code{plot.means}.
  }
  \item{include.full.length}{
    Logical. If FALSE, omit the full-length read length from the x-axis of the plot.
  }
  \item{cumulative}{
    Logical. If TRUE, plot shows cumulative rates.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
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
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
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
  x-axis: Read Length
  
  y-axis: Percentage of reads with length equal to the given length. If cumulative == TRUE, then it is the percentage of reads with length less than or equal to the given length.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.readLengthDist(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}