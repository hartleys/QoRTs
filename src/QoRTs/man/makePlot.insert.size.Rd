\name{makePlot.insert.size}
\docType{methods}
\alias{makePlot.insert.size}
\title{
   Plot Insert Size Distribution.
}
\description{
   Plots the distribution of the size of the region covered by the
   two paired reads.
}
\usage{
  makePlot.insert.size(plotter, calc.rate = TRUE, pct.cutoff = 0.98,
                   plot.medians = TRUE, plot.means = NULL,
                   debugMode, singleEndMode, xlim = NULL,
                   rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                   \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{calc.rate}{
    A logical value indicating whether or not to calculate and
    plot the rate at which each insert size occurs as the y-axis.
    If this is set to false, it will instead plot the total number
    of times each insert size occurs.
  }
  \item{pct.cutoff}{
    A numeric value between 0 and 1, indicating the quantile 
    within which to limit the x-axis. Generally the default value
    of 0.98 is perfectly usable.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average for each bamfile at the bottom of the plot.
  }
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average for each bamfile at the bottom of the plot. 
    Overrides \code{plot.means}.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{xlim}{
      Numeric vector of length 2. Sets the x-axis limits. By default
      QoRTs will attempt to autodetect reasonable values for this, but
      there are always cases where you want to set this explicitly.
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
  
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: The insert size. This is the genomic distance from the 
  start of one read to the end of the other. In other words, the
  size of the full region covered by the paired reads.
  
  y-axis: The rate at which the given insert size occurs, for each
  bamfile.
  
  Note: The insert size is calculated by using the alignment  
  as well as the provided gene/splice-junction annotation.
  \enumerate{
    \item If the paired reads align to overlapping regions of the
    reference genome, then the insert size can be calculated
    exactly. This is NOT dependant on the annotation. If the two 
    reads align inconsistantly (for example, one read showing a 
    splice junction and the other not), then the read is ignored.
    \item If the paired reads do NOT overlap, then the annotation
    information is required. Using the full set of all known splice
    junctions that lie between the inner alignment endpoints for 
    the two reads, the shortest possible path from the two 
    endpoints is found. 
    In some rare cases this may underestimate the insert size, if
    the actual insert does not take the minimal path, but this is
    rare.
    In other (somewhat more common) cases this may overestimate 
    the insert size, when the region between the paired reads 
    bridges an unannotated splice junction.
  }
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.insert.size(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}