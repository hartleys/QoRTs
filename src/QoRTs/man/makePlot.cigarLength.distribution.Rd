\name{makePlot.cigarLength.distribution}
\docType{methods}
\alias{makePlot.cigarLength.distribution}
\title{
   Plot the length distribution of a given cigar operation
}
\description{
   Plots the length distribution of a given cigar operation
}
\usage{
  makePlot.cigarLength.distribution(plotter, op, 
                                r2.buffer = NULL,
                                perMillion = TRUE,
                                log.x = FALSE, 
                                log.y = FALSE, 
                                debugMode, singleEndMode,
                                rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                plot = TRUE,
                                \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{op}{
    A character string naming which cigar op to plot.
    Must be one of the following:
      \itemize{
        \item \code{"SoftClip"} Soft Clip (Cigar Op "S")
        \item \code{"HardClip"} Hard Clip (Cigar Op "H")
        \item \code{"Del"} Deletion from reference (Cigar Op "D")
        \item \code{"Ins"} Insertion from reference (Cigar Op "I")
        \item \code{"Pad"} Padding (Cigar Op "P")
        \item \code{"Splice"} Splice Junction (Cigar Op "N")
        \item \code{"Aln"} Aligned to reference (Cigar Op "M")
      }
    Note that makePlot.cigarOp.byCycle(plotter,"SoftClip") gives
    identical results as makePlot.clipping, although the default
    value for the \code{rate.per.million} argument is different.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1 and read 2.
    By default this will choose a reasonable value.
  }
  \item{perMillion}{
    A logical value indicating whether or not to scale the y 
    axis to rate-per-million-reads, rather than rate-per-read.
    Some people may find the results more readable this way, even
    though the plots themselves will appear the same.
  }
  \item{log.x}{
    A logical value indicating whether or not to log-scale the x 
    axis.
  }
  \item{log.y}{
    A logical value indicating whether or not to log-scale the y 
    axis.
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
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: Length of the cigar operation.
  
  y-axis: Frequency of the given length.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.cigarLength.distribution(plotter, op = "Del");
  makePlot.cigarLength.distribution(plotter, op = "Ins");
}
\seealso{
  \code{\link{build.plotter}}
}