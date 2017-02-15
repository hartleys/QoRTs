\name{makePlot.NVC.clip.matchByClipPosition}
\docType{methods}
\alias{makePlot.NVC.clip.matchByClipPosition}
\alias{makePlot.NVC.lead.clip.matchByClipPosition}
\alias{makePlot.NVC.tail.clip.matchByClipPosition}
\title{
   Plot clipped nucleotide rates, aligned by the distance 
   from the point of clipping.
}
\description{
   WARNING: THESE FUNCTIONS ARE BETA, AND ARE NOT FULLY TESTED OR 
   READY FOR GENERAL USE.
}
\usage{
makePlot.NVC.lead.clip.matchByClipPosition(plotter, 
            clip.amt,  r2.buffer, 
            label.majority.bases = TRUE, 
            label.majority.bases.threshold = 0.5, 
            label.majority.bases.cex = 1, 
            rasterize.plotting.area = FALSE, 
            raster.height = 1000, 
            raster.width = 1000,
            show.base.legend = TRUE, 
            load.results = TRUE, 
            debugMode, singleEndMode, 
            plot = TRUE,
            ...)

makePlot.NVC.tail.clip.matchByClipPosition(plotter, 
            clip.amt, r2.buffer,
            label.majority.bases = TRUE, 
            label.majority.bases.threshold = 0.5, 
            label.majority.bases.cex = 1, 
            rasterize.plotting.area = FALSE, 
            raster.height = 1000, 
            raster.width = 1000,
            show.base.legend = TRUE, 
            debugMode, singleEndMode, 
            plot = TRUE,
            ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{clip.amt}{
    UNDOCUMENTED
  }
  \item{r2.buffer}{
    UNDOCUMENTED
  }
  \item{label.majority.bases}{
    UNDOCUMENTED
  }
  \item{label.majority.bases.threshold}{
    UNDOCUMENTED
  }
  \item{label.majority.bases.cex}{
    UNDOCUMENTED
  }
  \item{show.base.legend}{
    UNDOCUMENTED
  }
  \item{rasterize.plotting.area}{
    Logical. If TRUE, the plotting area will be written to a temp
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
  \item{load.results}{
    UNDOCUMENTED
  }
  \item{debugMode}{
    UNDOCUMENTED
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
   This function is currently BETA, and is not intended for general use. Documentation and testing is still pending.
}

\seealso{
  \code{\link{build.plotter}}, \code{\link{makePlot.NVC}}
}