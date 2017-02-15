\name{makePlot.NVC.lead.clip}
\docType{methods}
\alias{makePlot.NVC.lead.clip}
\alias{makePlot.NVC.tail.clip}
\title{
   Clipped sequence nucleotide-by-position plot
}
\description{
   Plots the nucleotide rates for clipped segments
}
\usage{
makePlot.NVC.lead.clip(plotter, clip.amt,  r2.buffer, 
                   points.highlighted = TRUE,
                   label.majority.bases = TRUE, 
                   label.majority.bases.threshold = 0.5, 
                   label.majority.bases.cex = 1, 
                   rasterize.plotting.area = FALSE, 
                   raster.height = 1000, 
                   raster.width = 1000,
                   show.base.legend = TRUE, 
                   debugMode, singleEndMode,
                   plot = TRUE,
                   \dots)


makePlot.NVC.tail.clip(plotter, clip.amt,  r2.buffer, 
                       points.highlighted = TRUE, 
                       label.majority.bases = TRUE, 
                       label.majority.bases.threshold = 0.5, 
                       label.majority.bases.cex = 1, 
                       rasterize.plotting.area = FALSE, 
                       raster.height = 1000, 
                       raster.width = 1000,
                       show.base.legend = TRUE, 
                       debugMode, singleEndMode,
                       plot = TRUE,
                       \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{clip.amt}{
    Numeric value. The number of bases clipped. These methods
    only plot the sequence for a single specific clip.amt. 
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
  }
  \item{points.highlighted}{
    A logical value. Determines whether to ever plot points on top of the lines. This
    can be useful for identifying outliers. If TRUE, then all highlighted lane-bams
    will be overlaid with points using their designated pch symbol. If the plotter does
    not highlight any lanebams, then points will be overlaid on ALL lines.
  }
  \item{label.majority.bases}{
    A logical value indicating whether to label the genotypes of
    read cycles in which the most common base has a frequency 
    exceeding label.majority.bases.threshold (see below).
  }
  \item{label.majority.bases.threshold}{
    A numeric value indicating the cutoff above which the most 
    common base will be labelled, if label.majority.bases is set
    TRUE (see above).
  }
  \item{label.majority.bases.cex}{
    The cex value fed to text() that is used to determine how
    big the labels are to be, if label.majority.bases is TRUE.
    (see \code{\link{par}} for information on cex).
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
  \item{show.base.legend}{
    A logical value indicating whether to print the base-color
    legend along the right edge of the plot.
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
   For general information on reading NVC plots, see \code{\link{makePlot.NVC}}.
   
   This function plots the nucleotide rates for the sections of reads that were soft-clipped by the aligner. Note that these will not function on reads that have been aligned using an aligner that does not practice soft clipping (such as TopHat2).
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.NVC.lead.clip(plotter, clip.amt = 12);
  makePlot.NVC.tail.clip(plotter, clip.amt = 12);
}
\seealso{
  \code{\link{build.plotter}}, \code{\link{makePlot.NVC}}
}