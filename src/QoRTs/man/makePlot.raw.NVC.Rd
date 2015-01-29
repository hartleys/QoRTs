\name{makePlot.raw.NVC}
\docType{methods}
\alias{makePlot.NVC}
\alias{makePlot.raw.NVC}
\alias{makePlot.minus.clipping.NVC}
\title{
   Nucleotide-rate by read cycle plot
}
\description{
   Plots the nucleotide rate by read cycle
}
\usage{
makePlot.raw.NVC(plotter,  r2.buffer = NULL, 
                         points.highlighted = TRUE,
                         label.majority.bases = FALSE, 
                         label.majority.bases.threshold = 0.5, 
                         label.majority.bases.cex = 0.5, 
                         rasterize.plotting.area = FALSE, 
                         raster.height = 1000, 
                         raster.width = 2000,
                         show.base.legend = TRUE, 
                         debugMode, singleEndMode,
                         \dots)
                         
makePlot.minus.clipping.NVC(plotter,  r2.buffer = NULL, 
                         points.highlighted = TRUE,
                         label.majority.bases = FALSE, 
                         label.majority.bases.threshold = 0.5, 
                         label.majority.bases.cex = 0.5, 
                         rasterize.plotting.area = FALSE, 
                         raster.height = 1000, 
                         raster.width = 2000,
                         show.base.legend = TRUE, 
                         debugMode, singleEndMode,
                         \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
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
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   This plot displays the nucleotide rates by read cycle. The color scheme for NVC plots is different from all other plots. Rather than being used for emphasis or to allow cross-comparisons by sample, biological-condition, or lane, the colors are used to indicate the four nucleotides: A (green), T (red), G (orange), or C (blue). Depending on the type of QoRT_Plotter being used, sample-runs will be marked and differentiated by marking the lines with shapes (R points). In many cases the points will be unreadable due to overplotting, but clear outliers that stray from the general trends can be readily identified.
   
   When used with a "sample.highlight" type QoRT_Plotter (see \code{\link{build.plotter}}), "highlighted" samples will be drawn with a deeper shade of the given color.
   
   \code{\link{makePlot.raw.NVC}} plots the nucleotides of all cycles of all aligned reads. \code{\link{makePlot.minus.clipping.NVC}} plots the nucleotides for all cycles that are NOT soft-clipped by the aligner. Note that for reads aligned with an aligner that does not perform soft-clipping (such as TopHat2), the two plots will be identical.
}

\seealso{
  \code{\link{build.plotter}}
}