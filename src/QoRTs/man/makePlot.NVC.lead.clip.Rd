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
                   show.base.legend = TRUE, 
                   debugMode,
                   \dots)


makePlot.NVC.tail.clip(plotter, clip.amt,  r2.buffer, 
                       points.highlighted = TRUE, 
                       label.majority.bases = TRUE, 
                       label.majority.bases.threshold = 0.5, 
                       label.majority.bases.cex = 1, 
                       show.base.legend = TRUE, 
                       debugMode, \dots)
}
\arguments{
  \item{plotter}{
    A QoRT_Plotter object.
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
  \item{show.base.legend}{
    A logical value indicating whether to print the base-color
    legend along the right edge of the plot.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   TODO!
}

\seealso{

}