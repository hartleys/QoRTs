\name{makePlot.cigarOp.byCycle}
\docType{methods}
\alias{makePlot.cigarOp.byCycle}
\title{
   Plot Cigar Operator Rate
}
\description{
   Plots the rate at which the given CIGAR operator occurs, by
   read cycle.
}
\usage{
  makePlot.cigarOp.byCycle(plotter, op, 
                       r2.buffer = NULL, 
                       rate.per.million = TRUE, 
                       debugMode,
                       ...)
}
\arguments{
  \item{plotter}{
    A QoRTs_Plotter object.
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
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
  }
  \item{rate.per.million}{
    A logical value indicating whether or not to scale the y 
    axis to rate-per-million-reads, rather than rate-per-read.
    Some people may find the results more readable this way, even
    though the plots themselves will appear the same.
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
  x-axis: The read cycle (ie. the base-pair position in the read).
  
  y-axis: The rate at which the bases at the given read-cycle is
  clipped off.
}

\seealso{

}