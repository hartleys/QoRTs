\name{makePlot.overlap}
\docType{methods}
\alias{makePlot.overlap}
\alias{makePlot.overlap.coverage}
\alias{makePlot.overlapMismatch.perRead}
\alias{makePlot.overlapMismatch.size}
\alias{makePlot.overlapMismatch.byQual.min}
\alias{makePlot.overlapMismatch.byQual.read}
\alias{makePlot.overlapMismatch.byBase}
\alias{makePlot.overlapMismatch.byBase.atScore}
\alias{makePlot.overlapMismatch.byCycle}
\title{
   Overlap Mismatch Rates
}
\description{
   For paired-end data only, plots the rate at which the two reads have point-mismatch on overlapping regions.
}
\usage{
  makePlot.overlap.coverage(plotter,
                            plot.rates = TRUE, 
                            singleEndMode,
                            rasterize.plotting.area = FALSE, 
                            raster.height = 1000, 
                            raster.width = 1000,
                            debugMode,
                            r2.buffer=NULL,
                            plot = TRUE,
                            ...)

  makePlot.overlapMismatch.size(plotter,
                                plot.rates = TRUE,
                                log.y = TRUE,
                                noIndel = TRUE,
                                draw.decade.lines = c(TRUE,TRUE),
                                xlim = NULL,
                                pct.cutoff=0.95,
                                singleEndMode,
                                debugMode,
                                rasterize.plotting.area = FALSE, 
                                raster.height = 1000, 
                                raster.width = 1000,
                                plot = TRUE,
                                ...)
                                
  makePlot.overlapMismatch.byQual.min(plotter,
                                      plot.rates = TRUE,
                                      log.y = TRUE,
                                      noIndel = TRUE,
                                      draw.decade.lines = TRUE,
                                      singleEndMode,
                                      debugMode,
                                      rasterize.plotting.area = FALSE, 
                                      raster.height = 1000, 
                                      raster.width = 1000,
                                      plot = TRUE,
                                      ...)

  makePlot.overlapMismatch.byQual.read(plotter,
                                       plot.rates = TRUE,
                                       noIndel = TRUE,
                                       log.y = TRUE,
                                       draw.decade.lines = TRUE,
                                       singleEndMode,
                                       debugMode,
                                       rasterize.plotting.area = FALSE, 
                                       raster.height = 1000, 
                                       raster.width = 1000,
                                       r2.buffer = NULL,
                                       plot = TRUE,
                                       ...)

  makePlot.overlapMismatch.byBase(plotter, 
                                  noIndel = TRUE, 
                                  plot.rates = TRUE,
                                  log.y = FALSE, 
                                  y.rate.per.kb = FALSE, 
                                  debugMode, 
                                  draw.vert.dotted.lines = TRUE,
                                  singleEndMode, 
                                  plot = TRUE,
                                  ...)
  
  makePlot.overlapMismatch.byBase.atScore(plotter,
                                     atScore = 41,
                                     overlapScoreMethod = c("pair","min"),
                                     plot.rates = TRUE,
                                     noIndel = TRUE,
                                     log.y = FALSE,
                                     singleEndMode,
                                     rasterize.plotting.area = FALSE,
                                     raster.height = 1000,
                                     raster.width = 1000,
                                     debugMode = DEFAULTDEBUGMODE,
                                     plot = TRUE,
                                     ...)
  
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{atScore}{
    Numeric integer. For makePlot.overlapMismatch.byBase.atScore, plot the overlap mismatch rates for each base-pair swap at the given phred score.
    
    If \code{overlapScoreMethod} is "min", then this will plot the rate of overlap mismatches where the minimum of the two reads' quality score is \code{atScore}. If \code{overlapScoreMethod} is "pair", then 
    this will plot the rate of overlap mismatches where both reads have the same quality score, equal to \code{atScore}.
  }
  \item{overlapScoreMethod}{
    Character string. See documentation for atScore, above.
  }
  \item{plot.rates}{
    A logical value indicating whether or not to plot mismatch rates or mismatch counts.
  }
  \item{log.y}{
    Logical. If TRUE, the y-axis will be log-scaled.
  }
  \item{noIndel}{
    Logical. If TRUE, only count overlapping reads if both reads have no aligned indels.
  }
  \item{draw.decade.lines}{
    Logical. If TRUE, draw mini tick lines at decade points on the y-axis when plotting in log-scale.
  }
  \item{xlim}{
    Numeric. The x-axis limits. If NULL, the x-limits will be autodetected using the pct.cutoff value.
  }
  \item{pct.cutoff}{
    Numeric. The percentile cutoff value for the x-axis upper limit.
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
  \item{r2.buffer}{
      Numeric. The desired distance between the read1/read2 sub-plots.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
 
  
  \item{y.rate.per.kb}{
    Logical. If TRUE, the y axis should be labelled in overlap rate per kilobase.
  }

  \item{draw.vert.dotted.lines}{
    Logical. If TRUE, draw dotted lines at decade y axis points.
  }
  
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  These plotting functions create plots that summarize the overlap-mismatch between paired end reads. 
  When paired-end reads overlap, the sequence on the two reads can be compared in order to estimate 
  the sequencer error rate. This error rate can be visualized in a variety of ways to detect sequencer artifacts
  or errors.
  
  makePlot.overlap.coverage: The rates at which each position in read 1 and read 2 are covered by the other read.
     x-axis: Read position.
     y-axis: How often the given position on the given read is covered by the other read.
  
  makePlot.overlap.mismatchCountPerRead: A histogram of the amount of mismatch between the two reads.
     x-axis: Number of mismatches.
     y-axis: How often a read is observed with the given number of mismatches.
  
  makePlot.overlap.mismatch.byMinQual: The mismatch rate as a function of the minimum of the two base quality scores.
     x-axis: Minimum quality score
     y-axis: Mismatch rate.
  
  makePlot.overlap.mismatch.byAvgQual: The mismatch rate as a function of the average of the two base quality scores.
     x-axis: Average quality score
     y-axis: Mismatch rate.
  
  makePlot.overlap.mismatch.byReadQual: The mismatch rate as a function of each read's quality score.
     x-axis: Quality score for read 1 and read 2.
     y-axis: Mismatch rate.
  
  makePlot.overlap.mismatch.byBase: The rate at which overlap-mismatches occur for each possible base-pair swap.
     x-axis: Read1/read2 base-pair values (note: read2 base is complemented).
     y-axis: Mismatch rate.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.overlap.coverage(plotter)
  makePlot.overlapMismatch.size(plotter)
  makePlot.overlapMismatch.byQual.min(plotter)
  makePlot.overlapMismatch.byQual.read(plotter)
  makePlot.overlapMismatch.byBase(plotter)
  makePlot.overlapMismatch.byBase.atScore(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}