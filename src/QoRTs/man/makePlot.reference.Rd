\name{makePlot.reference}
\docType{methods}

\alias{makePlot.reference}
\alias{makePlot.referenceMismatch}
\alias{makePlot.referenceMismatch.byScore}
\alias{makePlot.referenceMismatch.perRead}
\alias{makePlot.referenceMismatch.byBase}
\alias{makePlot.referenceMismatch.byBase.atScore}
\alias{makePlot.referenceMismatch.byCycle}

\title{
   Plot Reference Mismatch Rates
}
\description{
   Plots various rates of single-base mismatches against the reference genome.
}
\usage{
makePlot.referenceMismatch.byScore(plotter,
             plot.rates = TRUE,
             draw.decade.lines = c(FALSE,TRUE),
             log.y = TRUE,
             singleEndMode = plotter$res@singleEnd,
             debugMode = DEFAULTDEBUGMODE,
             rasterize.plotting.area = FALSE, 
             raster.height = 1000, 
             raster.width = 1000,
             r2.buffer = NULL,
             plot = TRUE,
             ...)

makePlot.referenceMismatch.byBase(plotter, 
             y.rate.per.kb = TRUE, draw.vert.dotted.lines = TRUE,
             debugMode = DEFAULTDEBUGMODE, 
             singleEndMode = plotter$res@singleEnd, 
             plot = TRUE,
             ...)

makePlot.referenceMismatch.byBase.atScore(plotter, 
             atScore = 41,
             forRead = c("BOTH","R1","R2"),
             plot.rates = TRUE,
             log.y = FALSE,
             singleEndMode = plotter$res@singleEnd,
             rasterize.plotting.area = FALSE, 
             raster.height = 1000, 
             raster.width = 1000,
             debugMode = DEFAULTDEBUGMODE,
             plot = TRUE,...)

makePlot.referenceMismatch.byCycle(plotter,
             plot.rates = TRUE,
             log.y = TRUE,
             ylim = NULL,
             singleEndMode = plotter$res@singleEnd,
             debugMode = DEFAULTDEBUGMODE,
             rasterize.plotting.area = FALSE, 
             raster.height = 1000, 
             raster.width = 1000,
             r2.buffer = NULL,
             plot = TRUE,
             ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.rates}{
    A logical value indicating whether or not to plot mismatch rates or mismatch counts.
  }
  \item{log.y}{
    Logical. If TRUE, the y-axis will be log-scaled.
  }
  \item{draw.decade.lines}{
    Logical. If TRUE, draw mini tick lines at decade points on the y-axis when plotting in log-scale.
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
  \item{y.rate.per.kb}{
    Logical. If TRUE, then y-axis should be labelled as the rate per kilobase.
  }
  \item{draw.vert.dotted.lines}{
    Logical. If TRUE, then draw dotted guide-lines at reasonable intervals.
  }
  \item{atScore}{
    Numeric integer. For \code{makePlot.referenceMismatch.byBase.atScore}, this sets the PHRED score.
  }
  \item{forRead}{
    Character string. For \code{makePlot.referenceMismatch.byBase.atScore}, this sets whether the plot should be for read 1, read 2, or both side by side.
  }
  \item{ylim}{
    Numeric vector of length 2. This sets the y-axis limits.
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
  These plotting functions create plots that summarize the rate of single-base-pair reference-mismatches of reads against the reference genome. A "reference mismatch"
  is defined as an aligned read base that does not match the genomic base.

  These reference mismatches may be caused by a number of different factors and may not necessarily indicate a serious quality issue. To begin with, many
  reference mismatches may be caused by real single nucleotide polymorphisms present in the subject. QoRTs produces several different plots that allow these 
  mismatch rates to be measured in a variety of ways.
  
  \emph{\code{makePlot.referenceMismatch.byCycle}}: Plots the overall reference mismatch rate as a function of read cycle.

  \emph{\code{makePlot.referenceMismatch.byBase}}: Plots the overall reference mismatch rate for each of the 12 possible X -> Y base-pair substitutions.
  
  \emph{\code{makePlot.referenceMismatch.byScore}}: Plots the reference mismatch rate as a function of the PHRED score.
  
  \emph{\code{makePlot.referenceMismatch.byBase.atScore}}: Plots the reference mismatch rate for each of the 12 possible X -> Y base-pair substitutions, for the
  subset of bases with PHRED score equal to \code{atScore}.

  \emph{NOTE:} Creation of the reference mismatch plots requires that the QoRTs QC run be executed with additional optional parameters. The BAM file must be sorted
  by coordinate, and a genome fasta file must be specified via the \code{--genomeFA} parameter.

}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  
}
\seealso{
  \code{\link{build.plotter}}
}