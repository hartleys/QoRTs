\name{makePlot.splice.junction.event.rates}
\docType{methods}
\alias{makePlot.splice.junction.event}
\alias{makePlot.splice.junction.event.ratesPerRead}
\alias{makePlot.splice.junction.event.counts}
\alias{makePlot.splice.junction.event.proportions}
\alias{makePlot.splice.junction.event.proportionsByType}
\title{
   Plot Splice Junction Event Rates
}
\description{
   Plots the rates at which splice junctions occur.
}
\usage{
   makePlot.splice.junction.event.counts(plotter,
                                    high.low.cutoff = 4, 
                                    debugMode, 
                                    ...)
   makePlot.splice.junction.event.ratesPerRead(plotter,
                                    high.low.cutoff = 4, 
                                    debugMode, 
                                    ...)   
   makePlot.splice.junction.event.proportions(plotter,
                                    high.low.cutoff = 4, 
                                    debugMode, 
                                    ...) 
   makePlot.splice.junction.event.proportionsByType(plotter, 
                                    high.low.cutoff = 4, 
                                    debugMode, 
                                    ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{high.low.cutoff}{
    For Advanced users only! The cutoff between high and low expression splice junctions. Note that in order to function, this same cutoff MUST be used by the QoRTs jar utility that generates these counts.
  }
  \item{debugMode}{
    Activates debug mode, which causes more verbose reporting.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   These functions plot various metrics for the rate at which splice junction "events" occur. A splice junction "event" is an occurrance of a mapped read-pair bridging a splice junction. Some read-pairs may contain multiple splice junction events, and some read-pairs may contain none.
   
   Splice junctions are characterized into six categories:
      \itemize{
        \item Known: The splice junction locus is found in the supplied transcript annotation gtf file.
        \item Novel: The splice junction locus is NOT found in the supplied transcript annotation gtf file.
        \item Known, 1-3 reads: The locus is known, and is only covered by 1-3 read-pairs.
        \item Known, 4+ reads: The locus is known, and is covered by 4 or more read-pairs.
        \item Novel, 1-3 reads: The locus is novel, and is only covered by 1-3 read-pairs.
        \item Novel, 4+ reads: The locus is novel, and is covered by 4 or more read-pairs.
      }
    
   \code{makePlot.splice.junction.event.counts} plots the number (y-axis) of all splice junction events falling into each of six categories. Note that because different samples/runs may have different total read counts and/or library sizes, this function is generally not the best for comparing between samples. For most purposes, \code{makePlot.splice.junction.event.ratesPerRead} will be preferable.
   
   \code{makePlot.splice.junction.event.ratesPerRead} plots the rate (y-axis) at which each type of splice junction events appear, per read-pair. 
   
   \code{makePlot.splice.junction.event.proportions} plots the proportion (y-axis) of all splice junction events falling into the six categories.
   
   \code{makePlot.splice.junction.event.proportionsByType} plots the proportion (y-axis) at which splice junction events appear on known vs novel splice junction loci, the proportion of known splice junction events that occur on low-coverage junctions (1-3 read-pairs) vs high-coverage junctions (4 or more read-pairs), and the proportion of novel splice junction events that occur on low vs high coverage junctions.
   
   All of these plots are generally used to detect whether sample-specific or batch effects have a substantial or biased effect on splice junction appearance, either due to differences in the original RNA, or due to artifacts that alter the rate at which the aligner maps across splice junctions.
}

\seealso{
  \code{\link{build.plotter}}, \code{\link{makePlot.splice.junction.loci.counts}}
}