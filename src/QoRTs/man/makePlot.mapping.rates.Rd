\name{makePlot.mapping.rates}
\docType{methods}
\alias{makePlot.mapping.rates}
\title{
   Plot mapping rates
}
\description{
   Plots the rates at which reads map to the genome.
}
\usage{
makePlot.mapping.rates(plotter, plot.mm = NULL, 
                   y.counts.in.millions = TRUE, 
                   debugMode, singleEndMode,
                   \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.mm}{
    Plot multi-mapping rates. By default, it autodetects whether multimapping data was included in the original decoder.
  }
  \item{y.counts.in.millions}{
    Label/scale the y-axis in millions of reads.
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
For each sample-run, this function plots the number of mapped reads and the rate at which reads map (if the total reads is provided).
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.mapping.rates(plotter);
}
\seealso{
  \code{\link{build.plotter}}
  \code{\link{read.qc.results.data}}
}