\name{makePlot.runTimePerformance}
\docType{methods}
\alias{makePlot.runTimePerformance}
\title{
   Chromosome type rate plot
}
\description{
   Plots the number or percent of read-pairs falling on each type of chromosome.
}
\usage{
  makePlot.runTimePerformance(plotter,
                 debugMode, singleEndMode,
                 plot = TRUE,
                 ...)

}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
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
   For each sample run, indicates the amount of time spent running the QoRTs QC data processing tool
}
\value{
   By default, this function returns nothing. If the return.table parameter is TRUE, then it returns a data.frame with the data that was plotted.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.runTimePerformance(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}
