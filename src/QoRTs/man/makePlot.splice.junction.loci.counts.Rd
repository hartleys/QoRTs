\name{makePlot.splice.junction.loci.counts}
\docType{methods}
\alias{makePlot.splice.junction.loci.counts}
\title{
   Splice Junction Loci Count Plot
}
\description{
   Plots the rate at which splice junction loci
   fall into various categories.
}
\usage{
makePlot.splice.junction.loci.counts(plotter, 
                                 calc.rate = FALSE, 
                                 high.low.cutoff = 4,
                                 debugMode, \dots)
}
\arguments{
  \item{plotter}{
    A QoRT_Plotter object.
  }
  \item{calc.rate}{
    Logical. If TRUE, the proportion of loci in each category will be calculated and plotted, rather than the raw number.
  }
  \item{high.low.cutoff}{
    (For advanced users only!) The cutoff between "high" and "low" expression junction loci. Note that this must match the cutoff used by the jarfile QC execution.
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