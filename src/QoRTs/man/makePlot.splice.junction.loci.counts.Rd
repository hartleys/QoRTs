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
                                 debugMode, singleEndMode, 
                                 plot = TRUE, \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
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
   This function plots the number (y-axis) of splice junction \emph{loci} of each type that appear in the bam-file's reads. Splice junctions are split into 4 groups, first by whether the splice junction appears in the transcript annotation gtf ("known" vs "novel"), and then by whether the splice junction has 4 or more reads covering it, or 1-3 reads ("Hi" vs "Lo").
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.splice.junction.loci.counts(plotter);
  #Add a legend:
  makePlot.legend.over("topright", plotter)
}
\seealso{
  \code{\link{build.plotter}}, \code{\link{makePlot.splice.junction.event}}
}