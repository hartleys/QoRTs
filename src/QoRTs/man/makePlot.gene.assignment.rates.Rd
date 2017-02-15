\name{makePlot.gene.assignment.rates}
\docType{methods}
\alias{makePlot.gene.assignment.rates}
\title{
   Gene assignment rate plot
}
\description{
   Plots the rate at which reads are assigned into various 
   categories.
}
\usage{
makePlot.gene.assignment.rates(plotter,  debugMode, singleEndMode, 
                               plot = TRUE, \dots)
                         
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
  For each bam-file, this function plots the rate (y-axis) for which the bam-file's read-pairs are assigned to the given categories.
  
The categories are:
  \itemize{
    \item \emph{Unique Gene}: The read-pair overlaps with the exonic segments of one and only one gene. For many downstream analyses tools, such as DESeq, DESeq2, and EdgeR, only read-pairs in this category are used.
    \item \emph{Ambig Gene}: The read-pair overlaps with the exons of more than one gene.
    \item \emph{No Gene}: The read-pair does not overlap with the exons of any annotated gene.
    \item \emph{No Gene, Intronic}: The read-pair does not overlap with the exons of any annotated gene, but appears in a region that is bridged by an annotated splice junction.
    \item \emph{No Gene, 1kb from gene}: The read-pair does not overlap with the exons of any annotated gene, but is within 1 kilobase from the nearest annotated gene.
    \item \emph{No Gene, 10kb from gene}: The read-pair does not overlap with the exons of any annotated gene, but is within 10 kilobases from the nearest annotated gene.
    \item \emph{No Gene, middle of nowhere}: The read-pair does not overlap with the exons of any annotated gene, and is \emph{more} than 10 kilobases from the nearest annotated gene.
  }
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.gene.assignment.rates(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}