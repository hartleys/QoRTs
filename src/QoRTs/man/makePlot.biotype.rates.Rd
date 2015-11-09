\name{makePlot.biotype.rates}
\docType{methods}
\alias{makePlot.biotype.rates}
\title{
   Plot Biotype rates
}
\description{
   Plots counts for each gene biotype. This plot is only useful and informative when QoRTs is run on a GTF file that contains "gene_biotype" tags.
}
\usage{
  makePlot.biotype.rates(plotter, 
              plot.rates = TRUE,
              count.type = c("all","unambigOnly"),
              log.y = TRUE,
              return.table = FALSE, 
              debugMode = DEFAULTDEBUGMODE,  
              singleEndMode = plotter$res@singleEnd, ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.rates}{
    A logical value. If \code{TRUE}, then reads/read-pairs 
    per million should be plotted, rather than raw counts. 
    Default is \code{TRUE}.
  }
  \item{count.type}{
    A logical value. If \code{TRUE}, then both ambiguous and unambiguous
    reads will be counted. Otherwise only unambiguous reads will be counted.
  }
  \item{log.y}{
    A logical value indicating that the y-axis should be log-scaled.
  }
  \item{return.table}{
    Logical. If \code{TRUE}, then return a data table.
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
  x-axis: Gene "BioTypes", from the annotation GTF.
  
  y-axis: Total read counts across all genes in the BioType, 
  measured either in reads or reads-per-million.
}

\seealso{
  \code{\link{build.plotter}}
}