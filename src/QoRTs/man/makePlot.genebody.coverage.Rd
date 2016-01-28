\name{makePlot.genebody.coverage}
\docType{methods}
\alias{makePlot.genebody}
\alias{makePlot.genebody.coverage}
\alias{makePlot.genebody.coverage.UMQuartile}
\alias{makePlot.genebody.coverage.lowExpress}
\title{
   Plot Gene-Body coverage distribution
}
\description{
   Plots Gene-Body coverage distribution
}
\usage{
  makePlot.genebody(plotter, 
                  geneset = c("Overall","90-100","75-90","50-75","0-50"),
                  avgMethod = c("TotalCounts", "AvgPercentile"), 
                  plot.medians, 
                  plot.means = TRUE, 
                  debugMode, 
                  singleEndMode, 
                  ... )

  #NOTE: The preferred method to access the below functions
  #  is to use makeplot.genebody and set:
  #  avgMethod = "TotalCounts"
  #DEPRECIATED:
  makePlot.genebody.coverage(plotter, plot.medians,
                       plot.means = TRUE,  
                       debugMode, singleEndMode, ...)
  #DEPRECIATED:
  makePlot.genebody.coverage.UMQuartile(plotter, plot.medians, 
                       plot.means = TRUE, 
                       debugMode, singleEndMode, ...)
  #DEPRECIATED:
  makePlot.genebody.coverage.lowExpress(plotter, plot.medians, 
                       plot.means = TRUE, 
                       debugMode, singleEndMode, ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  
  \item{geneset}{
    The set of genes to plot the gene-body coverage over. Genes are grouped
    into four quantiles by their by their total read counts: the 90-100 quantile,
    the 75-90 quantile, the 50-75 quantile (or "upper-middle quartile"), and the
    0-50 quantile. By default it will plot the overall gene-body coverage across
    all genes.
  }
  \item{avgMethod}{
    The method used to generate average gene body coverage. 
    
    The "TotalCounts" (default) method simply takes the sum of each bin across all genes in the geneset,
    then for each replicate normalizes the coverage into frequencies.
    
    The "AvgPercentile" calculates the coverage distribution across each percentile bin 
    and normalizes the counts so they add up to 1, FOR EACH GENE, and 
    THEN averages those values across all genes in the geneset. 
    
    The latter variant is experimental. the idea was for it to reduce the impact of 
    highly-expressed genes and gain a better estimate of the degree of 5'/3' bias.
  }

  
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average for each bamfile at the bottom of the plot. 
    Overrides \code{plot.means}.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average for each bamfile at the bottom of the plot.
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
  x-axis: Gene-body quantile. By default this is broken up into
  40 quantiles containing 2.5% each.
  
  y-axis: Rate at which reads falls into the given quantile of the
  genes' bodies.
  
  \code{makePlot.genebody.coverage} plots the gene body coverage across all genes.
  \code{makePlot.genebody.coverage.UMQuartile} plots the gene body coverage across the genes that fall in the upper-middle quartile of total expression for each sample-run (excluding genes with zero reads).
  \code{makePlot.genebody.coverage.lowExpress} plots the gene body coverage across the genes that fall in the lower two quartiles of total expression for each sample-run (excluding genes with zero reads).
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.genebody(plotter, geneset = "Overall");
  makePlot.genebody(plotter, geneset = "90-100");
  makePlot.genebody(plotter, geneset = "50-75");
}
\seealso{
  \code{\link{build.plotter}}
}