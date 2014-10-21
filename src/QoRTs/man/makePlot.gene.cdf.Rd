\name{makePlot.gene.cdf}
\docType{methods}
\alias{makePlot.gene.cdf}
\title{
   Plot the cumulative gene diversity curve
}
\description{
   Plots the cumulative gene diversity curve
}
\usage{
  makePlot.gene.cdf(plotter, 
                sampleWise = FALSE,
                plot.intercepts = TRUE,
                label.intercepts = FALSE,
                debugMode, 
                rasterize.plotting.area = FALSE, 
                raster.height = 1000, 
                raster.width = 1000,
                ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{sampleWise}{
    A logical value indicating whether to compile each sample
    together across all runs.
  }
  \item{plot.intercepts}{
    A logical value indicating whether or not to plot the
    intercepts with the round numbers on the x-axis, in 
    dotted lines.
  }
  \item{label.intercepts}{
    A logical value indicating whether or not to label the
    intercepts. No effect if plot.intercepts is FALSE.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{rasterize.plotting.area}{
    Logical. If TRUE, the plotting area will be written to a temp
    png file then drawn to the current device as a raster image.
    This option is generally preferred for vector devices, because
    NVC plots can be very large when drawn in vector format. 
    \emph{Note: this requires the \code{png} package!}
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  For each bam-file, this function plots the cumulative gene diversity. For each bam-file, the genes are sorted by read-count. Then, a cumulative function is calculated for the percent of the total proportion of reads as a function of the number of genes. Intercepts are plotted as well, showing the cumulative percent for 1 gene, 10 genes, 100 genes, 1000 genes, and 10000 genes. 

  So, for example, across all the bam-files, around 50 to 55 percent of the read-pairs were found to map to the top 1000 genes. Around 20 percent of the reads were found in the top 100 genes. And so on.

  This can be used as an indicator of whether a large proportion of the reads stem from of a small number of genes.

  Note that this is restricted to only the reads that map to a single unique gene. Reads that map to more than one gene, or that map to intronic or intergenic areas are ignored.
}

\seealso{
  \code{\link{build.plotter}}
}