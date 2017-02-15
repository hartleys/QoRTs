\name{makeMultiPlot.highlightSample.all}
\docType{methods}
\alias{makeMultiPlot.highlightSample.all}
\alias{makeMultiPlot.highlightSample.colorByLane.all}
\title{
   Generating sample highlight multi-plots
}
\description{
   Generates multiple sample highlight summary plots, 
   one for every sample.
}
\usage{
  makeMultiPlot.highlightSample.all(res, outfile.dir = "./",
                plotter.params = list(), 
                plot.device.name = "png", 
                plotting.device.params = list(), 
                verbose = TRUE,
                debugMode , 
                rasterize.large.plots, 
                rasterize.medium.plots,
                raster.height = 1050, 
                raster.width = 1050,
                exclude.autosomes.chrom.rate.plot = TRUE,
                chromosome.name.style = "UCSC",
                fig.res = 150, fig.base.height.inches = 7,
                insertSize.plot.xlim,
                sequencing.type = c("RNA","Exome","Genome"),
                maxColumns,
                maxRows,
                plotList,
                labelPlots=TRUE,
                ...)

makeMultiPlot.highlightSample.colorByLane.all(res, 
                outfile.dir = "./",
                plotter.params = list(), 
                plot.device.name = "png", 
                plotting.device.params = list(), 
                verbose = TRUE, 
                debugMode ,
                rasterize.large.plots, 
                rasterize.medium.plots,
                raster.height = 1050, 
                raster.width = 1050,
                exclude.autosomes.chrom.rate.plot = TRUE,
                chromosome.name.style = "UCSC",
                fig.res = 150, fig.base.height.inches = 7,
                insertSize.plot.xlim,
                sequencing.type = c("RNA","Exome","Genome"),
                maxColumns,
                maxRows,
                plotList,
                labelPlots=TRUE,
                ...)
}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  \item{outfile.dir}{
    A file prefix, used for all output files. Usually the directory to
    which you want all files to be written.
  }
  \item{plotter.params}{
    Additional parameters (advanced) used in creation of the Plotter 
    objects. See \link{build.plotter}. 
  }
  \item{plot.device.name}{ 
    The method to use to save plots. Can be one of:
      \itemize{
        \item \code{"png"} for standard png compression,
        \item \code{"CairoPNG"} for png compression using package Cairo.
          Note that this requires the package Cairo.
      }
  }
  \item{plotting.device.params}{ 
    A named list of parameters to be passed to the graphics device.
    For example:
      \itemize{
        \item \code{"width = 2000"}
      }
    Reasonable values for height, width, and pointsize will be chosen
    by default.
  }
  \item{verbose}{
    Logical. If TRUE, more info will be printed to the console.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{rasterize.large.plots}{
    Logical. If TRUE, then if the currently selected plotting device
    is a vector device (or the "curr" device), then certain large plots
    will have their plotting areas rasterized. The axis labels, titles,
    and text will all remain vectorized, only the plotting areas will
    be flattened. \emph{Note that this requires the png package.}
    
    By default, this parameter will be set to TRUE when a vector
    device is selected.
  }
  \item{rasterize.medium.plots}{
    Logical. as rasterize.large.plots, but applies to moderately-large plots. By default, this parameter will be set to TRUE for pdf/CairoPDF output only.
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels. Double-pane plots
    will be twice this width.
  }
  \item{exclude.autosomes.chrom.rate.plot}{
    A logical value indicating whether to exclude autosomes from the plot.
    See \code{\link{makePlot.chrom.type.rates}}
  }
  \item{chromosome.name.style}{
    A string value indicating the style of the chromosome names, 
    and also how to split up the categories. See 
    \code{\link{makePlot.chrom.type.rates}}
  }
  \item{fig.res}{
    Numeric value. The number of pixels per "inch" (for raster devices only).
    For some plotting devices the figure height will be in pixels not inches, and
    will equal this value multiplied by the fig.base.height.inches value.
  }
  \item{fig.base.height.inches}{
    Numeric value. The base height, in inches, of each sub-figure in the plot. This
    will be equal to the height for vector devices, or used to calculate the height
    in pixels using the fig.res value (see above).
  }
  \item{insertSize.plot.xlim}{
    A numeric vector of length 2. The x-axis limits for the insert size plot. By default QoRTs will attempt to find reasonable values for this, 
    but there are always situations where the default behavior is not ideal. Using this parameter you can set
    it explicitly.
  }
  \item{sequencing.type}{
      The type of sequencing data being analyzed. This only changes the default plot set, which can be overriden with the \code{plotList} parameter.
  }
  \item{maxColumns}{
    If set, QoRTs will attempt to create a multiplot that has (at most) maxColumns columns. Extra rows will be added to fit all the plots, as needed.
  }
  \item{maxRows}{
    If set, QoRTs will attempt to create a multiplot that has (at most) maxRows rows. Extra columns will be added to fit all the plots as needed. To set the number of rows and columns manually, you can set both maxColumns and maxRows.
  }
  \item{plotList}{
    The list of desired plots.
  }
  \item{labelPlots}{
    Logical. If TRUE, then label each plot with a letter.
  }
  
  
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  Generates a sample highlight plot for each sample in the dataset.
}

\examples{
  \dontrun{
  data(res,package="QoRTsExampleData");
    makeMultiPlot.highlightSample.all(res);
    makeMultiPlot.highlightSample.colorByLane.all(res);
  }
}
\seealso{
  \code{\link{build.plotter}}, \code{\link{makeMultiPlot}}
}