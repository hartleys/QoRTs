\name{makePlot.all.std}
\docType{methods}
\alias{makeMultiPlot.all}
\title{
   Generating all default plots
}
\description{
   Saves MANY compiled QC plots for the given dataset.
}
\usage{
  makeMultiPlot.all(res, outfile.dir = "./", 
                     plotter.params = list(), 
                     plot.device.name = "png", 
                     plotting.device.params = list(), 
                     debugMode,
                     rasterize.large.plots,
                     raster.height = 1000, 
                     raster.width = 1000,
                     exclude.autosomes.chrom.rate.plot = TRUE,
                     chromosome.name.style = "UCSC",
                     fig.res = 150, fig.base.height.inches = 7,
                     \dots)
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
    Additional (advanced) parameters used in creation of the Plotter 
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
  
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  Saves MANY compiled QC plots for the given dataset, calling all 
  standard variants of \code{\link{makeMultiPlot}}.
}

\seealso{
  \code{\link{read.qc.results.data}}, \code{\link{build.plotter}}, \code{\link{makeMultiPlot}}
}