\name{makePlot.summary}
\docType{methods}
\alias{makePlot.summary}
\alias{makePlot.summary.basic}
\alias{makePlot.summary.colorByGroup}
\alias{makePlot.summary.colorByLane}
\alias{makePlot.summary.sample.highlight}
\alias{makePlot.summary.sample.highlight.andColorByLane}
\title{
   Generating summary plots
}
\description{
   Creates a large multi-frame summary plot, or a report PDF file.
}
\usage{
  makePlot.summary.basic(res, outfile,
                      outfile.dir = "./", 
                      outfile.prefix = "plot.summary.basic",
                      outfile.ext,
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      rasterize.large.plots,
                      raster.height = 1000, 
                      raster.width = 1000,
                      separatePlots = FALSE,
                      \dots)
                      
  makePlot.summary.colorByGroup(res, outfile,
                      outfile.dir = "./", 
                      outfile.prefix = "plot.summary.colorByGroup",
                      outfile.ext,
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      rasterize.large.plots,
                      raster.height = 1000, 
                      raster.width = 1000,
                      separatePlots = FALSE,
                      \dots)
                      
  makePlot.summary.colorByLane(res, outfile,
                      outfile.dir = "./", 
                      outfile.prefix = "plot.summary.colorByLane",
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      rasterize.large.plots,
                      raster.height = 1000, 
                      raster.width = 1000,
                      separatePlots = FALSE,
                      \dots)                

  makePlot.summary.sample.highlight(res,
                      curr.sample,
                      outfile,
                      outfile.dir = "./", 
                      outfile.prefix,
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      rasterize.large.plots,
                      raster.height = 1000, 
                      raster.width = 1000,
                      separatePlots = FALSE,
                      \dots)   

  makePlot.summary.sample.highlight.andColorByLane(res,
                      curr.sample,
                      outfile,
                      outfile.dir = "./", 
                      outfile.prefix,
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      rasterize.large.plots,
                      raster.height = 1000, 
                      raster.width = 1000,
                      separatePlots = FALSE,
                      \dots)   

}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by 
    \code{\link{read.qc.results.data}}.
  }
  \item{curr.sample}{
    A character string. For the sample highlight summary plots,
    this should be the sample.ID of the sample that is to be
    highlighted.
  }
  \item{outfile}{
    The output file can be specified in two ways: first, by 
    setting the outfile parameter, which should be the full path
    for the output file. Alternatively, the file path will
    be the concatenation of the three parameters:
    outfile.dir, outfile.prefix, and outfile.ext.
    
    If plot.device.name == "curr", then all the outfile 
    parameters will be nonfunctional.
  }
  \item{outfile.dir}{
    A prefix to be appended to the output filename. Usually 
    the path to the destination directory. By default, this
    will be the current working directory.
  }
  \item{outfile.prefix}{
    A second prefix to be appended to the output filename, after
    outfile.dir. This is usually the name of the file, without
    the file extension. By default, a reasonable file name will be
    selected.
  }
  \item{outfile.ext}{
    The file extension. By default, this will be set to match
    the selected graphics device.
  }
  \item{plotter.params}{
    Additional parameters (advanced) used in creation of the Plotter 
    objects. See \link{build.plotter}.
  }
  \item{plot.device.name}{ 
    The method to use to save plots. Several formats and devices are supported:
      \itemize{
        \item \code{"curr"} (default) Do not create output files. Plot directly to the current or default device.
        \item \code{"png"} for standard png compression.
        \item \code{"CairoPNG"} for png compression using package Cairo. Note that this requires the package Cairo.
        \item \code{"tiff"} for the tiff device.
        \item \code{"jpeg"} for the jpeg device (\emph{not recommended!})
        \item \code{"svg"} for the vector svg device. Note that for all vector devices, by default, certain large plots will be rasterized if the \code{png} package is found. If not, the function will still work, but output files may be very large and may have trouble rendering and printing.
        \item \code{"pdf"} for a multi-page pdf report. 
        \item \code{"CairoPDF"} for a multi-page pdf report using package Cairo. Note that this requires the package Cairo. 
      }
  }
  \item{plotting.device.params}{ 
    A named list of parameters to be passed to the graphics device.
    For example:
      \itemize{
        \item \code{width = 2000}
        \item \code{height = 2000}
        \item \code{units = "px"}
      }
    Reasonable values for height, width, and pointsize will be chosen
    by default. Generally all raster plots will be set to 1000 by 1000 
    pixels, and all vector plots will be set to 7 inches by 7 inches.
  }
  \item{verbose}{
    Logical. If TRUE, more information will be printed to the console.
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
  \item{separatePlots}{
    Logical. If TRUE, then separate image files will be saved for 
    each individual plot, rather than saving one large multi-pane 
    plot. Note that this is not compatible with the "curr" device.
    Also note: if this parameter is set to TRUE, then the output
    files will be saved using the outfile.dir, outfile.prefix
    and outfile.ext parameters. The "outfile" parameter cannot be
    set. The files will be saved to the contatenation of 
    outfile.dir, outfile.prefix, the name of the plot type, 
    and then outfile.ext.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
    Produces large, multi-frame summary plots, with a large number of different
    plots and graphs.
    
    For devices pdf and CairoPDF, this function will produce a multi-page document
    suitable for printing, with 6 panels on each page.
    
    Note that for all vector devices (svg, pdf, and CairoPDF), by default, 
    QoRTs will attempt to load the \code{png} package. If this package is 
    found, then certain large plots (the NVC plots and the cumulative 
    diversity plots) will be rasterized within the plotting area. The labels
    and text will still be vectorized. If the \code{png} package is not installed, 
    the function will still work, but output files may be very large and may 
    have trouble rendering and printing.
    
    Rasteration of large plots can be turned off via the rasterize.large.plots option.
}

\seealso{
  \code{\link{build.plotter}}
}