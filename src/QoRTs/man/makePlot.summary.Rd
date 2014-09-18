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
  makePlot.summary.basic(res, outfile.prefix = NULL, 
                      outfile.suffix = "plot.full.summary.png",
                      plotter.params = list(), 
                      plot.device.name, 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      \dots)
                      
  makePlot.summary.colorByGroup(res, outfile.prefix = NULL, 
                      outfile.suffix = "plot.colorByGroup.png",
                      plotter.params = list(), 
                      plot.device.name, 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      \dots)
                      
  makePlot.summary.colorByLane(res, outfile.prefix = NULL, 
                      outfile.suffix = "plot.colorByLane.png",
                      plotter.params = list(), 
                      plot.device.name, 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      \dots)                

  makePlot.summary.sample.highlight(res, 
                      curr.sample,
                      outfile.prefix = NULL, 
                      outfile.suffix,
                      plotter.params = list(), 
                      plot.device.name, 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      \dots)   

  makePlot.summary.sample.highlight.andColorByLane(res, 
                      curr.sample,
                      outfile.prefix = NULL, 
                      outfile.suffix,
                      plotter.params = list(), 
                      plot.device.name, 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode,
                      \dots)   

}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  \item{curr.sample}{
    A character string. For the sample highlight summary plots,
    this should be the sample.ID of the sample that is to be
    highlighted.
  }
  \item{outfile.prefix}{
    A file prefix, used for all output files. Usually the file 
    directory to which you want all files to be written. If 
    left as the default value, then it will attempt to plot to 
    the default device rather than to file (for most setups, when 
    running R in interactive mode this will make the plot pop 
    up as a seperate window).
  }
  \item{outfile.suffix}{
    A file suffix. This will be set to a descriptively-named string ending 
    with ".png", by default.
  }
  \item{plotter.params}{
    Additional parameters (advanced) used in creation of the Plotter 
    objects. See \link{build.plotter}.
  }
  \item{plot.device.name}{ 
    The method to use to save plots. Can be one of:
      \itemize{
        \item \code{"png"} for standard png compression,
        \item \code{"CairoPNG"} for png compression using package Cairo. Note that this requires the package Cairo.
        \item \code{"pdf"} for a multi-page pdf report.
        \item \code{"CairoPDF"} for a multi-page pdf report using package Cairo. Note that this requires the package Cairo.
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
    Logical. If TRUE, more information will be printed to the console.
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
  WIP, Documentation Incomplete!
}

\seealso{

}