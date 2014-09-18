\name{makePlot.summary.sample.highlight.all}
\docType{methods}
\alias{makePlot.summary.sample.highlight.all}
\alias{makePlot.summary.sample.highlight.andColorByLane.all}
\title{
   Generating sample highlight plots
}
\description{
   Generates multiple sample highlight summary plots, 
   one for every sample.
}
\usage{
  makePlot.summary.sample.highlight.all(res, 
                     outfile.prefix = "./", 
                     plotter.params = list(), 
                     plot.device.name = "png", 
                     plotting.device.params = list(),
                     verbose = TRUE, debugMode,
                     \dots)

  makePlot.summary.sample.highlight.andColorByLane.all(res, 
                     outfile.prefix = "./", 
                     plotter.params = list(), 
                     plot.device.name = "png", 
                     plotting.device.params = list(),
                     verbose = TRUE, debugMode,
                     \dots)
}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  \item{outfile.prefix}{
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