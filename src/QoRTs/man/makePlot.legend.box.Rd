\name{makePlot.legend.box}
\docType{methods}
\alias{makePlot.legend.box}
\alias{makePlot.legend.over}
\title{
   Plot legend
}
\description{
   Plots the universal legend for a given plotter object.
}
\usage{
makePlot.legend.box(plotter, debugMode, singleEndMode, \dots)
makePlot.legend.over(position, plotter, debugMode, singleEndMode, \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{position}{
    For makePlot.legend.over, a character string indicating the location 
    where you want the legend to appear. This is passed to 
    \code{\link{legend}}, and can be any keyword allowed by 
    \code{\link{xy.coords}}. For example: "top", 
    "topleft","bottomright", etc.
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
   makePlot.legend.box creates a new plot (opening the next graphics
   frame), and writes a small description of the given plotter 
   type along with plotting a color-coded legend (if applicable).
   
   makePlot.legend.over adds a legend to the current plotting frame.
}

\seealso{
  \code{\link{build.plotter}}
}