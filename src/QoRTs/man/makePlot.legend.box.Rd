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
makePlot.legend.box(plotter, debugMode, singleEndMode, cex.legend, ncol, \dots)
makePlot.legend.over(position, plotter, debugMode, singleEndMode, ncol, \dots)
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
  \item{cex.legend}{
      Numeric. The cex expansion factor passed to \code{legend}. By default
      the cex is autofitted to fill the available space.
  }
  \item{ncol}{
      Integer value. The number of columns in the legend. By default this will be
      selected automatically depending on the number of items in the legend.
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
\examples{
data(res,package="QoRTsExampleData");
plotter <- build.plotter.colorByGroup(res);
#Add a legend to an existing plot:
makePlot.strandedness.test(plotter);
makePlot.legend.over("topright", plotter)

#Or make a legend as a separate plot:
makePlot.legend.box(plotter);
}


\seealso{
  \code{\link{build.plotter}}
}