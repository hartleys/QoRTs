\name{build.plotter}
\docType{methods}
\alias{build.plotter}
\alias{build.plotter}
\alias{build.plotter.highlightSample.colorByLane}
\alias{build.plotter.highlightSample}
\alias{build.plotter.colorByLane}
\alias{build.plotter.colorByGroup}
\alias{build.plotter.colorBySample}
\alias{build.plotter.basic}
\alias{build.plotter.colorByX}
\alias{QoRTs.default.plotting.params}
\alias{plotter}
\title{
   Generating plotters
}
\description{
   Generating QC_Plotter objects, which can be used in many of the QoRT utilities to organize samples in various ways to allow for easy comparison and detection of consistent biases and artifacts. 
}
\usage{
  build.plotter.basic(res, plotter.params = list())

  build.plotter.colorByGroup(res, plotter.params = list())
  
  build.plotter.colorByLane(res, plotter.params = list())
  
  build.plotter.colorBySample(res, plotter.params = list())
  
  build.plotter.highlightSample(curr.sample,
                                  res,
                                  plotter.params = list(),
                                  merge.offset.outgroup = TRUE)
                                  
  build.plotter.highlightSample.colorByLane(curr.sample,
                                  res,
                                  plotter.params = list(),
                                  merge.offset.outgroup = TRUE)

  build.plotter.colorByX(res, color.by.name,
                         color.by.title.name = color.by.name, 
                         plotter.params = list())

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
  \item{merge.offset.outgroup}{
    (For advanced users) A logical value. For the sample highlight plots, determines whether the all lanes that do not include the current sample should be treated as a single "outgroup".
  }
  \item{plotter.params}{
    (For advanced users) A named list. Allows you to specify 
    colors, offsets, and other similar patterns. By default 
    these will all be set to reasonable values, however, if 
    you want more control over colors, line-transparency, point 
    plotting characters, or similar, then you can specify a 
    named list.
    
    Any parameters that are not specified in the 
    \code{plotter.params} list will be left as default.
    
    Legal parameters are:
      \itemize{
        \item \code{"contrasting.colors"}: colors to use for contrast. By default these are set to a series of reasonably-contrasting colors. However, if you have too many different categories then it may be hard to tell some colors apart.
        \item \code{"contrasting.pch"}: point types to use for contrast (see \code{pch} in \link{graphical parameters}). By default this is set to the basic point types, then following through the upper and lower case letters.

        \item \code{"std.color"}: Color to use for the "highlighted" replicates.
        \item \code{"std.lines.lty"}: Line type to use for the "highlighted" replicates. (see \code{lty} in \link{graphical parameters})
        \item \code{"std.lines.lwd"}: Line width to use for the "highlighted" replicates. (see \code{lwd} in \link{graphical parameters})
        \item \code{"std.lines.alpha"}: Alpha transparency value to use on lines for the "highlighted" replicates. Numeric value between 0 and 255.
        \item \code{"std.points.pch"}: Character to use for points for the "highlighted" replicates. (see \code{pch} in \link{graphical parameters})
        \item \code{"std.points.alpha"}: Alpha transparency value to use on points for the "highlighted" replicates. Numeric value between 0 and 255.
        \item \code{"std.points.color"}: Color to use for the "highlighted" replicates.
        \item \code{"std.NVC.colors"}: A named list with elements named "A", "T","C", and "G", with each element specifying a color. The colors used to indicate each base for the "highlighted" replicates in the nucleotide-rate-by-position plots.

        \item \code{"alt.color"}: Color to use for the "non-highlighted" replicates.
        \item \code{"alt.lines.lty"}: Line type to use for the "non-highlighted" replicates. (see \code{lty} in \link{graphical parameters})
        \item \code{"alt.lines.lwd"}: Line width to use for the "non-highlighted" replicates. (see \code{lwd} in \link{graphical parameters})
        \item \code{"alt.lines.alpha"}: Alpha transparency value to use on lines for the "non-highlighted" replicates. Numeric value between 0 and 255.
        \item \code{"alt.points.pch"}: Character to use for points for the "non-highlighted" replicates. (see \code{pch} in \link{graphical parameters})
        \item \code{"alt.points.alpha"}: Alpha transparency value to use on points for the "non-highlighted" replicates. Numeric value between 0 and 255.
        \item \code{"alt.NVC.colors"}: A named list with elements named "A", "T","C", and "G", with each element specifying a color. The colors used to indicate each base for the "non-highlighted" replicates in the nucleotide-rate-by-position plots.
        
        \item \code{"show.legend"}: DEPRECIATED. Currently nonfunctional.

      }
  }
  \item{color.by.name}{
    (For advanced users) A character string. (TODO: document functionality)
  }
  \item{color.by.title.name}{
    (For advanced users) A character string. (TODO: document functionality)
  }
}
\value{
  A QoRT_Plotter reference object used to create QC summary plots. Depending on which plotter is used, samples/lane-bams can be organized by group, sample, lane, or any arbitrary variable found in the decoder.
}

\seealso{
  \code{\link{read.qc.results.data}} 
}