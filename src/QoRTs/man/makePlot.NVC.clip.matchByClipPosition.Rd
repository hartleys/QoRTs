\name{makePlot.NVC.clip.matchByClipPosition}
\docType{methods}
\alias{makePlot.NVC.clip.matchByClipPosition}
\alias{makePlot.NVC.lead.clip.matchByClipPosition}
\alias{makePlot.NVC.tail.clip.matchByClipPosition}
\title{
   Plot clipped nucleotide rates, aligned by the distance 
   from the point of clipping.
}
\description{
   WARNING: THESE FUNCTIONS ARE BETA, AND ARE NOT FULLY TESTED OR 
   READY FOR GENERAL USE.
}
\usage{
makePlot.NVC.lead.clip.matchByClipPosition(plotter, 
            clip.amt,  r2.buffer, 
            label.majority.bases = TRUE, 
            label.majority.bases.threshold = 0.5, 
            label.majority.bases.cex = 1, 
            show.base.legend = TRUE, 
            load.results = TRUE, 
            debugMode, ...)

makePlot.NVC.tail.clip.matchByClipPosition(plotter, 
            clip.amt, r2.buffer,
            label.majority.bases = TRUE, 
            label.majority.bases.threshold = 0.5, 
            label.majority.bases.cex = 1, 
            show.base.legend = TRUE, 
            debugMode, ...)
}
\arguments{
  \item{plotter}{
    A QoRT_Plotter object.
  }
  \item{clip.amt}{
    UNDOCUMENTED
  }
  \item{r2.buffer}{
    UNDOCUMENTED
  }
  \item{label.majority.bases}{
    UNDOCUMENTED
  }
  \item{label.majority.bases.threshold}{
    UNDOCUMENTED
  }
  \item{label.majority.bases.cex}{
    UNDOCUMENTED
  }
  \item{show.base.legend}{
    UNDOCUMENTED
  }
  \item{load.results}{
    UNDOCUMENTED
  }
  \item{debugMode}{
    UNDOCUMENTED
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   TODO!
}

\seealso{

}