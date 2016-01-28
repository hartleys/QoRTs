\name{makePlot.chrom.type.rates}
\docType{methods}
\alias{makePlot.chrom.type.rates}
\title{
   Chromosome type rate plot
}
\description{
   Plots the number or percent of read-pairs falling on each type of chromosome.
}
\usage{
makePlot.chrom.type.rates(plotter, 
                      plot.rates = TRUE,
                      chromosome.name.style = "UCSC",
                      exclude.autosomes = FALSE,
                      chrom.norm.factors = NULL,
                      custom.chromosome.style.def.function = NULL,
                      return.table = FALSE, 
                      debugMode, singleEndMode,
                      ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.rates}{
    A logical value indicating whether to plot percent of the total 
    (for each bam file), rather than read-counts.
  }
  \item{chromosome.name.style}{
    A string value indicating the style of the chromosome names, 
    and also how to split up the categories. There are 4 legal
    options:
    \itemize{
      \item "UCSC": The default. Chromosomes are named: "chr1,chr2,...,chrX,chrY,chrXY,chrM". There are 6 categories: autosome, X, Y, XY, mitochondrial, and other.
      \item "ENSEMBL": Chromosomes are named: "1,2,...,X,Y,XY,MT". There are 6 categories: autosome, X, Y, XY, mitochondrial, and other.
      \item "UCSC_WITH_ERCC": As UCSC, but there is an additional category, which contains all chromosomes that begin with the text "ERCC".
      \item "ENSEMBL_WITH_ERCC": As ENSEMBL, but there is an additional category, which contains all chromosomes that begin with the text "ERCC". 
    }
  }
  \item{chrom.norm.factors}{
    (Advanced users)
  }
  \item{exclude.autosomes}{
    A logical value indicating whether to exclude autosomes from the plot.
  }
  \item{custom.chromosome.style.def.function}{
    (For advanced users) If your chromosomes do not match any of the above styles, then you can set your own chromosome style by handing this option a function. The function must take one argument. When handed NULL, it must return a list of all legal categories. When handed a single chromosome name, it must return one of those categories. 
  }
  \item{return.table}{
    A Logical value. If TRUE, the function will return a table containing the plotted data.
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
   For each sample-run, this function plots the number of read-pairs mapping to each category of chromosome. 
}
\value{
   By default, this function returns nothing. If the return.table parameter is TRUE, then it returns a data.frame with the data that was plotted.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.chrom.type.rates(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}