\name{read.qc.results.data}
\docType{methods}
\alias{read.qc.results.data}
\title{
   Reading QC results data
}
\description{
   Creates a QoRT_QC_Results object using a set of QC result 
   data files.
}
\usage{
  read.qc.results.data(infile.dir, 
                       decoder.files, 
                       decoder.data.frame, 
                       calc.DESeq2 = FALSE, 
                       calc.edgeR = FALSE, 
                       debugMode)
}
\arguments{
  \item{infile.dir}{
    REQUIRED. The base file directory where all the QC results data
    is stored.
  }
  \item{decoder.files}{
    Either a character string or a character vector, containing the decoder. 
    
    Either decoder.files OR decoder.data.frame must be set. Not both.
  }
  \item{decoder.data.frame}{
    A data.frame containing the decoder information.
  }
  \item{calc.DESeq2}{
    Logical. If TRUE, this function will attempt to load the DESeq2 package. If the DESeq2 package is found, it will then calculate DESeq2's geometric normalization factors (also known as "size factors") for each replicate and for each sample.
  }
  \item{calc.edgeR}{
    Logical. If TRUE, this function will attempt to load the edgeR package. If the edgeR package is found, it will then calculate all of edgeR's  normalization factors (also known as "size factors") for each replicate and for each sample.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  
}
\details{
  Single Decoder: The "decoder" is a file or data frame that describes each replicate. It must have one row per replicate, with the following columns:
  \itemize{
    \item lanebam.ID: The base identifier for the individual replicate. Must be unique!
    \item lane.ID: The identifier for the lane, run, or batch.
    \item group.ID: The identifier for the biological condition for the given replicate.
    \item sample.ID: The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "lanebam.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis.
    \item qc.data.dir: (OPTIONAL) This column indicates the subdirectory in which the replicate's QC data was written. If this column is missing, it is assumed to be equal to the lanebam.ID.
    \item input.read.pair.count: (OPTIONAL) This column contains the number of input reads, before alignment. This is used later to calculate mapping rate.
    \item multi.mapped.read.pair.count: (OPTIONAL) This column contains the number of reads that were multi-mapped. This must be included for multi-mapping rate to be calculated.
    \item cycle.CT: (OPTIONAL) This column contains the length of the reads for the replicate, in base-pairs. By default this is calculated automatically, so this column never needs to be set.
  }
  
  Dual Decoder: Alternatively, two decoders can be set. In this case the first decoder should be the lanebam decoder, and the second should be the sample decoder.
  The lanebam decoder should have one row per replicate, with the following columns:
  \itemize{
    \item lanebam.ID: The base identifier for the individual replicate. Must be unique!
    \item lane.ID: The identifier for the lane, run, or batch.
    \item sample.ID: The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "lanebam.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis.
    \item qc.data.dir: (OPTIONAL) This column indicates the subdirectory in which the replicate's QC data was written. If this column is missing, it is assumed to be equal to the lanebam.ID. Must be unique!
    \item input.read.pair.count: (OPTIONAL) This column contains the number of input reads, before alignment. This is used later to calculate mapping rate.
    \item multi.mapped.read.pair.count: (OPTIONAL) This column contains the number of reads that were multi-mapped. This must be included for multi-mapping rate to be calculated.
    \item cycle.CT: (OPTIONAL) This column contains the length of the reads for the replicate, in base-pairs. By default this is calculated automatically, so this column never needs to be set.
  }
  
  The sample decoder should have one row per replicate, with the following columns:
  \itemize{
    \item group.ID: The identifier for the biological condition for the given replicate.
    \item sample.ID: The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "lanebam.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis. Must be unique!
  }
  
  All decoders are allowed to contain other columns in addition to the ones listed here, so long as their names are distinct. Columns do not need to appear in any particular order, so long as they are named according to the specifications above.
}

\seealso{
  \code{\link{build.plotter}}
}