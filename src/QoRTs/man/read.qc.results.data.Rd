\name{read.qc.results.data}
\docType{methods}
\alias{read.qc.results.data}
\alias{completeAndCheckDecoder}
\title{
   Reading QC results data
}
\description{
   Creates a QoRT_QC_Results object using a set of QC result 
   data files.
}
\usage{
  read.qc.results.data(infile.dir, 
                       decoder, 
                       decoder.files, 
                       calc.DESeq2 = FALSE, 
                       calc.edgeR = FALSE, 
                       debugMode)
                       
  completeAndCheckDecoder(decoder, decoder.files)
}
\arguments{
  \item{infile.dir}{
    REQUIRED. The base file directory where all the QC results data
    is stored.
  }
  \item{decoder}{
    A character vector or data.frame containing the decoder information. See details below.
  }
  \item{decoder.files}{
    Character vector. Either one or two character strings. Either decoder.files OR decoder must be set, never both. See details below.
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
  read.qc.results.data reads in a full QoRTs dataset of multiple QoRTs QC runs and compiles them into a QoRTs_Results object.
  
  completeAndCheckDecoder simply reads a decoder and "fills in" all missing parameters, returning a data.frame.
  
  The "decoder" is used to describe each replicate/sample. The standard decoder is a data frame that has one row per replicate, with the following columns:
  \itemize{
    \item unique.ID: The base identifier for the individual replicate. Must be unique. lanebam.ID is a synonym.
    \item lane.ID: (OPTIONAL) The identifier for the lane, run, or batch. The default is "UNKNOWN".
    \item group.ID: (OPTIONAL) The identifier for the biological condition for the given replicate. The default is "UNKNOWN".
    \item sample.ID: (OPTIONAL) The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "lanebam.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis. By default, it is assumed that each replicate comes from a different sample, and sample.ID is set to equal unique.ID.
    \item qc.data.dir: (OPTIONAL) This column indicates the subdirectory in which the replicate's QC data was written. If this column is missing, it is assumed to be equal to the unique.ID.
    \item input.read.pair.count: (OPTIONAL) This column contains the number of read-pairs (or just reads, for single-end data), before alignment. This is used later to calculate mapping rate.
    \item multi.mapped.read.pair.count: (OPTIONAL) This column contains the number of read-pairs (or just reads, for single-end data) that were multi-mapped. This must be included for multi-mapping rate to be calculated.
  }
  All the parameters except for unique.ID are optional. The decoder can even be supplied as a simple character vector, which is assumed to be the unique.ID's. All the other variables will be set to their default values.
  
  Alternatively, the decoder can be supplied as a file given by the decoder.files parameter.
  
  Dual Decoder: Optionally, two decoders can be supplied. In this case the first decoder should be the technical-replicate decoder and the second should be the biological-replicate decoder.
  The technical-replicate decoder should have one row per unique.ID, with the following columns:
  \itemize{
    \item unique.ID: The base identifier for the individual replicate. Must be unique!
    \item lane.ID: (OPTIONAL) The identifier for the lane, run, or batch.
    \item sample.ID: (OPTIONAL) The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "lanebam.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis.
    \item qc.data.dir: (OPTIONAL) This column indicates the subdirectory in which the replicate's QC data was written. If this column is missing, it is assumed to be equal to the lanebam.ID. Must be unique!
    \item input.read.pair.count: (OPTIONAL) This column contains the number of input reads, before alignment. This is used later to calculate mapping rate.
    \item multi.mapped.read.pair.count: (OPTIONAL) This column contains the number of reads that were multi-mapped. This must be included for multi-mapping rate to be calculated.
  }
  
  The biological-replicate decoder should have one row per sample.ID, with the following columns:
  \itemize{
    \item sample.ID: The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "unique.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis.
    \item group.ID (OPTIONAL): The identifier for the biological condition for the given replicate.
  }
  
  All decoders are allowed to contain other columns in addition to the ones listed here, so long as their names are distinct. Columns do not need to appear in any particular order, so long as they are named according to the specifications above.
}

\examples{
  #Load the decoder from the example dataset:
  directory <- paste0(system.file("extdata/", 
                      package="QoRTsExampleData", 
                      mustWork=TRUE),"/");
  decoder.file <- system.file("extdata/decoder.txt", 
                              package="QoRTsExampleData", 
                              mustWork=TRUE);
  decoder.data <- read.table(decoder.file,
                             header=TRUE,
                             stringsAsFactors=FALSE);
  print(decoder.data);

  #This command produces the example dataset used in all the other
  #  examples:
  res <- read.qc.results.data(directory, 
                              decoder = decoder.data,
                              calc.DESeq2 = TRUE, 
                              calc.edgeR = TRUE);
  #Note that DESeq2 and edgeR are required in order to 
  #  calculate the size factors using the options above.
  
  #You can also specify incomplete decoders, and use
  #  the following command to fill in the defaults:
  completeAndCheckDecoder(c("SAMP1","SAMP2",
                            "SAMP3","SAMP4",
                            "SAMP5","SAMP6"))

  #You don't actually have to use completeAndCheckDecoder,
  #You can just pass the incomplete decoder directly to QoRTs.
  #For example, to load a small subset of the example data 
  #(without phenotype data):
  res <- read.qc.results.data(paste0(directory,"/ex/"), 
                              decoder = c("SAMP1_RG1","SAMP2_RG1",
                                          "SAMP3_RG1","SAMP4_RG1"));

}

\seealso{
  \code{\link{build.plotter}}
}