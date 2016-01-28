

manual <- 
"
#################################################################### \n
# \n
# This R script is made to take a set of QC output from the QoRTs \n
#    QC Utility, run across multiple replicates, and generate \n
#    a battery of cross-comparison reports, organizing replicates \n
#    by sample, lane, and group.ID.  \n
# \n
# This script assumes that the QC statistics have already been \n
#    calculated for each replicate by the (java .jar file) QC \n
#    utility. See the documentation for more details at \n
#    https://github.com/hartleys/QoRTs\n
# \n
# The syntax is: \n
# Rscript qortsGenMultiQC.R infile/dir/ decoderFile.txt outfile/dir/ \n
# \n
# The decoder must be a tab-delimited text file with one row per  \n
#    sample (plus a title row). \n
# \n
# It must contain a column titled \"unique.ID\". This defines a  \n
#    name for each replicate. ALL OTHER COLUMNS ARE OPTIONAL. \n
# If there are multiple (technical) replicates per sample, there  \n
#    must also be a column titled \"sample.ID\", which lists the \n
#    sample ID of each replicate. \n
# By default, each replicate's QC data is assumed to be in a  \n
#    sub-directory of the supplied qcdata base dir with the same \n
#    name as the replicate's unique.ID. If this is not the case, \n
#    then the path to the replicate's qc data must be included in \n
#    a column titled \"qc.data.dir\" \n
# \n
# Optional additional columns: \n
# group.ID: the biological group for the replicate. \n
# lane.ID: the sequencer lane/run for the replicate.  \n
# input.read.pair.count: If the --seqReadCt or --rawfastq parameters  \n
#    were  not used in the QC step, then QoRTs needs some way to \n
#    know how many reads were generated for each replicate, before \n
#    alignment. You can set  this value using this optional column.  \n
#    If this is left out, then mapping rates will be skipped. \n
#\n
####################################################################\n
"


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 3){
  base.dir <- args[1];
  decoder <- args[2];
  output.dir <- args[3];
} else {
  message(manual);
  stop("FATAL ERROR: Syntax must be: Rscript qcdata/base/dir/ decoderFile.txt output/dir/");
}

#print input:
message("base directory: ", base.dir);
message("decoder file: ", decoder);
message("output directory: ", output.dir);


#check dependencies:
if(! require(QoRTs)){
  stop("The QoRTs companion R package is not installed. Install the QoRTs R package.");
}

pngReq <- require(png);
edgeRReq <- require(edgeR);
deseqReq <- require(DESeq2);
if(! pngReq){
  message("Note: the \"png\" package isn't installed. 
  While this is OPTIONAL, it is recommended for
  generating printable pdf files.
  Recommend installing package png with the command: 
  packages.install(\"png\");");
}
if(! edgeRReq){
  message("Note: the \"edgeR\" package isn't installed. 
  While this is OPTIONAL, it will add additional functionality
  if installed.
  Recommend installing package edgeR with the command: 
  source(\"http://bioconductor.org/biocLite.R\");
  biocLite(\"DESeq2\");");
}
if(! deseqReq){
  message("Note: the \"DESeq2\" package isn't installed. 
  While this is OPTIONAL, it will add additional functionality
  if installed.
  Recommend installing package DESeq2 with the command: 
  source(\"http://bioconductor.org/biocLite.R\");
  biocLite(\"DESeq2\");");
}

message("READING QC FILES...");

res <- read.qc.results.data(base.dir, decoder.files = decoder,
       calc.DESeq2 = deseqReq, calc.edgeR = edgeRReq);

message("DONE WITH FILE READING");

message("generating summary table...");
get.summary.table(res, outfile = paste0(output.dir,"summaryTable.txt"));

message("generating png QC multiplots...");
makeMultiPlot.all(res, outfile.dir = output.dir);
message("done.");


message("generating pdf QC reports...");
makeMultiPlot.all(res, outfile.dir = output.dir, plot.device.name = "pdf");
message("done.");

message("ALL DONE!");



