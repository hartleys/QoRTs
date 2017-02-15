
require(methods);
DEFAULTDEBUGMODE <- TRUE;

##############################
#Res builders:

#Input File requirements:
#Input can come in two forms:
#(1) Single Decoder
#    This requires a single "decoder" file, which MUST have the following column names:
#    unique.ID	lane.ID	group.ID	sample.ID	qc.data.dir
#       Other than these 4 required columns, it can have as many additional columns as desired. Column names must be unique.
#       OPTIONAL FIELDS:
#          input.read.pair.count: the # of input reads. this must be included for mapping rate to be calculated.
#          multi.mapped.read.pair.count: the # of reads that were multi mapped by the aligner. this must be included for multi-mapping rate to be calculated.
#       RESERVED FIELDS: Do not name any field this:
#          cycle.CT
#          lanebam.ID (a synonym for unique.ID)

#(2) Dual Decoder:
#    This requires two "decoder" files.
#    Replicate Decoder:
#        The Replicate Decoder lists the technical replicates of each sample, and must have the following column names:
#        unique.ID	lane.ID	sample.ID
#           Other than these 3 required columns, it can have any other columns desired, as long as the column names are unique.
#        OPTIONAL FIELDS:
#           input.read.pair.count: the # of input reads. this must be included for mapping rate to be calculated.
#           multi.mapped.read.pair.count: the # of reads that were multi mapped by the aligner. 
#                                         Depending on how your aligner implements multi-mapping, this might be required
#                                         for multi-mapping rate to be calculated.
#        RESERVED FIELDS: Do not name any field this:
#           cycle.CT
#           lanebam.ID (a synonym for unique.ID)
#    Sample Decoder:
#        The Sample Decoder must have the following column names:
#        sample.ID	group.ID
#            other than these 2 required columns, it can have any other columns desired, as long as the column names are unique.
#

read.qc.results.data <- function(infile.dir, 
                                 decoder = NULL, 
                                 decoder.files = NULL, 
                                 calc.DESeq2 = FALSE, 
                                 calc.edgeR = FALSE, 
                                 debugMode = DEFAULTDEBUGMODE,
                                 autodetectMissingSamples = FALSE,
                                 skip.files = c()
                                 ) {
   decoder.final <- completeAndCheckDecoder(decoder = decoder, decoder.files = decoder.files);
   
   return(read.in.results.data.with.decoder(decoder = decoder.final, 
                                            infile.dir = infile.dir, 
                                            calc.DESeq2 = calc.DESeq2 , 
                                            calc.edgeR = calc.edgeR, 
                                            debugMode = debugMode,
                                            autodetectMissingSamples=autodetectMissingSamples,
                                            skip.files = skip.files) );
}

completeAndCheckDecoder <- function(decoder = NULL, decoder.files = NULL){
   if(! is.null(decoder)){
     #Do nothing!
   } else if(! is.null(decoder.files)){
     decoder.file <- as.character(decoder.files);
     if(length(decoder.files) == 1){
       decoder <- get.decoder.from.single.file(decoder.files);
     } else if(length(decoder.files) == 2){
       decoder <- get.decoder.from.dual.file( lanebam.decoder = decoder.files[1], sample.decoder = decoder.files[2]);
     } else {
       stop("Fatal Error: QoRTs.read.results.data(): decoder.files must have length of 1 or 2!");
     }
   } else {
     stop("Fatal Error: QoRTs.read.results.data(): either decoder or decoder.files must be set!");
   }
   
   decoder.final <- expandAndCheckDecoder(decoder);
   
   return(decoder.final);
}

expandAndCheckDecoder <- function(decoder) {
  #Simplest decoder possible: just a list of sample ID's:
  if(class(decoder) == "character"){
    decoder <- data.frame(unique.ID = decoder, stringsAsFactors = FALSE);
    message("Simple decoder found, list of sample/unique ID's. Building complete decoder.");
  }
  if(class(decoder) == "matrix"){
    decoder <- data.frame(decoder, stringsAsFactors = FALSE);
  }
  if(class(decoder) != "data.frame"){
    message("Decoder must be either a matrix, a character vector (of unique ID's), or a data.frame. Decoder follows:");
    print(decoder);
    stop("ERROR: Decoder must be either a matrix, a character vector (of unique ID's), or a data.frame.");
  }
  #Convert strings to characters:
     decoder <- data.frame(decoder,stringsAsFactors=F);
     for(i in 1:length(names(decoder))){
        if(class(decoder[[i]]) == "factor"){
          message("Converting column ",names(decoder)[i]," to character vector.");
          decoder[[i]] <- as.character(decoder[[i]]);
        }
     }
  
  if(is.null(names(decoder))){
    message("Decoder must have column names! Decoder follows:");
    print(decoder);
    stop("ERROR: decoder must have column names!");
  }
  if(any(duplicated(names(decoder)))) stop("Decoder error: names(decoder) contains duplicates.");
  
  if((! ("unique.ID" %in% names(decoder))) & (! ("lanebam.ID" %in% names(decoder)))  & (! ("sample.ID" %in% names(decoder)))  ){
    err <- "Decoder must contain a unique.ID or sample.ID row.";
    message(err, " Decoder follows:");
    print(decoder);
    stop("ERROR: ",err);
  }
  
  #If unique.ID AND lanebam.ID is unspecified, define it as the sample.ID 
  #   (note: this assumes 1-1 sample to readGroup)
  if((! ("unique.ID" %in% names(decoder))) & (! ("lanebam.ID" %in% names(decoder)))){
    decoder$unique.ID = decoder$sample.ID;
  }
  #If unique.ID is specified but lanebam.ID exists, then use lanebam.ID as a synonym for unique.ID
  if(! ("unique.ID" %in% names(decoder))){
    decoder$unique.ID = decoder$lanebam.ID;
  }
  #If they decided to specify by sample.ID, then use the sample.ID to define the unique id 
  #   (note: this assumes 1-1 sample to readGroup)
  if(! ("sample.ID" %in% names(decoder))){
    decoder$sample.ID = decoder$unique.ID;
  }
  # if lane.ID is unspecified, mark them UNKNOWN.
  if(! ("lane.ID" %in% names(decoder))){
    decoder$lane.ID = "UNKNOWN";
  }
  # if group.ID is unspecified, mark them UNKNOWN.
  if(! ("group.ID" %in% names(decoder))){ 
    decoder$group.ID = "UNKNOWN";
  }
  # if qc.data.dir is unspecified, assume they used the unique.ID's as the qc.data.dir.
  #   This one might cause confusion, so print a minor warning just to be safe:
  if(! ("qc.data.dir" %in% names(decoder))){ 
    message("column 'qc.data.dir' not found in the decoder, assuming qc.data.dir = unique.ID");
    decoder$qc.data.dir = decoder$unique.ID;
  }
  
  decoder$qc.data.dir <- as.character(decoder$qc.data.dir);
  decoder$group.ID <- as.character(decoder$group.ID);
  decoder$lane.ID <- as.character(decoder$lane.ID);
  decoder$sample.ID <- as.character(decoder$sample.ID);
  decoder$unique.ID <- as.character(decoder$unique.ID);
  
  #Now do checks for validity:
  
  #These should never be triggered, since these variables were all set above. But just in case:
  if(is.null(decoder$unique.ID)) stop("Decoder formatting error: no column labelled unique.ID");
  if(is.null(decoder$lane.ID)) stop("Decoder formatting error: no column labelled lane.ID");
  if(is.null(decoder$group.ID)) stop("Decoder formatting error: no column labelled group.ID");
  if(is.null(decoder$sample.ID)) stop("Decoder formatting error: no column labelled sample.ID");
  
  if(any(duplicated( decoder$unique.ID ))) stop("Decoder error: unique.ID's must be unique.");
  if(any(duplicated(names(decoder)))) stop("Decoder error: names(decoder) contains duplicates.");
  
  for(samp in unique(decoder$sample.ID)){
    subsetDecoder <- decoder[ decoder$sample.ID == samp,];
    if(nrow(subsetDecoder) > 1){
      if(length(unique(subsetDecoder$group.ID)) != 1){
        stop("Decoder error: sample ",samp," has inconsistent group.ID depending on decoder row: (",paste0(unique,collapse=","),")");
      }
    }
  }
  
  return(decoder);
  
}

##########################################################################################

get.decoder.from.single.file <- function(decoder.file){
  decoder <- read.table(decoder.file,header=TRUE,stringsAsFactors=F);
  
  #if(! is.null(decoder$lanebam.ID)){
  #  decoder$unique.ID <- decoder$lanebam.ID;
  #}
  #if(is.null(decoder$unique.ID)) stop("Decoder formatting error: no column labelled unique.ID");
  #if(is.null(decoder$lane.ID)) stop("Decoder formatting error: no column labelled lane.ID");
  #if(is.null(decoder$group.ID)) stop("Decoder formatting error: no column labelled group.ID");
  #if(is.null(decoder$sample.ID)) stop("Decoder formatting error: no column labelled sample.ID");
  
  #if(is.null(decoder$qc.data.dir)){
  #   message("qc.data.dir not found, assuming qc.data.dir = unique.ID");
  #   decoder$qc.data.dir <- decoder$unique.ID;
  #}  
  
  #if(length(unique(decoder$unique.ID)) < length(decoder$unique.ID)) stop("Decoder formatting error: unique.ID must be unique!");
  ##more checks?
  #if( length(unique(names(decoder))) != length(names(decoder)) ) stop("Decoder formatting error: column names must be unique!");
  return(decoder);
}

get.decoder.from.dual.file <- function(lanebam.decoder, sample.decoder){
  sample.decoder <- read.table(sample.decoder,header=TRUE,stringsAsFactors=F);
  lanebam.decoder <- read.table(lanebam.decoder,header=TRUE,stringsAsFactors=F);
  
  if(! is.null(lanebam.decoder$lanebam.ID)){
    lanebam.decoder$unique.ID <- lanebam.decoder$lanebam.ID;
  }
  if(is.null(lanebam.decoder$unique.ID)) stop("Replicate Decoder (File: ",lanebam.decoder,") formatting error: no column labelled unique.ID");
  if(is.null(lanebam.decoder$sample.ID)) stop("Replicate Decoder (File: ",lanebam.decoder,") formatting error: no column labelled sample.ID");
  if(any(duplicated(lanebam.decoder$unique.ID))) stop("Replicate Decoder (File: ",lanebam.decoder,") formatting error: unique.ID must be unique!");
  if(any(duplicated(names(lanebam.decoder)))) stop("Replicate Decoder (File: ",lanebam.decoder,") formatting error: column names must be unique!");

  if(is.null(sample.decoder$sample.ID)) stop("Sample Decoder (File: ",sample.decoder,") formatting error: no column labelled sample.ID");
  if(any(duplicated(sample.decoder$sample.ID))) stop("Sample Decoder (File: ",sample.decoder,") formatting error: sample.ID must be unique!");
  if(any(duplicated(names(sample.decoder)))) stop("Sample Decoder (File: ",sample.decoder,") formatting error: column names must be unique!");

  if(! all(unique(sample.decoder$sample.ID) %in% unique(lanebam.decoder$sample.ID) )) {
    stop("Sample/Lanebam Decoders formatting error: no 1-to-1 matching between the two decoders sample.IDs!\n     Missing from replicate decoder: ", paste0(unique(sample.decoder$sample.ID)[ ! unique(sample.decoder$sample.ID) %in% unique(lanebam.decoder$sample.ID) ], collapse=",")  )
  }
  if(! all(unique(lanebam.decoder$sample.ID) %in% unique(sample.decoder$sample.ID) )) {
    stop("Sample/Lanebam Decoders formatting error: no 1-to-1 matching between the two decoders sample.IDs!\n     Missing from sample decoder: ", paste0(unique(lanebam.decoder$sample.ID)[ ! unique(lanebam.decoder$sample.ID) %in% unique(sample.decoder$sample.ID) ],collapse=",")  )
  }
  
  merged.decoder <- lanebam.decoder;

  for(pheno.var in names(sample.decoder)){
    if(! pheno.var %in% c("sample.ID")){
      merged.decoder[[pheno.var]] <- sapply(lanebam.decoder$sample.ID, FUN=function(samp){
        sample.decoder[[pheno.var]][sample.decoder$sample.ID == samp];
      });
    }
  }
  return(merged.decoder);
}


###############################################################

read.in.results.data.with.decoder <- function(decoder, infile.dir = "", 
                                              calc.DESeq2 = FALSE, 
                                              calc.edgeR = FALSE, 
                                              autodetectMissingSamples = FALSE, 
                                              skip.files = c(),
                                              debugMode = DEFAULTDEBUGMODE){
  if(! is.null(decoder$lanebam.ID)){
    decoder$unique.ID <- decoder$lanebam.ID;
  }
  if(is.null(decoder$unique.ID)) stop("Decoder formatting error: no column labelled unique.ID");
  if(is.null(decoder$lane.ID)) stop("Decoder formatting error: no column labelled lane.ID");
  if(is.null(decoder$group.ID)) stop("Decoder formatting error: no column labelled group.ID");
  if(is.null(decoder$sample.ID)) stop("Decoder formatting error: no column labelled sample.ID");
  if(is.null(decoder$qc.data.dir)) stop("Decoder formatting error: no column labelled qc.data.dir");
  if(length(unique(decoder$unique.ID)) < length(decoder$unique.ID)) stop("Decoder formatting error: unique.ID must be unique!");
  #more checks?
  if( length(unique(names(decoder))) != length(names(decoder)) ) stop("Decoder formatting error: column names must be unique!");

  if( is.null(decoder$input.read.pair.count)){
    message("Note: no input.read.pair.count column found. This column is optional, but without it mapping rates cannot be calculated.");
  }# else {
  #  message("Note: successfully found input.read.pair.count column. This will be used to calculate mapping rates.");
  #}
  if( is.null(decoder$multi.mapped.read.pair.count)){
    message("Note: no multi.mapped.read.pair.count column found. This column is optional, but without it (depending on how your aligner implements multi-mapping) multi-mapping rates might not be plotted.");
  }# else {
  #  message("Note: successfully found multi.mapped.read.pair.count: column. This will be used to plot multi-mapping rates.");
  #}
  

  dataDirs <- paste0(infile.dir,decoder$qc.data.dir);
  dataDirExists <- file.exists(dataDirs);
  if(any(! dataDirExists)){
    message("WARNING: QoRTs run absent! Dir not found: ",paste0(infile.dir,decoder$qc.data.dir[!dataDirExists],""),"!");
    if(autodetectMissingSamples){
      message("      Skipping missing samples!");
      decoder <- decoder[dataDirExists,,drop=FALSE];
    } else {
      stop("Fatal error: QoRTs run data not found! Use autodetectMissingSamples = TRUE to automatically skip these runs");
    }
  }

  
  compFiles <- paste0(infile.dir,decoder$qc.data.dir,"/QC.QORTS_COMPLETED_OK");
  compFileExists <- file.exists(compFiles)
  if(any(! compFileExists)){
    message("WARNING: QoRTs run may be incomplete! File not found: ",paste0(infile.dir,decoder$qc.data.dir[!compFileExists],"/QC.QORTS_COMPLETED_OK"),"!");
    if(autodetectMissingSamples){
      message("      Skipping missing samples!");
      decoder <- decoder[compFileExists,,drop=FALSE];
    } else {
      #stop("Fatal error: QoRTs run data not found! Use autodetectMissingSamples = TRUE to automatically skip these runs");
    }
  }
  
  res <- new("QoRTs_QC_Results");
  res@lanebam.list <- decoder$unique.ID;
  res@sample.list <- unique(decoder$sample.ID);
  res@lane.list <- unique(decoder$lane.ID);
  res@group.list <- unique(decoder$group.ID);
  
  res@decoder <- cbind.data.frame(decoder,cycle.CT = rep(-1,length(decoder$unique.ID)));

  lanebam.list <- as.list(decoder$unique.ID);
  names(lanebam.list) <- decoder$unique.ID
  
  qc.data.dir.list <- as.list(decoder$qc.data.dir);
  names(qc.data.dir.list) <- decoder$unique.ID

  if(debugMode) message("infile.dir = ",infile.dir);
  #if(debugMode) message("qc.data.dir.list = ",paste0(qc.data.dir.list,collapse=","));
  #if(debugMode) message("lanebam.list = ",paste0(lanebam.list,collapse=","));

  read.scalaqc.file.helper <- function(scalaqc_file, sep=""){
    if(debugMode) message(paste0("scalaqc_file = ",scalaqc_file), appendLF=FALSE);
    if(debugMode) ts <- timestamp();
    out <- read.in.scalaQC.files(infile.dir,lanebam.list, qc.data.dir.list,scalaqc_file,sep=sep);
    if(debugMode) reportTimeAndDiff(ts,prefix="   ");
    out;
  }
  
  res@qc.data <- list(summary = read.scalaqc.file.helper(QC_INTERNAL_SCALAQC_SUMMARY_FILE,sep='\t'));
  res@qc.data[["summary"]] <- lapply(res@qc.data[["summary"]],function(sumTable){sumTable[,1:2]})

  allSingleEnd <- any(sapply(res@qc.data[["summary"]], function(dl){
    dl$COUNT[dl$FIELD == "IS_SINGLE_END"] == 1;
  }));
  allPairedEnd <- any(sapply(res@qc.data[["summary"]], function(dl){
    dl$COUNT[dl$FIELD == "IS_SINGLE_END"] == 0;
  }));
  if(! (allSingleEnd | allPairedEnd)){
    stop("FATAL ERROR: Dataset is a mixture of single-end and paired-end data! You must analyse paired-end and single-end datasets separately.");
  }
  res@singleEnd <- allSingleEnd;

  if(res@singleEnd){
    if(debugMode) message("Autodetected Single-End mode.");
    USE.LIST <- QC_INTERNAL_SCALAQC_FILE_LIST_SINGLE_END
    #res@qc.data <- c(res@qc.data,lapply(QC_INTERNAL_SCALAQC_FILE_LIST_SINGLE_END, FUN=read.scalaqc.file.helper));
  } else {
    if(debugMode) message("Autodetected Paired-End mode.");
    USE.LIST <- QC_INTERNAL_SCALAQC_FILE_LIST
    #res@qc.data <- c(res@qc.data,lapply(QC_INTERNAL_SCALAQC_FILE_LIST, FUN=read.scalaqc.file.helper));
  }
  if(length(skip.files) > 0){
      USE.LIST <- USE.LIST[! USE.LIST %in% skip.files];
      if(debugMode) message("skipping ",length(skip.files)," files.");
      if(debugMode) message("   (\"",paste0(skip.files,collapse="\",\""),"\")");
  }
  
  
  read.scalaqc.file.helper <- function(scalaqc_file, sep=""){
    if(debugMode) message(paste0("(File ",which(USE.LIST == scalaqc_file)," of ",length(USE.LIST),"): ",scalaqc_file), appendLF=FALSE);
    if(debugMode) ts <- timestamp();
    out <- read.in.scalaQC.files(infile.dir,lanebam.list, qc.data.dir.list,scalaqc_file,sep=sep);
    if(debugMode) reportTimeAndDiff(ts,prefix="   ");
    out;
  }
  
  res@qc.data <- c(res@qc.data,lapply(USE.LIST, FUN=read.scalaqc.file.helper));

  res@calc.data <- list();
  
  message("calculating secondary data:");
  if(debugMode) ts <- timestamp();
  res <- calc.results.data(res, calc.DESeq2 = calc.DESeq2, calc.edgeR = calc.edgeR,debugMode=debugMode);
  message("done.");
  if(debugMode) reportTimeAndDiff(ts);
  
  return(res);
}

##################################################################################################################################
#LIST OF FILES:
##################################################################################################################################

QC_INTERNAL_SCALAQC_SUMMARY_FILE <- "QC.summary.txt";

QC_INTERNAL_SCALAQC_FILE_LIST <- list( gc.byPair = "QC.gc.byPair.txt.gz",
                                       gc.byRead = "QC.gc.byRead.txt.gz",
                                       gc.byRead.vsBaseCt = "QC.gc.byRead.vsBaseCt.txt.gz",
                                       quals.r1 = "QC.quals.r1.txt.gz", 
                                       quals.r2 = "QC.quals.r2.txt.gz",
                                       cigarOpDistribution.byReadCycle.R1 = "QC.cigarOpDistribution.byReadCycle.R1.txt.gz",
                                       cigarOpDistribution.byReadCycle.R2 = "QC.cigarOpDistribution.byReadCycle.R2.txt.gz",
                                       cigarOpLengths.byOp.R1 = "QC.cigarOpLengths.byOp.R1.txt.gz",
                                       cigarOpLengths.byOp.R2 = "QC.cigarOpLengths.byOp.R2.txt.gz",
                                       geneBodyCoverage.by.expression.level = "QC.geneBodyCoverage.by.expression.level.txt.gz",
                                       geneCounts = "QC.geneCounts.txt.gz",
                                       #NVC.raw = "scalaQC.ignoreClipping.NVC.raw.txt.gz",
                                       insert.size = "QC.insert.size.txt.gz",
                                       NVC.raw.R1 = "QC.NVC.raw.R1.txt.gz",
                                       NVC.raw.R2 = "QC.NVC.raw.R2.txt.gz",
                                       NVC.lead.clip.R1 = "QC.NVC.lead.clip.R1.txt.gz",
                                       NVC.lead.clip.R2 = "QC.NVC.lead.clip.R2.txt.gz",
                                       NVC.tail.clip.R1 = "QC.NVC.tail.clip.R1.txt.gz",
                                       NVC.tail.clip.R2 = "QC.NVC.tail.clip.R2.txt.gz",
                                       NVC.minus.clipping.R1 = "QC.NVC.minus.clipping.R1.txt.gz",
                                       NVC.minus.clipping.R2 = "QC.NVC.minus.clipping.R2.txt.gz",
                                       chrom.counts = "QC.chromCount.txt.gz",
                                       biotype.counts = "QC.biotypeCounts.txt.gz",
                                       geneBodyCoverage.pct = "QC.geneBodyCoverage.byExpr.avgPct.txt.gz",
                                       #onTarget = "QC.onTarget.txt.gz",
                                       overlapCoverage = "QC.overlapCoverage.txt.gz",
                                       #overlapMismatch = "QC.overlapMismatch.txt.gz",
                                       overlapMismatch.byRead = "QC.overlapMismatch.byRead.txt.gz",
                                       overlapMismatch.byScore = "QC.overlapMismatch.byScore.txt",
                                       overlapMismatchCombos = "QC.overlapMismatch.byBase.txt.gz",
                                       overlapMismatch.byScoreAndBP = "QC.overlapMismatch.byScoreAndBP.txt.gz",
                                       readLenDist = "QC.readLenDist.txt.gz",
                                       referenceMismatchCounts="QC.referenceMismatchCounts.txt.gz",
                                       referenceMismatchRaw.byReadStrand="QC.referenceMismatchRaw.byReadStrand.txt",
                                       referenceMismatch.byScore="QC.referenceMismatch.byScore.txt",
                                       referenceMismatch.byScoreAndBP="QC.referenceMismatch.byScoreAndBP.txt",
                                       mismatchSizeRates="QC.mismatchSizeRates.txt.gz",
                                       FQ.gc.byRead = "QC.FQ.gc.byRead.txt.gz",
                                       FQ.gc.byPair = "QC.FQ.gc.byPair.txt.gz",
                                       FQ.gc.r1 =     "QC.FQ.gc.R1.txt.gz",
                                       FQ.gc.r2 =     "QC.FQ.gc.R2.txt.gz",
                                       FQ.NVC.R1 =    "QC.FQ.NVC.R1.txt.gz",
                                       FQ.NVC.R2 =    "QC.FQ.NVC.R2.txt.gz",
                                       FQ.quals.r1 =  "QC.FQ.quals.r1.txt.gz",
                                       FQ.quals.r2 =  "QC.FQ.quals.r2.txt.gz"
                                       
                                       #,
                                       #spliceJunctionCounts.knownSplices = "scalaQC.spliceJunctionCounts.knownSplices.txt.gz"
                                       #spliceJunctionCounts.novelSplices = "scalaQC.spliceJunctionCounts.novelSplices.txt.gz"
                                       #
);

QC_INTERNAL_SCALAQC_FILE_LIST_SINGLE_END <- list(
                                       gc.byRead = "QC.gc.byRead.txt.gz",
                                       gc.byRead.vsBaseCt = "QC.gc.byRead.vsBaseCt.txt.gz",
                                       quals.r1 = "QC.quals.r1.txt.gz",
                                       cigarOpDistribution.byReadCycle.R1 = "QC.cigarOpDistribution.byReadCycle.R1.txt.gz",
                                       cigarOpLengths.byOp.R1 = "QC.cigarOpLengths.byOp.R1.txt.gz",
                                       geneBodyCoverage.by.expression.level = "QC.geneBodyCoverage.by.expression.level.txt.gz",
                                       geneCounts = "QC.geneCounts.txt.gz",
                                       NVC.raw.R1 =            "QC.NVC.raw.R1.txt.gz",
                                       NVC.lead.clip.R1 =      "QC.NVC.lead.clip.R1.txt.gz",
                                       NVC.tail.clip.R1 =      "QC.NVC.tail.clip.R1.txt.gz",
                                       NVC.minus.clipping.R1 = "QC.NVC.minus.clipping.R1.txt.gz",
                                       chrom.counts = "QC.chromCount.txt.gz",
                                       biotype.counts = "QC.biotypeCounts.txt.gz",
                                       geneBodyCoverage.pct = "QC.geneBodyCoverage.byExpr.avgPct.txt.gz",
                                       #onTarget = "QC.onTarget.txt.gz",
                                       readLenDist = "QC.readLenDist.txt.gz",
                                       referenceMismatchCounts="QC.referenceMismatchCounts.txt.gz",
                                       referenceMismatchRaw.byReadStrand="QC.referenceMismatchRaw.byReadStrand.txt",
                                       referenceMismatch.byScore="QC.referenceMismatch.byScore.txt",
                                       referenceMismatch.byScoreAndBP="QC.referenceMismatch.byScoreAndBP.txt",
                                       mismatchSizeRates="QC.mismatchSizeRates.txt.gz",
                                       FQ.gc.byRead = "QC.FQ.gc.byRead.txt.gz",
                                       FQ.NVC.R1 =    "QC.FQ.NVC.R1.txt.gz",
                                       FQ.quals.r1 =  "QC.FQ.quals.r1.txt.gz"
                                       #,
                                       #spliceJunctionCounts.knownSplices = "scalaQC.spliceJunctionCounts.knownSplices.txt.gz"
                                       #spliceJunctionCounts.novelSplices = "scalaQC.spliceJunctionCounts.novelSplices.txt.gz"
                                       #
);

##################################################################################################################################
##################################################################################################################################

find.compression.variant <- function(f){
  sapply(f,find.compression.variant.helper);
}

find.compression.variant.helper <- function(f){
  if(file.exists(f)){
    return(f);
  }
  if(substr(f,nchar(f)-2,nchar(f)) == ".gz"){
    f <- substr(f,1,nchar(f)-3);
  }
  if(substr(f,nchar(f)-3,nchar(f)) == ".zip"){
    f <- substr(f,1,nchar(f)-4);
  }
  
  if(file.exists(f)){
    return(f);
  } else if(file.exists(paste0(f,".gz"))){
    return(paste0(f,".gz"));
  } else if(file.exists(paste0(f,".zip"))){
    return(paste0(f,".zip"));
  } else {
    return(NA);
  }
}

read.in.scalaQC.files <- function(infile.prefix, lanebam.list, qc.data.dir.list, infile.suffix, sep = ""){
  #message(paste0("reading ",infile.suffix," files"),appendLF=FALSE);
  infiles <- find.compression.variant(paste0(infile.prefix,unlist(qc.data.dir.list),"/", infile.suffix));
  if(! is.na(infiles[1])){
    #print("!")
    for(i in 1:length(infiles)){
      if(is.na(infiles[i])){
        stop("File not found: ",infiles[i], "()");
      }
    }
    decade.ct <- 0;
    laneBamCt <- length(unlist(lanebam.list))
    #if(laneBamCt >= 10){ message(" [", appendLF=FALSE) }
    out <- lapply(lanebam.list,FUN=function(unique.ID){
      i <- which(unlist(lanebam.list) == unique.ID);
      if( floor(10 * i / laneBamCt) > decade.ct ){ 
          message(".", appendLF=FALSE); 
          decade.ct <<- floor(10 * i / laneBamCt) 
      }
      d <- tryCatch({
         read.table(infiles[i],header=T,stringsAsFactors=F,sep=sep,quote="");
      }, warning = function(w){
         message(paste0("Warning tripped on file: ",infiles[i]));
         #message(w);
      }, error = function(e){
         message(paste0("Error tripped on file: ",infiles[i]));
         stop(e);
      }, finally = {
         #do nothing.
      });
      return(d);
    })
    #if(laneBamCt >= 10){ message("] ", appendLF=FALSE) }
    message("done.");
    return(out);
  } else {
    message(paste0("Failed: Cannot find file: ",paste0(infile.prefix,unlist(qc.data.dir.list),"/", infile.suffix)[[1]],". Skipping tests that use this data."));
    return(NULL);
  }
}


#lanebam.list="character",
#sample.list="character",
#lane.list="character",
#group.list="character",
#decoder="data.frame", #decoder has columns: unique.ID	sample.ID	lane.ID	group.ID	cycle.CT	and then any number of user-defined columns (which are ignored internally)
#qc.data="list" #List of Lists. Each element corresponds to one qc test, and is composed of a list, one element for each lanebam.

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

check.all.completed.without.error <- function(res){
   summaries <- res@qc.data[["summary"]];
   if(is.null(summaries)){
     stop("FATAL ERROR: Samplerun summaries not found! Re-run QC jar utility!");
   }
   anyFailed <- c();
   
   for(lanebam in names(summaries)){
      currSummary <- summaries[[lanebam]];
      if(is.null(currSummary)){
         anyFailed <- c(anyFailed,lanebam);
      } else if(! ("COMPLETED_WITHOUT_ERROR" %in% currSummary$FIELD)){
         anyFailed <- c(anyFailed,lanebam);
      } else if(currSummary$COUNT[currSummary$FIELD == "COMPLETED_WITHOUT_ERROR"] != 1){
         anyFailed <- c(anyFailed,lanebam);
      }
   }
   
   if(length(anyFailed) > 0){
      for(lanebam in anyFailed){
        message("NOTE: SampleRun \"",lanebam,"\" marked incomplete!");
      }
      stop("FATAL ERROR: SampleRuns incomplete! Re-run QC jar utility for the noted samples above!");
   }
}

calc.results.data <- function(res, calc.DESeq2 = FALSE, calc.edgeR = FALSE, debugMode = FALSE){
   tryCatch({
     
     res <- fix.summary.to.numeric(res);
     check.all.completed.without.error(res);
     
     runTimedFunction(debugMode=debugMode,title="Calculating Quality Score Rates",expr={
       res <- calc.quals(res);
     })
     runTimedFunction(debugMode=debugMode,title="Calculating cumulative gene coverage, by replicate",expr={
       res <- calc.gene.CDF(res);
     })
     
     runTimedFunction(debugMode=debugMode,title="Calculating cumulative gene coverage, by sample",expr={
       res <- calc.gene.CDF.bySample(res);
     })
     #res <- calc.raw.NVC(res);
     runTimedFunction(debugMode=debugMode,title="Calculating Mapping Rates",expr={
       res <- calc.mapping.rates(res);
     })
     runTimedFunction(debugMode=debugMode,title="calculating normalization factors, by sample",expr={
       res <- calc.samplewise.norm.factors(res, calc.DESeq2,calc.edgeR);
     })
     runTimedFunction(debugMode=debugMode,title="calculating normalization factors, by replicate",expr={
       res <- calc.lanebamwise.norm.factors(res, calc.DESeq2,calc.edgeR);
     })
     runTimedFunction(debugMode=debugMode,title="calculating normalization factors, by sample/replicate",expr={
       res <- calc.lanebamwise.bysample.norm.factors(res);
     })
     runTimedFunction(debugMode=debugMode,title="Calculating summary stats",expr={
       res <- add.to.summary(res);
     })
     res <- calc.overlap.data(res,debugMode=debugMode);
     
     res <- calc.target.data(res,debugMode=debugMode);
     res <- calc.refmatch.data(res,debugMode=debugMode);
     
     runTimedFunction(debugMode=debugMode,title="Calculating summary table",expr={
       res@summaryTable <- get.summary.table(res);
     })
     
     runTimedFunction(debugMode=debugMode,title="Calculating overlap mismatch combos",expr={
       res <- add.overlap.mismatch.combos(res,debugMode=debugMode);
     })
     res <- calc.NVC.rates(res,debugMode=debugMode);
     
     return(res);
   },error = function(e){
        message("caught error: ",e);
        return(res);
   });
}


#overlapCoverage = "QC.overlapCoverage.txt.gz",
#overlapMismatch = "QC.overlapMismatch.txt.gz",
#overlapMismatch.byRead = "QC.overlapMismatch.byRead.txt.gz",
#readLenDist = "QC.readLenDist.txt.gz"


     mismatchCombos <- list(
c("A","C"),
c("T","G"),
c("C","A"),
c("G","T"),

c("A","G"),
c("T","C"),
c("G","A"),
c("C","T"),

c("A","T"),
c("T","A"),
c("C","G"),
c("G","C")
);

calc.refmatch.data <- function(res,debugMode = TRUE){
   tryCatch({
     #aa
     #referenceMismatchCounts
     if(! is.null(res@qc.data[["referenceMismatchCounts"]])){
       runTimedFunction(debugMode=debugMode,title="Calculating referenceMismatchCounts stats",expr={
       nameList <- as.list(names(res@qc.data[["referenceMismatchCounts"]]));
       names(nameList) <- nameList;
     
       res@calc.data[["referenceMismatchCounts"]] <- lapply(nameList,function(n){
         old <- res@qc.data[["referenceMismatchCounts"]][[n]];
         sumdata <- res@qc.data[["summary"]][[n]];
         #val bpcountR1 = if(! is.null(res@qc.data[["readLenDist"]])){
         #  sum( as.numeric(res@qc.data[["readLenDist"]][[n]]$CT_R1) * as.numeric(res@qc.data[["readLenDist"]][[n]]$LEN) )
         #} else {
         #  as.numeric(sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"] * sumdata$COUNT[sumdata$FIELD == "maxObservedReadLength"]);
         #}
         #val bpcountR1 = if(! is.null(res@qc.data[["readLenDist"]])){
         #  sum( as.numeric(res@qc.data[["readLenDist"]][[n]]$CT_R2) * as.numeric(res@qc.data[["readLenDist"]][[n]]$LEN) )
         #} else {
         #  as.numeric(sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"] * sumdata$COUNT[sumdata$FIELD == "maxObservedReadLength"]);
         #}
         
         readCt  <- sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
         
         dl <- old;
                            dl$RATE_R1 = 100 * dl$CT_R1 / dl$MAPPEDCT_R1;
         if(!res@singleEnd) dl$RATE_R2 = 100 * dl$CT_R2 / dl$MAPPEDCT_R2;
         
                            dl$CT_R1_LOG <- log10(dl$CT_R1);
         if(!res@singleEnd) dl$CT_R2_LOG <- log10(dl$CT_R2);
         
                            dl$RATE_R1_LOG <- log10(dl$RATE_R1)
         if(!res@singleEnd) dl$RATE_R2_LOG <- log10(dl$RATE_R2)
         
         dl;
       });
       });
     }
     
     if(! is.null(res@qc.data[["referenceMismatch.byScore"]])){
       runTimedFunction(debugMode=debugMode,title="Calculating referenceMismatch.byScore stats",expr={

       nameList <- as.list(names(res@qc.data[["referenceMismatch.byScore"]]));
       names(nameList) <- nameList;
       
       
       res@calc.data[["referenceMismatch.byScore"]] <- lapply(nameList,function(n){
         old <- res@qc.data[["referenceMismatch.byScore"]][[n]];
         
         dl <- old[old$QualCoverage_R1 > 0,];
         
         dl$RATE_R1 = 100 * dl$MismatchCt_R1 / dl$QualCoverage_R1
         
         if(!res@singleEnd) dl$RATE_R2 = 100 * dl$MismatchCt_R2 / dl$QualCoverage_R2;
         
                            dl$CT_R1_LOG <- log10(dl$MismatchCt_R1);
         if(!res@singleEnd) dl$CT_R2_LOG <- log10(dl$MismatchCt_R2);
         
                            dl$RATE_R1_LOG <- log10(dl$RATE_R1)
         if(!res@singleEnd) dl$RATE_R2_LOG <- log10(dl$RATE_R2)

         dl;
       });
       })
     }
     
   if(! is.null(res@qc.data[["referenceMismatchRaw.byReadStrand"]])){
     runTimedFunction(debugMode=debugMode,title="Calculating referenceMismatchRaw.byReadStrand stats",expr={

     nameList <- as.list(names(res@qc.data[["referenceMismatchRaw.byReadStrand"]]));
     names(nameList) <- nameList;
     
     res@calc.data[["referenceMismatchCombos"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["referenceMismatchRaw.byReadStrand"]][[n]];
       sumdata <- res@qc.data[["summary"]][[n]];
       readCt <- sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
       old$REFBASE <- as.character(old$REFBASE);
       old$READBASE <- as.character(old$READBASE);
       
       dl <- data.frame(
         refBase  = character(),
         readBase  = character(),
         combo = character(),
         CT_R1_FWD = numeric(),
         CT_R1_REV = numeric(),
         CT_R1 = numeric(),
         CT_R2_FWD = numeric(),
         CT_R2_REV = numeric(),
         CT_R2 = numeric(),stringsAsFactors=FALSE
       );
       for(mm in mismatchCombos){
         i <- nrow(dl)+1;
         A <- as.character(mm[1]);
         B <- as.character(mm[2]);
         dl[i,"refBase"] <- A
         dl[i,"readBase"] <- B
         dl[i,"combo"] <- paste0(B,"->",A);
         dl[i,"CT_R1_FWD"] <- sum(old$CT_R1_FWD[ old$REFBASE == A & old$READBASE == B]);
         dl[i,"CT_R1_REV"] <- sum(old$CT_R1_REV[ old$REFBASE == A & old$READBASE == B]);
         dl[i,"CT_R1"] <- sum(old$CT_R1[         old$REFBASE == A & old$READBASE == B]);
         dl[i,"CT_R2_FWD"] <- sum(old$CT_R2_FWD[ old$REFBASE == A & old$READBASE == B]);
         dl[i,"CT_R2_REV"] <- sum(old$CT_R2_REV[ old$REFBASE == A & old$READBASE == B]);
         dl[i,"CT_R2"] <- sum(old$CT_R2[         old$REFBASE == A & old$READBASE == B]);
         
         #dl[i,"RATE_R1"] <- dl[i,"CT_R1"] / readCt;
         #dl[i,"RATE_R2"] <- dl[i,"CT_R2"] / readCt;
       }
       dl$CT <- dl$CT_R1 + dl$CT_R2;
       dl$RATE_R1 <- dl$CT_R1 / readCt;
       dl$RATE_R2 <- dl$CT_R2 / readCt;
       dl$RATE <- dl$CT / (readCt * 2);
       dl;
     });
     })
   }
     
   if(! is.null(res@qc.data[["referenceMismatch.byScoreAndBP"]])){
     runTimedFunction(debugMode=debugMode,title="Calculating referenceMismatch.byScoreAndBP stats",expr={
     nameList <- as.list(names(res@qc.data[["referenceMismatch.byScoreAndBP"]]));
     names(nameList) <- nameList;
     
     res@calc.data[["referenceMismatch.byScoreAndBP"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["referenceMismatch.byScoreAndBP"]][[n]];
       sumdata <- res@qc.data[["summary"]][[n]];
       readCt <- sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
       
       #Score REFBASE READBASE  CT_R1_FWD CT_R1_REV CT_R2_FWD CT_R2_REV CT_R1 CT_R2
       names(old)[names(old) == "REFBASE"] <- "refBase";
       names(old)[names(old) == "READBASE"] <- "readBase";
       old$CT <- old$CT_R1 + old$CT_R2;
       old$combo <- paste0(as.character(old$refBase),"->",as.character(old$readBase))
       
       sumR1 <- sapply(1:nrow(old),function(i){
         sum(old$CT_R1[ old$Score == old$Score[[i]] ])
       })
       sumR2 <- sapply(1:nrow(old),function(i){
         sum(old$CT_R2[ old$Score == old$Score[[i]] ])
       })
       
       old$RATE_R1 <- 100 * old$CT_R1 / sumR1;
       old$RATE_R2 <- 100 * old$CT_R2 / sumR2;
       old$RATE    <- 100 * old$CT / (sumR1 + sumR2);
       
       old$CT_R1_LOG <- log10(old$CT_R1);
       old$CT_R2_LOG <- log10(old$CT_R2);
       old$CT_LOG <- log10(old$CT);
       old$RATE_R1_LOG <- log10(old$RATE_R1);
       old$RATE_R2_LOG <- log10(old$RATE_R2);
       old$RATE_LOG <- log10(old$RATE);
       old;
     });
     })
   }

     return(res)
   },error = function(e){
     message("caught error: ",e);
     return(res);
   });

}

calc.target.data <- function(res,debugMode = TRUE){
   nameList <- as.list(names(res@qc.data[["summary"]]));
   names(nameList) <- nameList;
   
   tryCatch({
   if(! is.null(res@qc.data[["onTarget"]])){
     
     res@calc.data[["onTargetCoverage"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["onTarget"]][[n]];
       sumdata <- res@qc.data[["summary"]][[n]];
       readCt  <- sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
       
       coverageDepth <- old$READPAIR_COVERAGE / (old$END - old$START)
       coverageDepth <- coverageDepth[order(coverageDepth)];
       dl <- data.frame(
         quantile = 100 * 1:length(coverageDepth) / length(coverageDepth),
         pairDepth = coverageDepth
       )
       coverageDepth <- old$READ_COVERAGE / (old$END - old$START)
       coverageDepth <- coverageDepth[order(coverageDepth)];
       dl$readDepth <- coverageDepth
       
       dl$logPairDepth <- log10(dl$pairDepth)
       dl$logReadDepth <- log10(dl$readDepth)
       
       dl;
     });
     
     maxPairCoverage <- max(sapply(res@calc.data[["onTargetCoverage"]],function(dl){ max(dl$pairDepth) }))
     maxReadCoverage <- max(sapply(res@calc.data[["onTargetCoverage"]],function(dl){ max(dl$readDepth) }))
     pairBreaks <- seq(0,maxPairCoverage+1);
     readBreaks <- seq(0,maxReadCoverage+1);
     res@calc.data[["maxPairCoverage"]] <- maxPairCoverage
     res@calc.data[["maxReadCoverage"]] <- maxReadCoverage
     
     res@calc.data[["onTargetCoverageBinsByRead"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["onTarget"]][[n]];
       sumdata <- res@qc.data[["summary"]][[n]];
       readCt  <- sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"];
       readDepth <- old$READ_COVERAGE / (old$END - old$START)
       readHist <- hist(readDepth,breaks=readBreaks,plot=F)
       dl <- data.frame(
         READBIN=readBreaks[-length(readBreaks)],
         readCt =readHist$counts,
         logReadCt = log10(readHist$counts)
       )
       dl;
     });
     res@calc.data[["onTargetCoverageBinsByPair"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["onTarget"]][[n]];
       sumdata <- res@qc.data[["summary"]][[n]];
       readCt  <- sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"];
       pairDepth <- old$READPAIR_COVERAGE / (old$END - old$START)
       pairHist <- hist(pairDepth,breaks=pairBreaks,plot=F)
       dl <- data.frame(
         PAIRBIN=pairBreaks[-length(pairBreaks)],
         pairCt =pairHist$counts,
         logPairCt = log10(pairHist$counts)
       )
       dl;
     });
     totalSpan <- sum(as.numeric(res@qc.data[["onTarget"]][[1]]$END - res@qc.data[["onTarget"]][[1]]$START));
     
     
     runTimedFunction(debugMode=debugMode,title="Calculating on-target summary info...",expr={
         res@qc.data[["summary"]] <- lapply(nameList,function(n){
           old <- res@qc.data[["onTarget"]][[n]];
           sumdata <- res@qc.data[["summary"]][[n]];
           readDepth <- sum(as.numeric(old$READ_COVERAGE)) / totalSpan;
           pairDepth <- sum(as.numeric(old$READPAIR_COVERAGE)) / totalSpan;
           readCt  <- sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"];

           if(! any(sumdata$FIELD == "meanOnTargetReadDepth")) sumdata[nrow(sumdata)+1,"FIELD"] <- "meanOnTargetReadDepth"
           if(! any(sumdata$FIELD == "meanOnTargetReadPairDepth")) sumdata[nrow(sumdata)+1,"FIELD"] <- "meanOnTargetReadPairDepth"
           if(! any(sumdata$FIELD == "readDepthPerMillion")) sumdata[nrow(sumdata)+1,"FIELD"] <- "readDepthPerMillion"
           if(! any(sumdata$FIELD == "pairDepthPerMillion")) sumdata[nrow(sumdata)+1,"FIELD"] <- "pairDepthPerMillion"
           
           sumdata[sumdata$FIELD == "meanOnTargetReadDepth",    "COUNT"]   <- readDepth
           sumdata[sumdata$FIELD == "meanOnTargetReadPairDepth", "COUNT"]   <- pairDepth
           
           sumdata[sumdata$FIELD == "readDepthPerMillion", "COUNT"]   <- readDepth / (readCt / 1000000);
           sumdata[sumdata$FIELD == "pairDepthPerMillion", "COUNT"]   <- pairDepth / (readCt / 1000000);

           sumdata
         });
     })
   }
     totalSpan <- res@qc.data[["summary"]][[1]][res@qc.data[["summary"]][[1]]$FIELD == "TotalTargetSpan", "COUNT"];
     
     if(all(c("fractionOfPairsOnTarget","fractionOfPairsOffTarget","readDepthPerMillion","pairDepthPerMillion","READ_PAIR_OK","OnTargetReadPairBases","OnTargetReadBases") %in% res@qc.data[["summary"]]$FIELD)) {
       runTimedFunction(debugMode=debugMode,title="Calculating more on-target summary info...",expr={
           res@qc.data[["summary"]] <- lapply(nameList,function(n){
             old <- res@qc.data[["onTarget"]][[n]];
             sumdata <- res@qc.data[["summary"]][[n]];

             readDepth <- sumdata[sumdata$FIELD == "OnTargetReadBases","COUNT"]     / totalSpan;
             pairDepth <- sumdata[sumdata$FIELD == "OnTargetReadPairBases","COUNT"] / totalSpan;
             readCt  <- sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"];

             if(! any(sumdata$FIELD == "fractionOfPairsOnTarget")) sumdata[nrow(sumdata)+1,"FIELD"] <- "fractionOfPairsOnTarget"
             if(! any(sumdata$FIELD == "fractionOfPairsOffTarget")) sumdata[nrow(sumdata)+1,"FIELD"] <- "fractionOfPairsOffTarget"
             if(! any(sumdata$FIELD == "readDepthPerMillion")) sumdata[nrow(sumdata)+1,"FIELD"] <- "readDepthPerMillion"
             if(! any(sumdata$FIELD == "pairDepthPerMillion")) sumdata[nrow(sumdata)+1,"FIELD"] <- "pairDepthPerMillion"

             sumdata[sumdata$FIELD == "fractionOfPairsOnTarget", "COUNT"]   <- sumdata[sumdata$FIELD == "OnTargetCount", "COUNT"] / (sumdata[sumdata$FIELD == "OnTargetCount", "COUNT"]+sumdata[sumdata$FIELD == "OffTargetCount", "COUNT"]);
             sumdata[sumdata$FIELD == "fractionOfPairsOffTarget", "COUNT"]   <- sumdata[sumdata$FIELD == "OffTargetCount", "COUNT"] / (sumdata[sumdata$FIELD == "OnTargetCount", "COUNT"]+sumdata[sumdata$FIELD == "OffTargetCount", "COUNT"]);
             sumdata[sumdata$FIELD == "readDepthPerMillion", "COUNT"]   <- readDepth / (readCt / 1000000);
             sumdata[sumdata$FIELD == "pairDepthPerMillion", "COUNT"]   <- pairDepth / (readCt / 1000000);

             sumdata
           });
       })
    }
   
   
   return(res);
   },error = function(e){
     message("caught error: ",e);
     return(res);
   });
}

calc.overlap.data <- function(res,debugMode = TRUE){

   tryCatch({
   if(! is.null(res@qc.data[["mismatchSizeRates"]])){
     nameList <- as.list(names(res@qc.data[["overlapCoverage"]]));
     names(nameList) <- nameList;
     runTimedFunction(debugMode=debugMode,title="Calculating overlap mismatch-size rates",expr={
     
     res@calc.data[["mismatchSizeRates"]] <- lapply(nameList,function(n){
            #OVERLAP_SIZE	BASES_MISMATCHED	CT_NOINDEL	CT
            dl <- res@qc.data[["mismatchSizeRates"]][[n]];
            sumdata <- res@qc.data[["summary"]][[n]];
            
            out <- data.frame(
              BASES_MISMATCHED = 1:sumdata$COUNT[sumdata$FIELD == "maxObservedReadLength"]
            );
            out$CT <- sapply(out$BASES_MISMATCHED,function(bmm){
              sum(dl$CT[dl$BASES_MISMATCHED == bmm]);
            })
            out$CT_NOINDEL <- sapply(out$BASES_MISMATCHED,function(bmm){
              sum(dl$CT[dl$BASES_MISMATCHED == bmm]);
            })
            
            out$RATE <- 100 * out$CT / sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
            out$RATE_NOINDEL <- 100 * out$CT_NOINDEL / sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
            
            out$RATE_LOG <- log10(out$RATE);
            out$RATE_NOINDEL_LOG <- log10(out$RATE_NOINDEL);
            out;
     });
     })
   }
   if(! is.null(res@qc.data[["mismatchSizeRates"]])){
     nameList <- as.list(names(res@qc.data[["overlapCoverage"]]));
     names(nameList) <- nameList;
     #if(debugMode) message("   Calculating cumulative overlap mismatch-size rates...");
     runTimedFunction(debugMode=debugMode,title="Calculating cumulative overlap mismatch-size rates",expr={

     res@calc.data[["mismatchSizeRatesCumulative"]] <- lapply(nameList,function(n){
            #OVERLAP_SIZE   BASES_MISMATCHED    CT_NOINDEL  CT
            dl <- res@qc.data[["mismatchSizeRates"]][[n]];
            sumdata <- res@qc.data[["summary"]][[n]];
            
            out <- data.frame(
              BASES_MISMATCHED = sumdata$COUNT[sumdata$FIELD == "maxObservedReadLength"]:1
            );
            out$CT <- sapply(out$BASES_MISMATCHED,function(bmm){
              sum(dl$CT[dl$BASES_MISMATCHED >= bmm]);
            })
            out$CT_NOINDEL <- sapply(out$BASES_MISMATCHED,function(bmm){
              sum(dl$CT[dl$BASES_MISMATCHED >= bmm]);
            })
            
            out$RATE <- 100 * out$CT / sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
            out$RATE_NOINDEL <- 100 * out$CT_NOINDEL / sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
            
            out$RATE_LOG <- log10(out$RATE);
            out$RATE_NOINDEL_LOG <- log10(out$RATE_NOINDEL);
            out;
     });
     })
   }
   
   if(! is.null(res@qc.data[["overlapCoverage"]])){
     nameList <- as.list(names(res@qc.data[["overlapCoverage"]]));
     names(nameList) <- nameList;
     #if(debugMode) message("   Calculating overlap coverage Rates...");
     runTimedFunction(debugMode=debugMode,title="Calculating overlap coverage Rates",expr={

     buf <- lapply(nameList,function(n){
       dl <- res@qc.data[["overlapCoverage"]][[n]];
       sumdata <- res@qc.data[["summary"]][[n]];
       dl$RATE_R1 = 100 * dl$CT_R1 / sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
       dl$RATE_R2 = 100 * dl$CT_R2 / sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"]
     });
     res@calc.data[["overlapCoverage"]] <- buf;
     
     res@calc.data[["overlapR1"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["overlapCoverage"]][[n]]
       sumdata <- res@qc.data[["summary"]][[n]];
       dl <- data.frame(
         POS  = old$POS,
         CT   = old$CT_R1,
         RATE = old$CT_R1 / (sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"] / 100)
       );
     });
     res@calc.data[["overlapR2"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["overlapCoverage"]][[n]];
       sumdata <- res@qc.data[["summary"]][[n]];
       dl <- data.frame(
         POS  = old$POS,
         CT   = old$CT_R2,
         RATE = old$CT_R2 / (sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"] / 100)
       );
     });
     });
   }
   if((! is.null(res@qc.data[["overlapMismatch.byRead"]])) && (! is.null(res@qc.data[["overlapCoverage"]]))){
     nameList <- as.list(names(res@qc.data[["overlapMismatch.byRead"]]));
     names(nameList) <- nameList;
     runTimedFunction(debugMode=debugMode,title="Calculating overlap coverage Rates By Read",expr={

     res@calc.data[["overlapMismatchR1"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["overlapMismatch.byRead"]][[n]]
       covdata <- res@qc.data[["overlapCoverage"]][[n]];

       dl <- data.frame(
         POS  = old$POS,
         CT   = old$CT_R1,
         RATE = 100* old$CT_R1 / covdata$CT_R1,
         CT_NOIND = old$CT_NOINDEL_R1,
         RATE_NOIND = 100* old$CT_NOINDEL_R1 / covdata$CT_NOINDEL_R1
       );
       dl$CT_LOG = log10(dl$CT);
       dl$RATE_LOG = log10(dl$RATE);
       dl$CT_LOG_NOIND = log10(dl$CT_NOIND);
       dl$RATE_LOG_NOIND = log10(dl$RATE_NOIND);
       dl;
     });
     res@calc.data[["overlapMismatchR2"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["overlapMismatch.byRead"]][[n]]
       covdata <- res@qc.data[["overlapCoverage"]][[n]];

       dl <- data.frame(
         POS  = old$POS,
         CT   = old$CT_R2,
         RATE = 100* old$CT_R2 / covdata$CT_R2,
         CT_NOIND = old$CT_NOINDEL_R2,
         RATE_NOIND = 100* old$CT_NOINDEL_R2 / covdata$CT_NOINDEL_R2
       );
       dl$CT_LOG = log10(dl$CT);
       dl$RATE_LOG = log10(dl$RATE);
       dl$CT_LOG_NOIND = log10(dl$CT_NOIND);
       dl$RATE_LOG_NOIND = log10(dl$RATE_NOIND);
       dl;
     });
     
     })
   }

   if(! is.null(res@qc.data[["readLenDist"]])){
     nameList <- as.list(names(res@qc.data[["readLenDist"]]));
     names(nameList) <- nameList;
     #if(debugMode) message("   Calculating read length distribution...");
     runTimedFunction(debugMode=debugMode,title="Calculating read length distribution",expr={

     res@calc.data[["readLenDistR1"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["readLenDist"]][[n]]
       sumdata <- res@qc.data[["summary"]][[n]];
       dl <- data.frame(
         LEN  = old$LEN,
         CT   = old$CT_R1,
         RATE = old$CT_R1 / (sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"] / 100)
       );
       dl$CUMRATE = cumsum(dl$RATE);
       dl$CUMCT = cumsum(dl$CT);
       dl;
     });
     if(! res@singleEnd){
         res@calc.data[["readLenDistR2"]] <- lapply(nameList,function(n){
           old <- res@qc.data[["readLenDist"]][[n]]
           sumdata <- res@qc.data[["summary"]][[n]];
           dl <- data.frame(
             LEN  = old$LEN,
             CT   = old$CT_R2,
             RATE = old$CT_R2 / (sumdata$COUNT[sumdata$FIELD == "READ_PAIR_OK"] / 100)
           );
           dl$CUMRATE = cumsum(dl$RATE);
           dl$CUMCT = cumsum(dl$CT);
           dl;
         });
     }
     })
   }

   if(! is.null(res@qc.data[["overlapMismatch.byScore"]])){
     
     nameList <- as.list(names(res@qc.data[["overlapMismatch.byScore"]]));
     names(nameList) <- nameList;
     maxQual <- max(res@qc.data[["overlapMismatch.byScore"]][[1]]$ScoreR1);
     avgQualSeq <- seq(from = 0, to = maxQual,by=0.5)
     qualSeq <- seq(from = 0, to = maxQual);
     
     runTimedFunction(debugMode=debugMode,title="Calculating overlap by AVG score",expr={

     res@calc.data[["overlapMismatch.byScore.avg"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["overlapMismatch.byScore"]][[n]]
       old <- old[old$OverlapCt > 0,]
       dl <- data.frame(
         avgQual  = avgQualSeq
       );
       avg <- rowMeans(old[,c("ScoreR1","ScoreR2"),drop=FALSE])
       dl$OverlapCt         <- sapply(avgQualSeq,function(q){ sum(as.numeric(old$OverlapCt[avg == q ])) })
       dl$OverlapCtNoIndel  <- sapply(avgQualSeq,function(q){ sum(as.numeric(old$OverlapCtNoIndel[avg == q ])) })
       dl$MismatchCt        <- sapply(avgQualSeq,function(q){ sum(as.numeric(old$MismatchCt[avg == q ])) })
       dl$MismatchCtNoIndel <- sapply(avgQualSeq,function(q){ sum(as.numeric(old$MismatchCtNoIndel[avg == q ])) })
       dl$MismatchRate        <- 100 * dl$MismatchCt/dl$OverlapCt;
       dl$MismatchRateNoIndel <- 100 * dl$MismatchCtNoIndel/dl$OverlapCtNoIndel;
       dl$MismatchRateLog <- log10(dl$MismatchRate)
       dl$MismatchRateLogNoIndel <- log10(dl$MismatchRateNoIndel)
       dl <- dl[dl$OverlapCt > 0,]
       dl;
     });
     })
     
     
          
     runTimedFunction(debugMode=debugMode,title="Calculating overlap by MIN score",expr={
     res@calc.data[["overlapMismatch.byScore.min"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["overlapMismatch.byScore"]][[n]]
       old <- old[old$OverlapCt > 0,]
       #avg <- rowMeans(old[,c("ScoreR1","ScoreR2"),drop=FALSE])
       minqual <- pmin(old$ScoreR1,old$ScoreR2);
       
       maxQual <- max(res@qc.data[["overlapMismatch.byScore"]][[1]]$ScoreR1);
       
       dl <- data.frame(
         minQual  = 0:maxQual
       );
       
       dl$OverlapCt         <- sapply(dl$minQual,function(q){ sum(as.numeric(old$OverlapCt[minqual == q ])) })
       dl$OverlapCtNoIndel  <- sapply(dl$minQual,function(q){ sum(as.numeric(old$OverlapCtNoIndel[minqual == q ])) })
       dl$MismatchCt        <- sapply(dl$minQual,function(q){ sum(as.numeric(old$MismatchCt[minqual == q ])) })
       dl$MismatchCtNoIndel <- sapply(dl$minQual,function(q){ sum(as.numeric(old$MismatchCtNoIndel[minqual == q ])) })
       dl$MismatchRate        <- 100 * dl$MismatchCt/dl$OverlapCt;
       dl$MismatchRateNoIndel <- 100 * dl$MismatchCtNoIndel/dl$OverlapCtNoIndel;
       dl$MismatchRateLog <- log10(dl$MismatchRate)
       dl$MismatchRateLogNoIndel <- log10(dl$MismatchRateNoIndel)
       dl <- dl[dl$OverlapCt > 0,]
       dl;
     });
     })
     
     runTimedFunction(debugMode=debugMode,title="Adding Min score error to summary tables",expr={
       res@qc.data[["summary"]] <- lapply(nameList,function(n){
         old <- res@qc.data[["summary"]][[n]]
         dl <- res@calc.data[["overlapMismatch.byScore.min"]][[n]];
         for(i in 1:length(dl$minQual)){
           x <- dl$minQual[i];
           fieldKey <- paste0("overlapMismatchRate.byMinQual.noIndel_",as.character(x));
           if(! any(fieldKey == old$FIELD)) old[nrow(old)+1,"FIELD"] <- fieldKey;
           old[old$FIELD == fieldKey,"COUNT"] <- dl$MismatchRateNoIndel[i];
         }
         old;
       })
     })

     
     
     runTimedFunction(debugMode=debugMode,title="Calculating overlap by R1 score",expr={

     res@calc.data[["overlapMismatch.byScore.R1"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["overlapMismatch.byScore"]][[n]]
       old <- old[old$OverlapCt > 0,]
       dl <- data.frame(
         qual  = qualSeq
       );
       key <- old$ScoreR1;
       dl$OverlapCt         <- sapply(qualSeq,function(q){ sum(as.numeric(old$OverlapCt[key == q ])) })
       dl$OverlapCtNoIndel  <- sapply(qualSeq,function(q){ sum(as.numeric(old$OverlapCtNoIndel[key == q ])) })
       dl$MismatchCt        <- sapply(qualSeq,function(q){ sum(as.numeric(old$MismatchCt[key == q ])) })
       dl$MismatchCtNoIndel <- sapply(qualSeq,function(q){ sum(as.numeric(old$MismatchCtNoIndel[key == q ])) })
       dl$MismatchRate        <- 100 * dl$MismatchCt/dl$OverlapCt;
       dl$MismatchRateNoIndel <- 100 * dl$MismatchCtNoIndel/dl$OverlapCtNoIndel;
       dl$MismatchRateLog <- log10(dl$MismatchRate)
       dl$MismatchRateLogNoIndel <- log10(dl$MismatchRateNoIndel)
       dl <- dl[dl$OverlapCt > 0,]
       dl;
     });
     })
     
     runTimedFunction(debugMode=debugMode,title="Calculating overlap by R2 score",expr={
     res@calc.data[["overlapMismatch.byScore.R2"]] <- lapply(nameList,function(n){
       old <- res@qc.data[["overlapMismatch.byScore"]][[n]]
       old <- old[old$OverlapCt > 0,]
       dl <- data.frame(
         qual  = qualSeq
       );
       key <- old$ScoreR2;
       dl$OverlapCt         <- sapply(qualSeq,function(q){ sum(as.numeric(old$OverlapCt[key == q ])) })
       dl$OverlapCtNoIndel  <- sapply(qualSeq,function(q){ sum(as.numeric(old$OverlapCtNoIndel[key == q ])) })
       dl$MismatchCt        <- sapply(qualSeq,function(q){ sum(as.numeric(old$MismatchCt[key == q ])) })
       dl$MismatchCtNoIndel <- sapply(qualSeq,function(q){ sum(as.numeric(old$MismatchCtNoIndel[key == q ])) })
       dl$MismatchRate        <- 100 * dl$MismatchCt/dl$OverlapCt;
       dl$MismatchRateNoIndel <- 100 * dl$MismatchCtNoIndel/dl$OverlapCtNoIndel;
       dl$MismatchRateLog <- log10(dl$MismatchRate)
       dl$MismatchRateLogNoIndel <- log10(dl$MismatchRateNoIndel)
       dl <- dl[dl$OverlapCt > 0,]
       dl;
     });
     })
     
     getPairedError <- function(phred1,phred2){
       getPhredFromError(1-((1-getErrorFromPhred(phred1)) * (1-getErrorFromPhred(phred2))))
     }
   }
   
   return(res);
   },error = function(e){
     message("caught error: ",e);
     return(res);
   });
}

add.overlap.mismatch.combos <- function(res,debugMode = TRUE){
  scalaqc_file = "QC.overlapMismatch.byBase.txt.gz";

  lanebam.list <- res@lanebam.list;
  names(lanebam.list) <- lanebam.list;
  qc.data.dir.list <- res@decoder$qc.data.dir
  names(qc.data.dir.list) <- res@decoder$unique.ID;
  
  #if(is.null(res@qc.data[["overlapMismatchCombos"]])){
  #  runTimedFunction(debugMode=debugMode,title="Reading mismatch combo data:",expr={
  #    res@qc.data[["overlapMismatchCombos"]] <- read.in.scalaQC.files(infile.prefix=infile.dir,lanebam.list = lanebam.list, qc.data.dir.list = qc.data.dir.list,infile.suffix=scalaqc_file)
  #  });
  #}
  #names(res@qc.data[["overlapMismatchCombos"]]) <- lanebam.list
  
  if(! is.null(res@qc.data[["overlapMismatchCombos"]])){
    runTimedFunction(debugMode=debugMode,title="Calculating mismatch combo rates:",expr={
      res@calc.data[["overlapMismatchCombos"]] <- lapply(lanebam.list,function(n){
        old <- res@qc.data[["overlapMismatchCombos"]][[n]];
        covdata <- res@qc.data[["overlapCoverage"]][[n]];
        totalOverlap <- sum(as.numeric(covdata$CT_R1));
        totalOverlap_NOINDEL <- sum(as.numeric(covdata$CT_NOINDEL_R1));

        data.frame(
          r1base  = old$baseA,
          r2base  = old$baseB,
          combo = paste0( old$baseA,"/", old$baseB),
          Ct = old$CT,
          rate = 100* old$CT / totalOverlap,
          Ct_NOINDEL = old$CT_NOINDEL,
          rate_NOINDEL = 100* old$CT_NOINDEL / totalOverlap_NOINDEL,
          Ct_LOG = log10(old$CT),
          rate_LOG = log10(100* old$CT / totalOverlap),
          Ct_NOINDEL_LOG = log10(old$CT_NOINDEL),
          rate_NOINDEL_LOG = log10(100* old$CT_NOINDEL / totalOverlap_NOINDEL),
          stringsAsFactors = FALSE
        );
      });

      names(res@calc.data[["overlapMismatchCombos"]]) <- lanebam.list
    });
  }
  
  if(! is.null(res@qc.data[["overlapMismatch.byScoreAndBP"]])){
    nameList <- as.list(names(res@qc.data[["overlapMismatch.byScoreAndBP"]]));
    names(nameList) <- nameList;
    
    runTimedFunction(debugMode=debugMode,title="Calculating overlapMismatch.byScoreAndBP stats",expr={
      res@calc.data[["overlapMismatch.byScoreAndBP"]] <- lapply(nameList,function(n){
        #ScoreR1  ScoreR2 baseR1  baseR2  CT_FWD  CT_REV  CT_NOINDEL_FWD  CT_NOINDEL_REV  CT  CT_NOINDEL
        
        old <- res@qc.data[["overlapMismatch.byScoreAndBP"]][[n]];
        covdata <- res@qc.data[["overlapCoverage"]][[n]];
        totalOverlap <- sum(as.numeric(covdata$CT_R1));
        totalOverlap_NOINDEL <- sum(as.numeric(covdata$CT_NOINDEL_R1));
        
        old$combos <- paste0( old$baseR1,"/", old$baseR2);
        old$pmin <- pmin(old$ScoreR1,old$ScoreR2);
        allCombos <- unique(old$combo);
        df <- data.frame(
          score = rep(unique(old$ScoreR1),each=length(allCombos)),
          combo = rep(allCombos,length(unique(old$ScoreR1))),
          stringsAsFactors = FALSE
        )
        
        df$CT_MIN <- sapply(1:nrow(df),function(i){
          sum(old$CT[ old$combos == df$combo[[i]] & old$pmin == df$score[[i]] ])
        })
        df$CT_MIN_NOINDEL <- sapply(1:nrow(df),function(i){
          sum(old$CT_NOINDEL[ old$combos == df$combo[[i]] & old$pmin == df$score[[i]] ])
        })
        df$CT_PAIR <- sapply(1:nrow(df),function(i){
          sum(old$CT[ old$combos == df$combo[[i]] & old$ScoreR1 == df$score[[i]] & old$ScoreR2 == df$score[[i]]])
        })
        df$CT_PAIR_NOINDEL <- sapply(1:nrow(df),function(i){
          sum(old$CT_NOINDEL[ old$combos == df$combo[[i]] & old$ScoreR1 == df$score[[i]] & old$ScoreR2 == df$score[[i]]])
        })
        
        df$RATE_MIN <- 100 * df$CT_MIN / sapply(1:nrow(df),function(i){
          sum(df$CT_MIN[ df$score == df$score[[i]] ])
        });
        df$RATE_MIN_NOINDEL <- 100 * df$CT_MIN_NOINDEL / sapply(1:nrow(df),function(i){
          sum(df$CT_MIN_NOINDEL[ df$score == df$score[[i]] ])
        });
        df$RATE_PAIR <- 100 * df$CT_PAIR / sapply(1:nrow(df),function(i){
          sum(df$CT_PAIR[ df$score == df$score[[i]] ])
        });
        df$RATE_PAIR_NOINDEL <- 100 * df$CT_PAIR_NOINDEL / sapply(1:nrow(df),function(i){
          sum(df$CT_PAIR_NOINDEL[ df$score == df$score[[i]] ])
        });
        
        df$CT_MIN_LOG <- log10(df$CT_MIN)
        df$CT_MIN_NOINDEL_LOG <- log10(df$CT_MIN_NOINDEL)
        df$CT_PAIR_LOG <- log10(df$CT_PAIR)
        df$CT_PAIR_NOINDEL_LOG <- log10(df$CT_PAIR_NOINDEL)
        df$RATE_MIN_LOG <- log10(df$RATE_MIN)
        df$RATE_MIN_NOINDEL_LOG <- log10(df$RATE_MIN_NOINDEL)
        df$RATE_PAIR_LOG <- log10(df$RATE_PAIR)
        df$RATE_PAIR_NOINDEL_LOG <- log10(df$RATE_PAIR_NOINDEL)
        
        return(df);
      });
    });
  }
  
  
  return(res);
  
  #res@calc.data[["mismatchCombos"]] <- lapply(tempList,function(tl){ tl[["mismatchCombos"]] })
  #res@calc.data[["mismatchCombosNoIndel"]] <- lapply(tempList,function(tl){ tl[["mismatchCombosNoIndel"]] })
}


add.overlap.mismatch.combos.OLD <- function(res,infile.dir,debugMode = TRUE){
    
  #overlapMismatch = "QC.overlapMismatch.txt.gz",
  scalaqc_file = "QC.overlapMismatch.txt.gz";
  
   #if(! is.null(res@qc.data[["overlapMismatch"]])){
     nameList <- as.list(res@lanebam.list);
     names(nameList) <- nameList;
     infiles <- paste0(infile.dir,"/",res@decoder$qc.data.dir,"/","QC.overlapMismatch.txt.gz");
     names(infiles) <- names(nameList);
     outdirs <- paste0(infile.dir,"/",res@decoder$qc.data.dir,"/");
     names(outdirs) <- names(nameList);
     
     outfile1 <- paste0(outdirs,"/QC.mismatchCombos.txt.gz");
     names(outfile1) <- names(nameList);
     outfile2 <- paste0(outdirs,"/QC.mismatchCombosNoIndel.txt.gz");
     names(outfile2) <- names(nameList);

     runTimedFunction(debugMode=debugMode,title="Calculating mismatch combos",expr={

     message("\n[", appendLF=FALSE);
     tempList <- lapply(1:length(nameList),function(i){

         n <- nameList[[i]];
         ibak <- i;

         if(file.exists(outfile1[[n]])){
           dl.withIndels <- read.table(outfile1[[n]],sep='\t',stringsAsFactors=F,header=TRUE);
           message("-", appendLF=FALSE);
         } else {
           message(".",appendLF=FALSE);
           old <- read.table(infiles[n],header=T,stringsAsFactors=F);
           covdata <- res@qc.data[["overlapCoverage"]][[n]];
           totalOverlap <- sum(as.numeric(covdata$CT_R1));
           totalOverlap_NOINDEL <- sum(as.numeric(covdata$CT_NOINDEL_R1));
           
           dl <- data.frame(
             r1base  = character(),
             r2base  = character(),
             combo = character(),
             #fwdCt = numeric(),
             #revCt = numeric(),
             Ct = numeric(),
             #fwdRate = numeric(),
             #revRate = numeric(),
             rate = numeric(),
             stringsAsFactors = FALSE
           );
           for(mm in mismatchCombos){
             i <- nrow(dl)+1;
             dl[i,"r1base"] <- mm[1];
             dl[i,"r2base"] <- mm[2];
             dl[i,"combo"] <- paste0(mm[1],"/",mm[2]);
             #dl[i,"fwdCt"] <- sum(old$CT_FWD[ old$baseA == mm[1] & old$baseB == mm[2]]);
             #dl[i,"revCt"] <- sum(old$CT_REV[ old$baseA == mm[1] & old$baseB == mm[2]]);
             dl[i,"Ct"] <- sum(old$CT[ old$baseA == mm[1] & old$baseB == mm[2]]);
             #dl[i,"fwdRate"] <- 100* dl[i,"fwdCt"] / totalOverlap;
             #dl[i,"revRate"] <- 100* dl[i,"revCt"] / totalOverlap;
             dl[i,"rate"] <- 100* dl[i,"Ct"] / totalOverlap;
           }
           dl.withIndels <- dl;
           gz1 <- gzfile(paste0(outdirs[[n]],"/QC.mismatchCombos.txt.gz"), "w")
           write.table(dl.withIndels,gz1,quote=F,row.names=F,sep='\t');
           close(gz1);
         }
         if(file.exists(outfile2[[n]])){
           dl <- read.table(outfile2[[n]],sep='\t',stringsAsFactors=F,header=TRUE);
         } else {
           dl <- data.frame(
             r1base  = character(),
             r2base  = character(),
             combo = character(),
             #fwdCt = numeric(),
             #revCt = numeric(),
             Ct = numeric(),
             #fwdRate = numeric(),
             #revRate = numeric(),
             rate = numeric(),
             stringsAsFactors = FALSE
           );
           for(mm in mismatchCombos){
             i <- nrow(dl)+1;
             dl[i,"r1base"] <- mm[1];
             dl[i,"r2base"] <- mm[2];
             dl[i,"combo"] <- paste0(mm[1],"/",mm[2]);
             #dl[i,"fwdCt"] <- sum(old$CT_NOINDEL_FWD[ old$baseA == mm[1] & old$baseB == mm[2]]);
             #dl[i,"revCt"] <- sum(old$CT_NOINDEL_REV[ old$baseA == mm[1] & old$baseB == mm[2]]);
             dl[i,"Ct"] <- sum(old$CT_NOINDEL[ old$baseA == mm[1] & old$baseB == mm[2]]);
             #dl[i,"fwdRate"] <- 100* dl[i,"fwdCt"] / totalOverlap;
             #dl[i,"revRate"] <- 100* dl[i,"revCt"] / totalOverlap;
             dl[i,"rate"] <- 100* dl[i,"Ct"] / totalOverlap;
           }
           dl;
           gz2 <- gzfile(paste0(outdirs[[n]],"/QC.mismatchCombosNoIndel.txt.gz"), "w")
           write.table(dl,gz2,quote=F,row.names=F,sep='\t');
           close(gz2);
         }
         
         if(ibak %% 100 == 0){ message(paste0("] (",ibak,") \n["), appendLF=FALSE)
         } else if(ibak %% 10 == 0){ message("] [", appendLF=FALSE);
         } else if(ibak %% 5 == 0) message(" ", appendLF=FALSE);
         
         list(mismatchCombos = dl.withIndels, mismatchCombosNoIndel = dl);
     })
     message("] done!");
     res@calc.data[["mismatchCombos"]] <- lapply(tempList,function(tl){ tl[["mismatchCombos"]] })
     res@calc.data[["mismatchCombosNoIndel"]] <- lapply(tempList,function(tl){ tl[["mismatchCombosNoIndel"]] })
     names(res@calc.data[["mismatchCombos"]]) <- names(res@qc.data[["summary"]]);
     names(res@calc.data[["mismatchCombosNoIndel"]]) <- names(res@qc.data[["summary"]]);
     
     
     })
     
     runTimedFunction(debugMode=debugMode,title="Calculating mismatch combo rates",expr={
     res@calc.data[["mismatchCombo.rate"]] <- lapply(nameList,function(n){
       old <- res@calc.data[["mismatchCombos"]][[n]];
       dl <- data.frame(
         FIELD = c(old$combo),
         COUNT = c(old$rate),stringsAsFactors = FALSE
       );
     });
     })
     
     return(res);
   #}
   return(res);
 
}



fix.summary.to.numeric <- function(res){
   buf <- res@qc.data[["summary"]];
   buf <- lapply(buf,function(dl){
     dl$COUNT = round(as.numeric(dl$COUNT),2);
     dl;
   });
   res@qc.data[["summary"]] <- buf;
   
   return(res);
}

add.to.summary <- function(res){
   #message("debug 0");
   summary <- res@qc.data[["summary"]];
   summary.names <- names(summary);
   #plotter$res@qc.data[["summary"]][[1]]$FIELD
   if(all( sapply(summary, function(s){ "SpliceLoci_Known_FewReads" %in% s$FIELD })) ){
     summary <- lapply(summary,function(dl){
        known.few <- dl$COUNT[dl$FIELD == "SpliceLoci_Known_FewReads"]
        known.many <- dl$COUNT[dl$FIELD == "SpliceLoci_Known_ManyReads"]
        novel <- dl$COUNT[dl$FIELD == "SpliceLoci_Novel"]

        known.any <- known.few + known.many;
        spliceLoci.any <- known.any + novel;

        dl[dim(dl)[1] + 1,] <- c(FIELD = "SpliceLoci_Known_Covered", COUNT = known.any);
        dl[dim(dl)[1] + 1,] <- c(FIELD = "SpliceLoci_Novel_Covered", COUNT = novel );
        dl[dim(dl)[1] + 1,] <- c(FIELD = "SpliceLoci_Covered", COUNT = spliceLoci.any);

        dl;
     })
     names(summary) <- summary.names;
     #message("debug 0B");
   }
   
 #message("names(res@qc.data) = ",paste0(names(res@qc.data),collapse=","));
 #message("debug 1");
   if(! res@singleEnd){
     if(! is.null(res@qc.data[["insert.size"]])){
       summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = res@qc.data[["insert.size"]], x.name = "InsertSize", y.name = "Ct", data.title = "InsertSize")
     }
     if(! is.null(res@qc.data[["gc.byPair"]])){
       summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = res@qc.data[["gc.byPair"]], x.name = "NUM_BASES_GC", y.name = "CT", data.title = "GC_byPair")
     }
     #message("debug 1B");
   }
 #message("debug 2");
 
   if(! is.null(res@qc.data[["gc.byRead"]])){
     summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = res@qc.data[["gc.byRead"]], x.name = "NUM_BASES_GC", y.name = "CT", data.title = "GC_byRead");
     #message("debug 2B");
   }
 #message("debug 3");
   if(! is.null(res@qc.data[["geneBodyCoverage.by.expression.level"]])){
     raw.data.list <- res@qc.data[["geneBodyCoverage.by.expression.level"]];
     step.size <- raw.data.list[[1]]$QUANTILE[2] - raw.data.list[[1]]$QUANTILE[1];
     genebody.data.list <- lapply(raw.data.list,function(df){ data.frame(Quantile = df$QUANTILE - (step.size / 2), GeneBodyCoverage = apply(df[,c("X1.bottomHalf","X2.upperMidQuartile","X3.75to90","X4.high")],1,sum)) });
     genebodyUM.data.list <- lapply(raw.data.list,function(df){ data.frame(Quantile = df$QUANTILE - (step.size / 2), GeneBodyCoverage = df[,"X2.upperMidQuartile"]) });
     genebodyLOW.data.list <- lapply(raw.data.list,function(df){ data.frame(Quantile = df$QUANTILE - (step.size / 2), GeneBodyCoverage = df[,"X1.bottomHalf"])});
   # message("debug 5");   
     summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = genebody.data.list, x.name = "Quantile", y.name = "GeneBodyCoverage", data.title = "GeneBodyCoverage_Overall");      
   # message("debug 6");
     summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = genebodyLOW.data.list, x.name = "Quantile", y.name = "GeneBodyCoverage", data.title = "GeneBodyCoverage_LowExpress");
   # message("debug 7");
     summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = genebodyUM.data.list, x.name = "Quantile", y.name = "GeneBodyCoverage", data.title = "GeneBodyCoverage_UMQuartile");
   # message("debug 8");
    # message("debug 3B");
   }
 
 #message("debug 5");

   if(! is.null(res@calc.data[["LANEBAM_GENE_CDF"]])){
    # message("debug 5B")
     summary <- lapply(names(res@calc.data[["LANEBAM_GENE_CDF"]]),function(n){
        genecdf <- res@calc.data[["LANEBAM_GENE_CDF"]][[n]];
        dl <- summary[[n]];

        dl[dim(dl)[1] + 1,] <- c(FIELD = "geneDiversityProfile_top1count", COUNT = genecdf[1]);
        dl[dim(dl)[1] + 1,] <- c(FIELD = "geneDiversityProfile_top10count", COUNT = genecdf[10]);
        dl[dim(dl)[1] + 1,] <- c(FIELD = "geneDiversityProfile_top100count", COUNT = genecdf[100]);
        dl[dim(dl)[1] + 1,] <- c(FIELD = "geneDiversityProfile_top1000count", COUNT = genecdf[1000]);
        readct <- genecdf[length(genecdf)];
        dl[dim(dl)[1] + 1,] <- c(FIELD = "geneDiversityProfile_top1pct", COUNT = genecdf[1] / readct);
        dl[dim(dl)[1] + 1,] <- c(FIELD = "geneDiversityProfile_top10pct", COUNT = genecdf[10] / readct);
        dl[dim(dl)[1] + 1,] <- c(FIELD = "geneDiversityProfile_top100pct", COUNT = genecdf[100] / readct);
        dl[dim(dl)[1] + 1,] <- c(FIELD = "geneDiversityProfile_top1000pct", COUNT = genecdf[1000] / readct);

        dl;
     })
     names(summary) <- summary.names;

   }

#message("debug 6");

   #message("1");
   if("input.read.pair.count" %in% names(res@decoder)){
     #message("2");
     #message("debug 6B")
     summary <- lapply(1:length(summary.names),function(i){
       #message("3");
       unique.ID <- summary.names[i];
       #message("4: unique.ID = ",unique.ID);
       dl <- summary[[i]];
       #message("5");
       irpc <- res@decoder$input.read.pair.count[res@decoder$unique.ID == unique.ID];
       totalReadPairs = dl$COUNT[dl$FIELD == "TOTAL_READ_PAIRS"];
       readPairsOk = dl$COUNT[dl$FIELD == "READ_PAIR_OK"];
       #message("6");
       dl[dim(dl)[1] + 1,] <- c(FIELD = "input.read.pair.count", COUNT = irpc);
       #message("6.1: irpc = ",irpc);
       #message("6.1: totalReadPairs = ",totalReadPairs);
       #message("6.1: readPairsOk = ",readPairsOk);
       dl[dim(dl)[1] + 1,] <- c(FIELD = "total.aligned.rate", COUNT = as.numeric(totalReadPairs) / as.numeric(irpc));
       #message("6.2");
       dl[dim(dl)[1] + 1,] <- c(FIELD = "aligned.and.pf.rate", COUNT = as.numeric(readPairsOk) / as.numeric(irpc));
       #message("7");
       ##READ_PAIR_OK
       ##TOTAL_READ_PAIRS
       
       if("multi.mapped.read.pair.count" %in% names(res@decoder)){
         mmrpc <- res@decoder$multi.mapped.read.pair.count[res@decoder$unique.ID == unique.ID];
         dl[dim(dl)[1] + 1,] <- c(FIELD = "multi.mapped.read.pair.count", COUNT = mmrpc);
         dl[dim(dl)[1] + 1,] <- c(FIELD = "multi.mapped.rate", COUNT = as.numeric(mmrpc) / as.numeric(irpc));

       }
       #message("8");
       dl;
     })
     
   }
#message("debug 7");
   #message("9");
   names(summary) <- summary.names;
   
   res@qc.data[["summary"]] <- summary;
   #message("10");
   
   res <- fix.summary.to.numeric(res);
   
   return(res);
}

INTERNAL.calcAndAdd.averages <- function(summary, summary.names, data.list, x.name, y.name, data.title){ 
   summary <- lapply(1:length(names(data.list)),function(i){
 #message(">   ",data.title, " - ",n," - ",length(data.list));
 #message(">      x.name contained? (", x.name %in% names(data.list[[1]]), ")");
 #message(">      y.name contained? (", y.name %in% names(data.list[[1]]), ")");
      dl <- summary[[i]];
      currdl <- data.list[[i]];
      dl[dim(dl)[1] + 1,] <- c(FIELD = paste0(data.title,"_Mean"), COUNT = INTERNAL.calc.mean(currdl,x.name,y.name));
      dl[dim(dl)[1] + 1,] <- c(FIELD = paste0(data.title,"_Median"), COUNT = INTERNAL.calc.median(currdl,x.name,y.name));
      dl;
   })
   names(summary) <- summary.names;
   summary;
}

INTERNAL.calc.median <- function(dl, x.name, y.name){
          curr.x <- dl[[x.name]];
          curr.y <- dl[[y.name]];
          curr.y <- curr.y / sum(curr.y);
          curr.cs <- cumsum(curr.y);
          index.median <- which(curr.cs > 0.5)[1];
          return(curr.x[index.median]);
}

INTERNAL.calc.mean <- function(dl, x.name, y.name){
   sum(
         as.numeric(dl[[x.name]]) *
         as.numeric(dl[[y.name]]) /
         sum(as.numeric(dl[[y.name]]))
      );
}

calc.NVC.rates <- function(res, debugMode = TRUE){
   runTimedFunction(debugMode=debugMode,title="Calculating NVC rates",expr={
     bases <- c("A","T","G","C");
     
     allBases <- c("A","T","G","C","N");
     
     tf <- function(df){
        dt <- sapply(bases, function(b){
          df$CT[ df$base == b ];
        });
        out <- data.frame(t(apply(dt, 1, function(r){ r / sum(r) })));
        return(out);
     }
     #if(!is.null(res@qc.data[["NVC.raw.R1"]]))
     if(!is.null(res@qc.data[["NVC.raw.R1"]])) res@calc.data[["NVC.raw.R1"]] <- lapply(res@qc.data[["NVC.raw.R1"]],tf);
     if(!is.null(res@qc.data[["NVC.raw.R2"]])) res@calc.data[["NVC.raw.R2"]] <- lapply(res@qc.data[["NVC.raw.R2"]],tf);

     if(!is.null(res@qc.data[["FQ.NVC.R1"]])) res@calc.data[["FQ.NVC.R1"]] <- lapply(res@qc.data[["FQ.NVC.R1"]],tf);
     if(!is.null(res@qc.data[["FQ.NVC.R2"]])) res@calc.data[["FQ.NVC.R2"]] <- lapply(res@qc.data[["FQ.NVC.R2"]],tf);
     if(!is.null(res@qc.data[["NVC.minus.clipping.R1"]])) res@calc.data[["NVC.minus.clipping.R1"]] <- lapply(res@qc.data[["NVC.minus.clipping.R1"]],tf);
     if(!is.null(res@qc.data[["NVC.minus.clipping.R2"]])) res@calc.data[["NVC.minus.clipping.R2"]] <- lapply(res@qc.data[["NVC.minus.clipping.R2"]],tf);
     
     tf <- function(df){
        dt <- sapply(bases, function(b){
          df$CT[ df$base == b ];
        });
        out <- data.frame(t(apply(dt, 1, function(r){ r / sum(r) })));
        out$readPos     <- df$readPos[ df$base == "A" ];
        out$leadClipLen <- df$leadClipLen[ df$base == "A" ];
        return(out);
     }
     if(!is.null(res@qc.data[["NVC.lead.clip.R1"]])) res@calc.data[["NVC.lead.clip.R1"]] <- lapply(res@qc.data[["NVC.lead.clip.R1"]],tf);
     if(!is.null(res@qc.data[["NVC.lead.clip.R2"]])) res@calc.data[["NVC.lead.clip.R2"]] <- lapply(res@qc.data[["NVC.lead.clip.R2"]],tf);
     
     tf <- function(df){
        dt <- sapply(bases, function(b){
          df$CT[ df$base == b ];
        });
        out <- data.frame(t(apply(dt, 1, function(r){ r / sum(r) })));
        out$readPos     <- df$readPos[ df$base == "A" ];
        out$tailClipLen <- df$tailClipLen[ df$base == "A" ];
        return(out);
     }
     if(!is.null(res@qc.data[["NVC.tail.clip.R1"]])) res@calc.data[["NVC.tail.clip.R1"]] <- lapply(res@qc.data[["NVC.tail.clip.R1"]],tf);
     if(!is.null(res@qc.data[["NVC.tail.clip.R2"]])) res@calc.data[["NVC.tail.clip.R2"]] <- lapply(res@qc.data[["NVC.tail.clip.R2"]],tf);
   });
   
   return(res);
}

calc.samplewise.norm.factors <- function(res, calc.DESeq2 , calc.edgeR ){
    samples <- res@sample.list;
   
    if(calc.DESeq2){
       if(! requireNamespace("DESeq2",quietly=TRUE)){
         message("Warning: DESeq2 installation not found. Skipping calculation of DESeq2 normalization factors.");
         calc.DESeq2 <- FALSE;
       }
       
    }
    if(calc.edgeR){
       if(! requireNamespace("edgeR",quietly=TRUE)){
         message("Warning: DESeq2 installation not found. Skipping calculation of DESeq2 normalization factors.");
         calc.edgeR <- FALSE;
       }
       
    }
    
    if(calc.edgeR | calc.DESeq2){
      mapped.reads <- lapply( res@qc.data[["summary"]], function(df){
         df$COUNT[df$FIELD == "READ_PAIR_OK"];
      });

      mapped.reads.by.sample <- sapply( samples, function(s){
        bamfiles.for.sample <- res@decoder$unique.ID[res@decoder$sample.ID == s];
        sum(unlist( mapped.reads[bamfiles.for.sample] ));
      })
      norm.factors <- data.frame(sample.ID = samples, M = mapped.reads.by.sample / 1000000, stringsAsFactors = FALSE);

      norm.factors$Norm_TC <- norm.factors$M / mean(norm.factors$M);

      if(calc.DESeq2 | calc.edgeR){
        count.matrix <- do.call(cbind.data.frame, res@calc.data[["SAMPLE_GENE_COUNTS"]]);
      }

      if(calc.DESeq2 & requireNamespace("DESeq2",quietly=TRUE)){
        message("Calculating DESeq2 Normalization Factors (Geometric normalization)...");
        #suppressPackageStartupMessages(requireNamespace("DESeq2"));

        tryCatch({
          norm.factors$Norm_Geo <- DESeq2::estimateSizeFactorsForMatrix(count.matrix);
        }, error = function(e){
          message("WARNING: DESeq2::estimateSizeFactorsForMatrix failed. Skipping DESeq2 normalization.",e);
          #norm.factors[,Norm_Geo:=NULL];
        });
      }

      if(calc.edgeR &  requireNamespace("edgeR",quietly=TRUE)){
        message("Calculating edgeR Normalization Factors (all edgeR normalizations)...");
        #suppressPackageStartupMessages(require("edgeR"));

        tryCatch({
          norm.factors$Norm_TMM <- edgeR::calcNormFactors(count.matrix, method="TMM");
        }, error = function(e){
          message("WARNING: edgeR::calcNormFactors(method=TMM) failed. Skipping edgeR TMM normalizations.",e);
          #norm.factors[,Norm_TMM:=NULL];
        });
         tryCatch({
           norm.factors$Norm_UQ <- edgeR::calcNormFactors(count.matrix, method="upperquartile");
         }, error = function(e){
           message("WARNING: edgeR::calcNormFactors(method=upperquartile) failed. Skipping edgeR upperquartile normalizations.",e);
           #norm.factors[,Norm_UQ:=NULL];
        });
        tryCatch({
          norm.factors$Norm_RLE <- edgeR::calcNormFactors(count.matrix, method="RLE");
        }, error = function(e){
          message("WARNING: edgeR::calcNormFactors(method=RLE) failed. Skipping edgeR RLE normalizations.",e);
          #norm.factors[,Norm_RLE:=NULL];
        });
      }

      res@calc.data[["norm.factors.bySample"]] <- norm.factors;
      return(res);
    } else {
      return(res);
    }
}

calc.lanebamwise.norm.factors <- function(res, calc.DESeq2 , calc.edgeR){
    if(calc.DESeq2){
       if(! requireNamespace("DESeq2",quietly=TRUE)){
         message("Warning: DESeq2 installation not found. Skipping calculation of DESeq2 normalization factors.");
         calc.DESeq2 <- FALSE;
       }
       
    }
    if(calc.edgeR){
       if(! requireNamespace("edgeR",quietly=TRUE)){
         message("Warning: edgeR installation not found. Skipping calculation of DESeq2 normalization factors.");
         calc.edgeR <- FALSE;
       }
       
    }
    
    if(calc.edgeR | calc.DESeq2){
      mapped.reads <- sapply( res@qc.data[["summary"]], function(df){
         as.numeric(df$COUNT[df$FIELD == "READ_PAIR_OK"]);
      });
      norm.factors <- data.frame(unique.ID = names(mapped.reads), M = unlist(mapped.reads) / 1000000, stringsAsFactors = FALSE);

      norm.factors$Norm_TC <- norm.factors$M / mean(norm.factors$M);



      if(calc.DESeq2 | calc.edgeR){
         df <- res@qc.data[["geneCounts"]][[1]];
         is.gene <- substr(df$GENEID,1,1) != "_";

         count.matrix <- do.call(cbind.data.frame, lapply(res@qc.data[["geneCounts"]], function(X){ X$COUNT } ));
         count.matrix <- count.matrix[is.gene,, drop = F];
      }

      if(calc.DESeq2 & requireNamespace("DESeq2",quietly=TRUE)){
        message("Calculating DESeq2 Normalization Factors (Geometric normalization)...");
        requireNamespace("DESeq2",quietly=TRUE);

        norm.factors$Norm_Geo <- DESeq2::estimateSizeFactorsForMatrix(count.matrix);
      }

      if(calc.edgeR & requireNamespace("edgeR",quietly=TRUE)){
        message("Calculating edgeR Normalization Factors (all edgeR normalizations)...");
        requireNamespace("edgeR",quietly=TRUE);
        norm.factors$Norm_TMM <- edgeR::calcNormFactors(count.matrix, method="TMM");
        norm.factors$Norm_UQ <- edgeR::calcNormFactors(count.matrix, method="upperquartile");
        norm.factors$Norm_RLE <- edgeR::calcNormFactors(count.matrix, method="RLE");
      }
      res@calc.data[["norm.factors.byLaneBam"]] <- norm.factors;
      return(res);
    } else {
      return(res);
    }
}


calc.lanebamwise.bysample.norm.factors <- function(res){


    if(is.null(res@calc.data[["norm.factors.bySample"]])){
      return(res);
    } else {
      samples <- res@sample.list;
      mapped.reads <- lapply( res@qc.data[["summary"]], function(df){
         df$COUNT[df$FIELD == "READ_PAIR_OK"];
      });
      norm.factors.bySample <- res@calc.data[["norm.factors.bySample"]];
      #norm.factors <- data.frame(unique.ID = res@decoder$unique.ID, sample.ID = res@decoder$sample.ID, lanebam.norm = rep(1,length(res@decoder$sample.ID)));

      bySample.list <- lapply(1:length(samples), function(i){
        s <- samples[i];
        in.df <- norm.factors.bySample[norm.factors.bySample$sample.ID == s,];
        bamfiles.for.sample <- res@decoder$unique.ID[res@decoder$sample.ID == s];
        curr.mapped.reads <- unlist(mapped.reads[bamfiles.for.sample]);
        lanebam.norm.factor <- curr.mapped.reads / mean(curr.mapped.reads);
        all.sample.norm.factors <- data.frame(t(sapply(1:length(bamfiles.for.sample),function(j){
          as.numeric(c(in.df[-1] * lanebam.norm.factor[j]));
        })));
        all.sample.norm.factors <- cbind(unique.ID = bamfiles.for.sample,all.sample.norm.factors);
        names(all.sample.norm.factors) <- c("unique.ID",names(in.df[-1]));
        all.sample.norm.factors;
      });

      norm.factors.byLaneBam <- do.call(rbind.data.frame, bySample.list);

      norm.factors.byLaneBam.reordering <- sapply(res@decoder$unique.ID, function(uid){
         which(norm.factors.byLaneBam$unique.ID == uid);
      });
      norm.factors.byLaneBam <- norm.factors.byLaneBam[norm.factors.byLaneBam.reordering,];

      res@calc.data[["norm.factors.bySample.splitToLaneBam"]] <- norm.factors.byLaneBam;
      return(res);
    }
}


calc.mapping.rates <- function(res){
  if(! is.null(res@decoder$input.read.pair.count)){
    total.reads <- res@decoder$input.read.pair.count;
  } else {
    total.reads <- sapply( res@qc.data[["summary"]], function(df){
       df$COUNT[df$FIELD == "PREALIGNMENT_READ_CT"];
    });
  }
  if(all(total.reads != -1)){
    mapped.reads <- sapply( res@qc.data[["summary"]], function(df){
       df$COUNT[df$FIELD == "READ_PAIR_OK"];
    });
    mapping.rate <- mapped.reads / total.reads;
    if(! is.null(res@decoder$multi.mapped.read.pair.count)){
      mm.reads <- res@decoder$multi.mapped.read.pair.count;
      mm.rate <- mm.reads / total.reads;
    } else {
      mm.reads <- sapply( res@qc.data[["summary"]], function(df){
        df$COUNT[df$FIELD == "KEPT_NOT_UNIQUE_ALIGNMENT"] + df$COUNT[df$FIELD == "DROPPED_NOT_UNIQUE_ALIGNMENT"];
      });
      mm.rate <- mm.reads / total.reads;
      
      if(all(mm.reads == 0)){
        message("NOTE: no multi-mapping found, and no multi.mapped.read.pair.count column found in decoder. All multi-mapping stats will be skipped.");
        mm.reads <- NULL;
      }
    }
    
    out.list <- lapply(1:length(res@decoder$unique.ID),function(i){
       if(is.null(mm.reads)){
         data.frame(FIELD = c("total.reads","mapped.reads","mapping.rate"), COUNT = c(total.reads[i], mapped.reads[i], mapping.rate[i]), stringsAsFactors = FALSE);
       } else {
         data.frame(FIELD = c("total.reads","mapped.reads","mapping.rate","mm.reads","mm.rate"), COUNT = c(total.reads[i], mapped.reads[i], mapping.rate[i],mm.reads[i], mm.rate[i]), stringsAsFactors = FALSE);
       }
    });
    names(out.list) <-  res@decoder$unique.ID;

    res@calc.data[["map.rates"]] <- out.list;
  }
  return(res);
}

calc.quals <- function(res){
   #if(! is.null(res@qc.data[["quals.r1"]])){
   #   res@decoder$cycle.CT <- get.cycleCt(res@qc.data[["quals.r1"]], res@decoder$unique.ID);
   #}
   res@decoder$cycle.CT <- sapply(res@qc.data[["summary"]], function(dl){
       dl$COUNT[dl$FIELD == "READ_LENGTH"];
   })
   return(res);
}

get.cycleCt <- function(data.list, unique.IDs){
   return(sapply(unique.IDs, FUN=function(lid){
      get.lanebam.cycleCt(data.list[[lid]]);
   }));
}
get.lanebam.cycleCt <- function(df){
   return(dim(df)[1]);
}

calc.gene.CDF <- function(res){

  temp <- res@qc.data[["geneCounts"]];
  if(! is.null(temp)){
    out <- lapply(temp, function(df){
      df <- df[substr(df$GENEID,1,1) != "_", c("GENEID","COUNT")];
      ord <- order(df$COUNT, decreasing=TRUE);
      sum <- sum(df$COUNT);
      cumsum(df$COUNT[ord]);
    });
    res@calc.data[["LANEBAM_GENE_CDF"]] <- out;
  }
  return(res);
}

calc.gene.CDF.TEMPVERSION <- function(res){
  temp <- res@qc.data[["geneCounts.deseq"]];
  out <- lapply(temp, function(df){
    df <- df[1:(dim(df)[1] - 5), c(1,2)];
    names(df) <- c("GENEID","COUNT");
    
    ord <- order(df$COUNT, decreasing=TRUE);
    sum <- sum(df$COUNT);
    cumsum(df$COUNT[ord]);
  });
  res@calc.data[["LANEBAM_GENE_CDF"]] <- out;
  return(res);
}

calc.gene.CDF.bySample <- function(res){
  #res@qc.data[["geneCounts"]];
  gene.ct.bySample <- lapply(res@sample.list, function(sample){
    lanebam.set <- res@decoder$unique.ID[res@decoder$sample.ID == sample];
    if(length(lanebam.set) == 1){
       df <-  res@qc.data[["geneCounts"]][[lanebam.set[1]]];
       is.gene <- substr(df$GENEID,1,1) != "_";
       ct.so.far <- df$COUNT[ is.gene ];
       return(ct.so.far);
       #return(res@qc.data[["geneCounts"]][[lanebam.set[1]]]);
    } else {
       df <- res@qc.data[["geneCounts"]][[lanebam.set[1]]];
       is.gene <- substr(df$GENEID,1,1) != "_";
       ct.so.far <- df$COUNT[ is.gene ];

       for(i in 2:length(lanebam.set)){
         df <- res@qc.data[["geneCounts"]][[lanebam.set[i]]];
         ct.so.far <- ct.so.far + df$COUNT[is.gene];
       }
       ct.so.far;
    }
  });
  names(gene.ct.bySample) <- res@sample.list;
  
  out <- lapply(res@sample.list, function(sample){
    lanebam.set <- res@decoder$unique.ID[res@decoder$sample.ID == sample];
    if(length(lanebam.set) == 1){
       return(res@calc.data[["LANEBAM_GENE_CDF"]][[lanebam.set[1]]]);
    } else {
       ct.so.far <- gene.ct.bySample[[sample]];
       ord <- order(ct.so.far, decreasing = TRUE);
       cumsum(ct.so.far[ord]);
    }
  });
  names(out) <- res@sample.list;
  res@calc.data[["SAMPLE_GENE_CDF"]] <- out;
  res@calc.data[["SAMPLE_GENE_COUNTS"]] <- gene.ct.bySample;
  return(res);
}

#DEPRECIATED:
calc.raw.NVC <- function(res){
  res@calc.data[["NVC.raw.r1"]] <- lapply(res@qc.data[["NVC.raw"]],function(df){
    return(df[df$readNum == 1,]);
  });
  res@calc.data[["NVC.raw.r2"]] <- lapply(res@qc.data[["NVC.raw"]],function(df){
    return(df[df$readNum == 2,]);
  });
  return(res);
}

#DEPRECIATED:
calc.gene.CDF.bySample.OLDVERSION <- function(res){
  #res@qc.data[["geneCounts"]];
  out <- lapply(res@sample.list, function(sample){
    lanebam.set <- res@decoder$unique.ID[res@decoder$sample.ID == sample];
    if(length(lanebam.set) == 1){
       return(res@calc.data[["LANEBAM_GENE_CDF"]][[lanebam.set[1]]]);
    } else {
       df <- res@qc.data[["geneCounts"]][[lanebam.set[1]]];
       is.gene <- substr(df$GENEID,1,1) != "_";
       ct.so.far <- df$COUNT[ is.gene ];

       for(i in 2:length(lanebam.set)){
         df <- res@qc.data[["geneCounts"]][[lanebam.set[i]]];
         ct.so.far <- ct.so.far + df$COUNT[is.gene];
       }
       ord <- order(ct.so.far, decreasing = TRUE);
       #sum <- sum(ct.so.far);
       cumsum(ct.so.far[ord]);
    }
  });

  names(out) <- res@sample.list;

  res@calc.data[["SAMPLE_GENE_CDF"]] <- out;
  return(res);
}

#DEPRECIATED:
read.in.results.data.with.single.decoder.file <- function(infile.dir,decoder.file){
  decoder <- read.table(decoder.file,header=TRUE,stringsAsFactors=F);
  if(is.null(decoder$unique.ID)) stop("Decoder formatting error: no column labelled unique.ID");
  if(is.null(decoder$lane.ID)) stop("Decoder formatting error: no column labelled lane.ID");
  if(is.null(decoder$group.ID)) stop("Decoder formatting error: no column labelled group.ID");
  if(is.null(decoder$sample.ID)) stop("Decoder formatting error: no column labelled sample.ID");
  if(length(unique(decoder$unique.ID)) < length(decoder$unique.ID)) stop("Decoder formatting error: unique.ID must be unique!");
  #more checks?
  if( length(unique(names(decoder))) != length(names(decoder)) ) stop("Decoder formatting error: column names must be unique!");
  read.in.results.data.with.decoder(decoder, infile.dir);
}

#DEPRECIATED:
read.in.results.data.with.two.decoder.files <- function(infile.dir, lanebam.decoder, sample.decoder){
  sample.decoder <- read.table(sample.decoder,header=TRUE,stringsAsFactors=F);
  lanebam.decoder <- read.table(lanebam.decoder,header=TRUE,stringsAsFactors=F);
  if(is.null(lanebam.decoder$unique.ID)) stop("Lanebam Decoder formatting error: no column labelled unique.ID");
  if(is.null(lanebam.decoder$lane.ID)) stop("Lanebam Decoder formatting error: no column labelled lane.ID");
  if(is.null(lanebam.decoder$sample.ID)) stop("Lanebam Decoder formatting error: no column labelled sample.ID");
  if(length(unique(lanebam.decoder$unique.ID)) < length(lanebam.decoder$unique.ID)) stop("Lanebam Decoder formatting error: unique.ID must be unique!");
  if( length(unique(names(lanebam.decoder))) != length(names(lanebam.decoder)) ) stop("Lanebam Decoder formatting error: column names must be unique!");

  if(is.null(sample.decoder$sample.ID)) stop("Sample Decoder formatting error: no column labelled sample.ID");
  if(is.null(sample.decoder$group.ID)) stop("Sample Decoder formatting error: no column labelled group.ID");
  if(length(unique(sample.decoder$sample.ID)) < length(sample.decoder$sample.ID)) stop("Sample Decoder formatting error: sample.ID must be unique!");

  if(! all(unique(sample.decoder$sample.ID) %in% unique(lanebam.decoder$sample.ID) )) {
    stop("Sample/Lanebam Decoders formatting error: no 1-to-1 matching between the two decoders sample.IDs!\n     Missing from lanebam.decoder: ", paste0(unique(sample.decoder$sample.ID)[ ! unique(sample.decoder$sample.ID) %in% unique(lanebam.decoder$sample.ID) ], collapse=",")  )
  }
   
  if(! all(unique(lanebam.decoder$sample.ID) %in% unique(sample.decoder$sample.ID) )) {
    stop("Sample/Lanebam Decoders formatting error: no 1-to-1 matching between the two decoders sample.IDs!\n     Missing from sample.decoder: ", paste0(unique(lanebam.decoder$sample.ID)[ ! unique(lanebam.decoder$sample.ID) %in% unique(sample.decoder$sample.ID) ],collapse=",")  )
  }

  if( length(unique(names(sample.decoder))) != length(names(sample.decoder)) ) stop("Sample Decoder formatting error: column names must be unique!");

  merged.decoder <- lanebam.decoder;
  
  merged.decoder$group.ID <- sapply(lanebam.decoder$sample.ID,FUN=function(samp){
    sample.decoder$group.ID[sample.decoder$sample.ID == samp];
  });

  for(pheno.var in names(sample.decoder)){
    if(! pheno.var %in% c("sample.ID","group.ID")){
      merged.decoder[[pheno.var]] <- sapply(lanebam.decoder$sample.ID, FUN=function(samp){
        sample.decoder[[pheno.var]][sample.decoder$sample.ID == samp];
      });
    }
  }
  read.in.results.data.with.decoder( merged.decoder, infile.dir );
}
