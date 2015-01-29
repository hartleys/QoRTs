
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
#           multi.mapped.read.pair.count: the # of reads that were multi mapped by the aligner. this must be included for multi-mapping rate to be calculated.
#        RESERVED FIELDS: Do not name any field this:
#           cycle.CT
#           lanebam.ID (a synonym for unique.ID)
#    Sample Decoder:
#        The Sample Decoder must have the following column names:
#        sample.ID	group.ID
#            other than these 2 required columns, it can have any other columns desired, as long as the column names are unique.
#

read.qc.results.data <- function(infile.dir, decoder = NULL, decoder.files = NULL, calc.DESeq2 = FALSE, calc.edgeR = FALSE, debugMode = DEFAULTDEBUGMODE ) {
   decoder.final <- completeAndCheckDecoder(decoder = decoder, decoder.files = decoder.files);
   
   return(read.in.results.data.with.decoder(decoder = decoder.final, infile.dir = infile.dir, calc.DESeq2 = calc.DESeq2 , calc.edgeR = calc.edgeR, debugMode = debugMode) );
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
    decoder <- data.frame(unique.ID = decoder);
    message("Simple decoder found, list of sample/unique ID's. Building complete decoder.");
  }
  if(class(decoder) == "matrix"){
    decoder <- data.frame(decoder);
  }
  if(class(decoder) != "data.frame"){
    message("Decoder must be either a matrix, a character vector (of unique ID's), or a data.frame. Decoder follows:");
    print(decoder);
    stop("ERROR: Decoder must be either a matrix, a character vector (of unique ID's), or a data.frame.");
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

read.in.results.data.with.decoder <- function(decoder, infile.dir = "", calc.DESeq2 = FALSE, calc.edgeR = FALSE , debugMode = DEFAULTDEBUGMODE){
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
    message("Note: no multi.mapped.read.pair.count column found. This column is optional, but without it multi-mapping rates cannot be plotted.");
  }# else {
  #  message("Note: successfully found multi.mapped.read.pair.count: column. This will be used to plot multi-mapping rates.");
  #}
  
  for(i in 1:length(decoder$unique.ID)){
    if(! file.exists(paste0(infile.dir,decoder$qc.data.dir[i]))){
       stop("Directory not found: ",paste0(infile.dir,decoder$qc.data.dir[i]));
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

  read.scalaqc.file.helper <- function(scalaqc_file){
    if(debugMode) ts <- timestamp();
    out <- read.in.scalaQC.files(infile.dir,lanebam.list, qc.data.dir.list,scalaqc_file);
    if(debugMode) reportTimeAndDiff(ts);
    out;
  }
  
  res@qc.data <- list(summary = read.scalaqc.file.helper(QC_INTERNAL_SCALAQC_SUMMARY_FILE));

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
    res@qc.data <- c(res@qc.data,lapply(QC_INTERNAL_SCALAQC_FILE_LIST_SINGLE_END, FUN=read.scalaqc.file.helper));
  } else {
    if(debugMode) message("Autodetected Paired-End mode.");
    res@qc.data <- c(res@qc.data,lapply(QC_INTERNAL_SCALAQC_FILE_LIST, FUN=read.scalaqc.file.helper));
  }
  res@calc.data <- list();

  
  message("calculating secondary data:");
  if(debugMode) ts <- timestamp();
  res <- calc.results.data(res, calc.DESeq2 = calc.DESeq2, calc.edgeR = calc.edgeR);
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
                                       chrom.counts = "QC.chromCount.txt.gz"
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
                                       chrom.counts = "QC.chromCount.txt.gz"
                                       #,
                                       #spliceJunctionCounts.knownSplices = "scalaQC.spliceJunctionCounts.knownSplices.txt.gz"
                                       #spliceJunctionCounts.novelSplices = "scalaQC.spliceJunctionCounts.novelSplices.txt.gz"
                                       #
);


##################################################################################################################################
##################################################################################################################################

read.in.scalaQC.files <- function(infile.prefix, lanebam.list, qc.data.dir.list, infile.suffix){
  message(paste0("reading ",infile.suffix," files..."));
  infiles <- paste0(infile.prefix,unlist(qc.data.dir.list),"/", infile.suffix);
  if(file.exists(infiles[1])){
    #print("!")
    for(i in 1:length(infiles)){
      if(! file.exists(infiles[i])){
        stop("File not found: ",infiles[i]);
      }
    }
    out <- lapply(lanebam.list,FUN=function(unique.ID){
      i <- which(unlist(lanebam.list) == unique.ID);
      d <- tryCatch({
         read.table(infiles[i],header=T,stringsAsFactors=F);
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
    message("done.");
    return(out);
  } else {
    message(paste0("Failed: Cannot find file: ",infiles[1],". Skipping tests that use this data."));
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

calc.results.data <- function(res, calc.DESeq2, calc.edgeR ){
   res <- fix.summary.to.numeric(res);
   check.all.completed.without.error(res);
   res <- calc.quals(res);
   res <- calc.gene.CDF(res);
   res <- calc.gene.CDF.bySample(res);
   #res <- calc.raw.NVC(res);
   res <- calc.mapping.rates(res);
   res <- calc.samplewise.norm.factors(res, calc.DESeq2,calc.edgeR);
   res <- calc.lanebamwise.norm.factors(res, calc.DESeq2,calc.edgeR);
   res <- calc.lanebamwise.bysample.norm.factors(res);
   res <- add.to.summary(res);
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
   summary <- res@qc.data[["summary"]];
   summary.names <- names(summary);
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
   
 #message("names(res@qc.data) = ",paste0(names(res@qc.data),collapse=","));
 #message("debug 1");
   if(! res@singleEnd){
     summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = res@qc.data[["insert.size"]], x.name = "InsertSize", y.name = "Ct", data.title = "InsertSize")
     summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = res@qc.data[["gc.byPair"]], x.name = "NUM_BASES_GC", y.name = "CT", data.title = "GC_byPair")
   }
 #message("debug 2");
   summary <- INTERNAL.calcAndAdd.averages(summary = summary, summary.names = summary.names, data.list = res@qc.data[["gc.byRead"]], x.name = "NUM_BASES_GC", y.name = "CT", data.title = "GC_byRead")
 #message("debug 3");
 #message("debug 4");
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
   #message("1");
   if("input.read.pair.count" %in% names(res@decoder)){
     #message("2");
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
   #message("9");
   names(summary) <- summary.names;
   
   res@qc.data[["summary"]] <- summary;
   #message("10");
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

calc.NVC.rates <- function(res){
   bases <- c("A","T","G","C");
   tf <- function(df){
      dt <- sapply(bases, function(b){
        df$CT[ df$base == b ];
      });
      out <- data.frame(t(apply(dt, 1, function(r){ r / sum(r) })));
      return(out);
   }
   res@calc.data[["rates.NVC.raw.R1"]] <- lapply(res@qc.data[["NVC.raw.R1"]],tf);
   res@calc.data[["rates.NVC.raw.R2"]] <- lapply(res@qc.data[["NVC.raw.R2"]],tf);

   res@calc.data[["rates.NVC.lead.clip.R1"]] <- lapply(res@qc.data[["NVC.lead.clip.R1"]],tf);
   res@calc.data[["rates.NVC.lead.clip.R2"]] <- lapply(res@qc.data[["NVC.lead.clip.R2"]],tf);
   res@calc.data[["rates.NVC.tail.clip.R1"]] <- lapply(res@qc.data[["NVC.tail.clip.R1"]],tf);
   res@calc.data[["rates.NVC.tail.clip.R2"]] <- lapply(res@qc.data[["NVC.tail.clip.R2"]],tf);
   res@calc.data[["rates.NVC.minus.clipping.R1"]] <- lapply(res@qc.data[["NVC.minus.clipping.R1"]],tf);
   res@calc.data[["rates.NVC.minus.clipping.R2"]] <- lapply(res@qc.data[["NVC.minus.clipping.R2"]],tf);
   return(res);
}

calc.samplewise.norm.factors <- function(res, calc.DESeq2 , calc.edgeR ){
    samples <- res@sample.list;
   
    if(calc.DESeq2){
       if(! "DESeq2" %in% rownames(installed.packages())){
         message("Warning: DESeq2 installation not found. Skipping calculation of DESeq2 normalization factors.");
       }
    }
    if(calc.edgeR){
       if(! "edgeR" %in% rownames(installed.packages())){
         message("Warning: DESeq2 installation not found. Skipping calculation of DESeq2 normalization factors.");
       }
    }
   
    mapped.reads <- lapply( res@qc.data[["summary"]], function(df){
       df$COUNT[df$FIELD == "READ_PAIR_OK"];
    });
    
    mapped.reads.by.sample <- sapply( samples, function(s){
      bamfiles.for.sample <- res@decoder$unique.ID[res@decoder$sample.ID == s];
      sum(unlist( mapped.reads[bamfiles.for.sample] ));
    })
    norm.factors <- data.frame(sample.ID = samples, M = mapped.reads.by.sample / 1000000);
    
    norm.factors$Norm_TC <- norm.factors$M / mean(norm.factors$M);
    
    if(calc.DESeq2 | calc.edgeR){
      count.matrix <- do.call(cbind.data.frame, res@calc.data[["SAMPLE_GENE_COUNTS"]]);
    }
    
    if(calc.DESeq2 & "DESeq2" %in% rownames(installed.packages())){
      message("Calculating DESeq2 Normalization Factors (Geometric normalization)...");
      suppressPackageStartupMessages(require("DESeq2"));
      
      tryCatch({
        norm.factors$Norm_Geo <- estimateSizeFactorsForMatrix(count.matrix);
      }, warning = function(w){
        message("WARNING: DESeq2::estimateSizeFactorsForMatrix threw warnings: ");
        message(w);
      }, error = function(e){
        message("WARNING: DESeq2::estimateSizeFactorsForMatrix failed. Skipping DESeq2 normalization.",e);
        #norm.factors[,Norm_Geo:=NULL];
      });
    }
    
    if(calc.edgeR & "edgeR" %in% rownames(installed.packages())){
      message("Calculating edgeR Normalization Factors (all edgeR normalizations)...");
      suppressPackageStartupMessages(require("edgeR"));
      
      tryCatch({
        norm.factors$Norm_TMM <- calcNormFactors(count.matrix, method="TMM");
      }, error = function(e){
        message("WARNING: edgeR::calcNormFactors(method=TMM) failed. Skipping edgeR TMM normalizations.",e);
        #norm.factors[,Norm_TMM:=NULL];
      });
       tryCatch({
         norm.factors$Norm_UQ <- calcNormFactors(count.matrix, method="upperquartile");
       }, error = function(e){
         message("WARNING: edgeR::calcNormFactors(method=upperquartile) failed. Skipping edgeR upperquartile normalizations.",e);
         #norm.factors[,Norm_UQ:=NULL];
      });
      tryCatch({
        norm.factors$Norm_RLE <- calcNormFactors(count.matrix, method="RLE");
      }, error = function(e){
        message("WARNING: edgeR::calcNormFactors(method=RLE) failed. Skipping edgeR RLE normalizations.",e);
        #norm.factors[,Norm_RLE:=NULL];
      });
    }
    
    res@calc.data[["norm.factors.bySample"]] <- norm.factors;
    return(res);

}

calc.lanebamwise.norm.factors <- function(res, calc.DESeq2 , calc.edgeR){
    if(calc.DESeq2){
       if(! "DESeq2" %in% rownames(installed.packages())){
         message("Warning: DESeq2 installation not found. Skipping calculation of DESeq2 normalization factors.");
       }
    }
    if(calc.edgeR){
       if(! "edgeR" %in% rownames(installed.packages())){
         message("Warning: DESeq2 installation not found. Skipping calculation of DESeq2 normalization factors.");
       }
    }
    mapped.reads <- sapply( res@qc.data[["summary"]], function(df){
       as.numeric(df$COUNT[df$FIELD == "READ_PAIR_OK"]);
    });
    norm.factors <- data.frame(unique.ID = names(mapped.reads), M = unlist(mapped.reads) / 1000000);
    
    norm.factors$Norm_TC <- norm.factors$M / mean(norm.factors$M);



    if(calc.DESeq2 | calc.edgeR){
       df <- res@qc.data[["geneCounts"]][[1]];
       is.gene <- substr(df$GENEID,1,1) != "_";

       count.matrix <- do.call(cbind.data.frame, lapply(res@qc.data[["geneCounts"]], function(X){ X$COUNT } ));
       count.matrix <- count.matrix[is.gene,, drop = F];
    }
    
    if(calc.DESeq2 & "DESeq2" %in% rownames(installed.packages())){
      message("Calculating DESeq2 Normalization Factors (Geometric normalization)...");
      suppressPackageStartupMessages(require("DESeq2"));
      
      norm.factors$Norm_Geo <- estimateSizeFactorsForMatrix(count.matrix);
    }
    
    if(calc.edgeR & "edgeR" %in% rownames(installed.packages())){
      message("Calculating edgeR Normalization Factors (all edgeR normalizations)...");
      suppressPackageStartupMessages(require("edgeR"));
      norm.factors$Norm_TMM <- calcNormFactors(count.matrix, method="TMM");
      norm.factors$Norm_UQ <- calcNormFactors(count.matrix, method="upperquartile");
      norm.factors$Norm_RLE <- calcNormFactors(count.matrix, method="RLE");
    }
    res@calc.data[["norm.factors.byLaneBam"]] <- norm.factors;
    return(res);
}


calc.lanebamwise.bysample.norm.factors <- function(res){
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


calc.mapping.rates <- function(res){
  if(! is.null(res@decoder$input.read.pair.count)){
    total.reads <- res@decoder$input.read.pair.count;
    mapped.reads <- sapply( res@qc.data[["summary"]], function(df){
       df$COUNT[df$FIELD == "READ_PAIR_OK"];
    });
    mapping.rate <- mapped.reads / total.reads;
    if(! is.null(res@decoder$multi.mapped.read.pair.count)){
      mm.reads <- res@decoder$multi.mapped.read.pair.count;
      mm.rate <- mm.reads / total.reads;
    }
    out.list <- lapply(1:length(res@decoder$unique.ID),function(i){
       if(is.null(res@decoder$multi.mapped.read.pair.count)){
         data.frame(FIELD = c("total.reads","mapped.reads","mapping.rate"), COUNT = c(total.reads[i], mapped.reads[i], mapping.rate[i]));
       } else {
         data.frame(FIELD = c("total.reads","mapped.reads","mapping.rate","mm.reads","mm.rate"), COUNT = c(total.reads[i], mapped.reads[i], mapping.rate[i],mm.reads[i], mm.rate[i]));
       }
    });
    names(out.list) <-  res@decoder$unique.ID;

    res@calc.data[["map.rates"]] <- out.list;
  }
  return(res);
}

calc.quals <- function(res){
   if(! is.null(res@qc.data[["quals.r1"]])){
      res@decoder$cycle.CT <- get.cycleCt(res@qc.data[["quals.r1"]], res@decoder$unique.ID);
   }
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
  out <- lapply(temp, function(df){
    df <- df[substr(df$GENEID,1,1) != "_", c("GENEID","COUNT")];
    ord <- order(df$COUNT, decreasing=TRUE);
    sum <- sum(df$COUNT);
    cumsum(df$COUNT[ord]);
  });
  res@calc.data[["LANEBAM_GENE_CDF"]] <- out;
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