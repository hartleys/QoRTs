
require(methods);
DEFAULTDEBUGMODE <- TRUE;

##############################
#Res builders:

#Input File requirements:
#Input can come in two forms:
#(1) Single Decoder
#    This requires a single "decoder" file, which MUST have the following column names:
#    lanebam.ID	lane.ID	group.ID	sample.ID	qc.data.dir
#       Other than these 4 required columns, it can have as many other columns as desired.
#       However: all column names must be UNIQUE.
#       OPTIONAL FIELDS:
#          input.read.pair.count: the # of input reads. this must be included for mapping rate to be calculated.
#          multi.mapped.read.pair.count: the # of reads that were multi mapped by the aligner. this must be included for multi-mapping rate to be calculated.
#       RESERVED FIELDS: Do not name any field this:
#          cycle.CT

#(2) Dual Decoder:
#    This requires two "decoder" files.
#    Lanebam Decoder:
#        The Lanebam Decoder must have the following column names:
#        lanebam.ID	lane.ID	sample.ID
#           Other than these 3 required columns, it can have any other columns desired, as long as the column names are unique.
#        OPTIONAL FIELDS:
#           input.read.pair.count: the # of input reads. this must be included for mapping rate to be calculated.
#           multi.mapped.read.pair.count: the # of reads that were multi mapped by the aligner. this must be included for multi-mapping rate to be calculated.
#        RESERVED FIELDS: Do not name any field this:
#           cycle.CT
#    Sample Decoder:
#        The Sample Decoder must have the following column names:
#        sample.ID	group.ID
#            other than these 2 required columns, it can have any other columns desired, as long as the column names are unique.

#OPTIONAL FIELDS

read.qc.results.data <- function(infile.dir, decoder.files = NULL, decoder.data.frame = NULL, calc.DESeq2 = FALSE, calc.edgeR = FALSE, debugMode = DEFAULTDEBUGMODE ) {
   decoder <- NULL;
   if(! is.null(decoder.data.frame)){
     decoder <- decoder.data.frame;
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
     stop("Fatal Error: QoRTs.read.results.data(): either decoder.files or decoder.data.frame must be set!");
   }
   if(is.null(decoder)) stop("FATAL ERROR: Null Decoder?");
   
   return(read.in.results.data.with.decoder(decoder = decoder, infile.dir = infile.dir, calc.DESeq2 = calc.DESeq2 , calc.edgeR = calc.edgeR, debugMode = debugMode) );
}

##########################################################################################

get.decoder.from.single.file <- function(decoder.file){
  decoder <- read.table(decoder.file,header=TRUE,stringsAsFactors=F);
  if(is.null(decoder$lanebam.ID)) stop("Decoder formatting error: no column labelled lanebam.ID");
  if(is.null(decoder$lane.ID)) stop("Decoder formatting error: no column labelled lane.ID");
  if(is.null(decoder$group.ID)) stop("Decoder formatting error: no column labelled group.ID");
  if(is.null(decoder$sample.ID)) stop("Decoder formatting error: no column labelled sample.ID");
  
  if(is.null(decoder$qc.data.dir)){
     message("qc.data.dir not found, assuming qc.data.dir = lanebam.ID");
     decoder$qc.data.dir <- decoder$lanebam.ID;
  }  
  
  if(length(unique(decoder$lanebam.ID)) < length(decoder$lanebam.ID)) stop("Decoder formatting error: lanebam.ID must be unique!");
  #more checks?
  if( length(unique(names(decoder))) != length(names(decoder)) ) stop("Decoder formatting error: column names must be unique!");
  return(decoder);
}

get.decoder.from.dual.file <- function(lanebam.decoder, sample.decoder){
  sample.decoder <- read.table(sample.decoder,header=TRUE,stringsAsFactors=F);
  lanebam.decoder <- read.table(lanebam.decoder,header=TRUE,stringsAsFactors=F);
  if(is.null(lanebam.decoder$lanebam.ID)) stop("Lanebam Decoder formatting error: no column labelled lanebam.ID");
  if(is.null(lanebam.decoder$lane.ID)) stop("Lanebam Decoder formatting error: no column labelled lane.ID");
  if(is.null(lanebam.decoder$sample.ID)) stop("Lanebam Decoder formatting error: no column labelled sample.ID");
  if(is.null(lanebam.decoder$qc.data.dir)){
     message("qc.data.dir not found, assuming qc.data.dir = lanebam.ID");
     lanebam.decoder$qc.data.dir <- lanebam.decoder$lanebam.ID;
  } 
  if(length(unique(lanebam.decoder$lanebam.ID)) < length(lanebam.decoder$lanebam.ID)) stop("Lanebam Decoder formatting error: lanebam.ID must be unique!");
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
  return(merged.decoder);
}


###############################################################

read.in.results.data.with.decoder <- function(decoder, infile.dir = "", calc.DESeq2 = FALSE, calc.edgeR = FALSE , debugMode){

  if(is.null(decoder$lanebam.ID)) stop("Decoder formatting error: no column labelled lanebam.ID");
  if(is.null(decoder$lane.ID)) stop("Decoder formatting error: no column labelled lane.ID");
  if(is.null(decoder$group.ID)) stop("Decoder formatting error: no column labelled group.ID");
  if(is.null(decoder$sample.ID)) stop("Decoder formatting error: no column labelled sample.ID");
  if(is.null(decoder$qc.data.dir)) stop("Decoder formatting error: no column labelled qc.data.dir");
  if(length(unique(decoder$lanebam.ID)) < length(decoder$lanebam.ID)) stop("Decoder formatting error: lanebam.ID must be unique!");
  #more checks?
  if( length(unique(names(decoder))) != length(names(decoder)) ) stop("Decoder formatting error: column names must be unique!");

  if( is.null(decoder$input.read.pair.count)){
    message("Note: no input.read.pair.count column found. This column is optional, but without it mapping rates cannot be calculated.");
  } else {
    message("Note: successfully found input.read.pair.count column. This will be used to calculate mapping rates.");
  }

  
  for(i in 1:length(decoder$lanebam.ID)){
    if(! file.exists(paste0(infile.dir,decoder$qc.data.dir[i]))){
       stop("Directory not found: ",paste0(infile.dir,decoder$qc.data.dir[i]));
    }
  }

  res <- new("QoRT_QC_Results");
  res@lanebam.list <- decoder$lanebam.ID;
  res@sample.list <- unique(decoder$sample.ID);
  res@lane.list <- unique(decoder$lane.ID);
  res@group.list <- unique(decoder$group.ID);
  
  res@decoder <- cbind.data.frame(decoder,cycle.CT = rep(-1,length(decoder$lanebam.ID)));

  lanebam.list <- as.list(decoder$lanebam.ID);
  names(lanebam.list) <- decoder$lanebam.ID
  
  qc.data.dir.list <- as.list(decoder$qc.data.dir);
  names(qc.data.dir.list) <- decoder$lanebam.ID

  res@qc.data <- lapply(QC_INTERNAL_SCALAQC_FILE_LIST, FUN=function(scalaqc_file){
    if(debugMode) ts <- timestamp();
    out <- read.in.scalaQC.files(infile.dir,lanebam.list, qc.data.dir.list,scalaqc_file);
    if(debugMode) reportTimeAndDiff(ts);
    out;
  });
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

QC_INTERNAL_SCALAQC_FILE_LIST <- list( summary = "QCsummary.txt",
                                       gc = "QC.gc.txt.gz",
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
                                       #
                                       #
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
    out <- lapply(lanebam.list,FUN=function(lanebam.ID){
      i <- which(unlist(lanebam.list) == lanebam.ID);
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
#decoder="data.frame", #decoder has columns: lanebam.ID	sample.ID	lane.ID	group.ID	cycle.CT	and then any number of user-defined columns (which are ignored internally)
#qc.data="list" #List of Lists. Each element corresponds to one qc test, and is composed of a list, one element for each lanebam.

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

calc.results.data <- function(res, calc.DESeq2, calc.edgeR ){
   res <- fix.summary.to.numeric(res);
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
   res@qc.data[["summary"]] <- summary;

   return(res);
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
      bamfiles.for.sample <- res@decoder$lanebam.ID[res@decoder$sample.ID == s];
      sum(unlist( mapped.reads[bamfiles.for.sample] ));
    })
    norm.factors <- data.frame(sample.ID = samples, M = mapped.reads.by.sample / 1000000);
    
    norm.factors$Norm_TC <- norm.factors$M / mean(norm.factors$M);
    
    if(calc.DESeq2 | calc.edgeR){
      count.matrix <- do.call(cbind.data.frame, res@calc.data[["SAMPLE_GENE_COUNTS"]]);
    }
    
    if(calc.DESeq2 & "DESeq2" %in% rownames(installed.packages())){
      message("Calculating DESeq2 Normalization Factors (Geometric normalization)...\n");
      require("DESeq2");
      
      tryCatch({
        norm.factors$Norm_Geo <- estimateSizeFactorsForMatrix(count.matrix);
      }, warning = function(w){
        message("WARNING: DESeq2::estimateSizeFactorsForMatrix threw warnings: ");
        message(w);
      }, error = function(e){
        message("WARNING: DESeq2::estimateSizeFactorsForMatrix failed. Skipping DESeq2 normalization.");
        #norm.factors[,Norm_Geo:=NULL];
      });
    }
    
    if(calc.edgeR & "edgeR" %in% rownames(installed.packages())){
      message("Calculating edgeR Normalization Factors (all edgeR normalizations)...\n");
      require("edgeR");
      
      tryCatch({
        norm.factors$Norm_TMM <- calcNormFactors(count.matrix, method="TMM");
      }, error = function(e){
        message("WARNING: edgeR::calcNormFactors(method=TMM) failed. Skipping edgeR TMM normalizations.");
        #norm.factors[,Norm_TMM:=NULL];
      });
       tryCatch({
         norm.factors$Norm_UQ <- calcNormFactors(count.matrix, method="upperquartile");
       }, error = function(e){
         message("WARNING: edgeR::calcNormFactors(method=upperquartile) failed. Skipping edgeR upperquartile normalizations.");
         #norm.factors[,Norm_UQ:=NULL];
      });
      tryCatch({
        norm.factors$Norm_RLE <- calcNormFactors(count.matrix, method="RLE");
      }, error = function(e){
        message("WARNING: edgeR::calcNormFactors(method=RLE) failed. Skipping edgeR RLE normalizations.");
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
    norm.factors <- data.frame(lanebam.ID = names(mapped.reads), M = unlist(mapped.reads) / 1000000);
    
    norm.factors$Norm_TC <- norm.factors$M / mean(norm.factors$M);



    if(calc.DESeq2 | calc.edgeR){
       df <- res@qc.data[["geneCounts"]][[1]];
       is.gene <- substr(df$GENEID,1,1) != "_";

       count.matrix <- do.call(cbind.data.frame, lapply(res@qc.data[["geneCounts"]], function(X){ X$COUNT } ));
       count.matrix <- count.matrix[is.gene,];
    }
    
    if(calc.DESeq2 & "DESeq2" %in% rownames(installed.packages())){
      message("Calculating DESeq2 Normalization Factors (Geometric normalization)...\n");
      require("DESeq2");
      
      norm.factors$Norm_Geo <- estimateSizeFactorsForMatrix(count.matrix);
    }
    
    if(calc.edgeR & "edgeR" %in% rownames(installed.packages())){
      message("Calculating edgeR Normalization Factors (all edgeR normalizations)...\n");
      require("edgeR");
      norm.factors$Norm_TMM <- calcNormFactors(count.matrix, method="TMM");
      norm.factors$Norm_UQ <- calcNormFactors(count.matrix, method="upperquartile");
      norm.factors$Norm_RLE <- calcNormFactors(count.matrix, method="RLE");
    }
    res@calc.data[["norm.factors.byLaneBam"]] <- norm.factors;
    return(res);
}


calc.lanebamwise.bysample.norm.factors <- function(res){
    #INCOMPLETE!
    samples <- res@sample.list;
    mapped.reads <- lapply( res@qc.data[["summary"]], function(df){
       df$COUNT[df$FIELD == "READ_PAIR_OK"];
    });

    norm.factors.bySample <- res@calc.data[["norm.factors.bySample"]];
    #norm.factors <- data.frame(lanebam.ID = res@decoder$lanebam.ID, sample.ID = res@decoder$sample.ID, lanebam.norm = rep(1,length(res@decoder$sample.ID)));
    
    
    bySample.list <- lapply(1:length(samples), function(i){
      s <- samples[i];
      in.df <- norm.factors.bySample[norm.factors.bySample$sample.ID == s,];
      bamfiles.for.sample <- res@decoder$lanebam.ID[res@decoder$sample.ID == s];
      curr.mapped.reads <- unlist(mapped.reads[bamfiles.for.sample]);
      lanebam.norm.factor <- curr.mapped.reads / mean(curr.mapped.reads);
      all.sample.norm.factors <- data.frame(t(sapply(1:length(bamfiles.for.sample),function(j){
        as.numeric(c(in.df[-1] * lanebam.norm.factor[j]));
      })));
      all.sample.norm.factors <- cbind(lanebam.ID = bamfiles.for.sample,all.sample.norm.factors);
      names(all.sample.norm.factors) <- c("lanebam.ID",names(in.df[-1]));
      all.sample.norm.factors;
    });

    norm.factors.byLaneBam <- do.call(rbind.data.frame, bySample.list);
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
    out.list <- lapply(1:length(res@decoder$lanebam.ID),function(i){
       if(is.null(res@decoder$multi.mapped.read.pair.count)){
         data.frame(FIELD = c("total.reads","mapped.reads","mapping.rate"), COUNT = c(total.reads[i], mapped.reads[i], mapping.rate[i]));
       } else {
         data.frame(FIELD = c("total.reads","mapped.reads","mapping.rate","mm.reads","mm.rate"), COUNT = c(total.reads[i], mapped.reads[i], mapping.rate[i],mm.reads[i], mm.rate[i]));
       }
    });
    names(out.list) <-  res@decoder$lanebam.ID;

    res@calc.data[["map.rates"]] <- out.list;
  }
  return(res);
}

calc.quals <- function(res){
   if(! is.null(res@qc.data[["quals.r1"]])){
      res@decoder$cycle.CT <- get.cycleCt(res@qc.data[["quals.r1"]], res@decoder$lanebam.ID);
   }
   return(res);
}

get.cycleCt <- function(data.list, lanebam.IDs){
   return(sapply(lanebam.IDs, FUN=function(lid){
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
    lanebam.set <- res@decoder$lanebam.ID[res@decoder$sample.ID == sample];
    if(length(lanebam.set) == 1){
       return(res@qc.data[["geneCounts"]][[lanebam.set[1]]]);
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
    lanebam.set <- res@decoder$lanebam.ID[res@decoder$sample.ID == sample];
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
    lanebam.set <- res@decoder$lanebam.ID[res@decoder$sample.ID == sample];
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
  if(is.null(decoder$lanebam.ID)) stop("Decoder formatting error: no column labelled lanebam.ID");
  if(is.null(decoder$lane.ID)) stop("Decoder formatting error: no column labelled lane.ID");
  if(is.null(decoder$group.ID)) stop("Decoder formatting error: no column labelled group.ID");
  if(is.null(decoder$sample.ID)) stop("Decoder formatting error: no column labelled sample.ID");
  if(length(unique(decoder$lanebam.ID)) < length(decoder$lanebam.ID)) stop("Decoder formatting error: lanebam.ID must be unique!");
  #more checks?
  if( length(unique(names(decoder))) != length(names(decoder)) ) stop("Decoder formatting error: column names must be unique!");
  read.in.results.data.with.decoder(decoder, infile.dir);
}

#DEPRECIATED:
read.in.results.data.with.two.decoder.files <- function(infile.dir, lanebam.decoder, sample.decoder){
  sample.decoder <- read.table(sample.decoder,header=TRUE,stringsAsFactors=F);
  lanebam.decoder <- read.table(lanebam.decoder,header=TRUE,stringsAsFactors=F);
  if(is.null(lanebam.decoder$lanebam.ID)) stop("Lanebam Decoder formatting error: no column labelled lanebam.ID");
  if(is.null(lanebam.decoder$lane.ID)) stop("Lanebam Decoder formatting error: no column labelled lane.ID");
  if(is.null(lanebam.decoder$sample.ID)) stop("Lanebam Decoder formatting error: no column labelled sample.ID");
  if(length(unique(lanebam.decoder$lanebam.ID)) < length(lanebam.decoder$lanebam.ID)) stop("Lanebam Decoder formatting error: lanebam.ID must be unique!");
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