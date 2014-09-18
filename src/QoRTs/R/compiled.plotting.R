
#######################################################################
#######################################################################
#######################################################################
# Plotting many plots:

DEFAULTDEBUGMODE = TRUE;

makePlot.all.std <- function(res, 
                         outfile.prefix = "./", 
                         plotter.params = list(), 
                         plot.device.name = "png", 
                         plotting.device.params = list(), 
                         debugMode = DEFAULTDEBUGMODE , ...){

  makePlot.summary.basic(res = res, outfile.prefix = outfile.prefix, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, ...);
  makePlot.summary.colorByGroup(res = res, outfile.prefix = outfile.prefix, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, ...);
  makePlot.summary.colorByLane(res = res, outfile.prefix = outfile.prefix, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, ...);

  makePlot.summary.sample.highlight.all(res = res,outfile.prefix = outfile.prefix,plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, ...);
  
  makePlot.summary.sample.highlight.andColorByLane.all(res = res,outfile.prefix = outfile.prefix,plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, ...);
  
}

#######################################################################
#######################################################################
#######################################################################

makePlot.summary.sample.highlight.all <- function(res, outfile.prefix = "./", plotter.params = list(), plot.device.name = "png", plotting.device.params = list(), verbose = TRUE, debugMode = DEFAULTDEBUGMODE , ...){
  for(curr.sample in unique(res@decoder$sample.ID)){
    makePlot.summary.sample.highlight(res = res, curr.sample = curr.sample,verbose = FALSE, outfile.prefix = outfile.prefix, plotter.params = plotter.params, plot.device.name = plot.device.name, ...);
    if(verbose) message(paste0(curr.sample," complete!"));
  }
}

makePlot.summary.sample.highlight.andColorByLane.all <- function(res, outfile.prefix = "./", plotter.params = list(), plot.device.name = "png", plotting.device.params = list(), verbose = TRUE, debugMode = DEFAULTDEBUGMODE , ...){
  for(curr.sample in unique(res@decoder$sample.ID)){
    makePlot.summary.sample.highlight.andColorByLane(res = res,curr.sample = curr.sample,verbose = FALSE, outfile.prefix = outfile.prefix, plotter.params = plotter.params, plot.device.name = plot.device.name, ...);
    if(verbose) message(paste0(curr.sample," complete!"));
  }
}


#######################################################################
#######################################################################
#######################################################################
# Summary Plot:

makePlot.summary.basic <- function(res,  outfile.prefix = NULL, outfile.suffix = "plot.full.summary.png", plotter.params = list(), plot.device.name = "png", plotting.device.params = list(),verbose = TRUE, debugMode = DEFAULTDEBUGMODE , ...){  

    build.plotter.function <- function(){
      build.plotter.summary(res, plotter.params = plotter.params);
    }
    
    plot.summary.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile.prefix = outfile.prefix, 
                         outfile.suffix = outfile.suffix, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = FALSE, ...);

}

#######################################################################
#######################################################################
#######################################################################

makePlot.summary.colorByGroup <- function(res, outfile.prefix = NULL, outfile.suffix = "plot.colorByGroup.png", plotter.params = list(), plot.device.name = "png", plotting.device.params = list(),verbose = TRUE, debugMode = DEFAULTDEBUGMODE , ...){
    build.plotter.function <- function(){
      build.plotter.colorByGroup(res, plotter.params = plotter.params);
    }
    
    plot.summary.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile.prefix = outfile.prefix, 
                         outfile.suffix = outfile.suffix, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = FALSE, ...);
}

#######################################################################
#######################################################################
#######################################################################

makePlot.summary.colorByLane <- function(res, outfile.prefix = NULL, outfile.suffix = "plot.colorByLane.png", plotter.params = list(), plot.device.name = "png", plotting.device.params = list(), verbose = TRUE, debugMode = DEFAULTDEBUGMODE , ...){
    build.plotter.function <- function(){
      build.plotter.colorByLane(res, plotter.params = plotter.params);
    }
    
    plot.summary.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile.prefix = outfile.prefix, 
                         outfile.suffix = outfile.suffix, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = FALSE, ...);
}


#######################################################################
#######################################################################
#######################################################################

makePlot.summary.sample.highlight <- function(res,curr.sample,
                                              outfile.prefix = NULL, 
                                              outfile.suffix = paste0("plot.summary.sample.",curr.sample,".png"), 
                                              plotter.params = list(), 
                                              plot.device.name = "png", 
                                              plotting.device.params = list(), 
                                              verbose = TRUE, debugMode = DEFAULTDEBUGMODE , ...){

    build.plotter.function <- function(){
      build.plotter.highlightBySample(curr.sample,res, merge.offset.outgroup = FALSE, plotter.params = plotter.params);
    }
    
    plot.summary.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile.prefix = outfile.prefix, 
                         outfile.suffix = outfile.suffix, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = TRUE, 
                         cdf.plotIntercepts = TRUE, ...);
}

#######################################################################
#######################################################################
#######################################################################

makePlot.summary.sample.highlight.andColorByLane <- function(res,curr.sample,  outfile.prefix = NULL, outfile.suffix = paste0("plot.summary.sample.coloredByLane.",curr.sample,".png"), plotter.params = list(),  plot.device.name = "png", plotting.device.params = list(), verbose = TRUE, debugMode = DEFAULTDEBUGMODE , ...){
    build.plotter.function <- function(){
      build.plotter.highlightBySample.colorByLane(curr.sample,res, merge.offset.outgroup = FALSE, plotter.params = plotter.params);
    }
    
    plot.summary.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile.prefix = outfile.prefix, 
                         outfile.suffix = outfile.suffix, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = TRUE, ...);
}

plot.summary.GENERIC <- function(res, build.plotter.function, outfile.prefix, outfile.suffix, plot.device.name = "png", plotting.device.params = list(), complete.extended = TRUE, debugMode, verbose, ... ){
  if(is.null(outfile.prefix)){
    if(debugMode) message("Plotting to the default device...");
    default.params <- list();
    dev.params <- overmerge.list(default.params,plotting.device.params);
    devCloseFunct <- QoRTs.open.plotting.device("",NULL,dev.params);
  } else {
    full.out.filepath <- paste0(outfile.prefix,outfile.suffix);
    if(debugMode) message(paste0("Plotting to file:",full.out.filepath, " using device: ",plot.device.name));
    default.params <- list(height = 4000, width = 5000, pointsize = 36);
    if(plot.device.name == "png") default.params[["units"]] <- "px";
    dev.params <- overmerge.list(default.params,plotting.device.params);
    devCloseFunct <- QoRTs.open.plotting.device(full.out.filepath,plot.device.name, dev.params);
  }
  
  tryCatch({
    plotter <- build.plotter.function();
    if(plot.device.name == "CairoPDF" | plot.device.name == "pdf"){
      if(debugMode) if(debugMode) message("Plotting pdf...");
      INTERNAL.plot.summaries.pdf(res = res, plotter = plotter, verbose = verbose, debugMode = debugMode, ...);
    } else if(! complete.extended) {
      if(debugMode) message("Plotting std...");
      INTERNAL.plot.summaries(res = res, plotter = plotter, verbose = verbose, debugMode = debugMode, ...);
    } else {
      if(debugMode) message("Plotting extended...");
      INTERNAL.plot.summaries.advanced(res = res, plotter = plotter, verbose = verbose, debugMode = debugMode, ...);
    }
  }, error = function(e){
    message(paste0("ERROR: ", e));
  }, finally = {
    #message("CRASHED!");
    if(debugMode)  message("Done with plot.");
    devCloseFunct();
  });
}




#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

QoRTs.open.plotting.device <- function(filename, plot.device.name = "png", device.params = list()){
   device.params.final <- device.params;
   device.params.final[["filename"]] <- filename;
   
   funct.do.nothing <- function(){
      #do nothing!
   }
   
   if(is.null(plot.device.name)){
     return(funct.do.nothing);
   } else if(plot.device.name == "png"){
     do.call(png,device.params.final);
     return(dev.off)
   } else if(plot.device.name == "CairoPNG"){
     require("Cairo");
     do.call(CairoPNG,device.params.final);
     return(dev.off)
   } else if(plot.device.name == "CairoPDF" | plot.device.name == "pdf"){
     #do nothing!
     return(funct.do.nothing);
   } else {
     stop(paste0("FATAL ERROR: QoRTs.open.plotting.device: Unrecognized device type: ",plot.device.name));
   }
}


#######################################################################
#######################################################################
#######################################################################
#######################################################################
######### INTERNAL FUNCTIONS:
#######################################################################
#######################################################################
#######################################################################
#######################################################################

INTERNAL.plot.summaries <- function(res, plotter,verbose = TRUE, cdf.bySample = TRUE, cdf.plotIntercepts = TRUE, debugMode, ...){

  layout(matrix(c(1:24,c(25,26,26,27,28,29)), 5, 6, byrow = TRUE));
  
  if(verbose) {message(paste0("Starting compiled plot..."));}
     ts <- timestamp();
  
  makePlot.legend.box(plotter, debugMode = debugMode, ...);
  makePlot.qual.pair(plotter,"min", debugMode = debugMode, ...);
  makePlot.qual.pair(plotter,"lowerQuartile", debugMode = debugMode, ...);
  makePlot.qual.pair(plotter,"median", debugMode = debugMode, ...);
  makePlot.qual.pair(plotter,"upperQuartile", debugMode = debugMode, ...);
  makePlot.qual.pair(plotter,"max", debugMode = debugMode, ...);

  makePlot.clipping(plotter, debugMode = debugMode, ...);
  makePlot.cigarOp.byCycle(plotter,"Del", debugMode = debugMode, ...);
  makePlot.cigarOp.byCycle(plotter,"Ins", debugMode = debugMode, ...);
  makePlot.cigarOp.byCycle(plotter,"Splice", debugMode = debugMode, ...);
  #makePlot.cigarLength.distribution(plotter,"Ins", log.y = TRUE, debugMode = debugMode, ...);
  #makePlot.cigarLength.distribution(plotter,"Del", log.y = TRUE, debugMode = debugMode, ...);
  makePlot.gc(plotter, debugMode = debugMode, ...);
  makePlot.missingness.rate(plotter, debugMode = debugMode, ...);
  
  makePlot.insert.size(plotter, debugMode = debugMode, ...);
  makePlot.gene.cdf(plotter, sampleWise = cdf.bySample, plot.intercepts = cdf.plotIntercepts, debugMode = debugMode, ...)
  makePlot.genebody.coverage(plotter, debugMode = debugMode, ...);
  makePlot.genebody.coverage.UMQuartile(plotter, debugMode = debugMode, ...);
  makePlot.genebody.coverage.lowExpress(plotter, debugMode = debugMode, ...);
  makePlot.gene.assignment.rates(plotter, debugMode = debugMode, ...);
  
  makePlot.splice.junction.loci.counts(plotter, debugMode = debugMode, ...);
  makePlot.splice.junction.event.rates.split.by.type(plotter, debugMode = debugMode, ...);
  makePlot.splice.junction.event.rates(plotter, debugMode = debugMode, ...);
  makePlot.strandedness.test(plotter, debugMode = debugMode, ...);
  makePlot.mapping.rates(plotter, debugMode = debugMode, ...);
  makePlot.chrom.type.rates(plotter, chromosome.name.style = "UCSC_WITH_ERCC", exclude.autosomes=TRUE, debugMode = debugMode, ...);
  
  makePlot.norm.factors(plotter, debugMode = debugMode, ...);
  makePlot.raw.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, ...);
  makePlot.NVC.lead.clip(plotter, clip.amt = 12, points.highlighted = TRUE, debugMode = debugMode, ...);
  makePlot.NVC.tail.clip(plotter, clip.amt = 12, points.highlighted = TRUE, debugMode = debugMode, ...);  
  #makePlot.minus.clipping.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, ...);
  makePlot.dropped.rates(plotter, debugMode = debugMode, ...);

  if(verbose) {message(paste0("Finished all plots"));}
}

INTERNAL.plot.summaries.advanced <- function(res, plotter,verbose = TRUE, cdf.bySample = TRUE, cdf.plotIntercepts = TRUE, debugMode, cex.corner.label = 2, ...){

#FIGURE OUT NEW LAYOUT!
  layout(matrix(c(1:30,31,31,32,32,33), 5, 7, byrow = TRUE));
  
  if(verbose) {message(paste0("Starting compiled plot..."));}
     ts <- timestamp();
     
  a.to.z <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
  corner.labels <- c(a.to.z, paste0(a.to.z[1],a.to.z), paste0(a.to.z[2],a.to.z), paste0(a.to.z[3],a.to.z));
  
  plot.corner.label <- function(i){
    devlim <- device.limits();
    text(devlim[1],devlim[4], corner.labels[i] , cex = cex.corner.label, adj = c(-0.1,1.1),  xpd=T);
  }
  
  makePlot.legend.box(plotter, debugMode = debugMode, ...); 
  makePlot.qual.pair(plotter,"min", debugMode = debugMode, ...); plot.corner.label(1);
  makePlot.qual.pair(plotter,"lowerQuartile", debugMode = debugMode, ...); plot.corner.label(2);
  makePlot.qual.pair(plotter,"median", debugMode = debugMode, ...); plot.corner.label(3);
  makePlot.qual.pair(plotter,"upperQuartile", debugMode = debugMode, ...); plot.corner.label(4);
  makePlot.qual.pair(plotter,"max", debugMode = debugMode, ...); plot.corner.label(5);
  makePlot.clipping(plotter, debugMode = debugMode, ...); plot.corner.label(6);
  
  makePlot.cigarOp.byCycle(plotter,"Del", debugMode = debugMode, ...); plot.corner.label(7);
  makePlot.cigarOp.byCycle(plotter,"Ins", debugMode = debugMode, ...); plot.corner.label(8);
  makePlot.cigarOp.byCycle(plotter,"Splice", debugMode = debugMode, ...); plot.corner.label(9);
  makePlot.cigarLength.distribution(plotter,"Ins", log.y = TRUE, debugMode = debugMode, ...); plot.corner.label(10);
  makePlot.cigarLength.distribution(plotter,"Del", log.y = TRUE, debugMode = debugMode, ...); plot.corner.label(11);
  makePlot.gc(plotter, debugMode = debugMode, ...); plot.corner.label(12);
  makePlot.missingness.rate(plotter, debugMode = debugMode, ...); plot.corner.label(13);
  
  makePlot.dropped.rates(plotter, debugMode = debugMode, ...); plot.corner.label(14);
  makePlot.insert.size(plotter, debugMode = debugMode, ...); plot.corner.label(15);
  makePlot.gene.cdf(plotter, sampleWise = cdf.bySample, plot.intercepts = cdf.plotIntercepts, debugMode = debugMode, ...); plot.corner.label(16);
  makePlot.genebody.coverage(plotter, debugMode = debugMode, ...); plot.corner.label(17);
  makePlot.genebody.coverage.UMQuartile(plotter, debugMode = debugMode, ...); plot.corner.label(18);
  makePlot.genebody.coverage.lowExpress(plotter, debugMode = debugMode, ...); plot.corner.label(19);
  makePlot.gene.assignment.rates(plotter, debugMode = debugMode, ...); plot.corner.label(20);
  
  makePlot.splice.junction.loci.counts(plotter, debugMode = debugMode, ...); plot.corner.label(21);
  makePlot.splice.junction.event.rates.split.by.type(plotter, debugMode = debugMode, ...); plot.corner.label(22);
  makePlot.splice.junction.event.rates(plotter, debugMode = debugMode, ...); plot.corner.label(23);
  makePlot.strandedness.test(plotter, debugMode = debugMode, ...); plot.corner.label(24);
  makePlot.mapping.rates(plotter, debugMode = debugMode, ...); plot.corner.label(25);
  makePlot.chrom.type.rates(plotter, chromosome.name.style = "UCSC_WITH_ERCC", exclude.autosomes=TRUE, debugMode = debugMode, ...); plot.corner.label(26);
  makePlot.norm.factors(plotter, debugMode = debugMode, ...); plot.corner.label(27);
  
  makePlot.NVC.lead.clip(plotter, clip.amt = 12, points.highlighted = TRUE, debugMode = debugMode, ...); plot.corner.label(28);
  makePlot.NVC.tail.clip(plotter, clip.amt = 12, points.highlighted = TRUE, debugMode = debugMode, ...); plot.corner.label(29);
  makePlot.raw.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, ...); plot.corner.label(30);
  makePlot.minus.clipping.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, ...); plot.corner.label(31);
  #makePlot.cigarMismatch(plotter, debugMode = debugMode, ...); plot.corner.label(33);
  makePlot.legend.box(plotter, debugMode = debugMode, ...); 

  if(verbose) {message(paste0("Finished all plots"));}
}

INTERNAL.plot.summaries.pdf <- function(pdf.outfile, pdf.device.name, res, plotter,verbose = TRUE, cdf.bySample = TRUE, cdf.plotIntercepts = TRUE, debugMode, cex.corner.label = 2, pdf.paper = "letter", pdf.pointsize = 10, pdf.width = 0, pdf.height = 0, ...){
  if(pdf.device.name == "pdf"){
     pdf(pdf.outfile, paper = pdf.paper, pointsize=pdf.pointsize, width=pdf.width, height=pdf.height);
  } else if(pdf.device.name == "CairoPDF"){
     require(Cairo);
     CairoPDF(pdf.outfile, paper = pdf.paper, pointsize=pdf.pointsize, width=pdf.width, height=pdf.height);
  } else {
     stop("unrecognized option for pdf.device.name: allowed values are pdf and CairoPDF");
  }

  a.to.z <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
  corner.labels <- c(a.to.z, paste0(a.to.z[1],a.to.z), paste0(a.to.z[2],a.to.z), paste0(a.to.z[3],a.to.z));
  
  plot.corner.label <- function(i){
    devlim <- device.limits();
    text(devlim[1],devlim[4], corner.labels[i] , cex = cex.corner.label, adj = c(-0.1,1.1), xpd=T);
  }
  
  if(verbose) {message(paste0("Starting compiled plot..."));}
     ts <- timestamp();

  layout(matrix(1:6,3,2,byrow=TRUE));
  makePlot.legend.box(plotter, debugMode = debugMode, ...); plot.corner.label(1);
  makePlot.qual.pair(plotter,"min", debugMode = debugMode, ...); plot.corner.label(2);
  makePlot.qual.pair(plotter,"lowerQuartile", debugMode = debugMode, ...); plot.corner.label(3);
  makePlot.qual.pair(plotter,"median", debugMode = debugMode, ...); plot.corner.label(4);
  makePlot.qual.pair(plotter,"upperQuartile", debugMode = debugMode, ...); plot.corner.label(5);
  makePlot.qual.pair(plotter,"max", debugMode = debugMode, ...); plot.corner.label(6);

  layout(matrix(1:6,3,2,byrow=TRUE));
  makePlot.clipping(plotter, debugMode = debugMode, ...); plot.corner.label(7);
  makePlot.cigarOp.byCycle(plotter,"Del", debugMode = debugMode, ...); plot.corner.label(8);
  makePlot.cigarOp.byCycle(plotter,"Ins", debugMode = debugMode, ...); plot.corner.label(9);
  makePlot.cigarOp.byCycle(plotter,"Splice", debugMode = debugMode, ...); plot.corner.label(10);
  makePlot.cigarLength.distribution(plotter,"Ins", log.y = TRUE, debugMode = debugMode, ...); plot.corner.label(11);
  makePlot.cigarLength.distribution(plotter,"Del", log.y = TRUE, debugMode = debugMode, ...); plot.corner.label(12);

  layout(matrix(1:6,3,2,byrow=TRUE));
  makePlot.gc(plotter, debugMode = debugMode, ...); plot.corner.label(13);
  makePlot.missingness.rate(plotter, debugMode = debugMode, ...); plot.corner.label(14);
  makePlot.dropped.rates(plotter, debugMode = debugMode, ...); plot.corner.label(15);
  makePlot.insert.size(plotter, debugMode = debugMode, ...); plot.corner.label(16);
  makePlot.gene.cdf(plotter, sampleWise = cdf.bySample, plot.intercepts = cdf.plotIntercepts, debugMode = debugMode, ...); plot.corner.label(17);
  makePlot.genebody.coverage(plotter, debugMode = debugMode, ...); plot.corner.label(18);
  
  layout(matrix(1:6,3,2,byrow=TRUE));
  makePlot.genebody.coverage.UMQuartile(plotter, debugMode = debugMode, ...); plot.corner.label(19);
  makePlot.genebody.coverage.lowExpress(plotter, debugMode = debugMode, ...); plot.corner.label(20);
  makePlot.gene.assignment.rates(plotter, debugMode = debugMode, ...); plot.corner.label(21);
  makePlot.splice.junction.loci.counts(plotter, debugMode = debugMode, ...); plot.corner.label(22);
  makePlot.splice.junction.event.rates.split.by.type(plotter, debugMode = debugMode, ...); plot.corner.label(23);
  makePlot.splice.junction.event.rates(plotter, debugMode = debugMode, ...); plot.corner.label(24);
  
  layout(matrix(1:6,3,2,byrow=TRUE));
  makePlot.strandedness.test(plotter, debugMode = debugMode, ...); plot.corner.label(25);
  makePlot.mapping.rates(plotter, debugMode = debugMode, ...); plot.corner.label(26);
  makePlot.chrom.type.rates(plotter, chromosome.name.style = "UCSC_WITH_ERCC", exclude.autosomes=TRUE, debugMode = debugMode, ...); plot.corner.label(27);
  makePlot.norm.factors(plotter, debugMode = debugMode, ...); plot.corner.label(28);
  makePlot.NVC.lead.clip(plotter, clip.amt = 12, points.highlighted = TRUE, debugMode = debugMode, ...); plot.corner.label(29);
  makePlot.NVC.tail.clip(plotter, clip.amt = 12, points.highlighted = TRUE, debugMode = debugMode, ...); plot.corner.label(30);
  
  layout(matrix(c(1,1,2,2),3,2,byrow=TRUE));
  makePlot.raw.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, ...); plot.corner.label(31);
  makePlot.minus.clipping.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, ...); plot.corner.label(32);
  #makePlot.cigarMismatch(plotter, debugMode = debugMode, ...); plot.corner.label(33);

  if(verbose) {message(paste0("Finished all plots"));}
  dev.off();
}