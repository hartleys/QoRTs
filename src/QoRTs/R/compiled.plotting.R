
#######################################################################
#######################################################################
#######################################################################
# Plotting many plots:

DEFAULTDEBUGMODE = TRUE;

makeMultiPlot.all <- function(res, outfile.dir = "./",
                         plotter.params = list(), 
                         plot.device.name = "png", 
                         plotting.device.params = list(), 
                         debugMode = DEFAULTDEBUGMODE , 
                         rasterize.large.plots = NULL, 
                         raster.height = 1000, 
                         raster.width = 1000,
                         exclude.autosomes.chrom.rate.plot = TRUE,
                         chromosome.name.style = "UCSC",
                         fig.res = 150, fig.base.height.inches = 7,
                         ...){

  get.summary.table(res, outfile = paste0(outfile.dir,"summary.table.txt"), debugMode = debugMode);

                                makeMultiPlot.basic(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,...);
                        makeMultiPlot.colorBySample(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches, ...);
                         makeMultiPlot.colorByGroup(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches, ...);
                          makeMultiPlot.colorByLane(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches, ...);
                  makeMultiPlot.highlightSample.all(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches, ...);
      makeMultiPlot.highlightSample.colorByLane.all(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches, ...);

}

#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.highlightSample.all <- function(res, outfile.dir = "./",
                                                  plotter.params = list(), 
                                                  plot.device.name = "png", 
                                                  plotting.device.params = list(), 
                                                  verbose = TRUE,
                                                  debugMode = DEFAULTDEBUGMODE , 
                                                  rasterize.large.plots = NULL, 
                                                  raster.height = 1000, 
                                                  raster.width = 1000,
                                                  exclude.autosomes.chrom.rate.plot = TRUE,
                                                  chromosome.name.style = "UCSC",
                                                  fig.res = 150, fig.base.height.inches = 7,
                                                  ...){
  for(curr.sample in unique(res@decoder$sample.ID)){
    makeMultiPlot.highlightSample(res = res, 
                                     curr.sample = curr.sample, 
                                     outfile.dir = outfile.dir,
                                     verbose = FALSE, 
                                     plotter.params = plotter.params,
                                     plot.device.name = plot.device.name,
                                     rasterize.large.plots = rasterize.large.plots, 
                                     raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                                     fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                                     ...);
    if(verbose) message(paste0(curr.sample," complete!"));
  }
}

makeMultiPlot.highlightSample.colorByLane.all <- function(res, outfile.dir = "./",
                                                                 plotter.params = list(), 
                                                                 plot.device.name = "png", 
                                                                 plotting.device.params = list(), 
                                                                 verbose = TRUE, 
                                                                 debugMode = DEFAULTDEBUGMODE ,
                                                                 rasterize.large.plots = NULL, 
                                                                 raster.height = 1000, 
                                                                 raster.width = 1000,
                                                                 exclude.autosomes.chrom.rate.plot = TRUE,
                                                                 chromosome.name.style = "UCSC",
                                                                 fig.res = 150, fig.base.height.inches = 7,
                                                                 ...){
  for(curr.sample in unique(res@decoder$sample.ID)){
    makeMultiPlot.highlightSample.colorByLane(res = res,
                    outfile.dir = outfile.dir,
                    curr.sample = curr.sample,
                    verbose = FALSE, 
                    plotter.params = plotter.params, 
                    plot.device.name = plot.device.name, 
                    rasterize.large.plots = rasterize.large.plots, 
                    raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                    fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                    ...);
    if(verbose) message(paste0(curr.sample," complete!"));
  }
}


#######################################################################
#######################################################################
#######################################################################
# Summary Plot:

makeMultiPlot.basic <- function(res,  outfile = NULL, 
                                   outfile.dir = "./",
                                   outfile.prefix = "plot-basic", 
                                   outfile.ext = NULL, 
                                   plotter.params = list(), 
                                   plot.device.name = "curr", 
                                   plotting.device.params = list(),
                                   verbose = TRUE, 
                                   debugMode = DEFAULTDEBUGMODE , 
                                   rasterize.large.plots = NULL, 
                                   raster.height = 1000, 
                                   raster.width = 1000,
                                   separatePlots = FALSE,
                                   exclude.autosomes.chrom.rate.plot = TRUE,
                                   chromosome.name.style = "UCSC",
                                   fig.res = 150, fig.base.height.inches = 7,
                                   ...){  

    build.plotter.function <- function(){
      build.plotter.basic(res, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile = outfile,
                         outfile.dir = outfile.dir,
                         outfile.prefix = outfile.prefix, 
                         outfile.ext = outfile.ext, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = FALSE, 
                         rasterize.large.plots = rasterize.large.plots,
                         raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                         separatePlots = separatePlots,
                         nvc.highlight.points = FALSE,
                         fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                         ...);

}

#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.colorByGroup <- function(res, outfile = NULL, 
                                          outfile.dir = "./",
                                          outfile.prefix = "plot-colorByGroup", 
                                          outfile.ext = NULL, 
                                          plotter.params = list(), 
                                          plot.device.name = "curr", 
                                          plotting.device.params = list(),
                                          verbose = TRUE, 
                                          debugMode = DEFAULTDEBUGMODE , 
                                          rasterize.large.plots = NULL,
                                          raster.height = 1000, 
                                          raster.width = 1000,
                                          separatePlots = FALSE,
                                          exclude.autosomes.chrom.rate.plot = TRUE,
                                          chromosome.name.style = "UCSC",
                                          fig.res = 150, fig.base.height.inches = 7,
                                          ...){
    build.plotter.function <- function(){
      build.plotter.colorByGroup(res, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile = outfile,
                         outfile.dir = outfile.dir,
                         outfile.prefix = outfile.prefix, 
                         outfile.ext = outfile.ext, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = FALSE, 
                         rasterize.large.plots = rasterize.large.plots, 
                         raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                         separatePlots = separatePlots, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                         ...);
}

#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.colorByLane <- function(res, outfile = NULL, 
                                         outfile.dir = "./",
                                         outfile.prefix = "plot-colorByLane", 
                                         outfile.ext = NULL, 
                                         plotter.params = list(), 
                                         plot.device.name = "curr",
                                         plotting.device.params = list(), 
                                         verbose = TRUE, 
                                         debugMode = DEFAULTDEBUGMODE , 
                                         rasterize.large.plots = NULL, 
                                         raster.height = 1000, 
                                         raster.width = 1000,
                                         separatePlots = FALSE,
                                         exclude.autosomes.chrom.rate.plot = TRUE,
                                         chromosome.name.style = "UCSC",
                                         fig.res = 150, fig.base.height.inches = 7,
                                         ...){
    build.plotter.function <- function(){
      build.plotter.colorByLane(res, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile = outfile,
                         outfile.dir = outfile.dir,
                         outfile.prefix = outfile.prefix, 
                         outfile.ext = outfile.ext, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = FALSE, 
                         rasterize.large.plots = rasterize.large.plots,
                         raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                         separatePlots = separatePlots, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                         ...);
}

#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.colorBySample <- function(res, outfile = NULL, 
                                         outfile.dir = "./",
                                         outfile.prefix = "plot-colorBySample", 
                                         outfile.ext = NULL, 
                                         plotter.params = list(), 
                                         plot.device.name = "curr",
                                         plotting.device.params = list(), 
                                         verbose = TRUE, 
                                         debugMode = DEFAULTDEBUGMODE , 
                                         rasterize.large.plots = NULL, 
                                         raster.height = 1000, 
                                         raster.width = 1000,
                                         separatePlots = FALSE,
                                         exclude.autosomes.chrom.rate.plot = TRUE,
                                         chromosome.name.style = "UCSC",
                                         fig.res = 150, fig.base.height.inches = 7,
                                         ...){
    build.plotter.function <- function(){
      build.plotter.colorBySample(res, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile = outfile,
                         outfile.dir = outfile.dir,
                         outfile.prefix = outfile.prefix, 
                         outfile.ext = outfile.ext, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = FALSE, 
                         rasterize.large.plots = rasterize.large.plots,
                         raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                         separatePlots = separatePlots, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                         ...);
}



#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.highlightSample <- function(res,curr.sample,outfile = NULL, 
                                              outfile.dir = "./",
                                              outfile.prefix = paste0("plot-sampleHL-",curr.sample), 
                                              outfile.ext = NULL,
                                              plotter.params = list(), 
                                              plot.device.name = "curr", 
                                              plotting.device.params = list(), 
                                              verbose = TRUE, debugMode = DEFAULTDEBUGMODE , 
                                              rasterize.large.plots = NULL, 
                                              raster.height = 1000,
                                              raster.width = 1000,
                                              separatePlots = FALSE,
                                              exclude.autosomes.chrom.rate.plot = TRUE,
                                              chromosome.name.style = "UCSC",
                                              fig.res = 150, fig.base.height.inches = 7,
                                              ...){

    build.plotter.function <- function(){
      build.plotter.highlightSample(curr.sample,res, merge.offset.outgroup = FALSE, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile = outfile,
                         outfile.dir = outfile.dir,
                         outfile.prefix = outfile.prefix, 
                         outfile.ext = outfile.ext, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = TRUE, 
                         cdf.plotIntercepts = TRUE, 
                         rasterize.large.plots = rasterize.large.plots,
                         raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                         separatePlots = separatePlots, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                         ...);
}

#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.highlightSample.colorByLane <- function(res,curr.sample, outfile = NULL, 
                                                             outfile.dir = "./",
                                                             outfile.prefix = paste0("plot-sampleHL-coloredByLane-",curr.sample), 
                                                             outfile.ext = NULL, 
                                                             plotter.params = list(),  
                                                             plot.device.name = "curr", 
                                                             plotting.device.params = list(), 
                                                             verbose = TRUE, 
                                                             debugMode = DEFAULTDEBUGMODE , 
                                                             rasterize.large.plots = NULL,
                                                             raster.height = 1000,
                                                             raster.width = 1000,
                                                             separatePlots = FALSE,
                                                             exclude.autosomes.chrom.rate.plot = TRUE,
                                                             chromosome.name.style = "UCSC",
                                                             fig.res = 150, fig.base.height.inches = 7,
                                                             ...){
    build.plotter.function <- function(){
      build.plotter.highlightSample.colorByLane(curr.sample,res, merge.offset.outgroup = FALSE, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC(res = res, 
                         build.plotter.function = build.plotter.function, 
                         outfile = outfile,
                         outfile.dir = outfile.dir,
                         outfile.prefix = outfile.prefix, 
                         outfile.ext = outfile.ext, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         debugMode = debugMode, 
                         verbose = verbose, 
                         cdf.bySample = FALSE, 
                         cdf.plotIntercepts = TRUE, 
                         rasterize.large.plots = rasterize.large.plots,
                         raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                         separatePlots = separatePlots, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                         ...);
}

supportedDevices <- c("png","CairoPNG","tiff","jpeg","CairoPDF","pdf","svg","curr");
supportedVectorDevices <- c("pdf","CairoPDF","svg");

makeMultiPlot.GENERIC <- function(res, 
                                 build.plotter.function, 
                                 outfile = NULL,
                                 outfile.dir, 
                                 outfile.prefix,
                                 outfile.ext = NULL, 
                                 plot.device.name, 
                                 plotting.device.params = list(), 
                                 complete.extended = TRUE, 
                                 debugMode, verbose, 
                                 rasterize.large.plots = NULL, 
                                 raster.height = 1000,
                                 raster.width = 1000,
                                 separatePlots = FALSE,
                                 nvc.highlight.points = TRUE,
                                 exclude.autosomes.chrom.rate.plot = TRUE,
                                 chromosome.name.style = "UCSC",
                                 fig.res = 150, fig.base.height.inches = 7,
                                 ... ){



  height.per.px <- fig.res * fig.base.height.inches;
  width.per.px <- fig.res * fig.base.height.inches;
  height.per.inches <- fig.base.height.inches;
  width.per.inches <- fig.base.height.inches;

  if(separatePlots){
    height.mult <- 1;
    width.mult <- 1;
  } else if(complete.extended){
    height.mult <- 5;
    width.mult <- 7;
  } else {
    height.mult <- 5;
    width.mult <- 6;
  }

  if(is.null(rasterize.large.plots)){
    if(plot.device.name %in% supportedVectorDevices){
      if(check.rasterize.or.warn("rasterize.large.plots")){
        rasterize.large.plots = TRUE;
      } else {
        rasterize.large.plots = FALSE;
      }
    } else {
      rasterize.large.plots = FALSE;
    }
  }
  if(rasterize.large.plots){
    check.rasterize.or.die("rasterize.large.plots");
    
    if((! plot.device.name %in% supportedVectorDevices) & (plot.device.name != "curr")){
      warning("rasterize.large.plots = TRUE should not be used with raster file formats (ie png, tiff, jpeg, etc). This will result in image degradation.");
    }
  }
  if(debugMode) message("Rasterize large plots: ", rasterize.large.plots);

  if(plot.device.name == "curr"){
    if(debugMode) message("Plotting to the currently-open device...");
    default.params <- list();
    dev.params <- list();
    devOpenFunct <- function(f,w){};
    devCloseFunct <- function(){};
  } else if(plot.device.name == "png"){
    if(is.null(outfile.ext)) outfile.ext = ".png";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(png,dev.params)
    };
    devCloseFunct <- function(){dev.off()};
  } else if(plot.device.name == "tiff"){
    if(is.null(outfile.ext)) outfile.ext = ".tiff";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(tiff,dev.params)
    };
    devCloseFunct <- function(){dev.off()}    
  } else if(plot.device.name == "jpeg"){
    if(is.null(outfile.ext)) outfile.ext = ".jpg";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(jpeg,dev.params)
    };
    devCloseFunct <- function(){dev.off()}        
  } else if(plot.device.name == "tiff"){
    if(is.null(outfile.ext)) outfile.ext = ".tiff";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(tiff,dev.params)
    };
    devCloseFunct <- function(){dev.off()}        
  } else if(plot.device.name == "CairoPNG"){
    if(! require(Cairo)) stop("Error: package Cairo not found. Install package Cairo or set plot.device.name to something other than CairoPNG or CairoPDF.");

    if(is.null(outfile.ext)) outfile.ext = ".png";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(CairoPNG,dev.params)
    };
    devCloseFunct <- function(){dev.off()}
  } else if(plot.device.name == "svg"){
    if(is.null(outfile.ext)) outfile.ext = ".svg";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w){
      default.params <- list(filename = f, height = height.per.inches * height.mult, width = width.per.inches * width.mult * w, pointsize = 24);
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(svg,dev.params)
    };
    devCloseFunct <- function(){dev.off()}
  } else if(plot.device.name == "pdf"){
    if(is.null(outfile.ext)) outfile.ext = ".pdf";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w){
      default.params <- list(file = f, height = 0, width = 0, pointsize = 10, paper = "letter");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(pdf,dev.params)
    };
    devCloseFunct <- function(){dev.off()}
  } else if(plot.device.name == "CairoPDF"){
    if(! require(Cairo)) stop("Error: package Cairo not found. Install package Cairo or set plot.device.name to something other than CairoPNG or CairoPDF.");
    if(is.null(outfile.ext)) outfile.ext = ".pdf";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w){
      default.params <- list(file = f, height = 0, width = 0, pointsize = 10, paper = "letter");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(CairoPDF,dev.params)
    };
    devCloseFunct <- function(){dev.off()}
  } else {
    stop(paste0("Error: graphics device \"",plot.device.name,"\" not supported! Legal options are: [",paste0(supportedDevices,collapse=","),"] Set plot.device.name to \"curr\" to plot to the currently-open and/or default device."));
  }
  
  if(! separatePlots){
    if(debugMode & plot.device.name != "curr") {
      message("Opening file: ",outfile);
      tempFunc <- devCloseFunct;
      devCloseFunct <- function(){
        message("Closing file: ",outfile);
        tempFunc();
      }
    }
    devOpenFunct(outfile,1);
    separatePlotFuncts <- list(separatePlots = FALSE, outfilePrefix = "", outfileExt = "", devOpenFunct = function(f,w){}, devCloseFunct = function(){});
  } else {
    separatePlotFuncts <- list(separatePlots = TRUE,
                               outfilePrefix = paste0(outfile.dir,outfile.prefix), 
                               outfileExt = outfile.ext,
                               devOpenFunct = devOpenFunct,
                               devCloseFunct = devCloseFunct);
  }
  tryCatch({
    plotter <- build.plotter.function();
    if(plot.device.name == "CairoPDF" | plot.device.name == "pdf"){
      if(debugMode) if(debugMode) message("Plotting pdf...");
      INTERNAL.plot.summaries.pdf( res = res, plotter = plotter, nvc.highlight.points = nvc.highlight.points, verbose = verbose, debugMode = debugMode, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, 
                                   exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                                   chromosome.name.style = chromosome.name.style,
                                   ...);
    } else if(! complete.extended) {
      if(debugMode) message("Plotting std...");
      INTERNAL.plot.summaries(res = res, plotter = plotter, nvc.highlight.points = nvc.highlight.points, verbose = verbose, debugMode = debugMode, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, separatePlotFuncts = separatePlotFuncts, ...);
    } else {
      if(debugMode) message("Plotting extended...");
      INTERNAL.plot.summaries.advanced(res = res, plotter = plotter, nvc.highlight.points = nvc.highlight.points, verbose = verbose, debugMode = debugMode, rasterize.large.plots = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, separatePlotFuncts = separatePlotFuncts,
                                       exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                                       chromosome.name.style = chromosome.name.style,
                                       ...);
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

#DEPRECIATED:
INTERNAL.plot.summaries <- function(res, plotter,verbose = TRUE, cdf.bySample = TRUE, 
                                             nvc.highlight.points = TRUE,
                                             cdf.plotIntercepts = TRUE, debugMode, 
                                             rasterize.large.plots = FALSE,
                                             raster.height,
                                             raster.width,
                                             separatePlotFuncts,
                                             ...){

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
  makePlot.gene.cdf(plotter, sampleWise = cdf.bySample, plot.intercepts = cdf.plotIntercepts, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...)
  makePlot.genebody.coverage(plotter, debugMode = debugMode, ...);
  makePlot.genebody.coverage.UMQuartile(plotter, debugMode = debugMode, ...);
  makePlot.genebody.coverage.lowExpress(plotter, debugMode = debugMode, ...);
  makePlot.gene.assignment.rates(plotter, debugMode = debugMode, ...);
  
  makePlot.splice.junction.loci.counts(plotter, debugMode = debugMode, ...);
  makePlot.splice.junction.event.proportionsByType(plotter, debugMode = debugMode, ...);
  makePlot.splice.junction.event.ratesPerRead(plotter, debugMode = debugMode, ...);
  makePlot.strandedness.test(plotter, debugMode = debugMode, ...);
  makePlot.mapping.rates(plotter, debugMode = debugMode, ...);
  makePlot.chrom.type.rates(plotter, chromosome.name.style = "UCSC", exclude.autosomes=TRUE, debugMode = debugMode, ...);
  
  makePlot.norm.factors(plotter, debugMode = debugMode, ...);
  makePlot.raw.NVC(plotter, points.highlighted = nvc.highlight.points, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = 2 * raster.width, ...);
  makePlot.NVC.lead.clip(plotter, clip.amt = 12, points.highlighted = nvc.highlight.points, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...);
  makePlot.NVC.tail.clip(plotter, clip.amt = 12, points.highlighted = nvc.highlight.points, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...);  
  #makePlot.minus.clipping.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, ...);
  makePlot.dropped.rates(plotter, debugMode = debugMode, ...);

  if(verbose) {message(paste0("Finished all plots"));}
}

INTERNAL.plot.summaries.advanced <- function(res, plotter,verbose = TRUE, cdf.bySample = TRUE, 
                                             nvc.highlight.points = TRUE,
                                             cdf.plotIntercepts = TRUE, debugMode, cex.corner.label = 2, 
                                             rasterize.large.plots = FALSE, 
                                             raster.height,
                                             raster.width,
                                             separatePlotFuncts,
                                             exclude.autosomes.chrom.rate.plot = TRUE,
                                             chromosome.name.style = "UCSC",
                                             ...){
  outfilePrefix <- separatePlotFuncts$outfilePrefix;
  outfileExt <- separatePlotFuncts$outfileExt;
  devOpenFunct <- separatePlotFuncts$devOpenFunct;
  devCloseFunct <- separatePlotFuncts$devCloseFunct;
  separatePlots <- separatePlotFuncts$separatePlots;
  
  openFunc <- function(n,w){
    devOpenFunct(paste0(outfilePrefix,".",n,outfileExt), w);
  }
  
#FIGURE OUT NEW LAYOUT!
  layout(matrix(c(1:31,32,32,33,33), 5, 7, byrow = TRUE));
  
  if(verbose) {message(paste0("Starting compiled plot..."));}
     ts <- timestamp();
     
  a.to.z <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
  corner.labels <- c(a.to.z, paste0(a.to.z[1],a.to.z), paste0(a.to.z[2],a.to.z), paste0(a.to.z[3],a.to.z));
  
  plot.corner.label <- function(i){
    if(separatePlots){
      #do nothing!
    } else {
      devlim <- device.limits();
      text(devlim[1],devlim[4], corner.labels[i] , cex = cex.corner.label, adj = c(-0.1,1.1),  xpd=T);
    }
  }
  
  openFunc("legend",1); makePlot.legend.box(plotter, debugMode = debugMode, ...); devCloseFunct();
  openFunc("qual.pair.min",1); makePlot.qual.pair(plotter,"min", debugMode = debugMode, ...); plot.corner.label(1); devCloseFunct();
  openFunc("qual.pair.lowerQuartile",1); makePlot.qual.pair(plotter,"lowerQuartile", debugMode = debugMode, ...); plot.corner.label(2); devCloseFunct();
  openFunc("qual.pair.median",1); makePlot.qual.pair(plotter,"median", debugMode = debugMode, ...); plot.corner.label(3); devCloseFunct();
  openFunc("qual.pair.upperQuartile",1); makePlot.qual.pair(plotter,"upperQuartile", debugMode = debugMode, ...); plot.corner.label(4); devCloseFunct();
  openFunc("qual.pair.max",1); makePlot.qual.pair(plotter,"max", debugMode = debugMode, ...); plot.corner.label(5); devCloseFunct();
  openFunc("clippingProfile",1); makePlot.clipping(plotter, debugMode = debugMode, ...); plot.corner.label(6); devCloseFunct();
  
  openFunc("DeletionProfile",1); makePlot.cigarOp.byCycle(plotter,"Del", debugMode = debugMode, ...); plot.corner.label(7); devCloseFunct();
  openFunc("InsertionProfile",1); makePlot.cigarOp.byCycle(plotter,"Ins", debugMode = debugMode, ...); plot.corner.label(8); devCloseFunct();
  openFunc("SpliceProfile",1); makePlot.cigarOp.byCycle(plotter,"Splice", debugMode = debugMode, ...); plot.corner.label(9); devCloseFunct();
  openFunc("InsertionLengthHisto",1); makePlot.cigarLength.distribution(plotter,"Ins", log.y = TRUE, debugMode = debugMode, ...); plot.corner.label(10); devCloseFunct();
  openFunc("DeletionLengthHisto",1); makePlot.cigarLength.distribution(plotter,"Del", log.y = TRUE, debugMode = debugMode, ...); plot.corner.label(11); devCloseFunct();
  openFunc("gc",1); makePlot.gc(plotter, debugMode = debugMode, ...); plot.corner.label(12); devCloseFunct();
  openFunc("missingness.rate",1); makePlot.missingness.rate(plotter, debugMode = debugMode, ...); plot.corner.label(13); devCloseFunct();
  
  openFunc("dropped.rate",1); makePlot.dropped.rates(plotter, debugMode = debugMode, ...); plot.corner.label(14); devCloseFunct();
  openFunc("insert.size",1); makePlot.insert.size(plotter, debugMode = debugMode, ...); plot.corner.label(15); devCloseFunct();
  openFunc("gene.diversity",1); makePlot.gene.cdf(plotter, sampleWise = cdf.bySample, plot.intercepts = cdf.plotIntercepts, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...); plot.corner.label(16); devCloseFunct();
  openFunc("genebody.coverage.allGenes",1); makePlot.genebody.coverage(plotter, debugMode = debugMode, ...); plot.corner.label(17); devCloseFunct();
  openFunc("genebody.coverage.umquartileExpressionGenes",1); makePlot.genebody.coverage.UMQuartile(plotter, debugMode = debugMode, ...); plot.corner.label(18); devCloseFunct();
  openFunc("genebody.coverage.lowExpressionGenes",1); makePlot.genebody.coverage.lowExpress(plotter, debugMode = debugMode, ...); plot.corner.label(19); devCloseFunct();
  openFunc("geneAssignmentRates",1); makePlot.gene.assignment.rates(plotter, debugMode = debugMode, ...); plot.corner.label(20); devCloseFunct();
  
  openFunc("sj.locus.ct",1); makePlot.splice.junction.loci.counts(plotter, debugMode = debugMode, ...); plot.corner.label(21); devCloseFunct();
  openFunc("sj.event.proportionByType",1); makePlot.splice.junction.event.proportionsByType(plotter, debugMode = debugMode, ...); plot.corner.label(22);  devCloseFunct();
  openFunc("sj.event.rate",1); makePlot.splice.junction.event.ratesPerRead(plotter, debugMode = debugMode, ...); plot.corner.label(23); devCloseFunct();
  openFunc("mapping.rates",1); makePlot.mapping.rates(plotter, debugMode = debugMode, ...); plot.corner.label(24); devCloseFunct();
  openFunc("chrom.rates",1); makePlot.chrom.type.rates(plotter, chromosome.name.style = chromosome.name.style, exclude.autosomes = exclude.autosomes.chrom.rate.plot, debugMode = debugMode, ...); plot.corner.label(25); devCloseFunct();
  openFunc("norm.factors",1); makePlot.norm.factors(plotter, debugMode = debugMode, ...); plot.corner.label(26); devCloseFunct();
  openFunc("norm.vs.TC",1); makePlot.norm.factors.vs.TC(plotter, debugMode = debugMode, ...); plot.corner.label(27); devCloseFunct();
  
  openFunc("strandedness.test",1); makePlot.strandedness.test(plotter, debugMode = debugMode, ...); plot.corner.label(28); devCloseFunct();
  openFunc("NVC.lead.clip",1); makePlot.NVC.lead.clip(plotter, clip.amt = 12, points.highlighted = nvc.highlight.points, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...); plot.corner.label(29); devCloseFunct();
  openFunc("NVC.tail.clip",1); makePlot.NVC.tail.clip(plotter, clip.amt = 12, points.highlighted = nvc.highlight.points, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...); plot.corner.label(30); devCloseFunct();
  openFunc("NVC.raw",2); makePlot.raw.NVC(plotter, points.highlighted = nvc.highlight.points, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = 2* raster.width, ...); plot.corner.label(31); devCloseFunct();
  openFunc("NVC.aligned",2); makePlot.minus.clipping.NVC(plotter, points.highlighted = nvc.highlight.points, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = 2* raster.width, ...); plot.corner.label(32); devCloseFunct();
  #makePlot.cigarMismatch(plotter, debugMode = debugMode, ...); plot.corner.label(33); devCloseFunct();

  #openFunc("legend",1); makePlot.legend.box(plotter, debugMode = debugMode, ...);  devCloseFunct();

  if(verbose) {message(paste0("Finished all plots"));}
}



INTERNAL.plot.summaries.pdf <- function(res, plotter,
                                        verbose = TRUE, cdf.bySample = TRUE, 
                                        nvc.highlight.points = TRUE,
                                        cdf.plotIntercepts = TRUE, debugMode, 
                                        cex.corner.label = 2, 
                                        rasterize.large.plots, 
                                        raster.height,
                                        raster.width,
                                        exclude.autosomes.chrom.rate.plot = TRUE,
                                        chromosome.name.style = "UCSC",
                                        ...){

  #if(pdf.device.name == "pdf"){
  #   pdf(pdf.outfile, paper = pdf.paper, pointsize=pdf.pointsize, width=pdf.width, height=pdf.height);
  #} else if(pdf.device.name == "CairoPDF"){
   #  require(Cairo);
   #  CairoPDF(pdf.outfile, paper = pdf.paper, pointsize=pdf.pointsize, width=pdf.width, height=pdf.height);
  #} else {
  #   stop("unrecognized option for pdf.device.name: allowed values are pdf and CairoPDF");
  #}
  

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
  makePlot.gene.cdf(plotter, sampleWise = cdf.bySample, plot.intercepts = cdf.plotIntercepts, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...); plot.corner.label(17);
  makePlot.genebody.coverage(plotter, debugMode = debugMode, ...); plot.corner.label(18);
  
  layout(matrix(1:6,3,2,byrow=TRUE));
  makePlot.genebody.coverage.UMQuartile(plotter, debugMode = debugMode, ...); plot.corner.label(19);
  makePlot.genebody.coverage.lowExpress(plotter, debugMode = debugMode, ...); plot.corner.label(20);
  makePlot.gene.assignment.rates(plotter, debugMode = debugMode, ...); plot.corner.label(21);
  makePlot.splice.junction.loci.counts(plotter, debugMode = debugMode, ...); plot.corner.label(22);
  makePlot.splice.junction.event.proportionsByType(plotter, debugMode = debugMode, ...); plot.corner.label(23);
  makePlot.splice.junction.event.ratesPerRead(plotter, debugMode = debugMode, ...); plot.corner.label(24);
  
  layout(matrix(1:6,3,2,byrow=TRUE));
  makePlot.strandedness.test(plotter, debugMode = debugMode, ...); plot.corner.label(25);
  makePlot.mapping.rates(plotter, debugMode = debugMode, ...); plot.corner.label(26);
  makePlot.chrom.type.rates(plotter, chromosome.name.style = chromosome.name.style, exclude.autosomes=exclude.autosomes.chrom.rate.plot, debugMode = debugMode, ...); plot.corner.label(27);
  makePlot.norm.factors(plotter, debugMode = debugMode, ...); plot.corner.label(28);
  makePlot.norm.factors.vs.TC(plotter, debugMode = debugMode, ...); plot.corner.label(29);
  plot.new();
  
  layout(matrix(c(1,1,2,2,3,4),3,2,byrow=TRUE));
  makePlot.raw.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = 2 * raster.width, ...); plot.corner.label(30);
  makePlot.minus.clipping.NVC(plotter, points.highlighted = TRUE, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = 2 * raster.width, ...); plot.corner.label(31);
  makePlot.NVC.lead.clip(plotter, clip.amt = 12, points.highlighted = TRUE, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...); plot.corner.label(32);
  makePlot.NVC.tail.clip(plotter, clip.amt = 12, points.highlighted = TRUE, debugMode = debugMode, rasterize.plotting.area = rasterize.large.plots, raster.height = raster.height, raster.width = raster.width, ...); plot.corner.label(33);

  #makePlot.cigarMismatch(plotter, debugMode = debugMode, ...); plot.corner.label(33);

  if(verbose) {message(paste0("Finished all plots"));}
  #dev.off();
}
