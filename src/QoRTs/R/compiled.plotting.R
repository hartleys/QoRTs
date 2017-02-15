
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
                         rasterize.medium.plots = NULL,
                         raster.height = 1050, 
                         raster.width = 1050,
                         exclude.autosomes.chrom.rate.plot = TRUE,
                         chromosome.name.style = "UCSC",
                         fig.res = 150, fig.base.height.inches = 7,
                         insertSize.plot.xlim = NULL,
                         sequencing.type = c("RNA","Exome","Genome"),
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                         labelPlots=TRUE,
                         ...){
  sequencing.type <- match.arg(sequencing.type);
  get.summary.table(res, outfile = paste0(outfile.dir,"summary.table.txt"), debugMode = debugMode);

                                makeMultiPlot.basic(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, rasterize.medium.plots=rasterize.medium.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,maxColumns = maxColumns,maxRows=maxRows,plotList=plotList,labelPlots=labelPlots,...);
                        makeMultiPlot.colorBySample(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, rasterize.medium.plots=rasterize.medium.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,maxColumns = maxColumns,maxRows=maxRows,plotList=plotList,labelPlots=labelPlots, ...);
                         makeMultiPlot.colorByGroup(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, rasterize.medium.plots=rasterize.medium.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,maxColumns = maxColumns,maxRows=maxRows,plotList=plotList,labelPlots=labelPlots, ...);
                          makeMultiPlot.colorByLane(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, rasterize.medium.plots=rasterize.medium.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,maxColumns = maxColumns,maxRows=maxRows,plotList=plotList,labelPlots=labelPlots, ...);
                  makeMultiPlot.highlightSample.all(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, rasterize.medium.plots=rasterize.medium.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,maxColumns = maxColumns,maxRows=maxRows,plotList=plotList,labelPlots=labelPlots, ...);
      makeMultiPlot.highlightSample.colorByLane.all(res = res, outfile.dir = outfile.dir, plotter.params = plotter.params, plot.device.name = plot.device.name, plotting.device.params = plotting.device.params, rasterize.large.plots = rasterize.large.plots, rasterize.medium.plots=rasterize.medium.plots, raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot, fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,maxColumns = maxColumns,maxRows=maxRows,plotList=plotList,labelPlots=labelPlots, ...);

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
                                                  rasterize.medium.plots = NULL,
                                                  raster.height = 1050, 
                                                  raster.width = 1050,
                                                  exclude.autosomes.chrom.rate.plot = TRUE,
                                                  chromosome.name.style = "UCSC",
                                                  fig.res = 150, fig.base.height.inches = 7,
                                                  insertSize.plot.xlim = NULL,
                                                  sequencing.type = c("RNA","Exome","Genome"),
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots=TRUE,
                                   ...){
  sequencing.type <- match.arg(sequencing.type);

  for(curr.sample in unique(res@decoder$sample.ID)){
    makeMultiPlot.highlightSample(res = res, 
                                     curr.sample = curr.sample, 
                                     outfile.dir = outfile.dir,
                                     verbose = FALSE, 
                                     plotter.params = plotter.params,
                                     plot.device.name = plot.device.name,
                                     rasterize.large.plots = rasterize.large.plots, 
                                     rasterize.medium.plots = rasterize.medium.plots,
                                     raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                                     fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                                     insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,
                                     maxColumns = maxColumns,maxRows=maxRows,plotList=plotList,labelPlots=labelPlots,
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
                                                                 rasterize.medium.plots = NULL,
                                                                 raster.height = 1050, 
                                                                 raster.width = 1050,
                                                                 exclude.autosomes.chrom.rate.plot = TRUE,
                                                                 chromosome.name.style = "UCSC",
                                                                 fig.res = 150, fig.base.height.inches = 7,
                                                                 insertSize.plot.xlim = NULL,
                                                                 sequencing.type = c("RNA","Exome","Genome"),
                                                                 maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots=TRUE,
                                                                 ...){

  sequencing.type <- match.arg(sequencing.type);

  for(curr.sample in unique(res@decoder$sample.ID)){
    makeMultiPlot.highlightSample.colorByLane(res = res,
                    outfile.dir = outfile.dir,
                    curr.sample = curr.sample,
                    verbose = FALSE, 
                    plotter.params = plotter.params, 
                    plot.device.name = plot.device.name, 
                    rasterize.large.plots = rasterize.large.plots, 
                    rasterize.medium.plots = rasterize.medium.plots,
                    raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                    fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                    insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,
                    maxColumns = maxColumns,maxRows=maxRows,plotList=plotList,labelPlots=labelPlots,
                    ...);
    if(verbose) message(paste0(curr.sample," complete!"));
  }
}


#######################################################################
#######################################################################
#######################################################################
# Summary Plot:

makeMultiPlot.withPlotter <- function(plotter, res = plotter$res, outfile = NULL,  
                                   outfile.dir = "./",
                                   outfile.prefix = "plot-custom", 
                                   outfile.ext = NULL, 
                                   plotter.params = list(), 
                                   plot.device.name = "curr", 
                                   plotting.device.params = list(),
                                   verbose = TRUE, 
                                   debugMode = DEFAULTDEBUGMODE , 
                                   rasterize.large.plots = NULL, 
                                   rasterize.medium.plots = NULL,
                                   raster.height = 1050, 
                                   raster.width = 1050,
                                   separatePlots = FALSE,
                                   exclude.autosomes.chrom.rate.plot = TRUE,
                                   chromosome.name.style = "UCSC",
                                   fig.res = 150, fig.base.height.inches = 7,
                                   insertSize.plot.xlim = NULL,
                                   sequencing.type = c("RNA","Exome","Genome"),
                                   nvc.mark.points = TRUE,
                                   
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots = TRUE,
                                   
                                   plot = TRUE,
                                   
                                   ...){
  sequencing.type <- match.arg(sequencing.type);

    build.plotter.function <- function(){
      plotter;
    }
    
    makeMultiPlot.GENERIC.v10(res=res, 
               build.plotter.function=build.plotter.function, 
               outfile = outfile,
               outfile.dir=outfile.dir, 
               outfile.prefix=outfile.prefix,
               outfile.ext = outfile.ext, 
               
               plot.device.name=plot.device.name, 
               plotting.device.params = plotting.device.params, 
               
               debugMode = debugMode, verbose = verbose, 

               fig.res = fig.res, 
               fig.base.height.inches = fig.base.height.inches,
               rasterize.large.plots = rasterize.large.plots, 
               rasterize.medium.plots = rasterize.medium.plots, 
               raster.height = raster.height,
               raster.width = raster.width,
               raster.res = fig.res,
               
               separatePlots = separatePlots,
               splitPlots = FALSE,
               
               nvc.highlight.points = nvc.mark.points,
               exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
               chromosome.name.style = chromosome.name.style,
               insertSize.plot.xlim=insertSize.plot.xlim,
               cdf.bySample = FALSE, 
               cdf.plotIntercepts = FALSE, 
                         
               makePlots = plotList,
               sequencing.type = sequencing.type,
               
               labelPlots = labelPlots,
               maxColumns = maxColumns,
               maxRows = maxRows,
               
               plot = plot,
               
               ... )
    
   # makeMultiPlot.GENERIC(res = res, 
   #                      build.plotter.function = build.plotter.function, 
   #                      outfile = outfile,
   #                      outfile.dir = outfile.dir,
   #                      outfile.prefix = outfile.prefix, 
   #                      outfile.ext = outfile.ext, 
   #                      plot.device.name = plot.device.name, 
   #                      plotting.device.params = plotting.device.params, 
   #                      debugMode = debugMode, 
   #                      verbose = verbose, 
   #                      cdf.bySample = FALSE, 
   #                      cdf.plotIntercepts = FALSE, 
   #                      rasterize.large.plots = rasterize.large.plots,
   #                      rasterize.medium.plots = rasterize.medium.plots,
   #                      raster.height = raster.height, raster.width = raster.width, chromosome.name.style = chromosome.name.style, exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
   #                      separatePlots = separatePlots,
   #                      nvc.highlight.points = nvc.mark.points,
   #                      fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
   #                      insertSize.plot.xlim=insertSize.plot.xlim,sequencing.type=sequencing.type,
   #                      maxColumns = maxColumns,
   #                      ...);

}



makeMultiPlot.basic <- function(res,
                                   outfile = NULL, 
                                   outfile.dir = "./",
                                   outfile.prefix = "plot-basic", 
                                   outfile.ext = NULL, 
                                   plotter.params = list(), 
                                   plot.device.name = "curr", 
                                   plotting.device.params = list(),
                                   verbose = TRUE, 
                                   debugMode = DEFAULTDEBUGMODE , 
                                   rasterize.large.plots = NULL, 
                                   rasterize.medium.plots = NULL,
                                   raster.height = 1050, 
                                   raster.width = 1050,
                                   separatePlots = FALSE,
                                   exclude.autosomes.chrom.rate.plot = TRUE,
                                   chromosome.name.style = "UCSC",
                                   fig.res = 150, fig.base.height.inches = 7,
                                   insertSize.plot.xlim = NULL,
                                   sequencing.type = c("RNA","Exome","Genome"),
                                   nvc.mark.points = TRUE,
                                   
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots = TRUE,
                                   
                                   plot = TRUE,
                                   ...){  
  sequencing.type <- match.arg(sequencing.type);

    build.plotter.function <- function(){
      build.plotter.basic(res, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC.v10(res=res, 
               build.plotter.function=build.plotter.function, 
               outfile = outfile,
               outfile.dir=outfile.dir, 
               outfile.prefix=outfile.prefix,
               outfile.ext = outfile.ext, 
               
               plot.device.name=plot.device.name, 
               plotting.device.params = plotting.device.params, 
               
               debugMode = debugMode, verbose = verbose, 

               fig.res = fig.res, 
               fig.base.height.inches = fig.base.height.inches,
               rasterize.large.plots = rasterize.large.plots, 
               rasterize.medium.plots = rasterize.medium.plots, 
               raster.height = raster.height,
               raster.width = raster.width,
               raster.res = fig.res,
               
               separatePlots = separatePlots,
               splitPlots = FALSE,
               
               nvc.highlight.points = nvc.mark.points,
               exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
               chromosome.name.style = chromosome.name.style,
               insertSize.plot.xlim=insertSize.plot.xlim,
               cdf.bySample = FALSE, 
               cdf.plotIntercepts = FALSE, 
                         
               makePlots = plotList,
               sequencing.type = sequencing.type,
               
               labelPlots = labelPlots,
               maxColumns = maxColumns,
               maxRows = maxRows,
               
               plot = plot,
               
               ... )

}



#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.colorByGroup <- function(res,
                                   outfile = NULL, 
                                   outfile.dir = "./",
                                   outfile.prefix = "plot-colorByGroup", 
                                   outfile.ext = NULL, 
                                   plotter.params = list(), 
                                   plot.device.name = "curr", 
                                   plotting.device.params = list(),
                                   verbose = TRUE, 
                                   debugMode = DEFAULTDEBUGMODE , 
                                   rasterize.large.plots = NULL, 
                                   rasterize.medium.plots = NULL,
                                   raster.height = 1050, 
                                   raster.width = 1050,
                                   separatePlots = FALSE,
                                   exclude.autosomes.chrom.rate.plot = TRUE,
                                   chromosome.name.style = "UCSC",
                                   fig.res = 150, fig.base.height.inches = 7,
                                   insertSize.plot.xlim = NULL,
                                   sequencing.type = c("RNA","Exome","Genome"),
                                   nvc.mark.points = TRUE,
                                   
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots = TRUE,
                                   
                                   plot = TRUE,
                                   ...){  
    sequencing.type <- match.arg(sequencing.type);

    build.plotter.function <- function(){
      build.plotter.colorByGroup(res, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC.v10(res=res, 
               build.plotter.function=build.plotter.function, 
               outfile = outfile,
               outfile.dir=outfile.dir, 
               outfile.prefix=outfile.prefix,
               outfile.ext = outfile.ext, 
               
               plot.device.name=plot.device.name, 
               plotting.device.params = plotting.device.params, 
               
               debugMode = debugMode, verbose = verbose, 

               fig.res = fig.res, 
               fig.base.height.inches = fig.base.height.inches,
               rasterize.large.plots = rasterize.large.plots, 
               rasterize.medium.plots = rasterize.medium.plots, 
               raster.height = raster.height,
               raster.width = raster.width,
               raster.res = fig.res,
               
               separatePlots = separatePlots,
               splitPlots = FALSE,
               
               nvc.highlight.points = nvc.mark.points,
               exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
               chromosome.name.style = chromosome.name.style,
               insertSize.plot.xlim=insertSize.plot.xlim,
               cdf.bySample = FALSE, 
               cdf.plotIntercepts = FALSE, 
                         
               makePlots = plotList,
               sequencing.type = sequencing.type,
               
               labelPlots = labelPlots,
               maxColumns = maxColumns,
               maxRows = maxRows,
               
               plot = plot,
               
               ... )
}

#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.colorByLane <- function(res,
                                   outfile = NULL, 
                                   outfile.dir = "./",
                                   outfile.prefix = "plot-colorByLane", 
                                   outfile.ext = NULL, 
                                   plotter.params = list(), 
                                   plot.device.name = "curr", 
                                   plotting.device.params = list(),
                                   verbose = TRUE, 
                                   debugMode = DEFAULTDEBUGMODE , 
                                   rasterize.large.plots = NULL, 
                                   rasterize.medium.plots = NULL,
                                   raster.height = 1050, 
                                   raster.width = 1050,
                                   separatePlots = FALSE,
                                   exclude.autosomes.chrom.rate.plot = TRUE,
                                   chromosome.name.style = "UCSC",
                                   fig.res = 150, fig.base.height.inches = 7,
                                   insertSize.plot.xlim = NULL,
                                   sequencing.type = c("RNA","Exome","Genome"),
                                   nvc.mark.points = TRUE,
                                   
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots = TRUE,
                                   
                                   plot = TRUE,
                                   ...){  
  sequencing.type <- match.arg(sequencing.type);

    build.plotter.function <- function(){
      build.plotter.colorByLane(res, plotter.params = plotter.params);
    }
    
  makeMultiPlot.GENERIC.v10(res=res, 
               build.plotter.function=build.plotter.function, 
               outfile = outfile,
               outfile.dir=outfile.dir, 
               outfile.prefix=outfile.prefix,
               outfile.ext = outfile.ext, 
               
               plot.device.name=plot.device.name, 
               plotting.device.params = plotting.device.params, 
               
               debugMode = debugMode, verbose = verbose, 

               fig.res = fig.res, 
               fig.base.height.inches = fig.base.height.inches,
               rasterize.large.plots = rasterize.large.plots, 
               rasterize.medium.plots = rasterize.medium.plots, 
               raster.height = raster.height,
               raster.width = raster.width,
               raster.res = fig.res,
               
               separatePlots = separatePlots,
               splitPlots = FALSE,
               
               nvc.highlight.points = nvc.mark.points,
               exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
               chromosome.name.style = chromosome.name.style,
               insertSize.plot.xlim=insertSize.plot.xlim,
               cdf.bySample = FALSE, 
               cdf.plotIntercepts = FALSE, 
                         
               makePlots = plotList,
               sequencing.type = sequencing.type,
               
               labelPlots = labelPlots,
               maxColumns = maxColumns,
               maxRows = maxRows,
               
               plot = plot,
               
               ... )
}

#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.colorBySample <- function(res,
                                   outfile = NULL, 
                                   outfile.dir = "./",
                                   outfile.prefix = "plot-colorByLane", 
                                   outfile.ext = NULL, 
                                   plotter.params = list(), 
                                   plot.device.name = "curr", 
                                   plotting.device.params = list(),
                                   verbose = TRUE, 
                                   debugMode = DEFAULTDEBUGMODE , 
                                   rasterize.large.plots = NULL, 
                                   rasterize.medium.plots = NULL,
                                   raster.height = 1050, 
                                   raster.width = 1050,
                                   separatePlots = FALSE,
                                   exclude.autosomes.chrom.rate.plot = TRUE,
                                   chromosome.name.style = "UCSC",
                                   fig.res = 150, fig.base.height.inches = 7,
                                   insertSize.plot.xlim = NULL,
                                   sequencing.type = c("RNA","Exome","Genome"),
                                   nvc.mark.points = TRUE,
                                   
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots = TRUE,
                                   
                                   plot = TRUE,
                                   ...){
  sequencing.type <- match.arg(sequencing.type);

    build.plotter.function <- function(){
      build.plotter.colorBySample(res, plotter.params = plotter.params);
    }
    
  makeMultiPlot.GENERIC.v10(res=res, 
               build.plotter.function=build.plotter.function, 
               outfile = outfile,
               outfile.dir=outfile.dir, 
               outfile.prefix=outfile.prefix,
               outfile.ext = outfile.ext, 
               
               plot.device.name=plot.device.name, 
               plotting.device.params = plotting.device.params, 
               
               debugMode = debugMode, verbose = verbose, 

               fig.res = fig.res, 
               fig.base.height.inches = fig.base.height.inches,
               rasterize.large.plots = rasterize.large.plots, 
               rasterize.medium.plots = rasterize.medium.plots, 
               raster.height = raster.height,
               raster.width = raster.width,
               raster.res = fig.res,
               
               separatePlots = separatePlots,
               splitPlots = FALSE,
               
               nvc.highlight.points = nvc.mark.points,
               exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
               chromosome.name.style = chromosome.name.style,
               insertSize.plot.xlim=insertSize.plot.xlim,
               cdf.bySample = FALSE, 
               cdf.plotIntercepts = FALSE, 
                         
               makePlots = plotList,
               sequencing.type = sequencing.type,
               
               labelPlots = labelPlots,
               maxColumns = maxColumns,
               maxRows = maxRows,
               
               plot = plot,
               
               ... )
}



#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.highlightSample <- function(res, curr.sample,
                                   outfile = NULL, 
                                   outfile.dir = "./",
                                   outfile.prefix = paste0("plot-sampleHL-",curr.sample), 
                                   outfile.ext = NULL, 
                                   plotter.params = list(), 
                                   plot.device.name = "curr", 
                                   plotting.device.params = list(),
                                   verbose = TRUE, 
                                   debugMode = DEFAULTDEBUGMODE , 
                                   rasterize.large.plots = NULL, 
                                   rasterize.medium.plots = NULL,
                                   raster.height = 1050, 
                                   raster.width = 1050,
                                   separatePlots = FALSE,
                                   exclude.autosomes.chrom.rate.plot = TRUE,
                                   chromosome.name.style = "UCSC",
                                   fig.res = 150, fig.base.height.inches = 7,
                                   insertSize.plot.xlim = NULL,
                                   sequencing.type = c("RNA","Exome","Genome"),
                                   nvc.mark.points = TRUE,
                                   
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots = TRUE,
                                   
                                   plot = TRUE,
                                   ...){
                                                
                                                
  sequencing.type <- match.arg(sequencing.type);

    build.plotter.function <- function(){
      build.plotter.highlightSample(curr.sample,res, merge.offset.outgroup = FALSE, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC.v10(res=res, 
                   build.plotter.function=build.plotter.function, 
                   outfile = outfile,
                   outfile.dir=outfile.dir, 
                   outfile.prefix=outfile.prefix,
                   outfile.ext = outfile.ext, 
                   
                   plot.device.name=plot.device.name, 
                   plotting.device.params = plotting.device.params, 
                   
                   debugMode = debugMode, verbose = verbose, 
    
                   fig.res = fig.res, 
                   fig.base.height.inches = fig.base.height.inches,
                   rasterize.large.plots = rasterize.large.plots, 
                   rasterize.medium.plots = rasterize.medium.plots, 
                   raster.height = raster.height,
                   raster.width = raster.width,
                   raster.res = fig.res,
                   
                   separatePlots = separatePlots,
                   splitPlots = FALSE,
                   
                   nvc.highlight.points = nvc.mark.points,
                   exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                   chromosome.name.style = chromosome.name.style,
                   insertSize.plot.xlim=insertSize.plot.xlim,
                   cdf.bySample = FALSE, 
                   cdf.plotIntercepts = FALSE, 
                             
                   makePlots = plotList,
                   sequencing.type = sequencing.type,
                   
                   labelPlots = labelPlots,
                   maxColumns = maxColumns,
                   maxRows = maxRows,
                   
                   plot = plot,
                   
               ... )
}

#######################################################################
#######################################################################
#######################################################################

makeMultiPlot.highlightSample.colorByLane <- function(res, curr.sample,
                                   outfile = NULL, 
                                   outfile.dir = "./",
                                   outfile.prefix = paste0("plot-sampleHL-coloredByLane-",curr.sample), 
                                   outfile.ext = NULL, 
                                   plotter.params = list(), 
                                   plot.device.name = "curr", 
                                   plotting.device.params = list(),
                                   verbose = TRUE, 
                                   debugMode = DEFAULTDEBUGMODE , 
                                   rasterize.large.plots = NULL, 
                                   rasterize.medium.plots = NULL,
                                   raster.height = 1050, 
                                   raster.width = 1050,
                                   separatePlots = FALSE,
                                   exclude.autosomes.chrom.rate.plot = TRUE,
                                   chromosome.name.style = "UCSC",
                                   fig.res = 150, fig.base.height.inches = 7,
                                   insertSize.plot.xlim = NULL,
                                   sequencing.type = c("RNA","Exome","Genome"),
                                   nvc.mark.points = TRUE,
                                   
                                   maxColumns = NULL,
                                   maxRows = NULL,
                                   plotList = NULL,
                                   labelPlots = TRUE,
                                   
                                   plot = TRUE,
                                   ...){
                                     
  sequencing.type <- match.arg(sequencing.type);

                                                           
    build.plotter.function <- function(){
      build.plotter.highlightSample.colorByLane(curr.sample,res, merge.offset.outgroup = FALSE, plotter.params = plotter.params);
    }
    
    makeMultiPlot.GENERIC.v10(res=res, 
                   build.plotter.function=build.plotter.function, 
                   outfile = outfile,
                   outfile.dir=outfile.dir, 
                   outfile.prefix=outfile.prefix,
                   outfile.ext = outfile.ext, 
                   
                   plot.device.name=plot.device.name, 
                   plotting.device.params = plotting.device.params, 
                   
                   debugMode = debugMode, verbose = verbose, 
    
                   fig.res = fig.res, 
                   fig.base.height.inches = fig.base.height.inches,
                   rasterize.large.plots = rasterize.large.plots, 
                   rasterize.medium.plots = rasterize.medium.plots, 
                   raster.height = raster.height,
                   raster.width = raster.width,
                   raster.res = fig.res,
                   
                   separatePlots = separatePlots,
                   splitPlots = FALSE,
                   
                   nvc.highlight.points = nvc.mark.points,
                   exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                   chromosome.name.style = chromosome.name.style,
                   insertSize.plot.xlim=insertSize.plot.xlim,
                   cdf.bySample = FALSE, 
                   cdf.plotIntercepts = FALSE, 
                             
                   makePlots = plotList,
                   sequencing.type = sequencing.type,
                   
                   labelPlots = labelPlots,
                   maxColumns = maxColumns,
                   maxRows = maxRows,
                   
                   plot = plot,
                   
               ... )
}

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

supportedDevices <- c("png","CairoPNG","tiff","jpeg","CairoPDF","pdf","svg","curr");
supportedVectorDevices <- c("pdf","CairoPDF","svg");
forPrintDevices <- c("pdf","CairoPDF");

getBestLayout <- function(plotWidths, maxColumns = NULL, maxRows = NULL, nrow = NULL){
  ct <- sum(plotWidths);
  plotCt <- length(plotWidths);
  
  #message("plotWidths = ",plotWidths,", maxColumns = ",maxColumns,", maxRows = ",maxRows,", nrow = ",nrow);
  
  if(is.null(maxColumns)){
      maxColumns <- ceiling(sqrt(sum(plotWidths)));
  }
  if(is.null(maxRows)){
      if(! is.null(nrow)){
          maxRows <- nrow;
      } else {
          maxRows <- Inf;
      }
  }

      soFar <- 1;
      mv <- c();
      while(soFar <= plotCt && (is.null(mv) || nrow(mv) < maxRows)){
         v <- rep(soFar,plotWidths[soFar]);
         soFar = soFar + 1;
         while(soFar <= plotCt && length(v) + plotWidths[soFar] <= maxColumns){
             v <- c(v,rep(soFar,plotWidths[soFar]));
             soFar = soFar + 1;
         }
         if(length(v) < maxColumns){
             v <- c(v, rep(0, maxColumns - length(v) ) )
         }
         mv <- rbind(mv,v);
      }
      ht <- nrow(mv);
      wd <- ncol(mv);
      mat <- mv;
      
      if( (! is.null(nrow)) ){
         while(nrow(mat) < nrow){
             mat <- rbind(mat,rep(0,ncol(mat)));
         }
      }
  
  initLayoutFunct <- function(){
      #message("initialized layout!");
      layout(mat);
  };
  
  #message("ht = ",ht,"wd = ",wd);
  #print(mat);
  return(list(ht = ht, wd = wd,initLayoutFunct=initLayoutFunct,mat=mat));
}


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
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
     do.call(grDevices::png,device.params.final);
     return(dev.off)
   } else if(plot.device.name == "CairoPNG"){
     requireNamespace("Cairo");
     do.call(Cairo::CairoPNG,device.params.final);
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

SKIP.FOR.RNA.DATA   <- c("onTarget.rates","onTarget.counts","overlap.mismatch.byAvgQual");

SKIP.FOR.EXOME.DATA <- c(
  "genebody.coverage.umquartileExpressionGenes","genebody.coverage.lowExpressionGenes",
  "sj.locus.ct","sj.event.proportionByType","sj.event.rate",
  "norm.factors","norm.vs.TC","SpliceProfile","overlapMismatch.byQual.avg"
);

PLOTTING.FUNCTION.COMMAND.LIST <- list(
legend = list(wd=1,FUN=function(plotter,debugMode,rast,params,...){ makePlot.legend.box(plotter, debugMode = debugMode, ...) }),

qual.pair.min             =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.qual.pair(plotter,"min",           debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["ht"]],plot=plot, ...)}),
qual.pair.lowerQuartile   =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.qual.pair(plotter,"lowerQuartile", debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
qual.pair.median          =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.qual.pair(plotter,"median",        debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot,...)}),
qual.pair.upperQuartile   =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.qual.pair(plotter,"upperQuartile", debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot,...)}),
qual.pair.max             =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.qual.pair(plotter,"max",           debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot,...)}),

clippingProfile      =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.clipping(plotter,                                     debugMode = debugMode,rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
DeletionProfile      =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.cigarOp.byCycle(plotter,"Del",                        debugMode = debugMode,rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
InsertionProfile     =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.cigarOp.byCycle(plotter,"Ins",                        debugMode = debugMode,rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
SpliceProfile        =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.cigarOp.byCycle(plotter,"Splice",                     debugMode = debugMode,rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
InsertionLengthHisto =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.cigarLength.distribution(plotter,"Ins", log.y = TRUE, debugMode = debugMode,rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
DeletionLengthHisto  =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.cigarLength.distribution(plotter,"Del", log.y = TRUE, debugMode = debugMode,rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
 
gc                           =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.gc(plotter,               debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
missingness.rate             =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.missingness.rate(plotter, debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot,...)}),
dropped.rate                 =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.dropped.rates(plotter,    debugMode = debugMode,plot=plot, ...)}),
insert.size                  =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.insert.size(plotter,      debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]], xlim = params[["insertSize.plot.xlim"]],plot=plot, ...)}),
overlap.coverage             =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.overlap.coverage(plotter, debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
readLengthDist               =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.readLengthDist(plotter,   debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),

overlapMismatch.byCycle       =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.overlapMismatch.byCycle(plotter,  debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
#overlapMismatch.byQual.avg    =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.overlapMismatch.byQual.avg(plotter, debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
overlapMismatch.byQual.min    =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.overlapMismatch.byQual.min( plotter,          debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
overlapMismatch.byQual.read   =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.overlapMismatch.byQual.read( plotter,         debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
overlapMismatch.byBase        =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.overlapMismatch.byBase(plotter,   debugMode = debugMode,plot=plot, ...)}),
overlapMismatch.size          =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.overlapMismatch.size( plotter,        debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
overlapMismatch.byBase.atScore=list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.overlapMismatch.byBase.atScore( plotter,atScore=41, debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),

referenceMismatch.byCycle    =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.referenceMismatch.byCycle( plotter,                   debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
referenceMismatch.byScore    =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.referenceMismatch.byScore( plotter,                   debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
referenceMismatch.byBase     =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.referenceMismatch.byBase( plotter,             debugMode = debugMode,plot=plot, ...)}),
referenceMismatch.byBase.atScore.R1=list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.referenceMismatch.byBase.atScore( plotter,atScore=41,forRead="R1", debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
referenceMismatch.byBase.atScore.R2=list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.referenceMismatch.byBase.atScore( plotter,atScore=41,forRead="R2", debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),

gene.diversity           =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.gene.cdf(plotter, sampleWise = params[["cdf.bySample"]], plot.intercepts = params[["cdf.plotIntercepts"]], 
                                    debugMode = debugMode, rasterize.plotting.area = rast[["big"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),

genebody.coverage.allGenes                  =list(wd=1,FUN=function(plotter,plot,debugMode,rast,params,...){makePlot.genebody(plotter, geneset="Overall",             debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
genebody.coverage.umquartileExpressionGenes =list(wd=1,FUN=function(plotter,plot,debugMode,rast,params,...){makePlot.genebody(plotter, geneset="50-75",               debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
genebody.coverage.lowExpressionGenes        =list(wd=1,FUN=function(plotter,plot,debugMode,rast,params,...){makePlot.genebody(plotter, geneset="0-50",                debugMode = debugMode, rasterize.plotting.area = rast[["med"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
geneAssignmentRates                         =list(wd=1,FUN=function(plotter,plot,debugMode,rast,params,...){makePlot.gene.assignment.rates(plotter,                   debugMode = debugMode,plot=plot, ...)}),
sj.locus.ct                                 =list(wd=1,FUN=function(plotter,plot,debugMode,rast,params,...){makePlot.splice.junction.loci.counts(plotter,             debugMode = debugMode,plot=plot, ...)}),
sj.event.proportionByType                   =list(wd=1,FUN=function(plotter,plot,debugMode,rast,params,...){makePlot.splice.junction.event.proportionsByType(plotter, debugMode = debugMode,plot=plot, ...)}),
sj.event.rate                               =list(wd=1,FUN=function(plotter,plot,debugMode,rast,params,...){makePlot.splice.junction.event.ratesPerRead(plotter,      debugMode = debugMode,plot=plot, ...)}),

mapping.rates       =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.mapping.rates(plotter,      debugMode = debugMode,plot=plot, ...)}),
chrom.rates         =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.chrom.type.rates(plotter, chromosome.name.style = params[["chromosome.name.style"]], exclude.autosomes = params[["exclude.autosomes.chrom.rate.plot"]], debugMode = debugMode,plot=plot, ...)}),
norm.factors        =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.norm.factors(plotter,       debugMode = debugMode,plot=plot, ...)}),
norm.vs.TC          =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.norm.factors.vs.TC(plotter, debugMode = debugMode,plot=plot, ...)}),
strandedness.test   =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.strandedness.test(plotter,  debugMode = debugMode,plot=plot, ...)}),
onTarget.counts     =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.onTarget.counts(plotter,    debugMode = debugMode,plot=plot, ...)}),
onTarget.rates      =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.onTarget.rates(plotter,     debugMode = debugMode,plot=plot, ...)}),
NVC.lead.clip       =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.NVC.lead.clip(plotter, clip.amt = params[["clip.amt"]], points.highlighted = params[["nvc.highlight.points"]], 
                                    debugMode = debugMode, rasterize.plotting.area = rast[["big"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),
NVC.tail.clip       =list(wd=1,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.NVC.tail.clip(plotter, clip.amt = 12, points.highlighted = params[["nvc.highlight.points"]], 
                                     debugMode = debugMode, rasterize.plotting.area = rast[["big"]], raster.height = rast[["ht"]], raster.width = rast[["wd"]],plot=plot, ...)}),

NVC.raw  =list(wd=2,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.raw.NVC(plotter, points.highlighted = params[["nvc.highlight.points"]], 
                                     debugMode = debugMode, rasterize.plotting.area = rast[["big"]], raster.height = rast[["ht"]], raster.width = 2* rast[["wd"]],plot=plot, ...)}),
NVC.aligned  =list(wd=2,FUN=function(plotter,debugMode,rast,params,plot,...){makePlot.minus.clipping.NVC(plotter, points.highlighted = params[["nvc.highlight.points"]], 
                                     debugMode = debugMode, rasterize.plotting.area = rast[["big"]], raster.height = rast[["ht"]], raster.width = 2* rast[["wd"]],plot=plot, ...)})

);



makeMultiPlot.GENERIC.v10 <- function(res, 
               build.plotter.function, 
               outfile = NULL,
               outfile.dir, 
               outfile.prefix,
               outfile.ext = NULL, 
               
               plot.device.name, 
               plotting.device.params = list(), 
               
               debugMode = DEFAULTDEBUGMODE, verbose = TRUE, 

               fig.res = 150, 
               fig.base.height.inches = 7,
               rasterize.large.plots = NULL, 
               rasterize.medium.plots = NULL, 
               raster.height = fig.res*fig.base.height.inches,
               raster.width = raster.height,
               raster.res = fig.res,
               
               separatePlots = FALSE,
               splitPlots = FALSE,
               
               nvc.highlight.points = TRUE,
               exclude.autosomes.chrom.rate.plot = TRUE,
               chromosome.name.style = "UCSC",
               insertSize.plot.xlim=NULL,
               cdf.bySample = TRUE,
               cdf.plotIntercepts = FALSE,
               clip.amt = 12,
               
               makePlots = NULL,
               skipPlots = NULL,
               sequencing.type = c("RNA","Exome","Genome"),
               skipMissingDataPlots = TRUE,
               
               labelPlots = TRUE,
               maxColumns = NULL,
               maxRows = NULL,
               
               plot = TRUE,
               ... ){
  
  sequencing.type <- match.arg(sequencing.type);
  #isExome <- sequencing.type == "Exome";
  
  plotter <- build.plotter.function();
  
  if(is.null(skipPlots)){
    if(sequencing.type == "Exome"){
      skipPlots <- SKIP.FOR.EXOME.DATA;
    } else if(sequencing.type == "RNA"){
      skipPlots <- SKIP.FOR.RNA.DATA;
    } else {
      #Not yet supported!
    }
    makePlots <- names(PLOTTING.FUNCTION.COMMAND.LIST)[! names(PLOTTING.FUNCTION.COMMAND.LIST) %in% skipPlots]
    
    if(debugMode) message("Skipping: \"",paste0(skipPlots,collapse="\",\""),"\"");
  }
  
  height.per.px <- fig.res * fig.base.height.inches;
  width.per.px <- fig.res * fig.base.height.inches;
  height.per.inches <- fig.base.height.inches;
  width.per.inches <- fig.base.height.inches;
  
  if(is.null(rasterize.large.plots)){
    if(plot.device.name %in% supportedVectorDevices){
      if(check.rasterize.or.warn("rasterize.large.plots")){
        message("Default: rasterizing large plots")
        rasterize.large.plots = TRUE;
      } else {
        rasterize.large.plots = FALSE;
      }
    } else {
      rasterize.large.plots = FALSE;
    }
  }
  if(is.null(rasterize.medium.plots)){
    if(plot.device.name %in% forPrintDevices){
      if(check.rasterize.or.warn("rasterize.medium.plots")){
        message("Default: rasterizing medium plots")
        rasterize.medium.plots = TRUE;
      } else {
        rasterize.medium.plots = FALSE;
      }
    } else {
      rasterize.medium.plots = FALSE;
    }
  }
  
  if(rasterize.large.plots || rasterize.medium.plots){
    check.rasterize.or.die("rasterize.large.plots");
    
    if((! plot.device.name %in% supportedVectorDevices) & (plot.device.name != "curr")){
      warning("rasterize.large.plots = TRUE should not be used with raster file formats (ie png, tiff, jpeg, etc). This will result in image degradation.");
    }
  }
  if(debugMode) message("Rasterize large plots: ", rasterize.large.plots);
  if(debugMode) message("Rasterize medium plots: ", rasterize.medium.plots);

  if(skipMissingDataPlots){
      runParams <- list(   nvc.highlight.points=nvc.highlight.points,
                    exclude.autosomes.chrom.rate.plot=exclude.autosomes.chrom.rate.plot,
                    rasterize.large.plots=rasterize.large.plots,
                    rasterize.medium.plots=rasterize.medium.plots,
                    chromosome.name.style=chromosome.name.style,
                    insertSize.plot.xlim=insertSize.plot.xlim,
                    cdf.plotIntercepts=cdf.plotIntercepts,
                    cdf.bySample=cdf.bySample,
                    clip.amt=clip.amt);
      rast <- list(big = rasterize.large.plots, med = rasterize.medium.plots,ht = raster.height, wd = raster.width, res = raster.res);
      
      dataMissingPlots <- sapply(makePlots,function(p){
        #message("p:",p);
        ! PLOTTING.FUNCTION.COMMAND.LIST[[p]]$FUN(plotter=plotter,debugMode=debugMode,rast=rast,params=runParams,plot=FALSE)
      });
      if(debugMode) message("Skipping due to missing data: \"",paste0(makePlots[dataMissingPlots],collapse="\",\""),"\"");

      makePlots <- makePlots[! dataMissingPlots]
      
  }

  if(plot.device.name == "curr"){
    if(debugMode) message("Plotting to the currently-open device...");
    default.params <- list();
    dev.params <- list();
    devOpenFunct <- function(f,w,height.mult,width.mult){};
    devCloseFunct <- function(){};
  } else if(plot.device.name == "png"){
    if(is.null(outfile.ext)) outfile.ext = ".png";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w,height.mult,width.mult){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(png,dev.params)
    };
    devCloseFunct <- function(){dev.off()};
  } else if(plot.device.name == "tiff"){
    if(is.null(outfile.ext)) outfile.ext = ".tiff";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w,height.mult,width.mult){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(tiff,dev.params)
    };
    devCloseFunct <- function(){dev.off()}    
  } else if(plot.device.name == "jpeg"){
    if(is.null(outfile.ext)) outfile.ext = ".jpg";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w,height.mult,width.mult){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(jpeg,dev.params)
    };
    devCloseFunct <- function(){dev.off()}        
  } else if(plot.device.name == "tiff"){
    if(is.null(outfile.ext)) outfile.ext = ".tiff";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w,height.mult,width.mult){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(tiff,dev.params)
    };
    devCloseFunct <- function(){dev.off()}        
  } else if(plot.device.name == "CairoPNG"){
    if(! requireNamespace("Cairo", quietly=TRUE)) stop("Error: package Cairo not found. Install package Cairo or set plot.device.name to something other than CairoPNG or CairoPDF.");

    if(is.null(outfile.ext)) outfile.ext = ".png";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w,height.mult,width.mult){
      default.params <- list(filename = f, height = height.per.px * height.mult, width = width.per.px * width.mult * w, pointsize = 36, units = "px");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      requireNamespace("Cairo", quietly=TRUE)
      do.call(Cairo::CairoPNG,dev.params)
    };
    devCloseFunct <- function(){dev.off()}
  } else if(plot.device.name == "svg"){
    if(is.null(outfile.ext)) outfile.ext = ".svg";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w,height.mult,width.mult){
      default.params <- list(filename = f, height = height.per.inches * height.mult, width = width.per.inches * width.mult * w, pointsize = 24);
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(svg,dev.params)
    };
    devCloseFunct <- function(){dev.off()}
  } else if(plot.device.name == "pdf"){
    if(is.null(outfile.ext)) outfile.ext = ".pdf";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w,height.mult,width.mult){
      default.params <- list(file = f, height = 0, width = 0, pointsize = 10, paper = "letter");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      do.call(pdf,dev.params)
    };
    devCloseFunct <- function(){dev.off()}
  } else if(plot.device.name == "CairoPDF"){
    if(! requireNamespace("Cairo",quietly=TRUE)) stop("Error: package Cairo not found. Install package Cairo or set plot.device.name to something other than CairoPNG or CairoPDF.");
    if(is.null(outfile.ext)) outfile.ext = ".pdf";
    if(is.null(outfile)) outfile <- paste0(outfile.dir, outfile.prefix,outfile.ext);
    
    devOpenFunct <- function(f,w,height.mult,width.mult){
      default.params <- list(file = f, height = 0, width = 0, pointsize = 10, paper = "letter");
      dev.params <- overmerge.list(default.params,plotting.device.params);
      requireNamespace("Cairo", quietly=TRUE)
      do.call(Cairo::CairoPDF,dev.params)
    };
    devCloseFunct <- function(){dev.off()}
  } else {
    stop(paste0("Error: graphics device \"",plot.device.name,"\" not supported! Legal options are: [",paste0(supportedDevices,collapse=","),"] Set plot.device.name to \"curr\" to plot to the currently-open and/or default device."));
  }
  
  #Default arrangement for multipage PDF files:
  multiPage <- (plot.device.name == "CairoPDF" || plot.device.name == "pdf");
  
  if(is.null(maxColumns) && multiPage){
      maxColumns <- 2;
  }
  if(is.null(maxRows) && multiPage){
      maxRows <- 3;
  }
  if(multiPage){
      nrow = maxRows;
  } else {
      nrow = NULL
  }
  
  tryCatch({
      #plotter <- build.plotter.function();

      if(debugMode) message("Plotting extended...");
      INTERNAL.plot.v10(               res = res, plotter = plotter, 
                                       
                                       devOpenFunct = devOpenFunct,
                                       devCloseFunct = devCloseFunct,
                                       
                                       outfile = outfile,
                                       outfile.dir=outfile.dir, 
                                       outfile.prefix=outfile.prefix,
                                       outfile.ext = outfile.ext,
                                       
                                       skipPlots = skipPlots,
                                       makePlots = makePlots,
                                       
                                       verbose = verbose, debugMode = debugMode, 
                                       
                                       separatePlots = separatePlots,
                                       labelPlots = labelPlots,
                                       maxColumns = maxColumns,
                                       maxRows = maxRows,
                                       nrow=nrow,
                                       multiPage = multiPage,

                                       cdf.bySample = cdf.bySample,
                                       cdf.plotIntercepts = cdf.plotIntercepts,
                                       nvc.highlight.points = nvc.highlight.points, 
                                       exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                                       chromosome.name.style = chromosome.name.style,
                                       insertSize.plot.xlim=insertSize.plot.xlim,
                                       clip.amt = clip.amt,

                                       rasterize.large.plots = rasterize.large.plots, 
                                       rasterize.medium.plots = rasterize.medium.plots, 
                                       raster.height = raster.height, 
                                       raster.width = raster.width, 
                                       raster.res = raster.res,
                                       
                                       plot=plot,
                                       
                                       ...);
  }, error = function(e){
    message(paste0("PLOTTING ERROR: ", e));
    devCloseFunct();
  }, finally = {
    #message("CRASHED!");
    if(debugMode)  message("Done with plot.");
  });
}

INTERNAL.plot.v10 <- function( res, 
                               plotter,
                               devOpenFunct,
                               devCloseFunct,
                               
                               outfile = NULL,
                               outfile.dir, 
                               outfile.prefix,
                               outfile.ext = NULL, 
                               
                               verbose = TRUE,
                               debugMode,
                               
                               separatePlots = FALSE,
                               labelPlots = TRUE,
                               maxColumns = NULL,
                               maxRows = NULL,
                               nrow = NULL,
                               
                               multiPage = FALSE,
                               
                               makePlots = names(PLOTTING.FUNCTION.COMMAND.LIST),
                               skipPlots = c(),
                               skipMissingDataPlots = FALSE,
                               
                               cex.corner.label = 2, 
                           #makeplot params:
                               cdf.bySample = FALSE, 
                               nvc.highlight.points = TRUE,
                               cdf.plotIntercepts = TRUE,  
                               exclude.autosomes.chrom.rate.plot = TRUE,
                               chromosome.name.style = "UCSC",
                               insertSize.plot.xlim=NULL,
                               clip.amt = 12,
                           #Rasterization params:
                               rasterize.large.plots = FALSE, 
                               rasterize.medium.plots = FALSE, 
                               raster.height = 1050,
                               raster.width  = 1050,
                               raster.res = 150,
                               
                               plot = TRUE,
                               
                           #graphics par:
                               ...){
  
  
  rast <- list(big = rasterize.large.plots, med = rasterize.medium.plots,ht = raster.height, wd = raster.width, res = raster.res);
  params <- list(   nvc.highlight.points=nvc.highlight.points,
                    exclude.autosomes.chrom.rate.plot=exclude.autosomes.chrom.rate.plot,
                    rasterize.large.plots=rasterize.large.plots,
                    rasterize.medium.plots=rasterize.medium.plots,
                    chromosome.name.style=chromosome.name.style,
                    insertSize.plot.xlim=insertSize.plot.xlim,
                    cdf.plotIntercepts=cdf.plotIntercepts,
                    cdf.bySample=cdf.bySample,
                    clip.amt=clip.amt);
  
  if(verbose) {message(paste0("Starting compiled plot..."));}
  
  ts <- timestamp();
     
  a.to.z <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
  corner.labels <- c(a.to.z, paste0(a.to.z[1],a.to.z), paste0(a.to.z[2],a.to.z), paste0(a.to.z[3],a.to.z));
  
  plot.corner.label <- function(i){
    if(! labelPlots){
      #do nothing!
    } else {
      devlim <- device.limits();
      text(devlim[1],devlim[4], corner.labels[i] , cex = cex.corner.label, adj = c(-0.1,1.1),  xpd=T);
    }
  }
  
  plotList <- PLOTTING.FUNCTION.COMMAND.LIST[makePlots];
  plotList <- plotList[! names(plotList) %in% skipPlots];
  
  if(skipMissingDataPlots){
    plotList <- plotList[ sapply(plotList,function(pf){
          pf$FUN(plotter=plotter,debugMode=debugMode,rast=rast,params=params,plot=FALSE)
        })];
  }
  plotWidths <- sapply(plotList,FUN=function(pf){ pf$wd })

  devOpenFunctSimple <- devOpenFunct;
  devOpenFunct <- function(...){
      #message("opening funct!");
      devOpenFunctSimple(...);
  }
  devCloseFunctSimple <- devCloseFunct;
  devCloseFunct <- function(...){
      #message("Closing funct!");
      devCloseFunctSimple(...);
  }
  
  
  if(separatePlots){
    #do nothing!
  } else if(is.null(maxRows)){
    blo <- getBestLayout(plotWidths,maxColumns=maxColumns);
    devOpenFunct(outfile,1,blo$ht,blo$wd);
    blo$initLayoutFunct();
  } else {
    blo <- getBestLayout(plotWidths,maxColumns=maxColumns,maxRows=maxRows,nrow=nrow);
    layoutCap <- max(blo$mat);
    layoutStart <- 1;
    if(multiPage){
        devOpenFunct(outfile,1,blo$ht,blo$wd);
    } else {
        devOpenFunct(paste0(outfile.dir,"/",outfile.prefix,".p",1,outfile.ext),1,blo$ht,blo$wd);
    }
    blo$initLayoutFunct();
  }
  
  plotCt = 0;
  for(i in seq_along(plotList)){
    plotFUN <- plotList[[i]]$FUN;
    plotName <- names(plotList)[[i]];
    plotWd <- plotWidths[[i]];
    
    if(separatePlots){
      devOpenFunct(paste0(outfile.dir,"/",outfile.prefix,".",plotName,outfile.ext),plotWd,1,1);
    } else if(is.null(maxRows)){
      #do nothing
    } else {
      if(layoutStart + layoutCap - 1 < i){
        blo <- getBestLayout(plotWidths[i:length(plotWidths)],maxColumns=maxColumns,maxRows=maxRows,nrow=nrow);
        
        layoutCap <- max(blo$mat);
        layoutStart <- i;
        
        message("  layoutCap = ",layoutCap);
        message("  layoutStart = ",layoutStart);
        message("  i = ",i);
        
        if(multiPage){  
            #devOpenFunct(outfile,1,blo$ht,blo$wd);
            #do nothing!
        } else       {
            devOpenFunct(paste0(outfile.dir,"/",outfile.prefix,".p",1,outfile.ext),1,blo$ht,blo$wd);
        }
        blo$initLayoutFunct();
      } else {
        #do nothing
      }
    }
    
    plotCt=plotCt+1;
    if(plot){
      plotFUN(plotter=plotter,debugMode=debugMode,rast=rast,params=params,plot=TRUE,...)
    } else {
      plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1));
      text(0.5,0.5,label=plotName);
    }
    plot.corner.label(i);
    
    
    if(separatePlots){
      devCloseFunct();
    }
  }
  if(! separatePlots){
      devCloseFunct();
  }
  
  if(debugMode) message("Finished Multiplot",getTimeAndDiff(ts))
}

