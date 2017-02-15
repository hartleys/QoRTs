
makePlot.outlierSlices <- function(
                         res,
                         outfile.dir = "./",
                         plotter.params = list(), 
                         plot.device.name = "png", 
                         plotting.device.params = list(), 
                         num.groups = length(use.colors),
                         use.colors = c("red","blue","green","orange","magenta","cyan","yellow3","deepskyblue4"),
                         color.titles = c("red","blue","green","orange","magenta","cyan","yellow","deepBlue"),
                         group.codes = c("A","B","C","D","E","F","G","H"),
                         debugMode = DEFAULTDEBUGMODE , 
                         rasterize.large.plots = NULL, 
                         rasterize.medium.plots = NULL,
                         raster.height = 1000, 
                         raster.width = 1000,
                         exclude.autosomes.chrom.rate.plot = TRUE,
                         chromosome.name.style = "UCSC",
                         fig.res = 150, fig.base.height.inches = 7,
                         insertSize.plot.xlim = NULL,
                         sequencing.type = c("RNA","Exome","Genome"),
                         generateHtmlIndex = TRUE,
                         generateSamplePlots = c("multiRun","all","none"),
                         maxColumns = 4,
                         verbose = TRUE,
                         plot = TRUE,
                         ...){
  
  if(num.groups > length(use.colors))   stop("Fatal Error: makePlot.outlierSlices: num.groups > length(use.colors)");
  if(num.groups > length(color.titles)) stop("Fatal Error: makePlot.outlierSlices: num.groups > length(color.titles)");
  if(num.groups > length(group.codes))  stop("Fatal Error: makePlot.outlierSlices: num.groups > length(group.codes)");

  use.colors <- use.colors[1:num.groups];
  color.titles <- use.colors[1:color.titles];
  group.codes <- use.colors[1:group.codes];
  
  sequencing.type <- match.arg(sequencing.type);
  generateSamplePlots <- match.arg(generateSamplePlots);
  
  splitPlotDir <- paste0(outfile.dir,"splitPlots/")
  htmlDir <- paste0(outfile.dir,"html/")
  htmlSplitDir <- paste0(outfile.dir,"html/split//")
  htmlSampPlotDir <- paste0(outfile.dir,"html/sampPlots//")
  htmlSampDir <- paste0(outfile.dir,"html/samps/")
  
  if(!dir.exists(outfile.dir)) dir.create(outfile.dir);
  if(!dir.exists(splitPlotDir)) dir.create(splitPlotDir);
  if(!dir.exists(htmlDir)) dir.create(htmlDir);
  if(!dir.exists(htmlSplitDir)) dir.create(htmlSplitDir);

  has.sample.plot <- c(); 
  if(generateSamplePlots == "none"){
     has.sample.plot <- c(); 
  } else if(generateSamplePlots == "all"){
     dir.create(htmlDir);
     dir.create(htmlSampPlotDir);
     has.sample.plot <- res@sample.list 
  } else {
     dir.create(htmlSampDir);
     dir.create(htmlSampPlotDir);
     has.sample.plot <- res@sample.list[sapply(res@sample.list,function(samp){
              sum(res@decoder$sample.ID == samp) > 1
     })]
  }
   
  message("STARTING SLICE PLOTS:");
   
  makePlot.outlierSlices.recursive(
                         super.group = NULL,
                         outfile.prefix = "set",
                         res=res,
                         outfile.dir = outfile.dir,
                         plotter.params = plotter.params, 
                         plot.device.name = plot.device.name, 
                         plotting.device.params = plotting.device.params, 
                         use.colors = use.colors,
                         color.titles = color.titles ,
                         group.codes = group.codes,
                         num.groups = num.groups,
                         debugMode = debugMode , 
                         rasterize.large.plots = rasterize.large.plots, 
                         rasterize.medium.plots = rasterize.medium.plots,
                         raster.height = raster.height, 
                         raster.width = raster.width,
                         exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                         chromosome.name.style = chromosome.name.style,
                         fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                         insertSize.plot.xlim = insertSize.plot.xlim,
                         sequencing.type = sequencing.type,
                         generateHtmlIndex = generateHtmlIndex,
                         maxColumns = maxColumns,
                         has.sample.plot = has.sample.plot,
                         plot = plot,
                         ...);
  
  message("FINISHED SLICE PLOTS:");
  
  if(length(has.sample.plot) > 0){
      message("STARTING SAMPLE PLOTS:");
      
      for(i in 1:length(has.sample.plot)){
        curr.sample <- has.sample.plot[[i]];
        cat(".");
        colorBy <- res@decoder$lane.ID;
        names(colorBy) <- res@decoder$unique.ID;
        
        plotter <- build.plotter.colorAndHighlightByVector(res, 
                     colorBy = colorBy,
                     color.by.title.name = "Other",
                     highlightBy = res@decoder$sample.ID, 
                     highlight = curr.sample, 
                     highlight.by.name = curr.sample,
                     highlight.by.title.name = "SAMPLE",
                     plotter.params = plotter.params);
        
        makeMultiPlot.withPlotter(plotter=plotter,
                                   outfile.dir = paste0(outfile.dir,"sampPlots/"),
                                   outfile.prefix = paste0(curr.sample), 
                                   plot.device.name = plot.device.name, 
                                   plotting.device.params = plotting.device.params,
                                   verbose = TRUE, 
                                   debugMode = FALSE , 
                                   rasterize.large.plots = rasterize.large.plots, 
                                   rasterize.medium.plots = rasterize.medium.plots,
                                   raster.height = raster.height , 
                                   raster.width = raster.width,
                                   exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                                   chromosome.name.style = chromosome.name.style,
                                   fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                                   insertSize.plot.xlim = insertSize.plot.xlim,
                                   sequencing.type = sequencing.type,
                                   nvc.mark.points = TRUE,
                                   maxColumns=maxColumns,
                                   plot = plot,...
        );
      }
      
      message("FINISHED SAMPLE PLOTS:");
  }
}

makePlot.outlierSlices.recursive <- function(
                         super.group = NULL,
                         outfile.prefix = "set",
                         res,
                         outfile.dir = "./",
                         plotter.params = list(), 
                         plot.device.name = "png", 
                         plotting.device.params = list(), 
                         use.colors = c("red","blue","green","orange","magenta","cyan","yellow3","deepskyblue4"),
                         color.titles = c("red","blue","green","orange","magenta","cyan","yellow","deepBlue"),
                         group.codes = c("A","B","C","D","E","F","G","H"),
                         num.groups = 8,
                         debugMode = DEFAULTDEBUGMODE , 
                         rasterize.large.plots = NULL, 
                         rasterize.medium.plots = NULL,
                         raster.height = 1000, 
                         raster.width = 1000,
                         exclude.autosomes.chrom.rate.plot = TRUE,
                         chromosome.name.style = "UCSC",
                         fig.res = 150, fig.base.height.inches = 7,
                         insertSize.plot.xlim = NULL,
                         sequencing.type = c("RNA","Exome","Genome"),
                         generateHtmlIndex = TRUE,
                         maxColumns = 4,
                         has.sample.plot = c(),
                         plot = TRUE,
                         ...){

  sequencing.type <- match.arg(sequencing.type);
  
  plot.dir <- paste0(outfile.dir,"/splitPlots/");
  html.dir <- paste0(outfile.dir,"/html/split/");
  
  #use.colors <- c("red","blue","green","orange","magenta","cyan","yellow3","deepskyblue4")
  #color.titles <- c("red","blue","green","orange","magenta","cyan","yellow","deepBlue");
  #group.codes <- c("A","B","C","D","E","F","G","H");
  #use.colors <- c("red","blue","green")
  #color.titles <- c("red","blue","green");
  #group.codes <- c("A","B","C");
  
  use.colors <- use.colors[1:num.groups];
  color.titles <- color.titles[1:num.groups];
  group.codes <- group.codes[1:num.groups];
  
  plotter.params <- list(
    contrasting.colors = use.colors,
    contrasting.pch    = strtoi(charToRaw(paste0(group.codes,collapse="")),16L)
  );
  group.titles <- paste0(group.codes,":",color.titles)
  outgroup.name <- "Z:other";
  
  if(debugMode) message("Starting makePlot.outlierSlices.recursive");

  if(is.null(super.group)){
      super.group <- res@sample.list;
  }
  inCurr <- which(res@sample.list %in% super.group);
  sample.group <- rep(outgroup.name,length(res@sample.list));
  sample.group[inCurr] <- rep(group.titles,length.out=length(inCurr));
  names(sample.group) <- res@sample.list;
  groupSize <- sum(sample.group == group.titles[1]);
  
  groupSizes <- sapply(1:length(group.titles),function(i){
    sum(sample.group == group.titles[i])
  })
  groupSizes <- range(groupSizes);
  if(groupSizes[1] == groupSizes[2]){
    groupRange <- paste0(groupSizes[1]);
  } else {
    groupRange <- paste0(groupSizes[1],"-",groupSizes[2]);
  }
  
  #message("debug AC")
  #message("use.colors = ",paste0(use.colors,collapse=","));
  
        if(generateHtmlIndex){
          generate.split.html.page(  html.dir=html.dir,
                                     plot.dir="../../split",
                                     sample.group=sample.group,
                                     outfile.prefix=paste0(outfile.prefix),
                                     use.colors=use.colors,
                                     color.titles=color.titles,
                                     group.codes=group.codes,
                                     group.titles=group.titles,
                                     outgroup.name=outgroup.name,
                                     has.sample.plot = has.sample.plot)
        }
  
  
  if(debugMode) message("(",outfile.prefix,")",": ",length(inCurr)," replicates split into ~",length(group.titles)," groups of ~",groupSize);
      #  plotter <- build.plotter.colorAndHighlightByVector(res,
      #             colorBy = sample.group,
      #             highlightList = res@sample.list[inCurr],
      #             color.by.title.name = paste0("Other"), 
      #             plotter.params = plotter.params)
      
      colorBy <- sample.group[res@decoder$sample.ID];
      highlightBy <- ifelse(colorBy == outgroup.name,"Other","Curr");
      
      plotter <- build.plotter.colorAndHighlightByVector(
                       res=res, 
                       colorBy = colorBy, 
                       highlightBy = highlightBy,
                       highlight = "Curr",
                       color.by.title.name = paste0(sum(groupSizes > 0)," groups of ",groupRange),
                       highlight.by.title.name = paste0(sum(groupSizes)," Samples"),
                       highlight.by.name = "Groups",
                       plotter.params = plotter.params);
                   
  if(debugMode) message("plotter generated");
        makeMultiPlot.withPlotter(plotter=plotter,
                                   outfile.dir = plot.dir,
                                   outfile.prefix = outfile.prefix, 
                                   plot.device.name = plot.device.name, 
                                   plotting.device.params = plotting.device.params,
                                   verbose = TRUE, 
                                   debugMode = FALSE , 
                                   rasterize.large.plots = rasterize.large.plots, 
                                   rasterize.medium.plots = rasterize.medium.plots,
                                   raster.height = raster.height , 
                                   raster.width = raster.width,
                                   exclude.autosomes.chrom.rate.plot = exclude.autosomes.chrom.rate.plot,
                                   chromosome.name.style = chromosome.name.style,
                                   fig.res = fig.res, fig.base.height.inches = fig.base.height.inches,
                                   insertSize.plot.xlim = insertSize.plot.xlim,
                                   sequencing.type = sequencing.type,
                                   nvc.mark.points = TRUE,
                                   maxColumns=maxColumns,
                                   plot = plot ,
                                   ...
        );
  if(debugMode) message("plot complete.");

    for(i in seq_along(group.titles)){
        numSamplesInGroup <- sum(sample.group == group.titles[i]);
        if(numSamplesInGroup > 1){
            makePlot.outlierSlices.recursive(super.group = res@sample.list[sample.group == group.titles[i]],
                                             outfile.prefix = paste0(outfile.prefix,".",group.codes[i]),
                                             res=res,
                                             outfile.dir=outfile.dir,
                                             plotter.params=plotter.params,
                                             plot.device.name=plot.device.name,debugMode=debugMode,
                                             rasterize.large.plots=rasterize.large.plots,rasterize.medium.plots=rasterize.medium.plots,
                                             raster.height=raster.height,raster.width=raster.width,
                                             exclude.autosomes.chrom.rate.plot=exclude.autosomes.chrom.rate.plot,
                                             chromosome.name.style=chromosome.name.style,
                                             fig.res=fig.res,fig.base.height.inches=fig.base.height.inches,insertSize.plot.xlim=insertSize.plot.xlim,
                                             sequencing.type=sequencing.type,
                                             maxColumns=maxColumns,
                                             plot = plot,
                                             ...);
        }
    }
}

generate.split.html.page <- function(html.dir,
                                     plot.dir,
                                     sample.group,
                                     has.sample.plot = c(),
                                     outfile.prefix,
                                     use.colors,
                                     color.titles,
                                     group.codes,
                                     group.titles,
                                     outgroup.name){
    #message("writing html...")
    #message("use.colors = ",paste0(use.colors,collapse=","));
    
    groups.used <- unique(sample.group);
    groups.used <- groups.used[groups.used != outgroup.name];
    file = paste0(html.dir,"/",outfile.prefix,".html");
    hex.colors <- color2hex(use.colors);
    
    message("writing html: ",file)
    
    #Start HTML:
    cat(paste0("<html><head>\n",
               "<title>", outfile.prefix, " </title>\n",
               "<link rel=\"stylesheet\" type=\"text/css\" href=\"styles.css\">",
               "</head>\n",
               "<body bgcolor=#D0D0D0>\n"
               ),
        file = file, append = FALSE);
    
    #start up table:
    cat(paste0(
          "<div class=\"linkTableContainer\">\n",
               "  <table>\n",
               "    <tr>"
        ),file = file, append = TRUE);
    
    for(i in seq_along(group.codes)){
        if(sum(sample.group == group.titles[i]) > 1){
          cat(paste0(
               "      <td>\n",
               "        <a href=\"","./",outfile.prefix,".",group.codes[i],".html","\" style=\"color:",hex.colors[i],"\">",outfile.prefix,".",group.codes[i],"</a>\n",
               "      </td>\n"
          ),file = file, append = TRUE);
        } else if(sum(sample.group == group.titles[i]) == 1){
          sample.ID <- names(sample.group)[sample.group == group.titles[i]];
          if(sample.ID %in% has.sample.plot){
            cat(paste0(
               "      <td>\n",
               "        <a href=\"","../samp/",sample.ID,".html","\" style=\"color:",hex.colors[i],"\">",sample.ID,"</a>\n",
               "      </td>\n"
            ),file = file, append = TRUE);
          } else {
            cat(paste0(
               "      <td>\n",
               "        <div style=\"color:",hex.colors[i],"\">",sample.ID,"</div>\n",
               "      </td>\n"
            ),file = file, append = TRUE);
          }
        }
    }
    
    #Finish up table:
    cat(paste0("    </tr>\n",
               "  </table>\n",
               "</div>\n"
        ),file = file, append = TRUE);

    #Add image:
    
    cat(paste0(
          "<div class=\"plotContainer\">\n",
          "<img src=\"../../splitPlots/",outfile.prefix,".png\"","width=1850>\n",
          "</div>"
        ),file = file, append = TRUE);
        
    #Finish up:
    cat(paste0("</body></html>"
        ),file = file, append = TRUE);

}


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

blank.plot <- function(label, newPlot = TRUE, cex = NULL, ...){
  if(newPlot){
    plot.new();
    plot.window(xlim=c(0,1),ylim=c(0,1));
  } else {
    plot.window(xlim=c(0,1),ylim=c(0,1));
    
  }
  #plot(0,0,col="transparent", xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="");
  label.text <- paste0(label,collapse="\n");
  if(is.null(cex)){
    cex <- fit.character.vector(label.text)
  }
  text(0.5,0.5, labels = label.text, cex = cex, ...);
  box();
}

plotter.error.wrapper <- function(plot.name, plotterFcn, ...){
   tryCatch({
     plotterFcn();
   }, error = function(e){
     message("    (Skipped ",plot.name," due to encountering errors)");
     message("     Error:\"",e,"\"");
   });
}

errorPlot <- function(plot.name, e, code = "UNKNOWN", newPlot = TRUE, cex = NULL, ...){
     message("WARNING: Error encountered (Code ",code,") while attempting ",plot.name,"\n",
             "    (Error Text: \"", removeLineBreaks(e),"\"");
     blank.plot(c("Error encountered attempting ",plot.name,"\nSkipping..."), newPlot = newPlot, ...);
}

########################################################################################################################

plot.miniticks <- function(ax = 2,log.y = TRUE, lwd = -1, lwd.ticks = 0.5, tcl = -0.25,...){
  
  if(log.y){
    miniticks <- log10(getAllOrdersOfTen(10 ^ par("usr")[3],10^par("usr")[4], skipBase=FALSE));
    axis(ax,at=miniticks,labels=FALSE,lwd = lwd, lwd.ticks = lwd.ticks,tcl=tcl, ...);
  } else {
    #do nothing!
  }
}

#makePlot.generic.pair <- function(plot.name, data.list.r1, data.list.r2, 
#                                                  plotter, x.name, y.name, norm.x, avg.y, xlim , ylim = NULL, vert.offset = 0, horiz.offset = 0, draw.horiz.lines = FALSE, plot.type = "lines", r2.buffer = 10, override.lty = -1, x.is.log = FALSE, y.is.log = FALSE, pre.plot.func = NULL, y.axis.las = 1, ...){


makePlot.generic <- function(plot.name, data.list, plotter, x.name, y.name, norm.x, avg.y, plot.type = "lines", xlim = NULL, ylim = NULL, 
                              plot.means = FALSE, plot.medians = FALSE, vert.offset = 0, pre.plot.func = NULL, override.lty = -1, 
                              draw.horiz.lines = FALSE, x.is.log = FALSE, y.is.log = FALSE, y.axis.las = 1, zeroBaseX = FALSE, 
                              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                              debugMode = FALSE, randomizePlotOrder = FALSE,
                              xaxis.pretty.n = 5,
                              ...){
  tryCatch({
    priorities <- sort(unique(plotter$lanebam.params$plot.priority));
    
    
    
    if(rasterize.plotting.area){
       rasterize.plotting.area <- check.rasterize.or.warn();
    }
    
    if(rasterize.plotting.area) rasterFunc <- make.mixed.raster.vector.function(raster.height = raster.height, raster.width = raster.width, debugMode = debugMode);
    
    tf.xy <- function(dl, i){
      curr.lanebam <- plotter$lanebam.params$unique.ID[i];
      curr.x <- dl[[curr.lanebam]][[x.name]];
      curr.y <- dl[[curr.lanebam]][[y.name]];
      if(zeroBaseX) curr.x <- curr.x + 1;
      if(norm.x) curr.x <- curr.x / max(curr.x) * 100;
      if(avg.y) {
        if(! all(curr.y == 0)){
          curr.y <- curr.y / sum(curr.y);
        }
      } 
      if(vert.offset != 0) curr.y <- curr.y + plotter$lanebam.params$vert.offsets[i] * vert.offset;
      df <- data.frame(curr.x,curr.y, stringsAsFactors = FALSE);
      names(df) <- c(x.name,y.name);
      df;
    }
    tf.list <- list();
    for(i in 1:length(plotter$lanebam.params$unique.ID)){
      curr.lanebam <- plotter$lanebam.params$unique.ID[i];
      tf.list[[curr.lanebam]] <- tf.xy(data.list,i);
    }
    if(is.null(xlim)){
      xlim.max <- max(sapply(tf.list,function(dl){ max(dl[[x.name]], na.rm = TRUE) }), na.rm = TRUE);
      xlim.min <- min(sapply(tf.list,function(dl){ min(dl[[x.name]], na.rm = TRUE) }), na.rm = TRUE);
      if(norm.x){
        xlim.max <- 100;
        xlim.min <- 0;
      }
      xlim <- c(xlim.min,xlim.max);
    }
    if(is.null(ylim)){
      ylim.max <- max(sapply(tf.list,function(dl){ max(dl[[y.name]], na.rm = TRUE) }), na.rm = TRUE);
      ylim.min <- min(sapply(tf.list,function(dl){ min(dl[[y.name]], na.rm = TRUE) }), na.rm = TRUE);
      if((! is.simple.number(ylim.max)) | (! is.simple.number(ylim.min) )){
        ylim.max <- max(sapply(tf.list,function(dl){ max(  nonsimple.replace(as.numeric(dl[[y.name]]),-Inf), na.rm = TRUE) }), na.rm = TRUE);
        ylim.min <- min(sapply(tf.list,function(dl){ min(  nonsimple.replace(as.numeric(dl[[y.name]]), Inf), na.rm = TRUE) }), na.rm = TRUE);
      }
      if(is.infinite(ylim.max) & is.infinite(ylim.min)){
        ylim.min <- 0;
        ylim.max <- 0;
      }
      if(ylim.max == ylim.min){
        ylim.max <- ylim.max + 0.0001;
      }
      ylim <- c(ylim.min,ylim.max);
    }
    xlim.max <- xlim[2];
    xlim.min <- xlim[1];
    ylim.max <- ylim[2];
    ylim.min <- ylim[1];

    if(plot.medians | plot.means){
      margin.sz <- (ylim.max - ylim.min) * 0.04;
      avg.plot.center <- ylim.min - margin.sz;
      if(plot.medians){
        avg <- lapply(tf.list, function(dl){
          curr.x <- dl[[x.name]];
          curr.y <- dl[[y.name]];
          curr.y <- curr.y / sum(curr.y);
          curr.cs <- cumsum(curr.y);
          index.median <- which(curr.cs > 0.5)[1];
          return(curr.x[index.median]);
        })
      } else {
        avg <- lapply(tf.list, function(dl){
          sum(dl[[x.name]] * dl[[y.name]]);
        });
      }
      ylim[1] <- ylim.min - 2 * margin.sz;
    }
  }, error = function(e){ errorPlot(plot.name,e, code = 1);  stop();});
  plot.new();
  tryCatch({
    plot.window(xlim=xlim,ylim=ylim, ...);
    #plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=xlim,ylim=ylim,axes=F, ...);
    if(! is.null(pre.plot.func)) pre.plot.func();
    if(draw.horiz.lines){
      horiz.lines.at <- floor(ylim[1]):floor(ylim[2] + 1) - 0.5
      #print(horiz.lines.at);
      abline(h = horiz.lines.at, lty=3,col="grey");
    }
    if(norm.x){
      abline(v=50,lty=3,col="grey");
    }
    if(plot.medians | plot.means){
      abline(h=ylim.min,lty=3,col="grey");
      abline(h=ylim.min - (2 * margin.sz),lty=3,col="grey");
      if(plot.medians) text(xlim.min,avg.plot.center,labels="Median:", adj=0);
      if(plot.means)   text(xlim.min,avg.plot.center,labels="Mean:", adj=0);
    }
    
    if(rasterize.plotting.area){
      rasterFunc$initRaster();
      plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2]),ylim=ylim,axes=F, bg = "transparent", ...);
    }
    
    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      for(j in 1:length(p.params$unique.ID)){
        #print(p.params$unique.ID[j]);
        curr.data <- tf.list[[ p.params$unique.ID[j] ]];
        curr.y <- curr.data[[y.name]]
        curr.x <- curr.data[[x.name]]
        if(plot.type == "lines"){
          lines(curr.x, curr.y, lty = ifelse(override.lty == -1, p.params$lines.lty[j], override.lty), 
                #col = color2transparent(p.params$lines.col[j],p.params$lines.alpha[j]),
                col  = p.params$lines.tcol[j],
                lwd = p.params$lines.lwd[j]);
        }
        if(plot.medians | plot.means){
          points(avg[[p.params$unique.ID[j]]], avg.plot.center + p.params$vert.offsets[j] * margin.sz * 2,pch = p.params$points.pch[j], 
                 col = color2transparent(p.params$points.col[j], p.params$points.alpha[j]),
                 lwd = p.params$lines.lwd[j]);
        }
      }
    }
    if(rasterize.plotting.area){
      rasterFunc$closeRaster();
      rasterFunc$printRaster();
    }
    
    #axis(1);
    #axis(2);
    box();
    
    pretty.axis.x <- limited.pretty(xlim,n = xaxis.pretty.n);
    pretty.axis.x.labels <- pretty.axis.x;
    if(x.is.log) {
      pretty.axis.x.is.int <- pretty.axis.x == round(pretty.axis.x);
      pretty.axis.x <- pretty.axis.x[pretty.axis.x.is.int];
      pretty.axis.x.labels <- 10 ^ pretty.axis.x;
    }
    axis(1,at=pretty.axis.x, labels=pretty.axis.x.labels, lwd = -1, lwd.ticks = par("lwd"));
    
    if(draw.horiz.lines){
      #axis(2, at = horiz.lines.at, labels = rep("",length(horiz.lines.at)));
      axis(2, at = c(pretty(ylim) - 0.5, pretty(ylim) + 0.5), rep("",length(pretty(ylim))*2) , lwd = -1, lwd.ticks = par("lwd"));
      axis(2, tcl=0, las=1, lwd = -1, lwd.ticks = par("lwd"));
    } else {
      #axis(2);
      pretty.axis.y <- pretty(ylim);
      pretty.axis.y.labels <- pretty.axis.y;
      if(y.is.log) {
        pretty.axis.y.is.int <- pretty.axis.y == round(pretty.axis.y);
        pretty.axis.y <- pretty.axis.y[pretty.axis.y.is.int];
        pretty.axis.y.labels <- 10 ^ pretty.axis.y ;
      }
      axis(2,at=pretty.axis.y, labels= pretty.axis.y.labels, las = y.axis.las, lwd = -1, lwd.ticks = par("lwd"));
      plot.miniticks(2,log.y = y.is.log,lwd = -1, lwd.ticks = par("lwd")/2);
    }
    box();
    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE); stop();});
}




#data.list.r1 <- get.clipping.from.dl(plotter$res@qc.data[["cigarOpDistribution.byReadCycle.R1"]])
#data.list.r2 <- get.clipping.from.dl(plotter$res@qc.data[["cigarOpDistribution.byReadCycle.R2"]])
#dl <- data.list.r1; x.name <- "CYCLE"; y.name <- "S_T"

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

makePlot.generic.pair <- function(plot.name, data.list.r1, data.list.r2, plotter, x.name, y.name, 
                                  y.names = NULL,
                                  norm.x, avg.y, xlim =NULL , ylim = NULL, vert.offset = 0, horiz.offset = 0, 
                                  draw.horiz.lines = FALSE, plot.type = "lines", r2.buffer = 10, override.lty = -1, 
                                  x.is.log = FALSE, y.is.log = FALSE, pre.plot.func = NULL, y.axis.las = 1, zeroBaseX = FALSE, 
                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                  axis.x.at = NULL, axis.y.at = NULL,
                                  debugMode = FALSE,
                                  ...){

  tryCatch({
    priorities <- sort(unique(plotter$lanebam.params$plot.priority));
    
    if(rasterize.plotting.area){
      rasterize.plotting.area <- check.rasterize.or.warn();
    }
    if(rasterize.plotting.area) rasterFunc <- make.mixed.raster.vector.function(raster.height = raster.height, raster.width = raster.width, debugMode = debugMode);
    
    if(is.null(y.names)){
      y.names <- c(y.name,y.name);
    }
    
    #print(data.list.r1[["CCSS_ALL_3006691_NA"]])
    
    tf.xy <- function(dl, i, r){
      curr.lanebam <- plotter$lanebam.params$unique.ID[i];
      curr.x <- dl[[curr.lanebam]][[x.name]];
      curr.y <- dl[[curr.lanebam]][[y.names[r]]];
      if(zeroBaseX) curr.x <- curr.x + 1;
      if(norm.x) curr.x <- curr.x / max(curr.x) * 100;
      if(avg.y)  curr.y <- curr.y / sum(curr.y);
      if(vert.offset != 0) curr.y <- curr.y + plotter$lanebam.params$vert.offsets[i] * vert.offset;
      df <- data.frame(curr.x,curr.y, stringsAsFactors = FALSE);
      names(df) <- c(x.name,y.name);
      df;
    }
    
    tf.list.r1 <- list();
    tf.list.r2 <- list();
    for(i in 1:length(plotter$lanebam.params$unique.ID)){
      curr.lanebam <- plotter$lanebam.params$unique.ID[i];
      tf.list.r1[[curr.lanebam]] <- tf.xy(data.list.r1,i,1);
      tf.list.r2[[curr.lanebam]] <- tf.xy(data.list.r2,i,2);
    }

    #print(tf.list.r1[["CCSS_ALL_3006691_NA"]])
    
    if(is.null(ylim)){
      ylim.min <- min(sapply(c(tf.list.r1,tf.list.r2),function(dl){ min(dl[[y.name]], na.rm = TRUE) }), na.rm = TRUE);
      ylim.max <- max(sapply(c(tf.list.r1,tf.list.r2),function(dl){ max(dl[[y.name]], na.rm = TRUE) }), na.rm = TRUE);
      if((! is.simple.number(ylim.max)) | (! is.simple.number(ylim.min) )){
        ylim.max <- max(sapply(c(tf.list.r1,tf.list.r2),function(dl){ max(  nonsimple.replace(as.numeric(dl[[y.name]]),-Inf), na.rm = TRUE) }), na.rm = TRUE);
        ylim.min <- min(sapply(c(tf.list.r1,tf.list.r2),function(dl){ min(  nonsimple.replace(as.numeric(dl[[y.name]]), Inf), na.rm = TRUE) }), na.rm = TRUE);
      }
      if(is.infinite(ylim.max) & is.infinite(ylim.min)){
        ylim.min <- 0;
        ylim.max <- 0;
      }
      if(ylim.max == ylim.min){
        ylim.max <- ylim.max + 0.0001;
      }
      ylim <- c(ylim.min, ylim.max);
    }
    if(is.null(xlim)){
      xlim.min <- min(sapply(c(tf.list.r1,tf.list.r2),function(dl){ min(dl[[x.name]], na.rm = TRUE) }), na.rm = TRUE);
      xlim.max <- max(sapply(c(tf.list.r1,tf.list.r2),function(dl){ max(dl[[x.name]], na.rm = TRUE) }), na.rm = TRUE);
      if((! is.simple.number(xlim.max)) | (! is.simple.number(xlim.min) )){
        xlim.max <- max(sapply(c(tf.list.r1,tf.list.r2),function(dl){ max(  nonsimple.replace(as.numeric(dl[[x.name]]),-Inf), na.rm = TRUE) }), na.rm = TRUE);
        xlim.min <- min(sapply(c(tf.list.r1,tf.list.r2),function(dl){ min(  nonsimple.replace(as.numeric(dl[[x.name]]), Inf), na.rm = TRUE) }), na.rm = TRUE);
      }
      if(is.infinite(xlim.max) & is.infinite(xlim.min)){
        xlim.min <- 0;
        xlim.max <- 0;
      }
      if(xlim.max == xlim.min){
        xlim.max <- xlim.max + 0.0001;
      }
      xlim <- c(xlim.min, xlim.max);
    }
    xlim.max <- xlim[2];
    xlim.min <- xlim[1];
    ylim.max <- ylim[2];
    ylim.min <- ylim[1];

    buffer <- r2.buffer;
    mini.buffer <- r2.buffer / 5;
    r2.offset <- xlim[2] + buffer ;
  
  }, error = function(e){ errorPlot(plot.name,e, code = 1);  stop();});

  plot.new();
  
  tryCatch({
  
    ####plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, ...);
    plot.window(xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim, ...);

    if(! is.null(pre.plot.func)) pre.plot.func();
    if(draw.horiz.lines){
      horiz.lines.at <- floor(ylim[1]):floor(ylim[2] + 1) - 0.5
      #print(horiz.lines.at);
      abline(h = horiz.lines.at, lty=3,col="grey");
    }
    
    if(rasterize.plotting.area){
      rasterFunc$initRaster();
      #plot.new();
      #plot.window(xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim, ...);
      plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, bg = "transparent", ...);
      
    }
    
    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      for(j in 1:length(p.params$unique.ID)){
        #print(p.params$unique.ID[j]);
        curr.data <- tf.list.r1[[ p.params$unique.ID[j] ]];
        curr.x <- curr.data[[x.name]];
        curr.y <- curr.data[[y.name]];
        if(plot.type == "lines"){
          #lines(curr.x, curr.y, lty = ifelse(override.lty == -1, p.params$lines.lty[j], override.lty), col = color2transparent(p.params$lines.col[j],p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          lines(curr.x, curr.y, lty = ifelse(override.lty == -1, p.params$lines.lty[j], override.lty), col = p.params$lines.tcol[j], lwd = p.params$lines.lwd[j]);
        }
        #message(min(curr.y), p.params$unique.ID[j]);
        curr.data <- tf.list.r2[[ p.params$unique.ID[j] ]];
        curr.x <- curr.data[[x.name]];
        curr.y <- curr.data[[y.name]];
        curr.x <- curr.x + r2.offset;
        if(plot.type == "lines"){
          lines(curr.x, curr.y, lty = ifelse(override.lty == -1, p.params$lines.lty[j], override.lty), col = p.params$lines.tcol[j], lwd = p.params$lines.lwd[j]);
        }
      }
    }
    
    if(rasterize.plotting.area){
      rasterFunc$closeRaster();
      rasterFunc$printRaster();
    }
    #axis.x.at = NULL, axis.y.at = NULL,
    if(is.null(axis.x.at)){
      pretty.axis.x <- limited.pretty(xlim);
    } else {
      pretty.axis.x <- axis.x.at;
    }
    
    pretty.axis.x.labels <- pretty.axis.x;
    if(x.is.log) {
      pretty.axis.x.is.int <- pretty.axis.x == round(pretty.axis.x);
      pretty.axis.x <- pretty.axis.x[pretty.axis.x.is.int];
      pretty.axis.x.labels <- 10 ^ pretty.axis.x;
    }
    axis(1,at=pretty.axis.x + r2.offset, labels = pretty.axis.x, lwd = -1, lwd.ticks = par("lwd"));
    axis(1,at=pretty.axis.x, labels=pretty.axis.x.labels, lwd = -1, lwd.ticks = par("lwd"));
    box();
    
    if(is.null(axis.y.at)){
      pretty.axis.y <- truncated.pretty(ylim);
    } else {
      pretty.axis.y <- axis.y.at;
    }
    
    if(draw.horiz.lines){
      pretty.axis.y <- pretty(ylim);
      #axis(2, at = horiz.lines.at, labels = rep("",length(horiz.lines.at)));
      axis(2, at = c(pretty.axis.y - 0.5, pretty.axis.y + 0.5), rep("",length(pretty.axis.y)*2) , lwd = -1, lwd.ticks = par("lwd"));
      axis(2, tcl=0, las=1, lwd = -1, lwd.ticks = par("lwd"));
    } else {
      #axis(2);
      pretty.axis.y.labels <- pretty.axis.y;
      if(y.is.log) {
        pretty.axis.y.is.int <- pretty.axis.y == round(pretty.axis.y);
        pretty.axis.y <- pretty.axis.y[pretty.axis.y.is.int];
        pretty.axis.y.labels <- 10 ^ pretty.axis.y ;
      }
      axis(2,at=pretty.axis.y, labels= pretty.axis.y.labels, las = y.axis.las, lwd = -1, lwd.ticks = par("lwd"));
      plot.miniticks(2,log.y = y.is.log,lwd = -1, lwd.ticks = par("lwd")/2);
    }
    usr <- par("usr");
    yrange <- abs(usr[4] - usr[3]);
    rect(xlim[2] + mini.buffer, usr[3] - yrange*0.02,  xlim[2]+buffer - mini.buffer, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
    rect(xlim[2] + mini.buffer, usr[3] - yrange     ,  xlim[2]+buffer - mini.buffer, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);
    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE);  stop();});
}

###################################################################################

#TO FINISH
makePlot.gene.cdf.helper <- function(plot.name, data.list, plotter,plot.intercepts = TRUE, label.intercepts = FALSE, override.lty = -1, pre.plot.func = NULL,
                                     rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                     debugMode = FALSE,
                                     ...){

  tryCatch({

    if(rasterize.plotting.area){
      rasterize.plotting.area <- check.rasterize.or.warn();
    }

   if(rasterize.plotting.area) rasterFunc <- make.mixed.raster.vector.function(raster.height = raster.height, raster.width = raster.width, debugMode = debugMode);
   
   priorities <- sort(unique(plotter$lanebam.params$plot.priority));
   
   gene.ct <- length(data.list[[1]])
   
   tf.y <- function(dl, i){
     curr.lanebam <- plotter$lanebam.params$unique.ID[i];
     read.ct <- dl[[curr.lanebam]][ gene.ct ];
     tfy <- dl[[curr.lanebam]] / read.ct;
     tfy;
   }
   tf.list <- list();
   for(i in 1:length(plotter$lanebam.params$unique.ID)){
     curr.lanebam <- plotter$lanebam.params$unique.ID[i];
     tf.list[[curr.lanebam]] <- tf.y(data.list,i);
   }

   curr.x <- 1:gene.ct;
   xlim <- c(1,gene.ct);
   if(label.intercepts) xlim[1] <- 0.65;

  }, error = function(e){ errorPlot(plot.name,e, code = 1);  stop();});

  plot.new();
  
  tryCatch({

   #plot(1,1, log="x",type='l',ylim=c(0,1),xlim=xlim,axes=F,xlab="",ylab="",col="white", ...);
   plot.window(xlim=xlim,ylim=c(0,1),log="x",...);
   if(! is.null(pre.plot.func)) pre.plot.func();
   
   if(rasterize.plotting.area){
     rasterFunc$initRaster();
     plot(1,1, log="x",type='l',ylim=c(0,1),xlim=xlim,axes=F,xlab="",ylab="",col="white", bg = "transparent", ...);
   }

   for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      for(j in 1:length(p.params$unique.ID)){
        curr.y <- tf.list[[ p.params$unique.ID[j] ]];
        lines(curr.x,curr.y, lty = ifelse(override.lty == -1, p.params$lines.lty[j], override.lty), col = p.params$lines.tcol[j], lwd = p.params$lines.lwd[j]);
      }
   }

   if(rasterize.plotting.area){
     rasterFunc$closeRaster();
     rasterFunc$printRaster();
   }

   box();
   title(ylab="Cumulative % of total reads",xlab="# Genes");
   axis(2,at=c(0,0.25,0.5,0.75,1.0),labels=c("0%","25%","50%","75%","100%"));
   axis(2,at=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),labels=rep("",9),tcl=-0.3);

   axis(1, at=c(1,10,100,1000,10000));
   axis(1, at=c(gene.ct));

   if(plot.intercepts){
     p <- max(priorities);
     p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];

     xs.2 <- 0.1;
     if(label.intercepts) xs.2 <- 0.9;

     internal.internal.plot.crosshairs <- function(i, j, curr.lanebam, curr.pct){
           ys.1 = -1;
           xs.1 = j;
           ye.1 = curr.pct;
           xe.1 = j;
           lines(c(xs.1,xe.1),c(ys.1,ye.1),lty = 3, col = p.params$lines.tcol[i], lwd = p.params$lines.lwd[j]);
           ys.2 = curr.pct;
           #xs.2 = 0.9;
           ye.2 = curr.pct;
           xe.2 = xe.1;
           lines(c(xs.2,xe.2),c(ys.2,ye.2),lty = 3, col = p.params$lines.tcol[i], lwd = p.params$lines.lwd[j]);
     }

     for(i in 1:length(p.params$unique.ID)){
        curr.lanebam <- p.params$unique.ID[i];
        for(j in c(1,10,100,1000,10000)){
           if(j < length( tf.list[[ curr.lanebam ]] )){
              curr.pct <- tf.list[[ curr.lanebam ]][j];
              if(curr.pct == 1){
                 break;
              }
              internal.internal.plot.crosshairs(i,j,curr.lanebam,curr.pct);
           }
        }
        j <- min(which(tf.list[[ curr.lanebam ]] == 1));
        curr.pct <- 1;
        internal.internal.plot.crosshairs(i,j,curr.lanebam,curr.pct);
     }
   }
   return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE);  stop();});

}

makePlot.gene.cdf.bySample.helper <- function(plot.name, data.list, curr.sample, plot.intercepts = TRUE, label.intercepts = TRUE, pre.plot.func = NULL, 
                                              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                              debugMode = FALSE,
                                              ...){
  tryCatch({
  
    if(rasterize.plotting.area){
      rasterize.plotting.area <- check.rasterize.or.warn();
    }
   if(rasterize.plotting.area) rasterFunc <- make.mixed.raster.vector.function(raster.height = raster.height, raster.width = raster.width, debugMode = debugMode);

   gene.ct <- length(data.list[[1]])
   
   tf.list <- list();
   for(i in 1:length(data.list)){
     tf.list[[i]] <- data.list[[i]] / data.list[[i]][ length(data.list[[i]]) ];
   }
   names(tf.list) <- names(data.list);

   curr.x <- 1:gene.ct;
   xlim <- c(1,gene.ct);
   if(label.intercepts) xlim[1] <- 0.65;

  }, error = function(e){ errorPlot(plot.name,e, code = 1); stop();});
  plot.new();
  
  tryCatch({
    plot.window(xlim=xlim,ylim=c(0,1),log="x",...);

    # plot(1,1, log="x",type='l',ylim=c(0,1),xlim=xlim,axes=F,xlab="",ylab="",col="white", ...);
     if(! is.null(pre.plot.func)) pre.plot.func();

     if(rasterize.plotting.area){
       rasterFunc$initRaster();
       plot(1,1, log="x",type='l',ylim=c(0,1),xlim=xlim,axes=F,xlab="",ylab="",col="white", ...);
     }

        for(j in 1:length(tf.list)){
          if(names(tf.list)[j] != curr.sample){
             curr.y <- tf.list[[j]];
             lines(curr.x,curr.y, lty = 1, col = "gray");
          }
        }
        j <- which(names(tf.list) == curr.sample);
        curr.y <- tf.list[[j]];
        lines(curr.x,curr.y, lty = 1, col = "red");

     if(rasterize.plotting.area){
       rasterFunc$closeRaster();
       rasterFunc$printRaster();
     }

     box();
     title(ylab="Cumulative % of total reads",xlab="# Genes");
     axis(2,at=c(0,0.25,0.5,0.75,1.0),labels=c("0%","25%","50%","75%","100%"), lwd = -1, lwd.ticks = par("lwd"));
     axis(2,at=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),labels=rep("",9),tcl=-0.3, lwd = -1, lwd.ticks = par("lwd"));

     axis(1, at=c(1,10,100,1000,10000), lwd = -1, lwd.ticks = par("lwd"));
     axis(1, at=c(gene.ct), lwd = -1, lwd.ticks = par("lwd"));

     if(plot.intercepts){
       xs.2 <- 0.1;
       if(label.intercepts) xs.2 <- 0.9;

       internal.internal.plot.crosshairs <- function( j, curr.pct){
             ys.1 = -1;
             xs.1 = j;
             ye.1 = curr.pct;
             xe.1 = j;
             lines(c(xs.1,xe.1),c(ys.1,ye.1),lty = 3, col = "red");
             ys.2 = curr.pct;
             #xs.2 = 0.9;
             ye.2 = curr.pct;
             xe.2 = xe.1;
             lines(c(xs.2,xe.2),c(ys.2,ye.2),lty = 3, col = "red");
       }


       i <- which(names(tf.list) == curr.sample);
       for(j in c(1,10,100,1000,10000)){
          if(j < length( tf.list[[ i ]] )){
             curr.pct <- tf.list[[ i ]][j];
            if(curr.pct == 1){
                break;
             }
             internal.internal.plot.crosshairs(j,curr.pct);
             if(label.intercepts) text(0.9,curr.pct,paste0(round(curr.pct * 100,1),""), adj=c(1,0.5),col="red", cex=0.75);
          }
        }
        j <- min(which(tf.list[[ i ]] == 1));
        curr.pct <- 1;
        xs.2 <- 0.1;
        internal.internal.plot.crosshairs(j,curr.pct);
        if(label.intercepts) text(j,0,j,adj=c(0.5,1),col="red");
      }
      return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE);  stop();});
}

################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################


makePlot.generic.NVC.single.OLD <- function(plot.name, data.list.r1, plotter, x.name, y.name, xlim, ylim,
                                  label.major = FALSE, label.cutoff = 0.5, blank.label.char = "", 
                                  points.highlighted = FALSE, color.unhighlighted.light = TRUE, 
                                  x.axis.labels  = NULL, count.x.from.end = FALSE, pre.plot.func = NULL, label.major.cex = 1, 
                                  show.base.legend = TRUE, nvc.colors = NULL, nvc.colors.light = NULL, 
                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 2000, debugMode = FALSE, ... ){

    if(rasterize.plotting.area){
      rasterize.plotting.area <- check.rasterize.or.warn();
    }

  if(rasterize.plotting.area){
    return(makePlot.generic.NVC.single.RASTERIZED(plot.name, data.list.r1 = data.list.r1, plotter = plotter, x.name = x.name, y.name = y.name, xlim = xlim, ylim = ylim,
                                      label.major = label.major, label.cutoff = label.cutoff, blank.label.char = blank.label.char, 
                                      points.highlighted = points.highlighted, color.unhighlighted.light = color.unhighlighted.light, 
                                      x.axis.labels  = x.axis.labels, count.x.from.end = count.x.from.end, pre.plot.func = pre.plot.func, label.major.cex = label.major.cex, 
                                      show.base.legend = show.base.legend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, 
                                      raster.height = raster.height, raster.width = raster.width, debugMode = debugMode, ... ))
  } else {
    return(makePlot.generic.NVC.single.DEFAULT(plot.name, data.list.r1 = data.list.r1, plotter = plotter, x.name = x.name, y.name = y.name, xlim = xlim, ylim = ylim,
                                      label.major = label.major, label.cutoff = label.cutoff, blank.label.char = blank.label.char, 
                                      points.highlighted = points.highlighted, color.unhighlighted.light = color.unhighlighted.light, 
                                      x.axis.labels  = x.axis.labels, count.x.from.end = count.x.from.end, pre.plot.func = pre.plot.func, label.major.cex = label.major.cex, 
                                      show.base.legend = show.base.legend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, debugMode = debugMode,
                                      ... ))
  }
}

makePlot.generic.NVC.single.DEFAULT <- function(plot.name, data.list.r1, plotter, x.name, y.name, xlim, ylim, 
                                  label.major = FALSE, label.cutoff = 0.5, blank.label.char = "", 
                                  points.highlighted = FALSE, color.unhighlighted.light = TRUE, 
                                  x.axis.labels  = NULL, count.x.from.end = FALSE, pre.plot.func = NULL, label.major.cex = 1, 
                                  show.base.legend = TRUE, nvc.colors = NULL, nvc.colors.light = NULL, debugMode = FALSE, ... ){
  
   tryCatch({
     bases <- names(plotter$nvc.colors);
     priorities <- sort(unique(plotter$lanebam.params$plot.priority));

     tf <- function(df){
        split <- lapply(1:length(bases), function(i){
          b <- bases[i];
          df$CT[ df$base == b ];
        });
        names(split) <- bases;
        sums <- split[[1]];
        for(i in 2:length(bases)){
          sums <- sums + split[[i]];
        }
        out <- lapply(split, function(x){
          x / sums;
        });
        return(out);
     }
     tf.list.r1 <- lapply(data.list.r1,tf);

    xlim.max <- xlim[2];
    xlim.min <- xlim[1];
    ylim.max <- ylim[2];
    ylim.min <- ylim[1];
  }, error = function(e){ errorPlot(plot.name,e, code = 1);  stop();});
  
  plot.new();
  
  tryCatch({
    #plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, ...);
    plot.window(xlim=c(xlim[1],xlim[2]),ylim=ylim, ...);

    if(! is.null(pre.plot.func)) pre.plot.func();

    abline(h=0.5,lty=3,col="gray");
    abline(h=0.25,lty=3,col="gray");
    if(label.major & label.cutoff != 0.5) abline(h = label.cutoff, lty=3,col="gray");

    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      curr.nvc.colors <- plotter$nvc.colors;
      if(p < 2) curr.nvc.colors <- plotter$nvc.colors.light;
      for(j in 1:length(p.params$unique.ID)){
        for(k in 1:length(bases)){
          curr.data <- tf.list.r1[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col =  color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
        }
      }
    }

    if(label.major){
      p <- 2;
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];

       #j <- 1;
       get.major.labels.by.lanebam <- function(dl,j){
         curr.lanebam <- p.params$unique.ID[j];
         ML <- rep(blank.label.char, xlim.max);
         maxRate <- rep(label.cutoff,xlim.max);
         for(k in 1:length(bases)){
           curr.y <- dl[[curr.lanebam]][[bases[k]]];
           for(m in 1:length(curr.y)){
             if(f.na(curr.y[m] > maxRate[m])){
               ML[m] <- bases[k];
               maxRate[m] <- curr.y[m];
             }
           }
         }
         return(ML);
       }

       collapse.major.labels <- function(ML.list){
         ML.final <- rep(blank.label.char, xlim.max);
         for(m in 1:xlim.max){
           ML.curr <- sapply(ML.list,function(ml){ ml[m] });
           if(f.na(all(ML.curr == ML.curr[1]))){
             ML.final[m] <- ML.curr[1];
           }
         }
         return(ML.final);
       }

       ML.list.r1 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r1, j) });

       ML.r1 <- collapse.major.labels(ML.list.r1);

       get.ML.col <- function(ML){
         sapply(ML, function(ml){
           if(ml == blank.label.char){ "gray";
           } else { plotter$nvc.colors[[ml]]; }
         });
       }

       #print("ML.list.r1:");
       #print(ML.list.r1);
       #print("ML.r1:");
       #print(ML.r1);
       #print("ML.r1.col:");
       #print(get.ML.col(ML.r1));

       text(1:xlim.max, rep(par("usr")[4],xlim.max), ML.r1, adj = c(0.5,1), col = get.ML.col(ML.r1) , cex = label.major.cex);
    }

    pretty.axis <- limited.pretty(xlim);
    pretty.axis <- pretty.axis[pretty.axis == round(pretty.axis)];

    if(count.x.from.end){ 
      max.cycle.ct <- max(plotter$res@decoder$cycle.CT);
      axis(1,at=pretty.axis,labels= max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
    } else if(! is.null(x.axis.labels)){
      x.labels <- x.axis.labels[pretty.axis];
      axis(1,at=pretty.axis,labels= x.labels, lwd = -1, lwd.ticks = par("lwd"));
    } else {
      axis(1,at=pretty.axis,labels=pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
    }

    if(show.base.legend){
      y.center <- (par("usr")[4] - par("usr")[3]) / 2 + par("usr")[3];
      bases.list <- names(nvc.colors);
      base.legend.adj = c(-0.2,-0.2);

      #test(par("usr")[2],par("usr")[3], );
      #mtext(c("A\nT\nG\nC"),side = 4);

      axis(4,at=y.center, labels = "A\n \n \n ",col="transparent",col.axis = nvc.colors[["A"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \nT\n \n ",col="transparent",col.axis = nvc.colors[["T"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \nG\n ",col="transparent",col.axis = nvc.colors[["G"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \n \nC",col="transparent",col.axis = nvc.colors[["C"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));

      #text(par("usr")[2],par("usr")[3],"A\n\n\n", adj=base.legend.adj, col = nvc.colors[["A"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\nT\n\n", adj=base.legend.adj, col = nvc.colors[["T"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\nG\n", adj=base.legend.adj, col = nvc.colors[["G"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\n\nC", adj=base.legend.adj, col = nvc.colors[["C"]], font=2, xpd=TRUE)

      #text(par("usr")[2],par("usr")[3],expression("A\n"*phantom("T\nG\nC")),             adj=base.legend.adj, col = nvc.colors[["A"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\n")*"T\n"*phantom("G\nC")), adj=base.legend.adj, col = nvc.colors[["T"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\n")*"G\n"*phantom("C")), adj=base.legend.adj, col = nvc.colors[["G"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\nG\n")*"C"),             adj=base.legend.adj, col = nvc.colors[["C"]])

    }


    box();
    axis(2, lwd = -1, lwd.ticks = par("lwd"));
    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE);  stop();});

}

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


makePlot.generic.NVC.single.RASTERIZED <- function(plot.name, data.list.r1,plotter, x.name, y.name, xlim, ylim, 
                                  label.major = FALSE, label.cutoff = 0.5, blank.label.char = "", 
                                  points.highlighted = FALSE, color.unhighlighted.light = TRUE, 
                                  x.axis.labels  = NULL, count.x.from.end = FALSE, pre.plot.func = NULL, label.major.cex = 1, 
                                  show.base.legend = TRUE, nvc.colors = NULL, nvc.colors.light = NULL, 
                                  raster.height = 1000, raster.width = 2000, debugMode = FALSE, ... ){

   tryCatch({

     rasterFunc <- make.mixed.raster.vector.function(raster.height = raster.height, raster.width = raster.width, debugMode = debugMode);

     bases <- names(plotter$nvc.colors);
     priorities <- sort(unique(plotter$lanebam.params$plot.priority));

     tf <- function(df){
        split <- lapply(1:length(bases), function(i){
          b <- bases[i];
          df$CT[ df$base == b ];
        });
        names(split) <- bases;
        sums <- split[[1]];
        for(i in 2:length(bases)){
          sums <- sums + split[[i]];
        }
        out <- lapply(split, function(x){
          x / sums;
        });
        return(out);
     }
     tf.list.r1 <- lapply(data.list.r1,tf);

    xlim.max <- xlim[2];
    xlim.min <- xlim[1];
    ylim.max <- ylim[2];
    ylim.min <- ylim[1];

    #mini.buffer <- (r2.buffer + 1) / 5;

  }, error = function(e){ errorPlot(plot.name,e, code = 1);  stop();});
  plot.new();
  
  tryCatch({
    #plot.window(xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim, ...);
    plot.window(xlim=c(xlim[1],xlim[2]),ylim=ylim, ...);
    #plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, ...);
    rasterPlot <- function(){
       plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2]),ylim=ylim,axes=F, bg = "transparent", ...);
    }

    if(! is.null(pre.plot.func)) pre.plot.func();

    abline(h=0.5,lty=3,col="gray");
    abline(h=0.25,lty=3,col="gray");
    if(label.major & label.cutoff != 0.5) abline(h = label.cutoff, lty=3,col="gray");

    rasterFunc$initRaster();
    rasterPlot();
    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      curr.nvc.colors <- plotter$nvc.colors;
      if(p < 2) curr.nvc.colors <- plotter$nvc.colors.light;
      for(j in 1:length(p.params$unique.ID)){
        for(k in 1:length(bases)){
          curr.data <- tf.list.r1[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col =  color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
        }
      }
    }
    rasterFunc$closeRaster();
    rasterFunc$printRaster();

    if(label.major){
      p <- 2;
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];

       #j <- 1;
       get.major.labels.by.lanebam <- function(dl,j){
         curr.lanebam <- p.params$unique.ID[j];
         ML <- rep(blank.label.char, xlim.max);
         maxRate <- rep(label.cutoff,xlim.max);
         for(k in 1:length(bases)){
           curr.y <- dl[[curr.lanebam]][[bases[k]]];
           for(m in 1:length(curr.y)){
             if(f.na(curr.y[m] > maxRate[m])){
               ML[m] <- bases[k];
               maxRate[m] <- curr.y[m];
             }
           }
         }
         return(ML);
       }

       collapse.major.labels <- function(ML.list){
         ML.final <- rep(blank.label.char, xlim.max);
         for(m in 1:xlim.max){
           ML.curr <- sapply(ML.list,function(ml){ ml[m] });
           if(f.na(all(ML.curr == ML.curr[1]))){
             ML.final[m] <- ML.curr[1];
           }
         }
         return(ML.final);
       }

       ML.list.r1 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r1, j) });

       ML.r1 <- collapse.major.labels(ML.list.r1);

       get.ML.col <- function(ML){
         sapply(ML, function(ml){
           if(ml == blank.label.char){ "gray";
           } else { plotter$nvc.colors[[ml]]; }
         });
       }

       #print("ML.list.r1:");
       #print(ML.list.r1);
       #print("ML.r1:");
       #print(ML.r1);
       #print("ML.r1.col:");
       #print(get.ML.col(ML.r1));

       text(1:xlim.max, rep(par("usr")[4],xlim.max), ML.r1, adj = c(0.5,1), col = get.ML.col(ML.r1) , cex = label.major.cex);
    }

    pretty.axis <- limited.pretty(xlim);
    pretty.axis <- pretty.axis[pretty.axis == round(pretty.axis)];

    if(count.x.from.end){ 
      max.cycle.ct <- max(plotter$res@decoder$cycle.CT);
      axis(1,at=pretty.axis,labels= max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
    } else if(! is.null(x.axis.labels)){
      x.labels <- x.axis.labels[pretty.axis];
      axis(1,at=pretty.axis,labels= x.labels, lwd = -1, lwd.ticks = par("lwd"));
    } else {
      axis(1,at=pretty.axis,labels=pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
    }

    if(show.base.legend){
      y.center <- (par("usr")[4] - par("usr")[3]) / 2 + par("usr")[3];
      bases.list <- names(nvc.colors);
      base.legend.adj = c(-0.2,-0.2);

      #test(par("usr")[2],par("usr")[3], );
      #mtext(c("A\nT\nG\nC"),side = 4);

      axis(4,at=y.center, labels = "A\n \n \n ",col="transparent",col.axis = nvc.colors[["A"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \nT\n \n ",col="transparent",col.axis = nvc.colors[["T"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \nG\n ",col="transparent",col.axis = nvc.colors[["G"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \n \nC",col="transparent",col.axis = nvc.colors[["C"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));

      #text(par("usr")[2],par("usr")[3],"A\n\n\n", adj=base.legend.adj, col = nvc.colors[["A"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\nT\n\n", adj=base.legend.adj, col = nvc.colors[["T"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\nG\n", adj=base.legend.adj, col = nvc.colors[["G"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\n\nC", adj=base.legend.adj, col = nvc.colors[["C"]], font=2, xpd=TRUE)

      #text(par("usr")[2],par("usr")[3],expression("A\n"*phantom("T\nG\nC")),             adj=base.legend.adj, col = nvc.colors[["A"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\n")*"T\n"*phantom("G\nC")), adj=base.legend.adj, col = nvc.colors[["T"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\n")*"G\n"*phantom("C")), adj=base.legend.adj, col = nvc.colors[["G"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\nG\n")*"C"),             adj=base.legend.adj, col = nvc.colors[["C"]])

    }


    box();
    axis(2, lwd = -1, lwd.ticks = par("lwd"));
    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE); stop();});
}


makePlot.generic.NVC.single <- function(plot.name, data.list.r1,plotter, x.name, y.name, xlim, ylim, 
                                  label.major = FALSE, label.cutoff = 0.5, blank.label.char = "", 
                                  points.highlighted = FALSE, color.unhighlighted.light = TRUE, 
                                  x.axis.labels  = NULL, count.x.from.end = FALSE, pre.plot.func = NULL, label.major.cex = 1, 
                                  show.base.legend = TRUE, nvc.colors = NULL, nvc.colors.light = NULL, 
                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 2000, debugMode = FALSE, ... ){
    if(rasterize.plotting.area){
      rasterize.plotting.area <- check.rasterize.or.warn();
    }

   tryCatch({

     if(rasterize.plotting.area) rasterFunc <- make.mixed.raster.vector.function(raster.height = raster.height, raster.width = raster.width, debugMode = debugMode);

     bases <- names(plotter$nvc.colors);
     priorities <- sort(unique(plotter$lanebam.params$plot.priority));
     tf.list.r1 <- data.list.r1;

     xlim.max <- xlim[2];
     xlim.min <- xlim[1];
     ylim.max <- ylim[2];
     ylim.min <- ylim[1];

    #mini.buffer <- (r2.buffer + 1) / 5;

  }, error = function(e){ errorPlot(plot.name,e, code = 1);  stop();});
  plot.new();
  
  tryCatch({
    plot.window(xlim=c(xlim[1],xlim[2]),ylim=ylim, ...);
    rasterPlot <- function(){
       plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2]),ylim=ylim,axes=F, bg = "transparent", ...);
    }
    if(! is.null(pre.plot.func)) pre.plot.func();
    abline(h=0.5,lty=3,col="gray");
    abline(h=0.25,lty=3,col="gray");
    if(label.major & label.cutoff != 0.5) abline(h = label.cutoff, lty=3,col="gray");
    if(rasterize.plotting.area) rasterFunc$initRaster();
    if(rasterize.plotting.area) rasterPlot();
    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      curr.nvc.colors <- plotter$nvc.colors;
      if(p < 2) curr.nvc.colors <- plotter$nvc.colors.light;
      for(j in 1:length(p.params$unique.ID)){
        for(k in 1:length(bases)){
          curr.data <- tf.list.r1[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col =  color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
        }
      }
    }
    if(rasterize.plotting.area) rasterFunc$closeRaster();
    if(rasterize.plotting.area) rasterFunc$printRaster();

    if(label.major){
      p <- 2;
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];

       get.major.labels.by.lanebam <- function(dl,j){
         curr.lanebam <- p.params$unique.ID[j];
         ML <- rep(blank.label.char, xlim.max);
         maxRate <- rep(label.cutoff,xlim.max);
         for(k in 1:length(bases)){
           curr.y <- dl[[curr.lanebam]][[bases[k]]];
           for(m in 1:length(curr.y)){
             if(f.na(curr.y[m] > maxRate[m])){
               ML[m] <- bases[k];
               maxRate[m] <- curr.y[m];
             }
           }
         }
         return(ML);
       }

       collapse.major.labels <- function(ML.list){
         ML.final <- rep(blank.label.char, xlim.max);
         for(m in 1:xlim.max){
           ML.curr <- sapply(ML.list,function(ml){ ml[m] });
           if(f.na(all(ML.curr == ML.curr[1]))){
             ML.final[m] <- ML.curr[1];
           }
         }
         return(ML.final);
       }

       ML.list.r1 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r1, j) });

       ML.r1 <- collapse.major.labels(ML.list.r1);

       get.ML.col <- function(ML){
         sapply(ML, function(ml){
           if(ml == blank.label.char){ "gray";
           } else { plotter$nvc.colors[[ml]]; }
         });
       }

       text(1:xlim.max, rep(par("usr")[4],xlim.max), ML.r1, adj = c(0.5,1), col = get.ML.col(ML.r1) , cex = label.major.cex);
    }

    pretty.axis <- limited.pretty(xlim);
    pretty.axis <- pretty.axis[pretty.axis == round(pretty.axis)];

    if(count.x.from.end){ 
      max.cycle.ct <- max(plotter$res@decoder$cycle.CT);
      axis(1,at=pretty.axis,labels= max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
    } else if(! is.null(x.axis.labels)){
      x.labels <- x.axis.labels[pretty.axis];
      axis(1,at=pretty.axis,labels= x.labels, lwd = -1, lwd.ticks = par("lwd"));
    } else {
      axis(1,at=pretty.axis,labels=pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
    }

    if(show.base.legend){
      y.center <- (par("usr")[4] - par("usr")[3]) / 2 + par("usr")[3];
      bases.list <- names(nvc.colors);
      base.legend.adj = c(-0.2,-0.2);
      axis(4,at=y.center, labels = "A\n \n \n ",col="transparent",col.axis = nvc.colors[["A"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \nT\n \n ",col="transparent",col.axis = nvc.colors[["T"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \nG\n ",col="transparent",col.axis = nvc.colors[["G"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \n \nC",col="transparent",col.axis = nvc.colors[["C"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
    }
    
    box();
    axis(2, lwd = -1, lwd.ticks = par("lwd"));
    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE); stop();});
}


################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################

makePlot.generic.NVC.pair.OLD <- function(plot.name, data.list.r1, data.list.r2, plotter, x.name, y.name, xlim, ylim, r2.buffer = 10, 
                                  label.major = FALSE, label.cutoff = 0.5, blank.label.char = "", 
                                  points.highlighted = FALSE, color.unhighlighted.light = TRUE, 
                                  x.axis.labels  = NULL, count.x.from.end = FALSE, pre.plot.func = NULL, label.major.cex = 1, 
                                  show.base.legend = TRUE, nvc.colors = NULL, nvc.colors.light = NULL, 
                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 2000, debugMode = FALSE, ... ){
    if(rasterize.plotting.area){
      rasterize.plotting.area <- check.rasterize.or.warn();
    }

  if(rasterize.plotting.area){
    return(makePlot.generic.NVC.pair.RASTERIZED(plot.name, data.list.r1 = data.list.r1, data.list.r2 = data.list.r2, plotter = plotter, x.name = x.name, y.name = y.name, xlim = xlim, ylim = ylim, r2.buffer = r2.buffer, 
                                      label.major = label.major, label.cutoff = label.cutoff, blank.label.char = blank.label.char, 
                                      points.highlighted = points.highlighted, color.unhighlighted.light = color.unhighlighted.light, 
                                      x.axis.labels  = x.axis.labels, count.x.from.end = count.x.from.end, pre.plot.func = pre.plot.func, label.major.cex = label.major.cex, 
                                      show.base.legend = show.base.legend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, 
                                      raster.height = raster.height, raster.width = raster.width, debugMode = debugMode, ... ))
  } else {
    return(makePlot.generic.NVC.pair.DEFAULT(plot.name, data.list.r1 = data.list.r1, data.list.r2 = data.list.r2, plotter = plotter, x.name = x.name, y.name = y.name, xlim = xlim, ylim = ylim, r2.buffer = r2.buffer, 
                                      label.major = label.major, label.cutoff = label.cutoff, blank.label.char = blank.label.char, 
                                      points.highlighted = points.highlighted, color.unhighlighted.light = color.unhighlighted.light, 
                                      x.axis.labels  = x.axis.labels, count.x.from.end = count.x.from.end, pre.plot.func = pre.plot.func, label.major.cex = label.major.cex, 
                                      show.base.legend = show.base.legend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, debugMode = debugMode,
                                      ... ))
  }
}

makePlot.generic.NVC.pair.DEFAULT <- function(plot.name, data.list.r1, data.list.r2, plotter, x.name, y.name, xlim, ylim, r2.buffer = 10, 
                                  label.major = FALSE, label.cutoff = 0.5, blank.label.char = "", 
                                  points.highlighted = FALSE, color.unhighlighted.light = TRUE, 
                                  x.axis.labels  = NULL, count.x.from.end = FALSE, pre.plot.func = NULL, label.major.cex = 1, 
                                  show.base.legend = TRUE, nvc.colors = NULL, nvc.colors.light = NULL, debugMode = FALSE, ... ){
  
   tryCatch({
     bases <- names(plotter$nvc.colors);
     priorities <- sort(unique(plotter$lanebam.params$plot.priority));

     tf <- function(df){
        split <- lapply(1:length(bases), function(i){
          b <- bases[i];
          df$CT[ df$base == b ];
        });
        names(split) <- bases;
        sums <- split[[1]];
        for(i in 2:length(bases)){
          sums <- sums + split[[i]];
        }
        out <- lapply(split, function(x){
          x / sums;
        });
        return(out);
     }
     tf.list.r1 <- lapply(data.list.r1,tf);
     tf.list.r2 <- lapply(data.list.r2,tf);

    xlim.max <- xlim[2];
    xlim.min <- xlim[1];
    ylim.max <- ylim[2];
    ylim.min <- ylim[1];

    mini.buffer <- (r2.buffer + 1) / 5;
    r2.offset <- xlim[2] + r2.buffer;
  }, error = function(e){ errorPlot(plot.name,e, code = 1); stop();});
  
  plot.new();
  
  tryCatch({
    #plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, ...);
    plot.window(xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim, ...);

    if(! is.null(pre.plot.func)) pre.plot.func();

    abline(h=0.5,lty=3,col="gray");
    abline(h=0.25,lty=3,col="gray");
    if(label.major & label.cutoff != 0.5) abline(h = label.cutoff, lty=3,col="gray");

    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      curr.nvc.colors <- plotter$nvc.colors;
      if(p < 2) curr.nvc.colors <- plotter$nvc.colors.light;
      for(j in 1:length(p.params$unique.ID)){
        for(k in 1:length(bases)){
          curr.data <- tf.list.r1[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col =  color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
          curr.data <- tf.list.r2[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          curr.x <- curr.x  + r2.offset;
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
        }
      }
    }

    if(label.major){
      p <- 2;
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];

       #j <- 1;
       get.major.labels.by.lanebam <- function(dl,j){
         curr.lanebam <- p.params$unique.ID[j];
         ML <- rep(blank.label.char, xlim.max);
         maxRate <- rep(label.cutoff,xlim.max);
         for(k in 1:length(bases)){
           curr.y <- dl[[curr.lanebam]][[bases[k]]];
           for(m in 1:length(curr.y)){
             if(f.na(curr.y[m] > maxRate[m])){
               ML[m] <- bases[k];
               maxRate[m] <- curr.y[m];
             }
           }
         }
         return(ML);
       }

       collapse.major.labels <- function(ML.list){
         ML.final <- rep(blank.label.char, xlim.max);
         for(m in 1:xlim.max){
           ML.curr <- sapply(ML.list,function(ml){ ml[m] });
           if(f.na(all(ML.curr == ML.curr[1]))){
             ML.final[m] <- ML.curr[1];
           }
         }
         return(ML.final);
       }

       ML.list.r1 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r1, j) });
       ML.list.r2 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r2, j) }); 

       ML.r1 <- collapse.major.labels(ML.list.r1);
       ML.r2 <- collapse.major.labels(ML.list.r2);

       get.ML.col <- function(ML){
         sapply(ML, function(ml){
           if(ml == blank.label.char){ "gray";
           } else { plotter$nvc.colors[[ml]]; }
         });
       }

       #print("ML.list.r1:");
       #print(ML.list.r1);
       #print("ML.r1:");
       #print(ML.r1);
       #print("ML.r1.col:");
       #print(get.ML.col(ML.r1));

       text(1:xlim.max, rep(par("usr")[4],xlim.max), ML.r1, adj = c(0.5,1), col = get.ML.col(ML.r1) , cex = label.major.cex);
       text(1:xlim.max + r2.offset, rep(par("usr")[4],xlim.max), ML.r2, adj = c(0.5,1), col = get.ML.col(ML.r2), cex = label.major.cex  );
    }

    pretty.axis <- limited.pretty(xlim);
    pretty.axis <- pretty.axis[pretty.axis == round(pretty.axis)];

    if(count.x.from.end){ 
      max.cycle.ct <- max(plotter$res@decoder$cycle.CT);
      axis(1,at=pretty.axis + r2.offset, labels = max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels= max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
    } else if(! is.null(x.axis.labels)){
      x.labels <- x.axis.labels[pretty.axis];
      axis(1,at=pretty.axis + r2.offset, labels = x.labels, lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels= x.labels, lwd = -1, lwd.ticks = par("lwd"));
    } else {
      axis(1,at=pretty.axis + r2.offset, labels = pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels=pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
    }

    if(show.base.legend){
      y.center <- (par("usr")[4] - par("usr")[3]) / 2 + par("usr")[3];
      bases.list <- names(nvc.colors);
      base.legend.adj = c(-0.2,-0.2);

      #test(par("usr")[2],par("usr")[3], );
      #mtext(c("A\nT\nG\nC"),side = 4);

      axis(4,at=y.center, labels = "A\n \n \n ",col="transparent",col.axis = nvc.colors[["A"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \nT\n \n ",col="transparent",col.axis = nvc.colors[["T"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \nG\n ",col="transparent",col.axis = nvc.colors[["G"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \n \nC",col="transparent",col.axis = nvc.colors[["C"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));

      #text(par("usr")[2],par("usr")[3],"A\n\n\n", adj=base.legend.adj, col = nvc.colors[["A"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\nT\n\n", adj=base.legend.adj, col = nvc.colors[["T"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\nG\n", adj=base.legend.adj, col = nvc.colors[["G"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\n\nC", adj=base.legend.adj, col = nvc.colors[["C"]], font=2, xpd=TRUE)

      #text(par("usr")[2],par("usr")[3],expression("A\n"*phantom("T\nG\nC")),             adj=base.legend.adj, col = nvc.colors[["A"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\n")*"T\n"*phantom("G\nC")), adj=base.legend.adj, col = nvc.colors[["T"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\n")*"G\n"*phantom("C")), adj=base.legend.adj, col = nvc.colors[["G"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\nG\n")*"C"),             adj=base.legend.adj, col = nvc.colors[["C"]])

    }


    box();
    axis(2, lwd = -1, lwd.ticks = par("lwd"));
    #rect(xlim[2] + mini.buffer, ylim[1] - ylim[2], xlim[2] + 1 + r2.buffer - mini.buffer, ylim[2] * 2, border="black",col="white");   
    usr <- par("usr");
    yrange <- abs(usr[4] - usr[3]);
    rect(xlim[2] + mini.buffer, usr[3] - yrange*0.02,  xlim[2]+1+r2.buffer - mini.buffer, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
    rect(xlim[2] + mini.buffer, usr[3] - yrange     ,  xlim[2]+1+r2.buffer - mini.buffer, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);
    
    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE); stop();});

}

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


makePlot.generic.NVC.pair.RASTERIZED <- function(plot.name, data.list.r1, data.list.r2, plotter, x.name, y.name, xlim, ylim, r2.buffer = 10, 
                                  label.major = FALSE, label.cutoff = 0.5, blank.label.char = "", 
                                  points.highlighted = FALSE, color.unhighlighted.light = TRUE, 
                                  x.axis.labels  = NULL, count.x.from.end = FALSE, pre.plot.func = NULL, label.major.cex = 1, 
                                  show.base.legend = TRUE, nvc.colors = NULL, nvc.colors.light = NULL, 
                                  raster.height = 1000, raster.width = 2000, debugMode = FALSE, ... ){

   tryCatch({

     rasterFunc <- make.mixed.raster.vector.function(raster.height = raster.height, raster.width = raster.width, debugMode = debugMode);

     bases <- names(plotter$nvc.colors);
     priorities <- sort(unique(plotter$lanebam.params$plot.priority));

     tf <- function(df){
        split <- lapply(1:length(bases), function(i){
          b <- bases[i];
          df$CT[ df$base == b ];
        });
        names(split) <- bases;
        sums <- split[[1]];
        for(i in 2:length(bases)){
          sums <- sums + split[[i]];
        }
        out <- lapply(split, function(x){
          x / sums;
        });
        return(out);
     }
     tf.list.r1 <- lapply(data.list.r1,tf);
     tf.list.r2 <- lapply(data.list.r2,tf);

    xlim.max <- xlim[2];
    xlim.min <- xlim[1];
    ylim.max <- ylim[2];
    ylim.min <- ylim[1];

    mini.buffer <- (r2.buffer + 1) / 5;
    r2.offset <- xlim[2] + r2.buffer;
  
  }, error = function(e){ errorPlot(plot.name,e, code = 1); stop();});
  plot.new();
  
  tryCatch({
    plot.window(xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim, ...);
    #plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, ...);
    rasterPlot <- function(){
       plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, bg = "transparent", ...);
    }

    if(! is.null(pre.plot.func)) pre.plot.func();

    abline(h=0.5,lty=3,col="gray");
    abline(h=0.25,lty=3,col="gray");
    if(label.major & label.cutoff != 0.5) abline(h = label.cutoff, lty=3,col="gray");

    rasterFunc$initRaster();
    rasterPlot();
    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      curr.nvc.colors <- plotter$nvc.colors;
      if(p < 2) curr.nvc.colors <- plotter$nvc.colors.light;
      for(j in 1:length(p.params$unique.ID)){
        for(k in 1:length(bases)){
          curr.data <- tf.list.r1[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col =  color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
          curr.data <- tf.list.r2[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          curr.x <- curr.x  + r2.offset;
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
        }
      }
    }
    rasterFunc$closeRaster();
    rasterFunc$printRaster();

    if(label.major){
      p <- 2;
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];

       #j <- 1;
       get.major.labels.by.lanebam <- function(dl,j){
         curr.lanebam <- p.params$unique.ID[j];
         ML <- rep(blank.label.char, xlim.max);
         maxRate <- rep(label.cutoff,xlim.max);
         for(k in 1:length(bases)){
           curr.y <- dl[[curr.lanebam]][[bases[k]]];
           for(m in 1:length(curr.y)){
             if(f.na(curr.y[m] > maxRate[m])){
               ML[m] <- bases[k];
               maxRate[m] <- curr.y[m];
             }
           }
         }
         return(ML);
       }

       collapse.major.labels <- function(ML.list){
         ML.final <- rep(blank.label.char, xlim.max);
         for(m in 1:xlim.max){
           ML.curr <- sapply(ML.list,function(ml){ ml[m] });
           if(f.na(all(ML.curr == ML.curr[1]))){
             ML.final[m] <- ML.curr[1];
           }
         }
         return(ML.final);
       }

       ML.list.r1 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r1, j) });
       ML.list.r2 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r2, j) }); 

       ML.r1 <- collapse.major.labels(ML.list.r1);
       ML.r2 <- collapse.major.labels(ML.list.r2);

       get.ML.col <- function(ML){
         sapply(ML, function(ml){
           if(ml == blank.label.char){ "gray";
           } else { plotter$nvc.colors[[ml]]; }
         });
       }

       #print("ML.list.r1:");
       #print(ML.list.r1);
       #print("ML.r1:");
       #print(ML.r1);
       #print("ML.r1.col:");
       #print(get.ML.col(ML.r1));

       text(1:xlim.max, rep(par("usr")[4],xlim.max), ML.r1, adj = c(0.5,1), col = get.ML.col(ML.r1) , cex = label.major.cex);
       text(1:xlim.max + r2.offset, rep(par("usr")[4],xlim.max), ML.r2, adj = c(0.5,1), col = get.ML.col(ML.r2), cex = label.major.cex  );
    }

    pretty.axis <- limited.pretty(xlim);
    pretty.axis <- pretty.axis[pretty.axis == round(pretty.axis)];

    if(count.x.from.end){ 
      max.cycle.ct <- max(plotter$res@decoder$cycle.CT);
      axis(1,at=pretty.axis + r2.offset, labels = max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels= max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
    } else if(! is.null(x.axis.labels)){
      x.labels <- x.axis.labels[pretty.axis];
      axis(1,at=pretty.axis + r2.offset, labels = x.labels, lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels= x.labels, lwd = -1, lwd.ticks = par("lwd"));
    } else {
      axis(1,at=pretty.axis + r2.offset, labels = pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels=pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
    }

    if(show.base.legend){
      y.center <- (par("usr")[4] - par("usr")[3]) / 2 + par("usr")[3];
      bases.list <- names(nvc.colors);
      base.legend.adj = c(-0.2,-0.2);

      #test(par("usr")[2],par("usr")[3], );
      #mtext(c("A\nT\nG\nC"),side = 4);

      axis(4,at=y.center, labels = "A\n \n \n ",col="transparent",col.axis = nvc.colors[["A"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \nT\n \n ",col="transparent",col.axis = nvc.colors[["T"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \nG\n ",col="transparent",col.axis = nvc.colors[["G"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \n \nC",col="transparent",col.axis = nvc.colors[["C"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));

      #text(par("usr")[2],par("usr")[3],"A\n\n\n", adj=base.legend.adj, col = nvc.colors[["A"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\nT\n\n", adj=base.legend.adj, col = nvc.colors[["T"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\nG\n", adj=base.legend.adj, col = nvc.colors[["G"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\n\nC", adj=base.legend.adj, col = nvc.colors[["C"]], font=2, xpd=TRUE)

      #text(par("usr")[2],par("usr")[3],expression("A\n"*phantom("T\nG\nC")),             adj=base.legend.adj, col = nvc.colors[["A"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\n")*"T\n"*phantom("G\nC")), adj=base.legend.adj, col = nvc.colors[["T"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\n")*"G\n"*phantom("C")), adj=base.legend.adj, col = nvc.colors[["G"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\nG\n")*"C"),             adj=base.legend.adj, col = nvc.colors[["C"]])

    }


    box();
    axis(2, lwd = -1, lwd.ticks = par("lwd"));
    #rect(xlim[2] + mini.buffer, ylim[1] - ylim[2], xlim[2] + 1 + r2.buffer - mini.buffer, ylim[2] * 2, border="black",col="white");
    #rect(xlim[2] + mini.buffer, ylim[1] - abs(ylim[2] - ylim[1])*0.02, xlim[2]+buffer - mini.buffer, (abs(ylim[1]) + abs(ylim[2])) * 2, border="white",col="white", xpd = FALSE,...);
    #rect(xlim[2] + mini.buffer, ylim[1] - abs(ylim[2] - ylim[1])*0.02, xlim[2]+buffer - mini.buffer, (abs(ylim[1]) + abs(ylim[2])) * 2, border="black",col="white", xpd = TRUE,...);
    usr <- par("usr");
    yrange <- abs(usr[4] - usr[3]);
    rect(xlim[2] + mini.buffer, usr[3] - yrange*0.02,  xlim[2]+1+r2.buffer - mini.buffer, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
    rect(xlim[2] + mini.buffer, usr[3] - yrange     ,  xlim[2]+1+r2.buffer - mini.buffer, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);

    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE); stop();});
}


makePlot.generic.NVC.pair <- function(plot.name, data.list.r1, data.list.r2, plotter, x.name, y.name, xlim, ylim, r2.buffer = 10, 
                                  label.major = FALSE, label.cutoff = 0.5, blank.label.char = "", 
                                  points.highlighted = FALSE, color.unhighlighted.light = TRUE, 
                                  x.axis.labels  = NULL, count.x.from.end = FALSE, pre.plot.func = NULL, label.major.cex = 1, 
                                  show.base.legend = TRUE, nvc.colors = NULL, nvc.colors.light = NULL, 
                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 2000, debugMode = FALSE, ... ){

    if(rasterize.plotting.area){
      rasterize.plotting.area <- check.rasterize.or.warn();
    }

   tryCatch({

     if(rasterize.plotting.area) rasterFunc <- make.mixed.raster.vector.function(raster.height = raster.height, raster.width = raster.width, debugMode = debugMode);

     bases <- names(plotter$nvc.colors);
     priorities <- sort(unique(plotter$lanebam.params$plot.priority));

     tf.list.r1 <- data.list.r1
     tf.list.r2 <- data.list.r2

     xlim.max <- xlim[2];
     xlim.min <- xlim[1];
     ylim.max <- ylim[2];
     ylim.min <- ylim[1];

     mini.buffer <- (r2.buffer + 1) / 5;
     r2.offset <- xlim[2] + r2.buffer;
  
  }, error = function(e){ errorPlot(plot.name,e, code = 1); stop();});
  plot.new();
  
  tryCatch({
    plot.window(xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim, ...);
    #plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, ...);
    rasterPlot <- function(){
       plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=c(xlim[1],xlim[2] * 2 + r2.buffer),ylim=ylim,axes=F, bg = "transparent", ...);
    }

    if(! is.null(pre.plot.func)) pre.plot.func();

    abline(h=0.5,lty=3,col="gray");
    abline(h=0.25,lty=3,col="gray");
    if(label.major & label.cutoff != 0.5) abline(h = label.cutoff, lty=3,col="gray");

    if(rasterize.plotting.area) rasterFunc$initRaster();
    if(rasterize.plotting.area) rasterPlot();
    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      curr.nvc.colors <- plotter$nvc.colors;
      if(p < 2) curr.nvc.colors <- plotter$nvc.colors.light;
      for(j in 1:length(p.params$unique.ID)){
        for(k in 1:length(bases)){
          curr.data <- tf.list.r1[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col =  color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
          curr.data <- tf.list.r2[[ p.params$unique.ID[j] ]][[bases[k]]];
          curr.x <- 1:length(curr.data)
          curr.y <- curr.data
          curr.x <- curr.x  + r2.offset;
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$lines.alpha[j]), lwd = p.params$lines.lwd[j]);
          if(p == 2 & points.highlighted){
            points(curr.x,curr.y, pch = p.params$points.pch[j], col = color2transparent(curr.nvc.colors[[bases[k]]], p.params$points.alpha[j]));
          }
        }
      }
    }
    if(rasterize.plotting.area) rasterFunc$closeRaster();
    if(rasterize.plotting.area) rasterFunc$printRaster();

    if(label.major){
      p <- 2;
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];

       #j <- 1;
       get.major.labels.by.lanebam <- function(dl,j){
         curr.lanebam <- p.params$unique.ID[j];
         ML <- rep(blank.label.char, xlim.max);
         maxRate <- rep(label.cutoff,xlim.max);
         for(k in 1:length(bases)){
           curr.y <- dl[[curr.lanebam]][[bases[k]]];
           for(m in 1:length(curr.y)){
             if(f.na(curr.y[m] > maxRate[m])){
               ML[m] <- bases[k];
               maxRate[m] <- curr.y[m];
             }
           }
         }
         return(ML);
       }

       collapse.major.labels <- function(ML.list){
         ML.final <- rep(blank.label.char, xlim.max);
         for(m in 1:xlim.max){
           ML.curr <- sapply(ML.list,function(ml){ ml[m] });
           if(f.na(all(ML.curr == ML.curr[1]))){
             ML.final[m] <- ML.curr[1];
           }
         }
         return(ML.final);
       }

       ML.list.r1 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r1, j) });
       ML.list.r2 <- lapply(1:length(p.params$unique.ID), function(j){ get.major.labels.by.lanebam(tf.list.r2, j) }); 

       ML.r1 <- collapse.major.labels(ML.list.r1);
       ML.r2 <- collapse.major.labels(ML.list.r2);

       get.ML.col <- function(ML){
         sapply(ML, function(ml){
           if(ml == blank.label.char){ "gray";
           } else { plotter$nvc.colors[[ml]]; }
         });
       }

       #print("ML.list.r1:");
       #print(ML.list.r1);
       #print("ML.r1:");
       #print(ML.r1);
       #print("ML.r1.col:");
       #print(get.ML.col(ML.r1));

       text(1:xlim.max, rep(par("usr")[4],xlim.max), ML.r1, adj = c(0.5,1), col = get.ML.col(ML.r1) , cex = label.major.cex);
       text(1:xlim.max + r2.offset, rep(par("usr")[4],xlim.max), ML.r2, adj = c(0.5,1), col = get.ML.col(ML.r2), cex = label.major.cex  );
    }

    pretty.axis <- limited.pretty(xlim);
    pretty.axis <- pretty.axis[pretty.axis == round(pretty.axis)];

    if(count.x.from.end){ 
      max.cycle.ct <- max(plotter$res@decoder$cycle.CT);
      axis(1,at=pretty.axis + r2.offset, labels = max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels= max.cycle.ct - (xlim.max - pretty.axis), lwd = -1, lwd.ticks = par("lwd"));
    } else if(! is.null(x.axis.labels)){
      x.labels <- x.axis.labels[pretty.axis];
      axis(1,at=pretty.axis + r2.offset, labels = x.labels, lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels= x.labels, lwd = -1, lwd.ticks = par("lwd"));
    } else {
      axis(1,at=pretty.axis + r2.offset, labels = pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
      axis(1,at=pretty.axis,labels=pretty.axis, lwd = -1, lwd.ticks = par("lwd"));
    }

    if(show.base.legend){
      y.center <- (par("usr")[4] - par("usr")[3]) / 2 + par("usr")[3];
      bases.list <- names(nvc.colors);
      base.legend.adj = c(-0.2,-0.2);

      #test(par("usr")[2],par("usr")[3], );
      #mtext(c("A\nT\nG\nC"),side = 4);

      axis(4,at=y.center, labels = "A\n \n \n ",col="transparent",col.axis = nvc.colors[["A"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \nT\n \n ",col="transparent",col.axis = nvc.colors[["T"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \nG\n ",col="transparent",col.axis = nvc.colors[["G"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));
      axis(4,at=y.center, labels = " \n \n \nC",col="transparent",col.axis = nvc.colors[["C"]], font = 2, las=1, hadj=0, tcl=0, line = -0.5, lwd = -1, lwd.ticks = par("lwd"));

      #text(par("usr")[2],par("usr")[3],"A\n\n\n", adj=base.legend.adj, col = nvc.colors[["A"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\nT\n\n", adj=base.legend.adj, col = nvc.colors[["T"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\nG\n", adj=base.legend.adj, col = nvc.colors[["G"]], font=2, xpd=TRUE)
      #text(par("usr")[2],par("usr")[3],"\n\n\nC", adj=base.legend.adj, col = nvc.colors[["C"]], font=2, xpd=TRUE)

      #text(par("usr")[2],par("usr")[3],expression("A\n"*phantom("T\nG\nC")),             adj=base.legend.adj, col = nvc.colors[["A"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\n")*"T\n"*phantom("G\nC")), adj=base.legend.adj, col = nvc.colors[["T"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\n")*"G\n"*phantom("C")), adj=base.legend.adj, col = nvc.colors[["G"]])
      #text(par("usr")[2],par("usr")[3],expression(phantom("A\nT\nG\n")*"C"),             adj=base.legend.adj, col = nvc.colors[["C"]])

    }


    box();
    axis(2, lwd = -1, lwd.ticks = par("lwd"));
    #rect(xlim[2] + mini.buffer, ylim[1] - ylim[2], xlim[2] + 1 + r2.buffer - mini.buffer, ylim[2] * 2, border="black",col="white");
    #rect(xlim[2] + mini.buffer, ylim[1] - abs(ylim[2] - ylim[1])*0.02, xlim[2]+buffer - mini.buffer, (abs(ylim[1]) + abs(ylim[2])) * 2, border="white",col="white", xpd = FALSE,...);
    #rect(xlim[2] + mini.buffer, ylim[1] - abs(ylim[2] - ylim[1])*0.02, xlim[2]+buffer - mini.buffer, (abs(ylim[1]) + abs(ylim[2])) * 2, border="black",col="white", xpd = TRUE,...);
    usr <- par("usr");
    yrange <- abs(usr[4] - usr[3]);
    rect(xlim[2] + mini.buffer, usr[3] - yrange*0.02,  xlim[2]+1+r2.buffer - mini.buffer, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
    rect(xlim[2] + mini.buffer, usr[3] - yrange     ,  xlim[2]+1+r2.buffer - mini.buffer, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);

    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE); stop();});
}

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


generic.points.tf <- function(data.list, x.names, x.titles, norm.by = NULL, field.name = "FIELD", 
                              count.name = "COUNT", offset.horiz = 0, horiz.offsets = NULL, x.start = 0,log.y = FALSE){
  tf.xy <- function(i){
    curr.lanebam <- names(data.list)[i];
    df <- data.list[[curr.lanebam]];
    curr.x <- 1:length(x.names) + x.start;
    curr.y <- as.numeric(sapply(x.names,function(n){ 
      df[[count.name]][df[[field.name]] == n];
    }))
    if(! is.null(norm.by)) {
      if(length(norm.by) == 1){
        norm.val <- as.numeric(df[[count.name]][df[[field.name]] == norm.by]);
      } else {
        norm.val <- sum(as.numeric(df[[count.name]][df[[field.name]] %in% norm.by]));
      }
      curr.y <- curr.y / norm.val;
    }
    if(offset.horiz != 0 ) curr.x <- curr.x + horiz.offsets[i] * offset.horiz;
    df <- data.frame(curr.x,curr.y,x.titles,stringsAsFactors=F);
    names(df) <- c("x","y","x.titles");
    if(log.y){
      df$y <- log10(df$y);
      df$y <- ifelse(is.finite(df$y),df$y,NA);
    }
    df;
  }
  tfed <- lapply(1:length(data.list), tf.xy);
  names(tfed) <- names(data.list);
  return( tfed );
}


generic.points.tf.from.df <- function(df, x.names, x.titles, offset.horiz = 0, horiz.offsets = NULL, skip.missing = FALSE){
  keep <- rep(TRUE,length(x.names));
  if(skip.missing){
    keep <- x.names %in% names(df);
  }
  x.names <- x.names[keep];
  x.titles <- x.titles[keep];
  
  tf.xy <- function(i){
    curr.lanebam <- df$unique.ID[i];
    
    curr.x <- 1:length(x.names);
    curr.y <- as.numeric(df[i, x.names]);
    
    if(offset.horiz != 0 ) curr.x <- curr.x + horiz.offsets[i] * offset.horiz;
    out.df <- cbind.data.frame(curr.x,curr.y,x.titles,stringsAsFactors=F);
    names(out.df) <- c("x","y","x.titles");
    out.df;
  }
  tfed <- lapply(1:length(df[,1]), tf.xy);
  names(tfed) <- df$unique.ID;
  return( tfed );
}

merge.points.tf <- function(tf.list.1,tf.list.2){
  return(lapply(1:length(tf.list.1), function(i){
    rbind(tf.list.1[[i]],tf.list.2[[i]]);
  }));
}

makePlot.generic.points <- function(plot.name, tf.list, plotter, plot.type = "points", 
                                    ylim = NULL, leave.blank.cols = 0, label.y = TRUE, 
                                    cex.x.axis = NULL, pre.plot.func = NULL, draw.lines = FALSE,
                                    y.is.log = FALSE,
                                     ...){

  tryCatch({
    priorities <- sort(unique(plotter$lanebam.params$plot.priority));

    x.titles <- tf.list[[1]]$x.titles

    #set X-axis limits:
    xlim.max <- length(x.titles);
    xlim.min <- 1;
    xlim <- c(xlim.min - 0.5, xlim.max + 0.5 + leave.blank.cols);

    #set Y-axis limits:
    if(is.null(ylim)){
        ylim.max <- max(sapply(tf.list,function(dl){ max(as.numeric(dl$y)) }));
        ylim.min <- min(sapply(tf.list,function(dl){ min(as.numeric(dl$y)) }));
        if((! is.simple.number(ylim.max)) | (! is.simple.number(ylim.min) )){
          ylim.max <- max(sapply(tf.list,function(dl){ max(  nonsimple.replace(as.numeric(dl$y),-Inf), na.rm = TRUE) }), na.rm = TRUE);
          ylim.min <- min(sapply(tf.list,function(dl){ min(  nonsimple.replace(as.numeric(dl$y), Inf), na.rm = TRUE) }), na.rm = TRUE);
        }
        if(is.infinite(ylim.max) & is.infinite(ylim.min)){
          ylim.min <- 0;
          ylim.max <- 0;
        }
        ylim <- c(ylim.min,ylim.max);
    } else {
      ylim.max <- ylim[2];
      ylim.min <- ylim[1];
    }
    if(ylim.max == ylim.min){
      ylim.max = ylim.min + 0.0001;
      ylim <- c(ylim.min,ylim.max);
    }
  }, error = function(e){ errorPlot(plot.name,e, code = 1); stop();});
  plot.new();
  
  tryCatch({
    #logxy <- if(ylog){"y"}else{""};
    plot.window(xlim=xlim,ylim=ylim,...);

    #plot(0,0,col="transparent",main="",xlab="",ylab="",xlim=xlim,ylim=ylim,axes=F, ...);
    if(! is.null(pre.plot.func)) pre.plot.func();


    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      for(j in 1:length(p.params$unique.ID)){
        curr.data <- tf.list[[ p.params$unique.ID[j] ]];
        curr.y <- curr.data$y
        curr.x <- curr.data$x
        if(plot.type == "points"){
          points(curr.x, curr.y, pch = p.params$points.pch[j], col = p.params$points.tcol[j], lwd =  p.params$lines.lwd[j], ...);
        }
        if(draw.lines){
          lines(curr.x, curr.y, lty = p.params$lines.lty[j], col = p.params$lines.tcol[j], lwd =  p.params$lines.lwd[j], ...)
        }
      }
    }
    #axis(1);
    #axis(1,at=1:xlim.max,labels=x.titles,tcl=0,adj=c(0.5,1));
    y.abs.min <- par("usr")[3];

    if(is.null(cex.x.axis)){
      cex.x.axis <- par("cex.axis");
      bottom.margin.size <- abs(device.limits()[3] - par("usr")[3]);
      max.x.axis.title.height <- bottom.margin.size;
      if(! any(grepl("\n",x.titles,fixed=TRUE))){
        max.x.axis.title.height <- bottom.margin.size * 0.5;
      }
      cex.x.axis <- fit.character.vector.helper(x.titles, cex.x.axis, min.width = 0.6, max.width = 0.95, max.width.per.char = 0.15, max.height = max.x.axis.title.height);
    }
    
    text(1:xlim.max,rep(y.abs.min,xlim.max),labels = x.titles, adj=c(0.5,1.05) , xpd=T, cex = cex.x.axis, ...);

    axis(1,at=(1:xlim.max - 0.5),labels=F, lwd = -1, lwd.ticks = par("lwd"), ...);
    axis(1,at=(1:xlim.max + 0.5),labels=F, lwd = -1, lwd.ticks = par("lwd"), ...);

    if(label.y) {
      if(y.is.log){
        #decades <- getAllDecades(10^ylim[1],10^ylim[2]);
        #axis(2,at=log10(decades),labels=decades,lwd=-1,lwd.ticks=par("lwd"),...);
        #plot.miniticks(2,...);
      } else {
        axis(2, lwd = -1, lwd.ticks = par("lwd"), ...);
      }
    }
    box();
    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE); stop();});
}



makePlot.generic.points.right <- function(plot.name, tf.list, plotter, plot.type = "points", section.offset = 0, 
                                          ylim.plot = NULL, ylim.data = NULL, label.y = TRUE, cex.x.axis = NULL, ...){
  tryCatch({

    priorities <- sort(unique(plotter$lanebam.params$plot.priority));

    x.titles <- tf.list[[1]]$x.titles

    if(is.null(ylim.plot)){
      ylim.plot.full <- par("usr")[c(3,4)];
      ylim.range.full <- ylim.plot.full[2] - ylim.plot.full[1];
      ylim.range <- ylim.range.full / 1.08;
      ylim.mar <- (ylim.range.full - ylim.range) / 2;
      ylim.plot <- c(ylim.plot.full[1] + ylim.mar, ylim.plot.full[2] - ylim.mar);
    }
    
    #set Y-axis limits:
    if(is.null(ylim.data)){
      ylim.max <- max(sapply(tf.list,function(dl){ max(as.numeric(dl$y)) }));
      ylim.min <- min(sapply(tf.list,function(dl){ min(as.numeric(dl$y)) }));
      if((! is.simple.number(ylim.max)) | (! is.simple.number(ylim.min) )){
        ylim.max <- max(sapply(tf.list,function(dl){ max(  nonsimple.replace(as.numeric(dl$y),-Inf), na.rm = TRUE) }), na.rm = TRUE);
        ylim.min <- min(sapply(tf.list,function(dl){ min(  nonsimple.replace(as.numeric(dl$y), Inf), na.rm = TRUE) }), na.rm = TRUE);
      }
      if(is.infinite(ylim.max) & is.infinite(ylim.min)){
        ylim.min <- 0;
        ylim.max <- 0;
      }
    } else {
      ylim.max <- ylim.data[2];
      ylim.min <- ylim.data[1];
    }
    if(ylim.max == ylim.min){
      ylim.max = ylim.min + 0.0001;
    }
    ylim.data <- c(ylim.min,ylim.max);
    ylim.plot.range <- ylim.plot[2] - ylim.plot[1];
    ylim.data.range <- ylim.data[2] - ylim.data[1];
    tf.y <- function(y){
      y.unit.scale <- (y - ylim.data[1]) / ylim.data.range;
      y.new.scale <- (y.unit.scale * ylim.plot.range) + ylim.plot[1];
      return(y.new.scale);
    }

    xlim.min <- section.offset;
    xlim.max <- section.offset + length(x.titles);

    for(p in priorities){
      p.params <- plotter$lanebam.params[plotter$lanebam.params$plot.priority == p,];
      for(j in 1:length(p.params$unique.ID)){
        curr.data <- tf.list[[ p.params$unique.ID[j] ]];
        curr.y <- tf.y(curr.data$y);
        curr.x <- curr.data$x + xlim.min;
        if(plot.type == "points"){
          points(curr.x, curr.y, pch = p.params$points.pch[j], col = p.params$points.tcol[j], lwd =  p.params$lines.lwd[j], ...);
        }
      }
    }
    #axis(1);
    #axis(1,at=1:xlim.max,labels=x.titles,tcl=0,adj=c(0.5,1));
    y.abs.min <- par("usr")[3];

    if(is.null(cex.x.axis)){
      cex.x.axis <- par("cex.axis");
      bottom.margin.size <- abs(device.limits()[3] - par("usr")[3]);
      max.x.axis.title.height <- bottom.margin.size;
      if(! any(grepl("\n",x.titles,fixed=TRUE))){
        max.x.axis.title.height <- bottom.margin.size * 0.5;
      }
      cex.x.axis <- fit.character.vector.helper(x.titles, cex.x.axis, min.width = 0.6, max.width = 0.95, max.width.per.char = 0.15, max.height = max.x.axis.title.height);
    }
    
    text((xlim.min + 1):xlim.max,rep(y.abs.min,length(x.titles)),labels = x.titles, adj=c(0.5,1.05) , xpd=T, cex = cex.x.axis, ...);
    axis(1,at=((xlim.min + 1):xlim.max - 0.5),labels=F, lwd = -1, lwd.ticks = par("lwd"), ...);
    axis(1,at=((xlim.min + 1):xlim.max + 0.5),labels=F, lwd = -1, lwd.ticks = par("lwd"), ...);

    if(label.y){
      pretty.ax <- pretty(ylim.data, n = 5);
      axis(4, at = tf.y(pretty.ax), labels= pretty.ax, lwd = -1, lwd.ticks = par("lwd"), ...);
    }
    box();
    return(TRUE);
  }, error = function(e){ errorPlot(plot.name, e, code = 2, newPlot = FALSE); stop();});

}






########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
################ Plotter internals:

internal.plot.legend <- function(plotter, legend.type, legend.pos, 
                                 autofit.limits = NULL, 
                                 ncol = NULL,
                                 lwd = NULL,
                                 ...){
  if(plotter$plot.type == "highlightSample.colorByLane"){
    
  } else if(plotter$plot.type == "highlightSample") {
    
  } else if(plotter$plot.type == "highlightSample.colorByLane.NVC") {
  
  } else if(plotter$plot.type == "colorByGroup"){
  
  } else if(plotter$plot.type == "colorByLane"){
  } else if(plotter$plot.type == "colorBySample"){ 
  } else if(plotter$plot.type == "colorByX"){ 
  } else if(plotter$plot.type == "colorByXhighlightByY"){ 

  } else if(plotter$plot.type == "summary"){
    #Do nothing. plot nothing.
    return("");
  } else {
    #Do nothing. Plot nothing.
    return("");
  }



  if(legend.type == "lnpt"){
    if(is.null(ncol)){
      ncol <- 1;
      if(length(plotter$legend.params$name) >= 6) ncol <- 2;
    }
    if(is.null(lwd)){
      lwd <- max(plotter$lanebam.params$lines.lwd) * 1.5
      pt.lwd <- max(plotter$lanebam.params$lines.lwd)
    } else {
      lwd <- lwd;
      pt.lwd <- lwd;
    }
    
    if(is.null(autofit.limits)){
      legend(legend.pos,
           legend=plotter$legend.params$name,
           lty=plotter$legend.params$lines.lty, 
           col=plotter$legend.params$lines.col, 
           pch=plotter$legend.params$points.pch,
           ncol = ncol,xpd=NA,
           lwd = lwd, pt.lwd = pt.lwd,
           ...);
    } else {
      x <- legend(legend.pos,
           legend=plotter$legend.params$name,
           lty=plotter$legend.params$lines.lty, 
           col=plotter$legend.params$lines.col, 
           pch=plotter$legend.params$points.pch,
           ncol = ncol,xpd=NA,
           plot=FALSE,
           lwd = lwd, pt.lwd = pt.lwd,
           ...);
      isTooBig <- function(x){
        x$rect$left < autofit.limits[1] || x$rect$left + x$rect$w > autofit.limits[2] || x$rect$top > autofit.limits[4] || x$rect$top - x$rect$h < autofit.limits[3];
      }
      #isTooBig <- function(x){
      #  x$rect$w > abs(autofit.limits[1] - autofit.limits[2]) || x$rect$h > abs(autofit.limits[3] - autofit.limits[4]);
      #}
      if(isTooBig(x)){
        new.cex <- 1;
        while(isTooBig(x)){
          new.cex <- new.cex * 0.9;
          x <- legend(legend.pos,
                 legend=plotter$legend.params$name,
                 lty=plotter$legend.params$lines.lty, 
                 col=plotter$legend.params$lines.col, 
                 pch=plotter$legend.params$points.pch,
                 ncol = ncol,
                 plot=FALSE,xpd=NA,
                 cex = new.cex,
                 lwd = lwd, pt.lwd = pt.lwd,
                 ...);
        }
      } else {
        test.cex <- 1;
        new.cex <- 1;
        while(! isTooBig(x)){
          new.cex <- test.cex;
          test.cex <- new.cex * 1.1;
          x <- legend(legend.pos,
                 legend=plotter$legend.params$name,
                 lty=plotter$legend.params$lines.lty, 
                 col=plotter$legend.params$lines.col, 
                 pch=plotter$legend.params$points.pch,
                 ncol = ncol,
                 plot=FALSE,xpd=NA,
                 cex = test.cex,
                 lwd = lwd, pt.lwd = pt.lwd,
                 ...);
        }
      }
      legend(legend.pos,
                 legend=plotter$legend.params$name,
                 lty=plotter$legend.params$lines.lty, 
                 col=plotter$legend.params$lines.col, 
                 pch=plotter$legend.params$points.pch,
                 ncol = ncol,xpd=NA,
                 cex = new.cex,
                 lwd = lwd, pt.lwd = pt.lwd,
                 ...);
    }
    
  } else if(legend.type == "lines"){
    legend(legend.pos,legend=plotter$legend.params$name,lty=plotter$legend.params$lines.lty, col=plotter$legend.params$lines.col,xpd=NA, ...);
  } else if(legend.type == "points"){
    legend(legend.pos,legend=plotter$legend.params$name,pch=plotter$legend.params$points.pch, col=plotter$legend.params$points.col,xpd=NA, ...);
  } else if(legend.type == "NVC.curr"){
    labels <- c(paste0(plotter$title.highlight.name,":"),names(plotter$nvc.colors));
    col <- c("transparent",unlist(plotter$nvc.colors));
    legend(legend.pos,legend=labels,lty=rep(plotter$legend.params$lines.lty[1], length(labels)), col=col, xpd=NA, lwd = rep(plotter$legend.params$lines.lwd[1], length(labels)), ...);
  } else if(legend.type == "NVC.other"){
    labels <- c(paste0("Other Samples",":"),names(plotter$nvc.colors.light));
    col <- c("transparent",unlist(plotter$nvc.colors.light));
    legend(legend.pos,legend=labels,lty=rep(plotter$legend.params$lines.lty[2], length(labels)), col=col,xpd=NA, ...);
  } else if(legend.type == "NVC.both"){
    labels.1 <- c(paste0(plotter$title.highlight.name,":"),names(plotter$nvc.colors));
    col.1 <- c("transparent",unlist(plotter$nvc.colors));
    lty.1 <- rep(plotter$legend.params$lines.lty[1], length(labels.1))
    lwd.1 <- rep(plotter$legend.params$lines.lwd[1], length(labels.1));
    labels.2 <- c(paste0("Other Samples",":"),names(plotter$nvc.colors.light));
    col.2 <- c("transparent",unlist(plotter$nvc.colors.light));
    lty.2 <- rep(plotter$legend.params$lines.lty[2], length(labels.2))
    lwd.2 <- rep(plotter$legend.params$lines.lwd[2], length(labels.2));
    
    legend(legend.pos,
           legend=c(labels.1,labels.2),
           lty=c(lty.1,lty.2), 
           col=c(col.1,col.2), 
           lwd = c(lwd.1,lwd.2),
           ncol=2, xpd=NA,
           ...);
  }
}

internal.plot.main.title <- function(plot.title, plotter, plot.type = "", print.title = TRUE, print.simplified.title = FALSE, cex.main = NULL, ...){


  if(print.title & print.simplified.title){
    title.text <- plot.title;
  } else if(print.title) {
    title.text <- paste0(plot.title,internal.get.main.title.fragment(plotter, plot.type))
  }
  
  if(print.title){
    if(! is.null(cex.main)){
      use.cex <- cex.main;
    } else {
      use.cex <- fit.title(title.text);
    }
    title(main = title.text, cex.main = use.cex, ...);
  }
  #else do nothing.
}

internal.get.main.title.fragment <- function(plotter, plot.type){
  if(plot.type == "NVC"){
    diff.text <- "Marked";
  } else {
    diff.text <- "Colored";
  }

  if(plotter$plot.type == "highlightSample.colorByLane" & plot.type == "NVC")        return( paste0("\nWith Sample ",plotter$title.highlight.name," in Bold"));
  if(plotter$plot.type == "highlightSample.colorByLane")        return( paste0("\nWith Sample ",plotter$title.highlight.name," ",diff.text," by lane"));
  if(plotter$plot.type == "highlightSample")                    return( paste0("\nWith Sample ",plotter$title.highlight.name," in ", plotter$legend.params$lines.col[1]));
  if(plotter$plot.type == "highlightSample.colorByLane.NVC")    return( paste0("\nWith Sample ",plotter$title.highlight.name," in Bold"));
  if(plotter$plot.type == "summary") return("");
  if(plotter$plot.type == "colorByLane") return(paste0("\n",diff.text," by Lane"));
  if(plotter$plot.type == "colorByGroup") return(paste0("\n",diff.text," by Group"));
  if(plotter$plot.type == "colorBySample") return(paste0("\n",diff.text," by Sample"));
  if(plotter$plot.type == "colorByX") return(paste0("\n",diff.text," by ",plotter$title.highlight.name));
  if(plotter$plot.type == "colorByXhighlightByY"){
    #return(paste0("\n",plotter$title.annotations[["highlight"]]," Highlighted, ",diff.text," by ",plotter$title.annotations[["color.by.title.name"]]));
    return(paste0("\n",plotter$title.annotations[["highlight.by.title.name"]],"=",plotter$title.annotations[["highlight"]]," Highlighted, ",diff.text," by ",plotter$title.annotations[["color.by.title.name"]]));
  }
  return("");
}

internal.plot.read.labels <- function(plotter, pos, r2.buffer, max.cycle.ct){
  if(pos == "top"){
    y.pos <- list(par("usr")[4], par("usr")[4]);
    adj <- list(c(0.5,1.1),c(0.5,1.1));
  } else if(pos == "bottom"){
    y.pos <- list(par("usr")[3], par("usr")[3]);
    adj <- list(c(0.5,-0.1),c(0.5,-0.1));
  } else {
    stop("Unrecognized position for read labels!");
  }
  x.pos <- list(max.cycle.ct / 2, r2.buffer + max.cycle.ct + (max.cycle.ct / 2));
  text(x.pos[[1]],y.pos[[1]],label="Read 1", adj = adj[[1]]);
  text(x.pos[[2]],y.pos[[2]],label="Read 2", adj = adj[[2]]);
}












