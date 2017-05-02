

DEFAULTDEBUGMODE <- TRUE;




######################################################################################################################################################################
########### New plots:

#FASTQ PLOTS:

#makePlot.overlapMismatch.perRead
######################################################################################################################################################################

#DEPRECATED?
makePlot.overlapMismatch.perRead <- function(plotter,plot.rates = TRUE,log.y = TRUE,noIndel = TRUE,
                                                  draw.decade.lines = c(TRUE,TRUE),
                                                  xlim = NULL,pct.cutoff=0.95,
                                                  singleEndMode = plotter$res@singleEnd,debugMode = DEFAULTDEBUGMODE,
                                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                                  plot = TRUE,
                                                  ...){
  res <- plotter$res
  plot.name <- "Num Mismatches per Read";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(res@calc.data[["overlapMismatch.byScore.min"]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        
        
        data.list = res@calc.data[["mismatchSizeRatesCumulative"]]
        x.name = "BASES_MISMATCHED"
        xlab = "# Mismatching Base-Pairs"
        #xlim <- c(max(plotter$res@decoder$cycle.CT),1);
        pre.plot.func <- function(){}
        if(plot.rates){
          y.name= "RATE"
          ylab = paste0("% ",readLabel);
        } else {
          y.name= "CT"
          ylab = paste0("# ",readLabel);
        }
        if(noIndel){
          y.name = paste0(y.name,"_NOINDEL");
        }
        if(log.y){
          y.name = paste0(y.name,"_LOG");
          
          pre.plot.func <- function(){
            decades <- seq(floor(par("usr")[3]),ceiling(par("usr")[4]),by=1)
            if(draw.decade.lines[2]) abline(h=decades,col="gray",lty=3);
            decades <- limited.pretty(c(par("usr")[1],1),n = 10);
            if(draw.decade.lines[1]) abline(v=decades,col="gray",lty=3);
          }
        } else {
          pre.plot.func <- function(){
          #  qual.decades <- seq(0,ceiling(par("usr")[4]),by=10)
          #  if(draw.decade.lines) abline(v=qual.decades,col="gray",lty=3);
          }
        }
        #message("y.name=",y.name);
        #message("x.name=",x.name);
        
        if(is.null(xlim)){
          xlim.min <- 1;
          xlim.max <- max(sapply(data.list,function(df){
            xvals <- rev(df[[x.name]]);
            cts <- rev(df[[y.name]]);
            keep <- is.finite(cts) & (! is.na(cts));
            cts <- cts[keep]
            cumsum <- cumsum(cts) / sum(cts);
            limit.val <- xvals[keep][which(cumsum >= pct.cutoff)[1]];
            limit.val
          }));
          xlim <- c(xlim.max,xlim.min);
        }
        
        makePlot.generic(
                 plot.name=plot.name,
                 data.list = data.list,
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,pre.plot.func=pre.plot.func,y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 xaxis.pretty.n = 10,
                 ...);
        title(xlab = xlab);
        title(ylab = ylab);
        
        if(log.y){
          miniticks <- log10(getAllOrdersOfTen(10 ^ par("usr")[3],10^par("usr")[4], skipBase=FALSE));
          axis(2,at=miniticks,labels=FALSE,lwd = -1, lwd.ticks = 1,tcl=-0.25, ...);
        }
        
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
}

makePlot.overlapMismatch.size <- function(plotter,plot.rates = TRUE,log.y = TRUE,noIndel = TRUE,
                                                  draw.decade.lines = c(TRUE,TRUE),
                                                  xlim = NULL,pct.cutoff=0.95,
                                                  singleEndMode = plotter$res@singleEnd,debugMode = DEFAULTDEBUGMODE,
                                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                                  plot = TRUE,
                                                  ...){
  res <- plotter$res
  plot.name <- "Overlap Mismatch Size Frequency";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(res@calc.data[["overlapMismatch.byScore.min"]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        
        
        data.list = res@calc.data[["mismatchSizeRates"]]
        x.name = "BASES_MISMATCHED"
        xlab = "# Mismatching Base-Pairs"
        #xlim <- c(1,max(plotter$res@decoder$cycle.CT));
        pre.plot.func <- function(){}
        if(plot.rates){
          y.name= "RATE"
          ylab = paste0("% ",readLabel);
        } else {
          y.name= "CT"
          ylab = paste0("# ",readLabel);
        }
        if(noIndel){
          y.name = paste0(y.name,"_NOINDEL");
        }
        if(log.y){
          y.name = paste0(y.name,"_LOG");
          
          pre.plot.func <- function(){
            decades <- seq(floor(par("usr")[3]),ceiling(par("usr")[4]),by=1)
            if(draw.decade.lines[2]) abline(h=decades,col="gray",lty=3);
            decades <- limited.pretty(c(1,par("usr")[2]),n = 10);
            if(draw.decade.lines[1]) abline(v=decades,col="gray",lty=3);
          }
        } else {
          pre.plot.func <- function(){
          #  qual.decades <- seq(0,ceiling(par("usr")[4]),by=10)
          #  if(draw.decade.lines) abline(v=qual.decades,col="gray",lty=3);
          }
        }
        #message("y.name=",y.name);
        #message("x.name=",x.name);
        
        if(is.null(xlim)){
          xlim.min <- 1;
          xlim.max <- max(sapply(data.list,function(df){
            xvals <- rev(df[[x.name]]);
            cts <- df[[y.name]];
            keep <- is.finite(cts) & (! is.na(cts))
            cts <- cts[keep]
            cumsum <- cumsum(cts) / sum(cts);
            xvals[keep][which(cumsum >= pct.cutoff)[1]];
          }));
          xlim <- c(xlim.min,xlim.max);
        }
        
        makePlot.generic(
                 plot.name=plot.name,
                 data.list = data.list,
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,pre.plot.func=pre.plot.func,y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 xaxis.pretty.n = 10,
                 ...);
        title(xlab = xlab);
        title(ylab = ylab);
        if(log.y){
          miniticks <- log10(getAllOrdersOfTen(10 ^ par("usr")[3],10^par("usr")[4], skipBase=FALSE));
          axis(2,at=miniticks,labels=FALSE,lwd = -1, lwd.ticks = 1,tcl=-0.25, ...);
        }
        
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
}

makePlot.overlapMismatch.byQual.min <- function(plotter,plot.rates = TRUE,log.y = TRUE,noIndel = TRUE,
                                                draw.decade.lines = TRUE,
                                                singleEndMode = plotter$res@singleEnd,debugMode = DEFAULTDEBUGMODE,
                                                rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                                plot = TRUE,
                                                ...){
  res <- plotter$res
  plot.name <- "Overlapping Mismatch by Min Qual";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(res@calc.data[["overlapMismatch.byScore.min"]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        maxQual <- max(res@calc.data[["overlapMismatch.byScore.min"]][[1]]$minQual)

        if(plot.rates){
          data.list = res@calc.data[["overlapMismatch.byScore.min"]]
          x.name = "minQual"
          y.name= "MismatchRate"
          ylab = "% Mismatch"
          xlab = "Minimum(R1,R2) Quality Score"
          xlim <- c(0,maxQual);
          pre.plot.func <- function(){
            qual.decades <- seq(0,maxQual,by=10)
            if(draw.decade.lines) abline(v=qual.decades,col="gray",lty=3);
          }
        } else {
          data.list = res@calc.data[["overlapMismatch.byScore.min"]]
          x.name = "minQual"
          y.name= "MismatchCt"
          ylab = "# Reads Mismatch"
          xlab = "Minimum(R1,R2) Quality Score"
          xlim <- c(0,maxQual);
          pre.plot.func <- function(){
            qual.decades <- seq(0,maxQual,by=10)
            if(draw.decade.lines) abline(v=qual.decades,col="gray",lty=3);
          }
        }
        if(log.y){
          y.name = paste0(y.name,"Log");
        }
        if(noIndel){
          y.name = paste0(y.name,"NoIndel");
        }
        #message("y.name=",y.name);
        #message("x.name=",x.name);
        
        
        makePlot.generic(
                 plot.name=plot.name,
                 data.list = data.list,
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,pre.plot.func=pre.plot.func,y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        title(xlab = xlab);
        title(ylab = ylab);
        if(log.y){
          miniticks <- log10(getAllOrdersOfTen(10 ^ par("usr")[3],10^par("usr")[4], skipBase=FALSE));
          axis(2,at=miniticks,labels=FALSE,lwd = -1, lwd.ticks = 1,tcl=-0.25, ...);
        }
        
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
}


#######################################################

makePlot.targetDistribution <- function(plotter,plot.rates = TRUE, plot.hist = TRUE, log.y = TRUE, 
                                      singleEndMode = plotter$res@singleEnd,
                                      byPair = !singleEndMode, 
                                      rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                      debugMode = DEFAULTDEBUGMODE,
                                      plot = TRUE,
                                      ...){
  res <- plotter$res
  plot.name <- "Target Coverage";
  readLabel <- if(singleEndMode || (! byPair)){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(res@calc.data[["onTargetCoverage"]]));
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        if(plot.hist){
          if(singleEndMode || (! byPair)){
            data.list = res@calc.data[["onTargetCoverageBinsByRead"]]
            x.name = "READBIN"
            if(log.y) y.name = "logReadCt"
            else y.name = "readCt"
            ylab = "# Target Regions"
            xlim <- c(0,res@calc.data[["maxPairCoverage"]])
          } else {
            data.list = res@calc.data[["onTargetCoverageBinsByPair"]]
            x.name = "PAIRBIN"
            if(log.y) y.name = "logPairCt"
            else y.name = "pairCt"
            ylab = "# Target Regions"
            xlim <- c(0,res@calc.data[["maxPairCoverage"]])
          }
          message("!")
          makePlot.generic(
                 plot.name=plot.name,
                 data.list = data.list,
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        } else {
          xlim <- c(0,100)
          data.list = res@calc.data[["onTargetCoverage"]]
          if(singleEndMode || (! byPair)){
            if(log.y) y.name = "logReadDepth"
            else y.name = "readDepth"
            ylab = "Read Coverage Depth"
            x.name = "quantile"
          } else {
            if(log.y) y.name = "logPairDepth"
            else y.name = "pairDepth"
            ylab = "Read-Pair Coverage Depth"
          }
          
          makePlot.generic(
                 plot.name=plot.name,
                 data.list = data.list,
                 plotter=plotter,
                 x.name = "quantile",y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        }
        
        
        

        title(xlab = "Target Percentile");
        
        title(ylab = ylab);
        
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
}

###########
makePlot.runTimePerformance <- function(plotter,  debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  plot.name <- "QoRTs Runtime Performance";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  
  isPlottable <- (! is.null(plotter$res@qc.data[["summary"]])) && all(c("BENCHMARK_MinutesOnSamIteration","BENCHMARK_MinutesPerMillionReads","BENCHMARK_MinutesPerMillionGoodReads") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD)
  if(! plot) return(isPlottable)
    
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if( ! isPlottable ){
      message("WARNING: Skipping ",plot.name," plotting. Data not found!\n");
      blank.plot(c(plot.name,"Data Not Found","Skipping..."));
    } else {
      oldmar <- par("mar");

         tryCatch({
           #par(mar = c(5,4,4,4) + 0.1);
           old.mar <- par("mar");
           curr.mar <- par("mar");
           curr.mar[4] <- curr.mar[2];
           par(mar = curr.mar);
           
           tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                        x.names = c("BENCHMARK_MinutesOnSamIteration"),
                                        x.titles = c("Runtime"),
                                        norm.by = NULL,
                                        offset.horiz = 0.5,
                                        horiz.offsets = plotter$lanebam.params$horiz.offsets);

           tf.list.2 <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                          x.names = c("BENCHMARK_MinutesPerMillionReads","BENCHMARK_MinutesPerMillionGoodReads"),
                                          x.titles = c("Runtime\nPer Million","Runtime\nPer Million\nPF"),
                                          norm.by = NULL,
                                          offset.horiz = 0.5,
                                          horiz.offsets = plotter$lanebam.params$horiz.offsets);
         }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
         
         makePlot.generic.points(plot.name,tf.list,plotter, leave.blank.cols = 3, label.y = T,...);
         makePlot.generic.points.right(plot.name,tf.list.2, plotter, section.offset = 2, ...);
         
         usr <- par("usr");
         yrange <- abs(usr[4] - usr[3]);
         rect(1.75,usr[3] - yrange*0.02,2.25, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
         rect(1.75,usr[3] - yrange     ,2.25, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);
         
         #text(2,usr[4],"On-Target Fraction",adj=c(0.5,1.1), ...);
         #text(6,usr[4],"Coverage Depth\nPer Million",adj=c(0.5,1.1), ...);
         
         axis(4, at = (yrange / 2) + usr[3], labels=c("Minutes per Million"),col="transparent",line=1, ...);
         axis(2, at = (yrange / 2) + usr[3], labels=c("Minutes"),col="transparent", line=1.5,...);
         
         tcl.wd <- par("cxy")[1] * 0.75
         hours.at <- seq(0,par("usr")[4],by=60);
         segments(x0=1.75,x1=tcl.wd+1.75,y0=hours.at,...)
         text(x=tcl.wd+1.75,y=hours.at, labels=paste0(" ",hours.at / 60,":00"), adj=c(0,0.5),...)
         segments(x0=1.75,x1=(tcl.wd/2)+1.75,y0=hours.at+30,...)
         
         par(mar = old.mar);
      
      #if(y.counts.in.millions){
      #  usr <- par("usr");
      #  pretty.y <- pretty(c(usr[3],usr[4]))
      #  pretty.y.labels <- pretty.y;
      #  pretty.y.labels <- paste0(pretty.y / 1000000, "M");
      #  axis(2,at=pretty.y,labels=pretty.y.labels,las=1, lwd = -1, lwd.ticks = 1, ...);
      #} else {
      #   axis(2, lwd = -1, lwd.ticks = 1, ...)
      #}
      lines(x=c(par("usr")[1],2.75),y=c(0,0),col="grey",lty=3, ...);
        
      internal.plot.main.title(plot.name, plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      par(mar = oldmar);
    }
  })
}

###########
makePlot.onTarget.rates <- function(plotter,  debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  plot.name <- "On-Target Stats";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(plotter$res@qc.data[["summary"]])) && all(c("OnTargetCount","OffTargetCount","OnTargetReadBases",
                                                                     "OnTargetReadPairBases","AvgTargetReadDepth",
                                                                     "AvgTargetReadPairDepth","TargetReadDepth","TargetReadPairDepth",
                                                                     "OnTargetFraction") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD)
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if( ! isPlottable ){
      message("WARNING: Skipping ",plot.name," plotting. Data not found!\n");
      blank.plot(c(plot.name,"Data Not Found","Skipping..."));
    } else {
      oldmar <- par("mar");

         tryCatch({
           #par(mar = c(5,4,4,4) + 0.1);
           old.mar <- par("mar");
           curr.mar <- par("mar");
           curr.mar[4] <- curr.mar[2];
           par(mar = curr.mar);
           
           tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                        x.names = c("fractionOfPairsOnTarget","fractionOfPairsOffTarget"),
                                        x.titles = c("% Read Pairs\nOn-Target","% Read Pairs\nOff-Target"),
                                        norm.by = NULL,
                                        offset.horiz = 0.5,
                                        horiz.offsets = plotter$lanebam.params$horiz.offsets);

           tf.list.2 <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                          x.names = c("readDepthPerMillion","pairDepthPerMillion"),
                                          x.titles = c("Read Depth\nPer Million","Read-pair\nDepth\nPer Million"),
                                          norm.by = NULL,
                                          offset.horiz = 0.5,
                                          horiz.offsets = plotter$lanebam.params$horiz.offsets);
         }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
         
         makePlot.generic.points(plot.name,tf.list,plotter, leave.blank.cols = 3, label.y = T,ylim=c(0,1),...);
         makePlot.generic.points.right(plot.name,tf.list.2, plotter, section.offset = 3, ...);
         
         usr <- par("usr");
         yrange <- abs(usr[4] - usr[3]);
         rect(2.75,usr[3] - yrange*0.02,3.25, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
         rect(2.75,usr[3] - yrange     ,3.25, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);
         
         #text(2,usr[4],"On-Target Fraction",adj=c(0.5,1.1), ...);
         #text(6,usr[4],"Coverage Depth\nPer Million",adj=c(0.5,1.1), ...);
         
         axis(4, at = (yrange / 2) + usr[3], labels=c("Depth"),col="transparent", line=1, ...);
         axis(2, at = (yrange / 2) + usr[3], labels=c("Read Counts"),col="transparent", line=1.5, ...);
         par(mar = old.mar);
      
      
      #if(y.counts.in.millions){
      #  usr <- par("usr");
      #  pretty.y <- pretty(c(usr[3],usr[4]))
      #  pretty.y.labels <- pretty.y;
      #  pretty.y.labels <- paste0(pretty.y / 1000000, "M");
      #  axis(2,at=pretty.y,labels=pretty.y.labels,las=1, lwd = -1, lwd.ticks = 1, ...);
      #} else {
      #   axis(2, lwd = -1, lwd.ticks = 1, ...)
      #}
      lines(x=c(par("usr")[1],2.75),y=c(0,0),col="grey",lty=3, ...);
        
      internal.plot.main.title(plot.name, plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      par(mar = oldmar);
    }
  })
}



makePlot.onTarget.counts <- function(plotter, y.counts.in.millions = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  plot.name <- "On-Target Count Stats";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (!is.null(plotter$res@qc.data[["summary"]])) && all(c("OnTargetCount","OffTargetCount","OnTargetReadBases",
                                                           "OnTargetReadPairBases","AvgTargetReadDepth",
                                                           "AvgTargetReadPairDepth","TargetReadDepth","TargetReadPairDepth",
                                                           "OnTargetFraction") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD)
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  

  
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if( ! isPlottable ){
      message("WARNING: Skipping ",plot.name," plotting. Data not found!\n");
      blank.plot(c(plot.name,"Data Not Found","Skipping..."));
    } else {
      oldmar <- par("mar");

         tryCatch({
           #par(mar = c(5,4,4,4) + 0.1);
           old.mar <- par("mar");
           curr.mar <- par("mar");
           curr.mar[4] <- curr.mar[2];
           par(mar = curr.mar);
           
           tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                        x.names = c("READ_PAIR_OK","OnTargetCount","OffTargetCount"),
                                        x.titles = c("Mapped\nRead\nPairs","Read-Pairs\nOn Target","Read-Pairs\nOff Target"),
                                        norm.by = NULL,
                                        offset.horiz = 0.5,
                                        horiz.offsets = plotter$lanebam.params$horiz.offsets);

           tf.list.2 <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                          x.names = c("TargetReadDepth","TargetReadPairDepth"),
                                          x.titles = c("Read\nDepth","Read-pair\nDepth"),
                                          norm.by = NULL,
                                          offset.horiz = 0.5,
                                          horiz.offsets = plotter$lanebam.params$horiz.offsets);
         }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
         makePlot.generic.points(plot.name,tf.list,plotter, leave.blank.cols = 4, label.y = F,...);
         makePlot.generic.points.right(plot.name,tf.list.2, plotter, section.offset = 4.5, ...);
         usr <- par("usr");
         yrange <- abs(usr[4] - usr[3]);
         rect(4,usr[3] - yrange*0.02,4.5, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
         rect(4,usr[3] - yrange     ,4.5, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);
         
         text(2,usr[4],"Read Counts",adj=c(0.5,1.1), ...);
         text(6,usr[4],"Coverage Depth",adj=c(0.5,1.1), ...);
         
         axis(4, at = (yrange / 2) + usr[3], labels=c("Depth"),col="transparent", line=1, ...);
         axis(2, at = (yrange / 2) + usr[3], labels=c("Read Counts"),col="transparent", line=1.5, ...);
         par(mar = old.mar);
      
      if(y.counts.in.millions){
        usr <- par("usr");
        pretty.y <- pretty(c(usr[3],usr[4]))
        pretty.y.labels <- pretty.y;
        pretty.y.labels <- paste0(pretty.y / 1000000, "M");
        axis(2,at=pretty.y,labels=pretty.y.labels,las=1, lwd = -1, lwd.ticks = 1, ...);
      } else {
         axis(2, lwd = -1, lwd.ticks = 1, ...)
      }
      abline(h=0,col="grey",lty=3, ...);
        
      internal.plot.main.title(plot.name, plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      par(mar = oldmar);
    }
  })
}

makePlot.overlapMismatch.byQual.avg <- function(plotter,plot.rates = TRUE,log.y = TRUE,noIndel = TRUE,
                                              draw.decade.lines = TRUE,
                                              singleEndMode = plotter$res@singleEnd,
                                              debugMode = DEFAULTDEBUGMODE,
                                              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                              plot = TRUE,
                                              ...){
  res <- plotter$res                                              
  byTheoreticalErrorRate = FALSE 
  plot.name <- "Overlapping Mismatch by Avg Qual";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(res@calc.data[["overlapMismatch.byScore.avg"]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        maxQual <- max(res@calc.data[["overlapMismatch.byScore.avg"]][[1]]$avgQual)

        if(byTheoreticalErrorRate){
          data.list = res@calc.data[["overlapMismatch.byScore.byErr"]]
          x.name = "errorRate"
          y.name= "ObservedPhred"
          ylab = "Estimated Phred Score (Observed)"
          xlab="Phred Score"
          xlim <- c(0,max(res@calc.data[["overlapMismatch.byScore.byErr"]][[1]]$errorRate));
          
          pre.plot.func <- function(){
            abline(a=0,b=1,col="gray",lty=3);
          }
          
        } else if(plot.rates){
          data.list = res@calc.data[["overlapMismatch.byScore.avg"]]
          x.name = "avgQual"
          y.name= "MismatchRate"
          ylab = "% Mismatch"
          xlab = "Avg Quality Score"
          xlim <- c(0,maxQual);
          pre.plot.func <- function(){
            qual.decades <- seq(0,maxQual,by=10)
            if(draw.decade.lines) abline(v=qual.decades,col="gray",lty=3);
          }
        } else {
          data.list = res@calc.data[["overlapMismatch.byScore.avg"]]
          x.name = "avgQual"
          y.name= "MismatchCt"
          ylab = "# Reads Mismatch"
          xlab = "Avg Quality Score"
          xlim <- c(0,maxQual);
          pre.plot.func <- function(){
            qual.decades <- seq(0,maxQual,by=10)
            if(draw.decade.lines) abline(v=qual.decades,col="gray",lty=3);
          }
        }
        if(log.y){
          y.name = paste0(y.name,"Log");
        }
        if(noIndel){
          y.name = paste0(y.name,"NoIndel");
        }
        
        
        makePlot.generic(
                 plot.name=plot.name,
                 data.list = data.list,
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,pre.plot.func=pre.plot.func,y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        title(xlab = xlab);
        title(ylab = ylab);
        if(log.y){
          miniticks <- log10(getAllOrdersOfTen(10 ^ par("usr")[3],10^par("usr")[4], skipBase=FALSE));
          axis(2,at=miniticks,labels=FALSE,lwd = -1, lwd.ticks = 1,tcl=-0.25, ...);
        }
        
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
  
}

makePlot.overlapMismatch.byQual.read <- function(plotter,plot.rates = TRUE,noIndel = TRUE,log.y = TRUE,
                                              draw.decade.lines = TRUE,
                                              singleEndMode = plotter$res@singleEnd,debugMode = DEFAULTDEBUGMODE,
                                              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                              r2.buffer = NULL,plot = TRUE,
                                              ...){
  res <- plotter$res
  plot.name <- "Overlap Mismatch by Read Qual";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(res@calc.data[["overlapMismatch.byScore.avg"]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        x.name = "qual"
        maxQual <- max(res@calc.data[["overlapMismatch.byScore.R1"]][[1]][[x.name]])
        if(plot.rates){
          y.name= "MismatchRate"
          ylab = "% Mismatch"
        } else {
          y.name= "MismatchCt"
          ylab = "# Reads Mismatch"
        }
        xlab = "Quality Score"
        xlim <- c(0,maxQual);
        pre.plot.func <- function(){
          qual.decades <- seq(0,maxQual,by=10)
          if(draw.decade.lines) abline(v=qual.decades,col="gray",lty=3);
        }
        
        if(log.y){
          y.name <- paste0(y.name,"Log");
        }
        if(noIndel){
          y.name = paste0(y.name,"NoIndel");
        }
        
        xlim <- c(0,maxQual);
        if(is.null(r2.buffer)) r2.buffer <- (xlim[2]-xlim[1])/10;
        
        if(singleEndMode){
        makePlot.generic(
                 plot.name=plot.name,
                 data.list = res@calc.data[["overlapMismatch.byScore.R1"]],
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        } else {
        makePlot.generic.pair(
                 plot.name=plot.name,
                 data.list.r1 = res@calc.data[["overlapMismatch.byScore.R1"]], 
                 data.list.r2 = res@calc.data[["overlapMismatch.byScore.R2"]],
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,r2.buffer=r2.buffer,y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        }
        if(! singleEndMode){
          internal.plot.read.labels(plotter,"top",r2.buffer,xlim[2]);
        }
        title(xlab = xlab);
        title(ylab = ylab);
        
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
  
}

makePlot.referenceMismatch.byScore <- function(plotter,plot.rates = TRUE,
                                                draw.decade.lines = c(FALSE,TRUE),
                                                log.y = TRUE,
                                                singleEndMode = plotter$res@singleEnd,debugMode = DEFAULTDEBUGMODE,
                                                rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                                r2.buffer = NULL,plot = TRUE,
                                                ...){
  res <- plotter$res
  plot.name <- "Reference Mismatch by Read Qual";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  data.list <- res@calc.data[["referenceMismatch.byScore"]];
  isPlottable <- (! is.null(data.list))
  
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        x.name = "Score"
        maxQual <- max(data.list[[1]][[x.name]])
        
        if(plot.rates){
          if(log.y){
            y.name = "RATE_R1_LOG";
            y.names = c("RATE_R1_LOG","RATE_R1_LOG");
          } else {
            y.name = "RATE_R1";
            y.names = c("RATE_R1","RATE_R2");
          }
          ylab = "% Mismatch"
        } else {
          if(log.y){
            y.name = "CT_R1_LOG";
            y.names = c("CT_R1_LOG","CT_R2_LOG");
          } else {
            y.name = "MismatchCt_R1";
            y.names = c("MismatchCt_R1","MismatchCt_R2");
          }
          ylab = "# Mismatch"
        }
        xlab = "Quality Score"
        xlim <- c(0,maxQual);
        #pre.plot.func <- function(){
        #    qual.decades <- seq(0,maxQual,by=10)
        #    if(draw.decade.lines) abline(v=qual.decades,col="gray",lty=3);
        #}
        
        xlim <- c(0,maxQual);
        if(is.null(r2.buffer)) r2.buffer <- (xlim[2]-xlim[1])/10;
        
        quals <- data.list[[1]][[x.name]];
        if(length(quals) > 10) quals <- NULL;
        
        pre.plot.func <- function(){
            decades <- seq(floor(par("usr")[3]),ceiling(par("usr")[4]),by=1)
            if(draw.decade.lines[2]) abline(h=decades,col="gray",lty=3);
            if(is.null(quals)) {
              decades <- limited.pretty(c(1,par("usr")[2]),n = 10);
            } else {
              decades <- quals;
            }
            if(draw.decade.lines[1]) abline(v=decades,col="gray",lty=3);
        }
        
        if(singleEndMode){
          makePlot.generic(
                 plot.name=plot.name,
                 data.list = data.list,
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim, y.is.log = log.y,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, pre.plot.func=pre.plot.func,
                 ...);
        } else {
          makePlot.generic.pair(
                 plot.name=plot.name,
                 data.list.r1 = data.list, 
                 data.list.r2 = data.list,
                 plotter=plotter,
                 x.name = x.name,y.name = y.name,y.names=y.names,
                 norm.x = FALSE,avg.y = FALSE, y.is.log = log.y, axis.x.at = quals,
                 plot.type = "lines", xlim=xlim,r2.buffer=r2.buffer,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, pre.plot.func=pre.plot.func,
                 ...);
        }
        
        if(! singleEndMode){
          internal.plot.read.labels(plotter,"top",r2.buffer,xlim[2]);
        }
        title(xlab = xlab);
        title(ylab = ylab);
        
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
  
}


makePlot.overlap.coverage <- function(plotter,plot.rates = TRUE, 
                                      singleEndMode = plotter$res@singleEnd,
                                      
                                      rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                      debugMode = DEFAULTDEBUGMODE,r2.buffer=NULL,
                                      plot = TRUE,
                                      ...){
  res <- plotter$res
  plot.name <- "Overlap Coverage";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(plotter$res@calc.data[["overlapR1"]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        if(plot.rates){
          y.name= "RATE"
          ylab = "% Reads w/ overlap"
        } else {
          y.name= "CT"
          ylab = "# Reads w/ overlap"
        }
        xlim <- c(0,max(plotter$res@decoder$cycle.CT))
        if(is.null(r2.buffer)) r2.buffer <- (xlim[2]-xlim[1])/10;
        makePlot.generic.pair(
                 plot.name=plot.name,
                 data.list.r1 = res@calc.data[["overlapR1"]], 
                 data.list.r2 = res@calc.data[["overlapR2"]],
                 plotter=plotter,
                 x.name = "POS",y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,r2.buffer=r2.buffer,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        if(! singleEndMode){
          internal.plot.read.labels(plotter,"top",r2.buffer,xlim[2]);
        }
        title(xlab = "Read Position");
        
        title(ylab = ylab);
        internal.plot.main.title("Overlapping Coverage", plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
}

makePlot.readLengthDist   <- function(plotter,plot.rates = TRUE, plot.means = TRUE, plot.medians = NULL,
                                      include.full.length = FALSE, cumulative = TRUE,
                                      singleEndMode = plotter$res@singleEnd,
                                      rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                      debugMode = DEFAULTDEBUGMODE,r2.buffer =NULL,
                                      plot = TRUE,
                                      ...){
  res <- plotter$res
  plot.name <- "Read Length Distribution";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- ! is.null(plotter$res@calc.data[["readLenDistR1"]]);
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        if(is.null(plot.medians)) {
            plot.medians <- FALSE;
        } else if(plot.means & plot.medians){
            plot.means <- FALSE;
        }
        if(plot.rates){
          if(cumulative){
            y.name <- "CUMRATE"
            ylab <- "% Reads <= length"
          } else {
            y.name <- "RATE"
            ylab <- "% Reads at given length"
          }
        } else {
          if(cumulative){
            y.name <- "CUMCT"
            ylab <- "# Reads <= length"
          } else {
            y.name <- "CT"
            ylab <- "# Reads at given length"
          }
        }
        if(include.full.length){
          xlim <- c(1,max(plotter$res@decoder$cycle.CT))
        } else {
          xlim <- c(1,max(plotter$res@decoder$cycle.CT)-1)
        }
        if(is.null(r2.buffer)) r2.buffer <- (xlim[2]-xlim[1])/10;
        
        if(singleEndMode){
          makePlot.generic(
                 plot.name=plot.name,
                 data.list = res@calc.data[["readLenDistR1"]], 
                 plotter=plotter, 
                 x.name = "LEN",y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        } else {
          #message("y.name = ",y.name);
          
          if(! include.full.length){
            tf.r1 <- res@calc.data[["readLenDistR1"]];
            tf.r2 <- res@calc.data[["readLenDistR2"]];
            tf.r1 <- lapply(tf.r1,function(tf){tf[tf$LEN <= xlim[2],,drop=FALSE]})
            tf.r2 <- lapply(tf.r2,function(tf){tf[tf$LEN <= xlim[2],,drop=FALSE]})
          } else {
            tf.r1 <- res@calc.data[["readLenDistR1"]];
            tf.r2 <- res@calc.data[["readLenDistR2"]];
          }
          
          makePlot.generic.pair(
                 plot.name=plot.name,
                 data.list.r1 = tf.r1, 
                 data.list.r2 = tf.r2,
                 plotter=plotter, 
                 x.name = "LEN",y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,r2.buffer=r2.buffer,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        }
        title(xlab = "Read Length");
        
        if(! singleEndMode){
          internal.plot.read.labels(plotter,"top",r2.buffer,xlim[2]);
        }
        
        title(ylab = ylab);
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
}

makePlot.referenceMismatch.perRead <- function(plotter,
                                             debugMode = DEFAULTDEBUGMODE, log.y = TRUE,
                                             singleEndMode = plotter$res@singleEnd, plot = TRUE,
                                             ...){
  res <- plotter$res
  plot.name <- "Reference Mismatch Summary";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  
  fieldList <- c("PAIR_CONTAINS_DEL","PAIR_CONTAINS_INS",
                 "READ_CONTAINS_DEL_R1","READ_CONTAINS_INS_R1",
                 "READ_CONTAINS_DEL_R2","READ_CONTAINS_INS_R2",
                 "HAS_REF_BASE_SWAP_R1","HAS_REF_BASE_SWAP_R2")
  fieldTitles <- c("Either\nDEL","Either\nINS","R1\nDEL","R1\nINS","R2\nDEL","R2\nINS","R1\nBase Swap","R2\nBase Swap");
  isPlottable <- (!is.null(plotter$res@qc.data[["summary"]])) && all(fieldList %in% plotter$res@qc.data[["summary"]][[1]]$FIELD)
  
  if(! plot) return(isPlottable)
  
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(! isPlottable){
      message("WARNING: Skipping ",plot.name," plotting. Data not found!\n");
      blank.plot(c(plot.name,"Data Not Found","Skipping..."));
    } else {
      oldmar <- par("mar");
         
         datatable = plotter$res@qc.data[["summary"]]
         count.name = "RATE"

         tryCatch({
           #par(mar = c(5,4,4,4) + 0.1);
           old.mar <- par("mar");
           #curr.mar <- par("mar");
           #curr.mar[4] <- curr.mar[2];
           #par(mar = curr.mar);
           
           tf.list <- generic.points.tf(datatable, 
                                        x.names = fieldList,
                                        x.titles = fieldTitles,
                                        field.name = "FIELD",
                                        count.name = "COUNT",
                                        norm.by = "READ_PAIR_OK",
                                        offset.horiz = 0.5,
                                        horiz.offsets = plotter$lanebam.params$horiz.offsets,
                                        log.y = TRUE);
                                        
         }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
         
         makePlot.generic.points(plot.name,tf.list,plotter, label.y = T,y.is.log = log.y, ...);

         usr <- par("usr");
         yrange <- abs(usr[4] - usr[3]);
         
        #usr <- par("usr");
        #pretty.y <- pretty(c(usr[3],usr[4]))
        #pretty.y.labels <- pretty.y;
        #pretty.y.labels <- paste0(pretty.y, "");
        #axis(2,at=pretty.y,labels=pretty.y.labels,lwd = -1, lwd.ticks = 1,las=1, ...);
        #axis(4,at=pretty.y,labels=FALSE,lwd = -1, lwd.ticks = 1,las=1, ...);
        
        title(ylab = "Proportion of Reads")
      
      #abline(h=0,col="grey",lty=3, ...);
      #title(xlab = "Read -> Reference Mismatches")
      
      internal.plot.main.title(plot.name, plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      par(mar = oldmar);
    }
  })
}

makePlot.referenceMismatch.byBase <- function(plotter, 
                                             y.rate.per.kb = TRUE, draw.vert.dotted.lines = TRUE,
                                             debugMode = DEFAULTDEBUGMODE, 
                                             singleEndMode = plotter$res@singleEnd, plot = TRUE,
                                             ...){
  res <- plotter$res
  plot.name <- "Reference Mismatch Combinations";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(plotter$res@calc.data[["referenceMismatchCombos"]]));
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }

  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(! isPlottable){
      message("WARNING: Skipping ",plot.name," plotting. Data not found!\n");
      blank.plot(c(plot.name,"Data Not Found","Skipping..."));
    } else {
      oldmar <- par("mar");
         
         tableName <- "referenceMismatchCombos"
         count.name = "RATE"

         tryCatch({
           #par(mar = c(5,4,4,4) + 0.1);
           old.mar <- par("mar");
           #curr.mar <- par("mar");
           #curr.mar[4] <- curr.mar[2];
           #par(mar = curr.mar);
           
           tf.list <- generic.points.tf(plotter$res@calc.data[[tableName]], 
                                        x.names = plotter$res@calc.data[[tableName]][[1]]$combo,
                                        x.titles = plotter$res@calc.data[[tableName]][[1]]$combo,
                                        field.name = "combo",
                                        count.name = count.name,
                                        norm.by = NULL,
                                        offset.horiz = 0.5,
                                        horiz.offsets = plotter$lanebam.params$horiz.offsets);
                                        
         }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
         
         makePlot.generic.points(plot.name,tf.list,plotter, label.y = F, ...);

         usr <- par("usr");
         yrange <- abs(usr[4] - usr[3]);
         
      #if(y.rate.per.kb){
        usr <- par("usr");
        pretty.y <- pretty(c(usr[3],usr[4]))
        pretty.y.labels <- pretty.y;
        pretty.y.labels <- paste0(pretty.y, "");
        axis(2,at=pretty.y,labels=pretty.y.labels,lwd = -1, lwd.ticks = 1,las=1, ...);
        axis(4,at=pretty.y,labels=FALSE,lwd = -1, lwd.ticks = 1,las=1, ...);
        
        title(ylab = "% Mismatch (for aligned bases)")
      #} else {
      #  axis(4, lwd = -1, lwd.ticks = 1, ...)
      #}
      
      if(draw.vert.dotted.lines){
        abline(v=c(4.5,8.5),col="gray",lty=3);
      }
      
      abline(h=0,col="grey",lty=3, ...);
      title(xlab = "Read -> Reference Mismatches")
      
      internal.plot.main.title(plot.name, plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      par(mar = oldmar);
    }
  })
}

#formerly makePlot.overlapMismatchCombos
makePlot.overlapMismatch.byBase <- function(plotter, noIndel = TRUE, plot.rates = TRUE,log.y = FALSE, 
                                           y.rate.per.kb = FALSE, debugMode = DEFAULTDEBUGMODE, 
                                           draw.vert.dotted.lines = TRUE,
                                           singleEndMode = plotter$res@singleEnd, plot = TRUE,...){
  res <- plotter$res
  plot.name <- "Overlap Mismatch Combinations";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  tableName <- "overlapMismatchCombos"
  isPlottable <- (! is.null(plotter$res@calc.data[[tableName]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(! isPlottable){
      message("WARNING: Skipping ",plot.name," plotting. Data not found!\n");
      blank.plot(c(plot.name,"Data Not Found","Skipping..."));
    } else {
      oldmar <- par("mar");
         if(plot.rates){
           count.name = "rate"
           ylab <- "% Mismatches (for overlapping bases)"
         } else {
           count.name = "Ct"
           ylab <- "# Mismatches (for overlapping bases)"
         }
         
         if(noIndel){
           count.name = paste0(count.name,"_NOINDEL");
         }
         if(log.y){
           count.name = paste0(count.name,"_LOG");
         }

         tryCatch({
           #par(mar = c(5,4,4,4) + 0.1);
           old.mar <- par("mar");
           #curr.mar <- par("mar");
           #curr.mar[4] <- curr.mar[2];
           #par(mar = curr.mar);
           
           tf.list <- generic.points.tf(plotter$res@calc.data[[tableName]], 
                                        x.names = plotter$res@calc.data[[tableName]][[1]]$combo,
                                        x.titles = plotter$res@calc.data[[tableName]][[1]]$combo,
                                        field.name = "combo",
                                        count.name = count.name,
                                        norm.by = NULL,
                                        offset.horiz = 0.5,
                                        horiz.offsets = plotter$lanebam.params$horiz.offsets);
                                        
         }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
         
         makePlot.generic.points(plot.name,tf.list,plotter, label.y = FALSE, y.is.log=log.y,...);

         usr <- par("usr");
         yrange <- abs(usr[4] - usr[3]);
      
      if(! log.y){
        if(y.rate.per.kb){
          usr <- par("usr");
          pretty.y <- pretty(c(usr[3],usr[4]))
          pretty.y.labels <- pretty.y;
          pretty.y.labels <- paste0(pretty.y * 10, "");
          axis(2,at=pretty.y,labels=pretty.y.labels,lwd = -1, lwd.ticks = 1,las=1, ...);
          axis(4,at=pretty.y,labels=FALSE,lwd = -1, lwd.ticks = 1,las=1, ...);

          title(ylab = paste0(ylab," per kBP overlap"))
        } else {
          pretty.y <- pretty(c(usr[3],usr[4]))
          axis(2, lwd = -1, lwd.ticks = 1,las=1, ...);
          axis(4, labels=FALSE, lwd = -1, lwd.ticks = 1, ...)

          title(ylab = ylab)
        }
      } else {
        title(ylab = ylab)
      }
      abline(h=0,col="grey",lty=3, ...);
      title(xlab = "Read1 -> Read2 Mismatching Bases")
      
      if(draw.vert.dotted.lines){
        abline(v=c(4.5,8.5),col="gray",lty=3);
      }
      
      internal.plot.main.title(plot.name, plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      par(mar = oldmar);
    }
  })
}

#TODO:
#makePlot.overlapMismatch.byBase.atScore 
#makePlot.referenceMismatch.byBase.atScore 

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
refMM.comboStrings <- sapply(mismatchCombos,function(mmc){ paste0(mmc[[1]],"->",mmc[[2]]) })
ovrMM.comboStrings <- sapply(mismatchCombos,function(mmc){ paste0(mmc[[1]],"/",mmc[[2]]) })

makePlot.overlapMismatch.byBase.atScore <- function(plotter,atScore = 41,overlapScoreMethod = c("pair","min"),
                                        plot.rates = TRUE, noIndel = TRUE, log.y = FALSE,
                                        singleEndMode = plotter$res@singleEnd,
                                        rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                        debugMode = DEFAULTDEBUGMODE,#r2.buffer = NULL,
                                        plot = TRUE,...){
  res <- plotter$res
  overlapScoreMethod <- match.arg(overlapScoreMethod);
  
  scoreMethodSymbol <- if(overlapScoreMethod == "pair") "==" else "<="
  
  plot.name <- paste0("Overlap Mismatch Pairs At Phred ",scoreMethodSymbol," ",atScore);
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  tableName <- "overlapMismatch.byScoreAndBP"
  
  isPlottable <- (! is.null(plotter$res@qc.data[[tableName]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(plot.rates){
          y.name= "RATE"
          ylab = "% of overlap mismatches"
        } else {
          y.name= "CT"
          ylab = "# of overlap mismatches"
        }
        y.name <- if(overlapScoreMethod == "pair") paste0(y.name,"_PAIR") else paste0(y.name,"_MIN");
        y.name <- if(noIndel) paste0(y.name,"_NOINDEL") else y.name;
        y.name <- if(log.y) paste0(y.name,"_LOG") else y.name;

        data.list <- lapply(res@calc.data[[tableName]],function(dl){
          dl[dl$score == atScore,];
        })
        tf.list <- generic.points.tf(data.list, 
                                     x.names = ovrMM.comboStrings,
                                     x.titles = ovrMM.comboStrings,
                                     field.name = "combo",
                                     count.name = y.name,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});

      makePlot.generic.points(plot.name=plot.name, tf.list, plotter, plot.type = "points", 
                                    ylim = NULL, leave.blank.cols = 0, label.y = TRUE, 
                                    cex.x.axis = NULL, pre.plot.func = NULL, draw.lines = FALSE,
                                    y.is.log = log.y,
                                    ...)
      if(log.y){
        logRangeY <- par("usr")[3:4];
        if(abs(logRangeY[[1]] - logRangeY[[2]]) < 2){
          miniticks <- log10(getAllOrdersOfTen(10 ^ par("usr")[3],10^par("usr")[4], skipBase=FALSE));
          axis(2,at=miniticks,labels=10 ^ miniticks,lwd = -1, lwd.ticks = 0.5,tcl=-0.25,las=1, ...);
        } else {
          plot.miniticks(2,log.y = log.y,...);
        }
      } else {
        plot.miniticks(2,log.y = log.y,...);
      }
      
      abline(v=c(4.5,8.5),col="gray",lty=3);
      title(ylab = ylab);
      title(xlab = "Read1/Read2 Base Swap")
      internal.plot.main.title(plot.name, plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.referenceMismatch.byBase.atScore <- function(plotter,atScore = 41,
                                        forRead = c("BOTH","R1","R2"),
                                        plot.rates = TRUE,log.y = FALSE,
                                        singleEndMode = plotter$res@singleEnd,
                                        rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                        debugMode = DEFAULTDEBUGMODE,
                                        plot = TRUE,...){
  res <- plotter$res
  forRead <- match.arg(forRead);
  
  scoreMethodSymbol <- "=="
  
  plot.name <- paste0(forRead," Ref Mismatches At Phred ",scoreMethodSymbol," ",atScore);
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  tableName <- "referenceMismatch.byScoreAndBP"
  data.table <- plotter$res@calc.data[[tableName]]
  
  isPlottable <- (! is.null(data.table));
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(plot.rates){
          y.name= "RATE"
          ylab = "% of ref mismatches"
        } else {
          y.name= "CT"
          ylab = "# of ref mismatches"
        }
        if(forRead != "BOTH") y.name <- paste0(y.name,"_",forRead);

        y.name <- if(log.y) paste0(y.name,"_LOG") else y.name;
        

        data.list <- lapply(res@calc.data[[tableName]],function(dl){
          dl <- dl[dl$Score == atScore,];
          dl$combo <- paste0(dl$refBase,"->",dl$readBase);
          dl
        })
        tf.list <- generic.points.tf(data.list, 
                                     x.names = refMM.comboStrings,
                                     x.titles = refMM.comboStrings,
                                     field.name = "combo",
                                     count.name = y.name,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});

      makePlot.generic.points(plot.name=plot.name, tf.list, plotter, plot.type = "points", 
                                    ylim = NULL, leave.blank.cols = 0, label.y = TRUE, 
                                    cex.x.axis = NULL, pre.plot.func = NULL, draw.lines = FALSE,
                                    y.is.log = log.y,
                                    ...)
      if(log.y){
        logRangeY <- par("usr")[3:4];
        if(abs(logRangeY[[1]] - logRangeY[[2]]) < 2){
          miniticks <- log10(getAllOrdersOfTen(10 ^ par("usr")[3],10^par("usr")[4], skipBase=FALSE));
          axis(2,at=miniticks,labels=10 ^ miniticks,lwd = -1, lwd.ticks = 0.5,tcl=-0.25,las=1, ...);
        } else {
          plot.miniticks(2,log.y = log.y,...);
        }
      } else {
        plot.miniticks(2,log.y = log.y,...);
      }
      abline(v=c(4.5,8.5),col="gray",lty=3);
      title(ylab = ylab);
      title(xlab = "Ref->Read Base Swap (Read Strand)")
      internal.plot.main.title(plot.name, plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

#formerly makePlot.overlapMismatch
makePlot.overlapMismatch.byCycle   <- function(plotter,
                                        plot.rates = TRUE, noIndel = TRUE, log.y = TRUE,
                                        singleEndMode = plotter$res@singleEnd,
                                        rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                        debugMode = DEFAULTDEBUGMODE,r2.buffer = NULL,
                                        plot = TRUE,
                                        ...){
  res <- plotter$res
  plot.name <- "Overlap Mismatch by Read Cycle";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  isPlottable <- (! is.null(res@calc.data[["overlapMismatchR1"]])) && (! singleEndMode);
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else {
        if(plot.rates){
          y.name= "RATE"
          ylab = "% Overlapping Reads Mismatch"
        } else {
          y.name= "CT"
          ylab = "# Overlapping Reads Mismatch"
        }
        if(log.y){
          y.name = paste0(y.name,"_LOG");
        }
        if(noIndel){
          y.name=paste0(y.name,"_NOIND");
        }

        xlim <- c(0,max(plotter$res@decoder$cycle.CT))
        if(is.null(r2.buffer)) r2.buffer <- (xlim[2]-xlim[1])/10;
        makePlot.generic.pair(
                 plot.name=plot.name,
                 data.list.r1 = res@calc.data[["overlapMismatchR1"]], 
                 data.list.r2 = res@calc.data[["overlapMismatchR2"]],
                 plotter=plotter, 
                 x.name = "POS",y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,y.is.log = log.y,r2.buffer=r2.buffer,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        title(xlab = "Read Position");
        
        if(! singleEndMode){
          internal.plot.read.labels(plotter,"top",r2.buffer,xlim[2]);
        }
        
        plot.miniticks(2,log.y = log.y,...);
        
        title(ylab = ylab);
        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      }
    })
}

makePlot.referenceMismatch.byCycle <- function(plotter,
                                 plot.rates = TRUE,log.y = TRUE, ylim = NULL,
                                 singleEndMode = plotter$res@singleEnd,debugMode = DEFAULTDEBUGMODE,
                                 rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                 r2.buffer = NULL,
                                 plot = TRUE,
                                 ...){
  res <- plotter$res
  plot.name <- "Reference Mismatch by Read Cycle";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  isPlottable <- (! is.null(res@calc.data[["referenceMismatchCounts"]]));
  if(! plot) return(isPlottable)
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(!isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        if(plot) blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else if(plot) {
        if(plot.rates){
          y.name = "RATE_R1";
          y.names = c("RATE_R1","RATE_R2");
          ylab = "% Mismatch"
        } else {
          y.name = "CT_R1";
          y.names = c("CT_R1","CT_R2");
          ylab = "# Mismatch"
        }
        if(singleEndMode) y.names=y.name;
        if(log.y){
          y.name <- paste0(y.name,"_LOG");
          y.names <- paste0(y.names,"_LOG");
        }
        xlim <- c(1,max(plotter$res@decoder$cycle.CT))
        pre.plot.func <- function(){}
        if(is.null(r2.buffer)) r2.buffer <- (xlim[2]-xlim[1])/25;
        
        yspan <- bad.to.na(unlist(lapply(res@calc.data[["referenceMismatchCounts"]],function(df){ df[y.names] })));
        if(is.null(ylim)){
          ylim <- c(quantile(yspan,probs=0.025,na.rm=TRUE),max(yspan,na.rm=TRUE))
        }
        
        if(singleEndMode){
          makePlot.generic(
                 plot.name=plot.name,
                 data.list = res@calc.data[["referenceMismatchCounts"]],
                 plotter=plotter,
                 x.name = "POS",y.name = y.name,norm.x = FALSE,avg.y = FALSE,
                 plot.type = "lines", xlim=xlim,pre.plot.func=pre.plot.func, y.is.log = log.y,ylim=ylim,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        } else {
          makePlot.generic.pair(
                 plot.name=plot.name,
                 data.list.r1 = res@calc.data[["referenceMismatchCounts"]],
                 data.list.r2 = res@calc.data[["referenceMismatchCounts"]],
                 plotter=plotter,
                 x.name = "POS",y.name = y.name,y.names = y.names,
                 norm.x = FALSE,avg.y = FALSE,r2.buffer=r2.buffer,
                 plot.type = "lines", xlim=xlim,pre.plot.func=pre.plot.func,y.is.log = log.y,ylim=ylim,
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        }
        
        if(! singleEndMode){
          internal.plot.read.labels(plotter,"top",r2.buffer,xlim[2]);
        }
        title(xlab = "");
        title(ylab = ylab);
        plot.miniticks(2,log.y = log.y,...);

        internal.plot.main.title(plot.name, plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      } else return(TRUE)
    })
  
}


######################################################################################################################################################################

#makePlot.generic.pair <- function(plot.name, data.list.r1, data.list.r2, plotter, x.name, y.name, norm.x, avg.y, xlim , ylim = NULL, vert.offset = 0, horiz.offset = 0, 
#                                  draw.horiz.lines = FALSE, plot.type = "lines", r2.buffer = 10, override.lty = -1, 
#                                  x.is.log = FALSE, y.is.log = FALSE, pre.plot.func = NULL, y.axis.las = 1, zeroBaseX = FALSE, 
#                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
#                                  debugMode = FALSE,
#                                  ...){

makePlot.biotype.rates <- function(plotter, 
                                   plot.rates = TRUE,
                                   count.type = c("all","unambigOnly"),
                                   log.y = TRUE,
                                   return.table = FALSE, 
                                   debugMode = DEFAULTDEBUGMODE,  
                                   singleEndMode = plotter$res@singleEnd,
                                   showTypes = NULL,
                                   plot = TRUE,
                                   ...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["biotype.counts"]]))
  if(! plot) return(isPlottable);

  plot.name <- "Biotype Rates";
  count.type <- match.arg(count.type);
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["biotype.counts"]])){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      if(plot){
        blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      }
    } else {
      tryCatch({
        bt.cts <- plotter$res@qc.data[["biotype.counts"]];
        pre.plot.func <- function(){
          #do nothing!
        }
        
        x.titles <- bt.cts[[1]]$BIOTYPE;
        names(x.titles) <- x.titles;
        
        if(! is.null(showTypes)){
          keepVector <- x.titles %in% showTypes;
          
          bt.cts <- lapply(bt.cts,function(td){
            if(! all(td$BIOTYPE == names(x.titles))){
              message("WARNING WARNING WARNING: Biotype list not identical across all samples!");
              warning("WARNING WARNING WARNING: Biotype list not identical across all samples!");
            }
            td[keepVector,];
          });
          x.titles <- x.titles[keepVector];
        }
        
        if(! plot){
          tf.list <- lapply(1:length(bt.cts), function(i){
            td <- bt.cts[[i]];
            cts <- if(count.type == "all"){ td$TOTAL } else { td$COUNT }
           
            if(! all(td$BIOTYPE == names(x.titles))){
              message("WARNING WARNING WARNING: Biotype list not identical across all samples!");
              warning("WARNING WARNING WARNING: Biotype list not identical across all samples!");
            }
            cts;
          });
          resMatrix <- do.call(cbind.data.frame, tf.list);
          rownames(resMatrix) <- x.titles;
          colnames(resMatrix) <- plotter$res@lanebam.list;
          return(resMatrix);
        }
        
        x.titles <- sapply(x.titles, wordwrap.string, width = 8)
        
        tf.list <- lapply(1:length(bt.cts), function(i){
           td <- bt.cts[[i]];
           curr.summary <- plotter$res@qc.data[["summary"]][[i]];
           
           cts <- if(count.type == "all"){ td$TOTAL } else { td$COUNT }
           if(plot.rates){
             normCt <- as.numeric(curr.summary$COUNT[curr.summary$FIELD == "READ_PAIR_OK"]);
             cts <- (cts / (normCt / 1000000));
           }
           
           if(! all(td$BIOTYPE == names(x.titles))){
             message("WARNING WARNING WARNING: Biotype list not identical across all samples!");
             warning("WARNING WARNING WARNING: Biotype list not identical across all samples!");
           }
           
           df <- data.frame( x = as.numeric(1:nrow(td)),y = as.numeric(cts), x.titles = x.titles, stringsAsFactors=F);
           if(log.y){
             df$y <- log10(df$y);
           }
           df$x <- df$x + plotter$lanebam.params$horiz.offsets[i] * 0.5;
           df;
        });
        names(tf.list) <- names(bt.cts);
        
        if(log.y){
          min.y <- min(sapply(tf.list,function(td){ ifelse(is.infinite(td$y), Inf,td$y) }),na.rm=T);
          max.y <- max(sapply(tf.list,function(td){ ifelse(is.infinite(td$y),-Inf,td$y) }),na.rm=T);
          NINF.VALUE <- min.y - abs(max.y - min.y)*0.03;
          tf.list <- lapply(tf.list, function(td){
            td$y <- ifelse(td$y == -Inf, NINF.VALUE, td$y);
            td;
          });
        }
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      #makePlot.generic.points(plot.name,tf.list,plotter,pre.plot.func=pre.plot.func, label.y = F, ylim = c(min.y,max.y), family="mono",...)
      
      if(plot){
        makePlot.generic.points(plot.name,tf.list,plotter,pre.plot.func=pre.plot.func, label.y = F, family="mono",...)

        if(log.y){
          abline(h=ceiling(min.y):ceiling(max.y), col="gray",lty=3,...);
          abline(h=min.y, col="gray",lty=3,...);
          abline(h=NINF.VALUE, col="gray",lty=3,...);
          draw.logyaxis.stdScalePlot(  ylim.truncate = c(min.y,max.y)  );
          qorts.axis.break(axis=2,breakpos = min.y - abs(min.y-NINF.VALUE)/2, fill = TRUE, cex = 0.75,...);
          qorts.axis.break(axis=4,breakpos = min.y - abs(min.y-NINF.VALUE)/2, fill = TRUE, cex = 0.75,...);
          axis(2,at=NINF.VALUE,labels=0, tcl = -0.5,lwd = -1, lwd.ticks = par("lwd"), las=2,...);
        } else {
          axis(2);
        }
        if(plot.rates){
          title(ylab=paste0(readLabel," per Million"));
          internal.plot.main.title("Coverage for each \"Biotype\"", plotter, ...);
        } else {
          title(ylab=paste0(readLabel));
          internal.plot.main.title("Coverage for each \"Biotype\"", plotter, ...);
        }
      }
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      
      #print("1:");
      #print(tf.list[[1]]);
      #print("2:");
      #print(tf.list[[2]]);

    }
  })
}


######################################################################################################################################################################
########### Specific Plots:






makePlot.qual.pair <- function(plotter, y.name, r2.buffer = NULL, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                               rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000, 
                               plot = TRUE,
                               ...) {
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["quals.r1"]]) )
  if(! plot) return(isPlottable);

  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
   if(y.name == "min"){ plot.name <- "Minimum";
   } else if(y.name == "max"){ plot.name <- "Maximum";
   } else if(y.name == "lowerQuartile"){ plot.name <- "Lower Quartile";
   } else if(y.name == "upperQuartile"){ plot.name <- "Upper Quartile";
   } else if(y.name == "median"){ plot.name <- "Median";
   } else { stop("FATAL ERROR: y.name must be one of: min,max,lowerQuartile,upperQuartile,median"); }

   if(debugMode) message("Starting: ",plot.name," plot."); 
   
   res <- plotter$res;
   defaultPhredMax <- 41;
   phredMax <- max(sapply(res@qc.data[["summary"]], FUN = function(td){
     idx <- which(td$FIELD == "maxLegalPhredScore");
     if(length(idx) == 1){
       as.numeric(as.character(td$COUNT[idx]));
     } else {
       defaultPhredMax;
     }
   }));
   
   plot.name <- paste0(plot.name," Phred Quality Score");
   plotter.error.wrapper(plot.name, plotterFcn = function(){
     if(debugMode){ ts <- timestamp() }
     if(is.null(plotter$res@qc.data[["quals.r1"]]) ){
       message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
       blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
     } else {
       xlim <- c(1,max(plotter$res@decoder$cycle.CT))
       if(is.null(r2.buffer)) r2.buffer <- xlim[2] / 10;
       
       if(singleEndMode){
         makePlot.generic(plot.name, plotter$res@qc.data[["quals.r1"]], 
                               plotter = plotter, 
                               x.name = "readLen",
                               y.name = y.name,
                               norm.x=FALSE,
                               avg.y=FALSE,plot.type = "lines",
                               xlim=xlim,
                               ylim=c(0,phredMax), 
                               vert.offset=1, 
                               draw.horiz.lines = TRUE,
                               override.lty = 1,
                               zeroBaseX = TRUE,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                               ...);
       } else {
         makePlot.generic.pair(plot.name, plotter$res@qc.data[["quals.r1"]], 
                               plotter$res@qc.data[["quals.r2"]], 
                               plotter = plotter, 
                               x.name = "readLen",
                               y.name = y.name,
                               norm.x=FALSE,
                               avg.y=FALSE,plot.type = "lines",
                               xlim=xlim,
                               ylim=c(0,phredMax), 
                               vert.offset=1, 
                               draw.horiz.lines = TRUE,
                               override.lty = 1,
                               r2.buffer =  r2.buffer,
                               zeroBaseX = TRUE,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                               ...);
       }
       
       #text(xlim[2] / 2,-0.5,adj=c(0.5,0), labels = "Read 1");
       #text(xlim[2] + r2.buffer + (xlim[2] / 2),-0.5, labels = "Read 2",adj=c(0.5,0));
       if(! singleEndMode){
         internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
       }
       title(xlab = "Read Cycle");
       title(ylab = paste0(plot.name));
   
       internal.plot.main.title(paste0(plot.name), plotter, ...);
   
       if(debugMode) message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); 
     }
  })
}


makePlot.gc <- function(plotter, plot.medians = NULL, plot.means = TRUE, plotRate = FALSE, 
                        byPair = FALSE, useFQ = FALSE,
                        debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                        rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000, 
                        plot = TRUE,...){
  res <- plotter$res
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  if(byPair & singleEndMode) stop("ERROR with makePlot.gc: byPair and singleEndMode cannot both be TRUE.");
  
  #if(! byPair){
  #  makePlot.gc.byRead(plotter = plotter, plot.medians = plot.medians, plot.means = plot.means, useFQ=useFQ,
  #                     plotRate = plotRate, debugMode = debugMode,rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
  #                     plot=plot, ...);
  #} else {
    data.list <- if(useFQ){ 
      if(byPair){ plotter$res@qc.data[["FQ.gc.byPair"]] 
      } else { plotter$res@qc.data[["FQ.gc.byRead"]] } 
    } else { 
      if(byPair){ plotter$res@qc.data[["gc.byPair"]] 
      } else { plotter$res@qc.data[["gc.byRead"]] }
    }
    
    isPlottable <- ! (is.null(data.list))
    if(! plot) return(isPlottable);
    
    plot.name <- "GC Content";
    if(byPair) plot.name <- paste0(plot.name," per Read-Pair");
    
    if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
      res <- plotter$res;
      if(debugMode){ ts <- timestamp() }
      if(! isPlottable){
        message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
        blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
      } else {
        if(is.null(plot.medians)){
            plot.medians <- FALSE;
        } else if(plot.means & plot.medians){
            plot.means <- FALSE;
        }
        makePlot.generic(plot.name,data.list, plotter, 
                 x.name = "NUM_BASES_GC",y.name = "CT",norm.x = plotRate,avg.y = TRUE,plot.type = "lines", plot.means = plot.means, plot.medians = plot.medians, 
                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                 ...);
        if(plotRate){ title(xlab = "% G/C");
        } else { title(xlab = "# G/C"); }
        
        title(ylab = "Frequency");
        #abline(v=50);
        internal.plot.main.title("GC Content", plotter, ...);
        if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
      }
    })
  #}
}

makePlot.gc.byRead <- function(plotter, plot.medians = NULL, plot.means = TRUE, plotRate = FALSE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                               rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                               plot = TRUE,...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["gc.byRead"]]) )
  if(! plot) return(isPlottable);

  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plot.name <- "GC Content";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["gc.byRead"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      if(is.null(plot.medians)){
        plot.medians <- FALSE;
      } else if(plot.means & plot.medians){
        plot.means <- FALSE;
      }

      makePlot.generic(plot.name,res@qc.data[["gc.byRead"]], plotter, 
               x.name = "NUM_BASES_GC",y.name = "CT",norm.x = plotRate,
               avg.y = TRUE,plot.type = "lines", 
               plot.means = plot.means, plot.medians = plot.medians, 
               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
               ...);
      if(plotRate){ title(xlab = "% G/C");
      } else { title(xlab = "# G/C"); }
      title(ylab = "Frequency");
      #abline(v=50);
      internal.plot.main.title("GC Content", plotter, ...);
  
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.clipping <- function(plotter, rate.per.million = FALSE, use.readLength.denominator = TRUE, 
                              r2.buffer = NULL , debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,
                              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                              plot = TRUE,...){
  res <- plotter$res
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  makePlot.cigarOp.byCycle(plotter,"SoftClip", r2.buffer = r2.buffer, rate.per.million = rate.per.million,use.readLength.denominator=use.readLength.denominator, singleEndMode = singleEndMode,debugMode=debugMode,
                           rasterize.plotting.area=rasterize.plotting.area,raster.height=raster.height,raster.width=raster.width,plot=plot, ... );
}

#tryCatch({
#}, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
makePlot.cigarOp.byCycle <- function(plotter,op = c("SoftClip","Del","Ins","HardClip","Pad","Splice","Aln"),
                                     r2.buffer = NULL, 
                                     rate.per.million = TRUE, 
                                     use.readLength.denominator = TRUE, 
                                     debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,
                                     rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                     plot = TRUE,...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["cigarOpDistribution.byReadCycle.R1"]]) )
  if(! plot) return(isPlottable);
  
  op <- match.arg(op);

  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  SUFFIX.FIELD.IDS <- c("_S","_M","_E","_B")
  H.FIELD.IDS <- paste0("H",SUFFIX.FIELD.IDS);
  S.FIELD.IDS <- paste0("S",SUFFIX.FIELD.IDS);
  I.FIELD.IDS <- paste0("I",SUFFIX.FIELD.IDS);
  D.FIELD.IDS <- paste0("D",SUFFIX.FIELD.IDS);
  P.FIELD.IDS <- paste0("P",SUFFIX.FIELD.IDS);
  N.FIELD.IDS <- paste0("N",SUFFIX.FIELD.IDS);
  M.FIELD.IDS <- paste0("M",SUFFIX.FIELD.IDS);
  if(op == "SoftClip"){
    op.title <- "Alignment Clipping Rate";
    op.field.list <- S.FIELD.IDS;
  } else if(op == "Del"){
    op.title <- "Deletion Rate";
    op.field.list <- D.FIELD.IDS
  } else if(op == "Ins"){
    op.title <- "Insertion Rate";
    op.field.list <- I.FIELD.IDS
  } else if(op == "HardClip"){
    op.title <- "Hard Clipping Rate";
    op.field.list <- H.FIELD.IDS
  } else if(op == "Pad"){
    op.title <- "Padding Rate";
    op.field.list <- P.FIELD.IDS
  } else if(op == "Splice"){
    op.title <- "Splice Junction Rate";
    op.field.list <- N.FIELD.IDS
  } else if(op == "Aln"){
    op.title <- "Aligned To Reference Rate";
    op.field.list <- M.FIELD.IDS
  } else {
    stop("Fatal error: cigarOp Code not recognized! Must be one of: SoftClip,HardClip,Del,Ins,Pad,Splice");
  }
  res <- plotter$res;
  plot.name <- paste0(op.title, ", by read cycle");
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["cigarOpDistribution.byReadCycle.R1"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(! singleEndMode){
      tryCatch({
        get.cigarOpRate.from.dl <- function(dl){
          lapply(dl, function(df){
              if(use.readLength.denominator){
                readCt <- rowSums(df[,c(M.FIELD.IDS,H.FIELD.IDS,S.FIELD.IDS,I.FIELD.IDS)]);
              } else {
                readCt <- sum(df[1,c(M.FIELD.IDS,H.FIELD.IDS,S.FIELD.IDS,I.FIELD.IDS)]);
              }
              if(rate.per.million){ cigarRate <- ((apply( df[,op.field.list], 1,sum)) * 1000000) / readCt;
              } else { cigarRate <- (apply( df[,op.field.list], 1,sum)) / readCt; }
              data.frame(CYCLE = df$CYCLE + 1, cigarRate = cigarRate );
          });
        }
        cigarOpRate.list.r1 <- get.cigarOpRate.from.dl(plotter$res@qc.data[["cigarOpDistribution.byReadCycle.R1"]]);
        cigarOpRate.list.r2 <- get.cigarOpRate.from.dl(plotter$res@qc.data[["cigarOpDistribution.byReadCycle.R2"]]);
        xlim <- c(1,max(plotter$res@decoder$cycle.CT))
        if(is.null(r2.buffer)) r2.buffer <- xlim[2] / 10;
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.pair(plot.name,cigarOpRate.list.r1, 
                             cigarOpRate.list.r2, 
                             plotter = plotter, 
                             x.name = "CYCLE",
                             y.name = "cigarRate",
                             norm.x = FALSE,
                             avg.y = FALSE,
                             plot.type = "lines",
                             xlim=xlim,
                             vert.offset=0, 
                             draw.horiz.lines = FALSE,
                             override.lty = 1,
                             r2.buffer = r2.buffer,
                             rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                             ...);
      
      if(! singleEndMode){
        internal.plot.read.labels(plotter,"top",r2.buffer,xlim[2]);
      }
      title(xlab="Read Cycle");
      if(rate.per.million){
         title(ylab=paste0(op.title," (per Million Reads)"));
      } else {
         title(ylab=op.title);
      }
   
      internal.plot.main.title(paste0(op.title, ", by read cycle"), plotter, ...);
   
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      tryCatch({
        get.cigarOpRate.from.dl <- function(dl){
          lapply(dl, function(df){
              readCt <- sum(df[1,c(M.FIELD.IDS,H.FIELD.IDS,S.FIELD.IDS,I.FIELD.IDS)]);
              if(rate.per.million){ cigarRate <- ((apply( df[,op.field.list], 1,sum)) * 1000000) / readCt;
              } else { cigarRate <- (apply( df[,op.field.list], 1,sum)) / readCt; }
              data.frame(CYCLE = df$CYCLE + 1, cigarRate = cigarRate );
          });
        }
        cigarOpRate.list.r1 <- get.cigarOpRate.from.dl(plotter$res@qc.data[["cigarOpDistribution.byReadCycle.R1"]]);
        xlim <- c(1,max(plotter$res@decoder$cycle.CT))
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic(plot.name,cigarOpRate.list.r1, 
                             plotter = plotter, 
                             x.name = "CYCLE",
                             y.name = "cigarRate",
                             norm.x = FALSE,
                             avg.y = FALSE,
                             plot.type = "lines",
                             xlim=xlim,
                             vert.offset=0, 
                             draw.horiz.lines = FALSE,
                             override.lty = 1,
                             rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                             ...);
      title(xlab="Read Cycle");
      if(rate.per.million){
         title(ylab=paste0(op.title," (per Million Reads)"));
      } else {
         title(ylab=op.title);
      }
   
      internal.plot.main.title(paste0(op.title, ", by read cycle"), plotter, ...);
   
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}



makePlot.insert.size <- function(plotter, calc.rate = TRUE, pct.cutoff = 0.98, plot.medians = TRUE, plot.means = NULL, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, xlim = NULL,
                                 rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                 plot = TRUE,...){
  res <- plotter$res
  isPlottable <- ! (singleEndMode || is.null(plotter$res@qc.data[["insert.size"]]))
  if(! plot) return(isPlottable);

  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plot.name <- "Insert Size";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(singleEndMode){
      message(paste0("Warning: Skipping ",plot.name," plotting.\n    (Insert Size Not Applicable in Single-End Mode.)"));
      blank.plot("Insert Size Not Applicable\nin Single-End Mode!\nSkipping...");
    } else if(is.null(plotter$res@qc.data[["insert.size"]])){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(is.null(plot.means )){
            plot.means <- FALSE;
        } else if(plot.means & plot.medians){
            plot.medians <- FALSE;
        }
        data.list <- lapply(res@qc.data[["insert.size"]],function(df){
          df[df$InsertSize >= 0,];
        });
        if(is.null(xlim)){
          xlim.min <- min(sapply(data.list,function(df){min(df$InsertSize)}));
          xlim.max <- max(sapply(data.list,function(df){
            cts <- df$Ct;
            total <- sum(cts);
            cumsum <- cumsum(cts) / sum(cts);
            limit.index <- which(cumsum >= pct.cutoff)[1];
            df$InsertSize[limit.index];
          }));
          xlim <- c(xlim.min,xlim.max);
        }
        
        #print(xlim);
        data.list <- lapply(data.list,function(df){
          rbind.data.frame(data.frame(InsertSize = c(-1), Ct = c(0)),df);
        });

        pre.plot.func <- function(){ 
          max.cycle.ct <- max(plotter$res@decoder$cycle.CT);
          abline(v=max.cycle.ct,col="gray",lty=3);
          abline(v=max.cycle.ct * 2,col="gray",lty=3);
        }
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic(plot.name,res@qc.data[["insert.size"]], plotter, 
               x.name = "InsertSize",y.name = "Ct",
               norm.x = FALSE, avg.y = calc.rate,
               plot.type = "lines", 
               xlim = xlim,
               plot.means = plot.means, plot.medians = plot.medians,pre.plot.func = pre.plot.func,
               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
               ...);
      internal.plot.main.title("Insert Size", plotter, ...);
      title(xlab = "Insert Size (bp)");
      title(ylab = "Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}



makePlot.genebody <- function(plotter, 
                  geneset = c("Overall","90-100","75-90","50-75","0-50"),
                  avgMethod = c("TotalCounts", "AvgPercentile"), 
                  plot.medians = NULL, 
                  plot.means = TRUE, 
                  debugMode = DEFAULTDEBUGMODE, 
                  singleEndMode = plotter$res@singleEnd, 
                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                  plot = TRUE,... ){
  res <- plotter$res
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  avgMethod <- match.arg(avgMethod);
  geneset <- match.arg(geneset);
  
  if(avgMethod == "TotalCounts"){
    if(geneset == "Overall"){
      makePlot.genebody.coverage(plotter = plotter, plot.medians = plot.medians, plot.means=plot.means, debugMode=debugMode, singleEndMode=singleEndMode, rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,plot=plot, ...);
    } else {
      makePlot.genebody.coverage.TC(plotter = plotter, geneset = geneset, plot.medians = plot.medians, plot.means = plot.means, debugMode = debugMode, singleEndMode=singleEndMode,rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,plot=plot, ...);
    }
  } else {
    makePlot.genebody.coverage.PCT(plotter = plotter, geneset = geneset, plot.medians = plot.medians, plot.means = plot.means, debugMode = debugMode, singleEndMode=singleEndMode,rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,plot=plot, ...);
  }
  
}

makePlot.genebody.coverage.PCT <- function(plotter, 
                    geneset = c("Overall","90-100","75-90","50-75","0-50"), 
                    plot.meanPercentile = TRUE, 
                    plot.medians = NULL, 
                    plot.means = TRUE, 
                    debugMode = DEFAULTDEBUGMODE, 
                    singleEndMode = plotter$res@singleEnd, 
                    rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000, 
                    plot = TRUE,
                    ... ){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["geneBodyCoverage.pct"]]) )
  if(! plot) return(isPlottable);
  
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  geneset <- match.arg(geneset);
  
  if(geneset == "Overall"){
    geneset.title <- "All";
    df.ID <- "TOTAL";
  } else if(geneset == "90-100"){
    geneset.title <- "90th Percentile";
    df.ID <- "X4.high";
  } else if(geneset == "75-90"){
    geneset.title <- "75-90 Percentile";
    df.ID <- "X3.75to90";
  } else if(geneset == "50-75"){
    geneset.title <- "Upper Middle Quartile";
    df.ID <- "X2.upperMidQuartile";
  } else if(geneset == "0-50"){
    geneset.title <- "Low Expression";
    df.ID <- "X1.bottomHalf";
  }
  
  plot.name <- paste0("Gene-Body Coverage, ",geneset.title," Genes");
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["geneBodyCoverage.pct"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(is.null(plot.medians)){
            plot.medians <- FALSE;
        } else if(plot.means & plot.medians){
            plot.means <- FALSE;
        }
        raw.data.list <- res@qc.data[["geneBodyCoverage.pct"]];
        step.size <- raw.data.list[[1]]$QUANTILE[2] - raw.data.list[[1]]$QUANTILE[1];
        #print(step.size);
        data.list <- lapply(raw.data.list,function(df){
          data.frame(Quantile = df$QUANTILE - (step.size / 2), GeneBodyCoverage = df[,df.ID]);
        });
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic(plot.name, data.list, plotter, 
                   x.name = "Quantile",y.name = "GeneBodyCoverage", 
                   norm.x = TRUE, avg.y = FALSE,
                   plot.type = "lines", 
                   plot.means = plot.means, plot.medians = plot.medians,
                   rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                   ...);

      internal.plot.main.title(plot.name, plotter, ...);
      title(xlab = "Percentile of Gene Body (5\'->3\')");
      title(ylab = "Avg Rate");
      
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.genebody.coverage.TC <- function(plotter, 
                    geneset = c("90-100","75-90","50-75","0-50"), 
                    plot.meanPercentile = TRUE, 
                    plot.medians = NULL, 
                    plot.means = TRUE, 
                    debugMode = DEFAULTDEBUGMODE, 
                    rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000, 
                    singleEndMode = plotter$res@singleEnd, 
                    plot = TRUE,... ){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["geneBodyCoverage.by.expression.level"]]) )
  if(! plot) return(isPlottable);

  geneset <- match.arg(geneset);
  
  if(geneset == "90-100"){
    geneset.title <- "Top 10 Percentile";
    df.ID <- "X4.high";
  } else if(geneset == "75-90"){
    geneset.title <- "75-90 Percentile";
    df.ID <- "X3.75to90";
  } else if(geneset == "50-75"){
    geneset.title <- "Upper Middle Quartile";
    df.ID <- "X2.upperMidQuartile";
  } else if(geneset == "0-50"){
    geneset.title <- "Low Expression";
    df.ID <- "X1.bottomHalf";
  }
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  
  plot.name <- paste0("Gene-Body Coverage, ",geneset.title," Genes");
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["geneBodyCoverage.by.expression.level"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(is.null(plot.medians)){
            plot.medians <- FALSE;
        } else if(plot.means & plot.medians){
            plot.means <- FALSE;
        }
        raw.data.list <- res@qc.data[["geneBodyCoverage.by.expression.level"]];
        step.size <- raw.data.list[[1]]$QUANTILE[2] - raw.data.list[[1]]$QUANTILE[1];
        #print(step.size);
        data.list <- lapply(raw.data.list,function(df){
          data.frame(Quantile = df$QUANTILE - (step.size / 2), GeneBodyCoverage = df[,df.ID]);
        });
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic(plot.name, data.list, plotter, 
                   x.name = "Quantile",y.name = "GeneBodyCoverage", 
                   norm.x = TRUE, avg.y = TRUE,
                   plot.type = "lines", 
                   plot.means = plot.means, plot.medians = plot.medians,
                   rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                   ...);

      internal.plot.main.title(plot.name, plotter, ...);
      title(xlab = "Percentile of Gene Body (5\'->3\')");
      title(ylab = paste0("Proportion of ",readLabel));

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}



makePlot.genebody.coverage <- function(plotter, plot.medians = NULL, plot.means = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                                       rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                       plot = TRUE,... ){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["geneBodyCoverage.by.expression.level"]]) )
  if(! plot) return(isPlottable);
  
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plot.name <- "Gene-Body Coverage";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["geneBodyCoverage.by.expression.level"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(is.null(plot.medians)){
          plot.medians <- FALSE;
        } else if(plot.means & plot.medians){
          plot.means <- FALSE;
        }
        raw.data.list <- res@qc.data[["geneBodyCoverage.by.expression.level"]];
        step.size <- raw.data.list[[1]]$QUANTILE[2] - raw.data.list[[1]]$QUANTILE[1];
        #print(step.size);
        data.list <- lapply(raw.data.list,function(df){
          data.frame(Quantile = df$QUANTILE - (step.size / 2), GeneBodyCoverage = apply(df[,c("X1.bottomHalf","X2.upperMidQuartile","X3.75to90","X4.high")],1,sum));
        });
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});

      makePlot.generic(plot.name, data.list, plotter, 
               x.name = "Quantile",y.name = "GeneBodyCoverage", 
               norm.x = TRUE, avg.y = TRUE,
               plot.type = "lines", 
               plot.means = plot.means, plot.medians = plot.medians, 
               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
               ...);

      internal.plot.main.title("Gene-Body Coverage", plotter, ...);
      title(xlab = "Percentile of Gene Body (5\'->3\')");
      title(ylab = paste0("Proportion of ",readLabel));
  
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}





makePlot.genebody.coverage.UMQuartile <- function(plotter, plot.medians = NULL, plot.means = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                                  plot = TRUE,... ){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["geneBodyCoverage.by.expression.level"]]) )
  if(! plot) return(isPlottable);

  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plot.name <- "Gene-Body Coverage, Upper Middle Quartile Genes";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
    plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["geneBodyCoverage.by.expression.level"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(is.null(plot.medians)){
            plot.medians <- FALSE;
        } else if(plot.means & plot.medians){
            plot.means <- FALSE;
        }
        raw.data.list <- res@qc.data[["geneBodyCoverage.by.expression.level"]];
        step.size <- raw.data.list[[1]]$QUANTILE[2] - raw.data.list[[1]]$QUANTILE[1];
        #print(step.size);
        data.list <- lapply(raw.data.list,function(df){
          data.frame(Quantile = df$QUANTILE - (step.size / 2), GeneBodyCoverage = df[,"X2.upperMidQuartile"]);
        });
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic(plot.name, data.list, plotter, 
                   x.name = "Quantile",y.name = "GeneBodyCoverage", 
                   norm.x = TRUE, avg.y = TRUE,
                   plot.type = "lines", 
                   plot.means = plot.means, plot.medians = plot.medians, 
                   rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                   ...);

      internal.plot.main.title("Gene-Body Coverage, Upper Middle Quartile Genes", plotter, ...);
      title(xlab = "Percentile of Gene Body (5\'->3\')");
      title(ylab = "Rate");

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.genebody.coverage.lowExpress <- function(plotter, plot.medians = NULL, plot.means = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                                                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                                  plot = TRUE,... ){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["geneBodyCoverage.by.expression.level"]]) )
  if(! plot) return(isPlottable);

  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plot.name <- "Gene-Body Coverage, Low Expression Genes";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["geneBodyCoverage.by.expression.level"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(is.null(plot.medians)){
            plot.medians <- FALSE;
        } else if(plot.means & plot.medians){
            plot.means <- FALSE;
        }
        raw.data.list <- res@qc.data[["geneBodyCoverage.by.expression.level"]];
        step.size <- raw.data.list[[1]]$QUANTILE[2] - raw.data.list[[1]]$QUANTILE[1];
        #print(step.size);
        data.list <- lapply(raw.data.list,function(df){
          data.frame(Quantile = df$QUANTILE - (step.size / 2), GeneBodyCoverage = df[,"X1.bottomHalf"]);
        });
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic(plot.name, data.list, plotter, 
                   x.name = "Quantile",y.name = "GeneBodyCoverage", 
                   norm.x = TRUE, avg.y = TRUE,
                   plot.type = "lines", 
                   plot.means = plot.means, plot.medians = plot.medians,
                   rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                   ...);

      internal.plot.main.title("Gene-Body Coverage, Low Expression Genes", plotter, ...);
      title(xlab = "Percentile of Gene Body (5\'->3\')");
      title(ylab = "Rate");

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}



makePlot.missingness.rate <- function(plotter,  r2.buffer = NULL, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                                      rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                      log.y = FALSE,ylim=NULL,
                                      plot = TRUE,...){

  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["NVC.raw.R1"]]) )
  if(! plot) return(isPlottable);
  
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plot.name <- "Missingness Rate, by read cycle";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["NVC.raw.R1"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(singleEndMode){
      tryCatch({
         raw.data.list.r1 <- res@qc.data[["NVC.raw.R1"]];
         data.list.r1 <- lapply(raw.data.list.r1,function(df){
           readCt <- sum(df$CT[df$readPos == 0]);
           subset <- df[df$base == "N",]
           data.frame(readPos = subset$readPos, CT = subset$CT, RATE = subset$CT / readCt, RATE_LOG = log10(subset$CT / readCt));
         });
         
         y.name <- "RATE";
         if(log.y){
           y.name <- paste0(y.name,"_LOG");
         }
         
         xlim <- c(1,max(plotter$res@decoder$cycle.CT))
       }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
       makePlot.generic(plot.name, data.list.r1, 
                                 plotter = plotter, 
                                 x.name = "readPos",
                                 y.name = y.name,
                                 norm.x = FALSE,
                                 avg.y = FALSE,
                                 plot.type = "lines",
                                 xlim=xlim,y.is.log = log.y,
                                 vert.offset=0, 
                                 draw.horiz.lines = FALSE,
                                 override.lty = 1,ylim=ylim,
                                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                                 ...);
      internal.plot.main.title("'N' Rate, by Read Cycle", plotter, ...);
      title(xlab = "Read Cycle");
      title(ylab = "Missing Nucleotide Rate");
      
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      tryCatch({
        raw.data.list.r1 <- res@qc.data[["NVC.raw.R1"]];
        raw.data.list.r2 <- res@qc.data[["NVC.raw.R2"]];  
        data.list.r1 <- lapply(raw.data.list.r1,function(df){
          readCt <- sum(df$CT[df$readPos == 0]);
          subset <- df[df$base == "N",]
          data.frame(readPos = subset$readPos, CT = subset$CT, RATE = subset$CT / readCt, RATE_LOG = log10(subset$CT / readCt));
        });
        data.list.r2 <- lapply(raw.data.list.r2,function(df){
          readCt <- sum(df$CT[df$readPos == 0]);
          subset <- df[df$base == "N",]
          data.frame(readPos = subset$readPos, CT = subset$CT, RATE = subset$CT / readCt, RATE_LOG = log10(subset$CT / readCt));
        });
         y.name <- "RATE";
         if(log.y){
           y.name <- paste0(y.name,"_LOG");
         }
         xlim <- c(1,max(plotter$res@decoder$cycle.CT))
         if(is.null(r2.buffer)) r2.buffer <- xlim[2] / 10;
       }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
       makePlot.generic.pair(plot.name,        data.list.r1, 
                                 data.list.r2, 
                                 plotter = plotter, 
                                 x.name = "readPos",
                                 y.name = y.name,
                                 norm.x = FALSE,
                                 avg.y = FALSE,
                                 plot.type = "lines",
                                 xlim=xlim,y.is.log = log.y,
                                 vert.offset=0, 
                                 draw.horiz.lines = FALSE,
                                 override.lty = 1,
                                 r2.buffer = r2.buffer,ylim=ylim,
                                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                                 ...);
      internal.plot.main.title("'N' Rate, by Read Cycle", plotter, ...);
      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab = "Read Cycle");
      title(ylab = "Missing Nucleotide Rate");

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

####################################################################################
####################################################################################
####################################################################################


makePlot.cigarLength.distribution <- function(plotter,op, r2.buffer = NULL, perMillion = TRUE, log.x = FALSE, log.y = FALSE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,
                                              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000, 
                                              plot = TRUE,...){

  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["cigarOpLengths.byOp.R1"]]) )
  if(! plot) return(isPlottable);
  
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(op == "SoftClip"){
    op.title <- "Alignment Clipping Length";
    op.field <- "S"
  } else if(op == "Del"){
    op.title <- "Deletion Length";
    op.field <- "D"
  } else if(op == "Ins"){
    op.title <- "Insertion Length";
    op.field <- "I"
  } else if(op == "HardClip"){
    op.title <- "Hard Clipping Length";
    op.field <- "H"
  } else if(op == "Pad"){
    op.title <- "Padding Length";
    op.field <- "P"
  } else if(op == "Splice"){
    op.title <- "Splice Junction Length";
    op.field <- "N"
  } else if(op == "Aln"){
    op.title <- "Exon Block Length";
    op.field <- "M"
  } else {
    stop("Fatal error: cigarOp Code not recognized! Must be one of: SoftClip,HardClip,Del,Ins,Pad,Splice");
  }
  plot.name <- paste0(op.title," Distribution", "");
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){ 
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["cigarOpLengths.byOp.R1"]]) ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(singleEndMode) {
      tryCatch({
        get.cigarLenData.from.dl <- function(dl){
          out.list <- lapply(1:length(dl), function(i){
            df <- dl[[i]];
            curr.lanebam <- names(dl)[i];
            readCt <- as.numeric(plotter$res@qc.data[["summary"]][[curr.lanebam]]$COUNT[plotter$res@qc.data[["summary"]][[curr.lanebam]]$FIELD == "READ_PAIR_OK"]);
            opDf <- df[df[["OP"]] == op.field,,drop = FALSE];
            opDf$RATE <- (opDf$CT ) / readCt;
            if(perMillion) opDf$RATE <- (opDf$CT * 1000000) / readCt;
            if(log.x) opDf$LEN <- log10(opDf$LEN);
            if(log.y) opDf$RATE <- log10(opDf$RATE);
            opDf[, c("LEN","RATE")];
          });
          names(out.list) <- names(dl);
          out.list;
        }
        cigarLenData.list.r1 <- get.cigarLenData.from.dl(plotter$res@qc.data[["cigarOpLengths.byOp.R1"]]);

        xlim.max <- max(sapply(cigarLenData.list.r1, function(df){max(df$LEN)}))
        xlim.min <- min(sapply(cigarLenData.list.r1, function(df){min(df$LEN)}))
        if(xlim.max == xlim.min){
          xlim.max <- xlim.max + 0.0001;
        }
        
        xlim <- c(xlim.min,xlim.max)
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      
      
      makePlot.generic(plot.name,cigarLenData.list.r1, 
                                 plotter = plotter, 
                                 x.name = "LEN",
                                 y.name = "RATE",
                                 norm.x = FALSE,
                                 avg.y = FALSE,
                                 plot.type = "lines",
                                 xlim=xlim,
                                 vert.offset=0, 
                                 draw.horiz.lines = FALSE,
                                 override.lty = 1,
                                 x.is.log = log.x,
                                 y.is.log = log.y,
                                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                                 ...);

      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab=op.title);
      if(perMillion){ title(ylab="Rate (per Million Reads)");
      } else { title(ylab="Rate (per Read)");}

      internal.plot.main.title(paste0(op.title," Distribution"), plotter, ...);

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      tryCatch({
        get.cigarLenData.from.dl <- function(dl){
          out.list <- lapply(1:length(dl), function(i){
            df <- dl[[i]];
            curr.lanebam <- names(dl)[i];
            readCt <- as.numeric(plotter$res@qc.data[["summary"]][[curr.lanebam]]$COUNT[plotter$res@qc.data[["summary"]][[curr.lanebam]]$FIELD == "READ_PAIR_OK"]);
            opDf <- df[df[["OP"]] == op.field,];
            opDf$RATE <- (opDf$CT ) / readCt;
            if(perMillion) opDf$RATE <- (opDf$CT * 1000000) / readCt;
            if(log.x) opDf$LEN <- log10(opDf$LEN);
            if(log.y) opDf$RATE <- log10(opDf$RATE);
            opDf[, c("LEN","RATE")];
          });
          names(out.list) <- names(dl);
          out.list;
        }

        cigarLenData.list.r1 <- get.cigarLenData.from.dl(plotter$res@qc.data[["cigarOpLengths.byOp.R1"]]);
        cigarLenData.list.r2 <- get.cigarLenData.from.dl(plotter$res@qc.data[["cigarOpLengths.byOp.R2"]]);

        xlim.max <- max(
                         max(sapply(cigarLenData.list.r1, function(df){max(df$LEN)})),
                         max(sapply(cigarLenData.list.r2, function(df){max(df$LEN)}))
                         );
        xlim.min <- max(
                         min(sapply(cigarLenData.list.r1, function(df){min(df$LEN)})),
                         min(sapply(cigarLenData.list.r2, function(df){min(df$LEN)}))
                         );

        xlim <- c(xlim.min,xlim.max)
        if(is.null(r2.buffer)) r2.buffer <- (xlim.max - xlim.min) / 10;
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.pair(plot.name,        cigarLenData.list.r1, 
                                 cigarLenData.list.r2, 
                                 plotter = plotter, 
                                 x.name = "LEN",
                                 y.name = "RATE",
                                 norm.x = FALSE,
                                 avg.y = FALSE,
                                 plot.type = "lines",
                                 xlim=xlim,
                                 vert.offset=0, 
                                 draw.horiz.lines = FALSE,
                                 override.lty = 1,
                                 r2.buffer = r2.buffer,
                                 x.is.log = log.x,
                                 y.is.log = log.y,
                                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                                 ...);

      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab=op.title);
      if(perMillion){ title(ylab="Rate (per Million Reads)");
      } else { title(ylab="Rate (per Read)");}

      internal.plot.main.title(paste0(op.title," Distribution"), plotter, ...);

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}




####################################################################################
####################################################################################
####################################################################################

makePlot.gene.cdf <- function(plotter, sampleWise = FALSE, plot.intercepts = TRUE, label.intercepts = FALSE, debugMode = DEFAULTDEBUGMODE,
                              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000, singleEndMode = plotter$res@singleEnd,
                              plot = TRUE,
                              ...){
  res <- plotter$res
  #plot.name <- "gene cdf";
  #if(debugMode){ message("Starting: ",plot.name," plot."); }
  #plotter.error.wrapper(plot.name, plotterFcn = function(){
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(rasterize.plotting.area){
    check.rasterize.or.die("rasterize.plotting.area");
  }
  if(sampleWise){
    makePlot.gene.cdf.sampleWise(plotter, plot.intercepts = plot.intercepts, label.intercepts = label.intercepts,
          rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode, singleEndMode = singleEndMode,plot=plot, ...);
  } else {
    makePlot.gene.cdf.bamWise(plotter, plot.intercepts = plot.intercepts, label.intercepts = label.intercepts, 
          rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode, singleEndMode = singleEndMode,plot=plot,...);
  } 
  #})
}

makePlot.gene.cdf.bamWise <- function(plotter, plot.intercepts = TRUE, label.intercepts = FALSE, debugMode = DEFAULTDEBUGMODE,
                                      rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000, singleEndMode = plotter$res@singleEnd,
                                      plot = TRUE,
                                      ...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@calc.data[["LANEBAM_GENE_CDF"]]) )
  if(! plot) return(isPlottable);

  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plot.name <- "Cumulative Gene Assignment Diversity";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@calc.data[["LANEBAM_GENE_CDF"]]) ){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      makePlot.gene.cdf.helper(plot.name,plotter$res@calc.data[["LANEBAM_GENE_CDF"]], plotter, plot.intercepts = plot.intercepts, label.intercepts = label.intercepts, rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode, ...);
      internal.plot.main.title("Cumulative Gene Assignment Diversity", plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.gene.cdf.sampleWise <- function(plotter, plot.intercepts = TRUE, label.intercepts = TRUE, debugMode = DEFAULTDEBUGMODE,
                                         rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000, singleEndMode = plotter$res@singleEnd,
                                         plot = TRUE,
                                         ...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@calc.data[["SAMPLE_GENE_CDF"]]) )
  if(! plot) return(isPlottable);
  
  plot.name <- "Cumulative Gene Assignment Diversity, sampleWise";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@calc.data[["SAMPLE_GENE_CDF"]]) ){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      makePlot.gene.cdf.bySample.helper(plot.name,plotter$res@calc.data[["SAMPLE_GENE_CDF"]], plotter$title.highlight.name, plot.intercepts = plot.intercepts, label.intercepts = label.intercepts, 
                                        rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                                        ...);
      internal.plot.main.title("Cumulative Gene Assignment Diversity", plotter, ...);

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

####################################################################################
####################################################################################
####################################################################################

makePlot.raw.NVC <- function(plotter,  r2.buffer = NULL,  points.highlighted = TRUE,
                         label.majority.bases = FALSE, label.majority.bases.threshold = 0.5, label.majority.bases.cex = 0.5, 
                         rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 2000,
                         show.base.legend = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                         useFQ = FALSE,
                         plot = TRUE,...){
  res <- plotter$res
  if(useFQ){
    data.list.r1 <- plotter$res@calc.data[["FQ.NVC.R1"]];
    data.list.r2 <- plotter$res@calc.data[["FQ.NVC.R2"]];
    plot.title <- "Raw Nucleotide Rate by Cycle (fastq)"
  } else {
    data.list.r1 <- plotter$res@calc.data[["NVC.raw.R1"]];
    data.list.r2 <- plotter$res@calc.data[["NVC.raw.R2"]];
    plot.title <- "Raw Nucleotide Rate by Cycle"
  }
  
  isPlottable <- ! (is.null(data.list.r1) )
  if(! plot) return(isPlottable);

  plot.name <- "NVC, Raw";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    if(rasterize.plotting.area){
      check.rasterize.or.die("rasterize.plotting.area");
    }

    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(! isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(singleEndMode){
      xlim <- c(1,max(plotter$res@decoder$cycle.CT))
      if(is.null(r2.buffer)) r2.buffer <- xlim[2] / 10;
      ylim <- c(0.1,0.55);

      blank.label.char = "";

      makePlot.generic.NVC.single(plot.name,data.list.r1, 
                               plotter = plotter, 
                               x.name = "readPos",
                               y.name = "CT",
                               xlim=xlim,
                               ylim=ylim,
                               points.highlighted = points.highlighted,
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = blank.label.char,label.major.cex = label.majority.bases.cex,
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                               ...);
      internal.plot.main.title(plot.title, plotter, plot.type = "NVC", ...);
      title(xlab = "Read Cycle");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      xlim <- c(1,max(plotter$res@decoder$cycle.CT))
      if(is.null(r2.buffer)) r2.buffer <- xlim[2] / 10;
      ylim <- c(0.1,0.55);

      blank.label.char = "";

      makePlot.generic.NVC.pair(plot.name,    
                               data.list.r1, 
                               data.list.r2, 
                               plotter = plotter, 
                               x.name = "readPos",
                               y.name = "CT",
                               xlim=xlim,
                               ylim=ylim,
                               r2.buffer = r2.buffer, points.highlighted = points.highlighted,
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = blank.label.char,label.major.cex = label.majority.bases.cex,
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                               ...);
      internal.plot.main.title(plot.title, plotter, plot.type = "NVC", ...);
      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab = "Read Cycle");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.minus.clipping.NVC <- function(plotter,  r2.buffer = NULL, points.highlighted = TRUE,
                                    label.majority.bases = FALSE, label.majority.bases.threshold = 0.5, label.majority.bases.cex = 0.5, 
                                    rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 2000,
                                    show.base.legend = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,
                                    plot = TRUE,...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["NVC.minus.clipping.R1"]]) )
  if(! plot) return(isPlottable);

  plot.name <- "NVC, Aligned bases only";
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    if(rasterize.plotting.area){
      check.rasterize.or.die("rasterize.plotting.area");
    }
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!  isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(singleEndMode){
      xlim <- c(1,max(plotter$res@decoder$cycle.CT))
      ylim <- c(0.1,0.55);
      data.list.r1 <- plotter$res@calc.data[["NVC.minus.clipping.R1"]];
      makePlot.generic.NVC.single(plot.name,   data.list.r1 , 
                                 plotter = plotter, 
                                 x.name = "readPos",
                                 y.name = "CT",
                                 xlim=xlim,
                                 ylim=ylim,
                                 points.highlighted = points.highlighted,
                                 label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = blank.label.char,label.major.cex = label.majority.bases.cex,
                                 show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                                 ...);
      internal.plot.main.title("Nucleotide Rate by Cycle, Aligned bases only", plotter, plot.type = "NVC", ...);
      title(xlab = "Read Cycle");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      xlim <- c(1,max(plotter$res@decoder$cycle.CT))
      if(is.null(r2.buffer)) r2.buffer <- xlim[2] / 10;
      #if(is.null(label.major.cutoff)){
      #   label.major = FALSE;
      #   label.cutoff = 0.5;
          blank.label.char = "";
      #   label.major.cex = label.major.cex;
      #} else {
      #   label.major = TRUE;
      #   label.cutoff = label.major.cutoff;
      #   blank.label.char = "";
      #   label.major.cex = label.major.cex;
      #}
      ylim <- c(0.1,0.55);
      data.list.r1 <- plotter$res@calc.data[["NVC.minus.clipping.R1"]];
      data.list.r2 <- plotter$res@calc.data[["NVC.minus.clipping.R2"]];

      makePlot.generic.NVC.pair(plot.name, data.list.r1, 
                                 data.list.r2,
                                 plotter = plotter,
                                 x.name = "readPos",
                                 y.name = "CT",
                                 xlim=xlim,
                                 ylim=ylim,
                                 r2.buffer = r2.buffer, points.highlighted = points.highlighted,
                                 label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = blank.label.char,label.major.cex = label.majority.bases.cex,
                                 show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                                 rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                                 ...);
      internal.plot.main.title("Nucleotide Rate by Cycle, Aligned bases only", plotter, plot.type = "NVC", ...);
      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab = "Read Cycle");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.NVC.lead.clip <- function(plotter, clip.amt = 10,  r2.buffer = clip.amt / 10, points.highlighted = TRUE,
                              label.majority.bases = TRUE, label.majority.bases.threshold = 0.5, label.majority.bases.cex = 1, 
                              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                              show.base.legend = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,
                              plot = TRUE,...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["NVC.lead.clip.R1"]]) )
  if(! plot) return(isPlottable);
  
  plot.name <- "Lead Clipping NVC";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  clip.amt <- min(c(clip.amt, plotter$res@decoder$cycle.CT));
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    if(rasterize.plotting.area){
      check.rasterize.or.die("rasterize.plotting.area");
    }
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["NVC.lead.clip.R1"]]) ){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(singleEndMode){
      xlim <- c(1,clip.amt)
      ylim <- c(0,1);
      data.list.r1 <- lapply(plotter$res@calc.data[["NVC.lead.clip.R1"]], function(df){
        df[df$leadClipLen == clip.amt,,drop = F];
      });
      makePlot.generic.NVC.single(plot.name,
                               data.list.r1,
                               plotter = plotter, 
                               x.name = "readPos",
                               y.name = "CT",
                               xlim=xlim,
                               ylim=ylim,
                               points.highlighted = points.highlighted,
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = "-",label.major.cex = label.majority.bases.cex,
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                               ...);
      internal.plot.main.title(paste0("Nucleotide Rate by Cycle, Leading Clipped bases (",clip.amt,")"), plotter, plot.type = "NVC", ...);
      title(xlab = "Read Cycle");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      xlim <- c(1,clip.amt)
      ylim <- c(0,1);
      data.list.r1 <- lapply(plotter$res@calc.data[["NVC.lead.clip.R1"]], function(df){
        df[df$leadClipLen == clip.amt,,drop = F];
      });
      data.list.r2 <- lapply(plotter$res@calc.data[["NVC.lead.clip.R2"]], function(df){
        df[df$leadClipLen == clip.amt,,drop = F];
      });
      makePlot.generic.NVC.pair(plot.name,    data.list.r1, 
                               data.list.r2, 
                               plotter = plotter, 
                               x.name = "readPos",
                               y.name = "CT",
                               xlim=xlim,
                               ylim=ylim,
                               r2.buffer = r2.buffer, points.highlighted = points.highlighted,
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = "-",label.major.cex = label.majority.bases.cex,
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                               ...);
      internal.plot.main.title(paste0("Nucleotide Rate by Cycle, Leading Clipped bases (",clip.amt,")"), plotter, plot.type = "NVC", ...);
      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab = "Read Cycle");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}


makePlot.NVC.tail.clip <- function(plotter, clip.amt = 10,  r2.buffer = clip.amt / 10, points.highlighted = TRUE, 
                               label.majority.bases = TRUE, label.majority.bases.threshold = 0.5, label.majority.bases.cex = 1, 
                               rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                               show.base.legend = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, 
                               plot = TRUE,...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["NVC.tail.clip.R1"]]) )
  if(! plot) return(isPlottable);

  plot.name <- "Tail Clipping NVC";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  clip.amt <- min(c(clip.amt, plotter$res@decoder$cycle.CT));
  
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    if(rasterize.plotting.area){
      check.rasterize.or.die("rasterize.plotting.area");
    }
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(is.null(plotter$res@qc.data[["NVC.tail.clip.R1"]]) ){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(singleEndMode){
      xlim <- c(1,clip.amt)
      ylim <- c(0,1);
      data.list.r1 <- lapply(plotter$res@calc.data[["NVC.tail.clip.R1"]], function(df){
        df[df$tailClipLen == clip.amt,,drop = F];
      })
      makePlot.generic.NVC.single(plot.name,    data.list.r1, 
                               plotter = plotter, 
                               x.name = "readPos", 
                               y.name = "CT", 
                               xlim=xlim, 
                               ylim=ylim, 
                               points.highlighted =  points.highlighted, 
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = "-",label.major.cex = label.majority.bases.cex,
                               count.x.from.end = TRUE, 
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                               ...);
      internal.plot.main.title(paste0("Nucleotide Rate by Cycle, Trailing Clipped bases (",clip.amt,")"), plotter, plot.type = "NVC", ...);
      title(xlab = "Read Cycle");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      xlim <- c(1,clip.amt)
      ylim <- c(0,1);
      data.list.r1 <- lapply(plotter$res@calc.data[["NVC.tail.clip.R1"]], function(df){
        df[df$tailClipLen == clip.amt,,drop = F];
      })
      data.list.r2 <- lapply(plotter$res@calc.data[["NVC.tail.clip.R2"]], function(df){
        df[df$tailClipLen == clip.amt,,drop = F];
      })
      makePlot.generic.NVC.pair(plot.name,    data.list.r1, 
                               data.list.r2, 
                               plotter = plotter, 
                               x.name = "readPos", 
                               y.name = "CT", 
                               xlim=xlim, 
                               ylim=ylim, 
                               r2.buffer = r2.buffer, 
                               points.highlighted =  points.highlighted, 
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = "-",label.major.cex = label.majority.bases.cex,
                               count.x.from.end = TRUE, 
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                               ...);
      internal.plot.main.title(paste0("Nucleotide Rate by Cycle, Trailing Clipped bases (",clip.amt,")"), plotter, plot.type = "NVC", ...);
      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab = "Read Cycle");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

####################################################################################
####################################################################################
####################################################################################

makePlot.NVC.lead.clip.matchByClipPosition <- function(plotter, clip.amt = 10,  r2.buffer = clip.amt / 10, 
                                                       label.majority.bases = TRUE, label.majority.bases.threshold = 0.5, label.majority.bases.cex = 1, 
                                                       rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                                       show.base.legend = TRUE, load.results = TRUE, debugMode = DEFAULTDEBUGMODE,  singleEndMode = plotter$res@singleEnd,
                                                       plot = TRUE,
                                                       ...){
  res <- plotter$res
  isPlottable <- !(is.null(plotter$res@qc.data[["NVC.lead.clip.R1"]]) )
  if(! plot) return(isPlottable);

  plot.name <- "Lead Clipping NVC";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    if(rasterize.plotting.area){
      check.rasterize.or.die("rasterize.plotting.area");
    }
    res <- plotter$res;
    #print(clip.amt)
    #print(r2.buffer)
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(singleEndMode){
      xlim <- c(1,clip.amt)
      ylim <- c(0,1);
      bases <- names(plotter$nvc.colors);

      get.data.list <- function(dl){
         lapply(dl, function(df){
            out <- data.frame(base = rep(bases,clip.amt), 
                              CT = unlist(lapply(1:clip.amt, function(i){
                                      sapply(bases, function(b){
                                         sum( df$CT[ df$base == b & df$leadClipLen <= clip.amt & (df$leadClipLen - df$readPos == i)] );
                                      });
                                   })),
                              stringsAsFactors = FALSE);
            #names(out) <- c("base","CT");
            return(out);
         })
      }

      list.r1.name <- paste0("NVC.lead.clip.matchByClipPosition.",clip.amt,".R1");

      if(! is.null(plotter$res@calc.data[[paste0("NVC.lead.clip.matchByClipPosition.",clip.amt,".R1")]])){
         data.list.r1 <- plotter$res@calc.data[[list.r1.name]];
      } else {
         message("generating data...");
         data.list.r1 <- get.data.list(plotter$res@qc.data[["NVC.lead.clip.R1"]]);
         plotter$res@calc.data[[list.r1.name]] <- data.list.r1;
      }
      makePlot.generic.NVC.single(plot.name, data.list.r1, 
                               plotter = plotter, 
                               x.name = "readPos",
                               y.name = "CT",
                               xlim=xlim,
                               ylim=ylim,
                               x.axis.labels = clip.amt:1,
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = "-",label.major.cex = label.majority.bases.cex,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light, debugMode = debugMode,
                               ...);
      internal.plot.main.title(paste0("Clipped Nucleotide Rate, Leading Clipped bases\n(1-",clip.amt,"), by distance to clip start"), plotter, plot.type = "NVC", ...);
      title(xlab = "Distance from Clip Start");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      xlim <- c(1,clip.amt)
      ylim <- c(0,1);
      bases <- names(plotter$nvc.colors);

      get.data.list <- function(dl){
         lapply(dl, function(df){
            out <- data.frame(base = rep(bases,clip.amt), 
                              CT = unlist(lapply(1:clip.amt, function(i){
                                      sapply(bases, function(b){
                                         sum( df$CT[ df$base == b & df$leadClipLen <= clip.amt & (df$leadClipLen - df$readPos == i)] );
                                      });
                                   })),
                              stringsAsFactors = FALSE);
            #names(out) <- c("base","CT");
            return(out);
         })
      }

      list.r1.name <- paste0("NVC.lead.clip.matchByClipPosition.",clip.amt,".R1");
      list.r2.name <- paste0("NVC.lead.clip.matchByClipPosition.",clip.amt,".R2");

      if(! is.null(plotter$res@calc.data[[paste0("NVC.lead.clip.matchByClipPosition.",clip.amt,".R1")]])){
         data.list.r1 <- plotter$res@calc.data[[list.r1.name]];
         data.list.r2 <- plotter$res@calc.data[[list.r2.name]];
      } else {
         message("generating data...");
         data.list.r1 <- get.data.list(plotter$res@qc.data[["NVC.lead.clip.R1"]]);
         data.list.r2 <- get.data.list(plotter$res@qc.data[["NVC.lead.clip.R2"]]);
         plotter$res@calc.data[[list.r1.name]] <- data.list.r1;
         plotter$res@calc.data[[list.r2.name]] <- data.list.r2;
      }
      makePlot.generic.NVC.pair(plot.name,    data.list.r1, 
                               data.list.r2, 
                               plotter = plotter, 
                               x.name = "readPos",
                               y.name = "CT",
                               xlim=xlim,
                               ylim=ylim,
                               r2.buffer = r2.buffer,
                               x.axis.labels = clip.amt:1,
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = "-",label.major.cex = label.majority.bases.cex,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width,
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light, debugMode = debugMode,
                               ...);
      internal.plot.main.title(paste0("Clipped Nucleotide Rate, Leading Clipped bases\n(1-",clip.amt,"), by distance to clip start"), plotter, plot.type = "NVC", ...);
      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab = "Distance from Clip Start");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

####################################################################################
####################################################################################
####################################################################################

makePlot.NVC.tail.clip.matchByClipPosition <- function(plotter, clip.amt = 10,  r2.buffer = clip.amt / 10, 
                                                       label.majority.bases = TRUE, label.majority.bases.threshold = 0.5, label.majority.bases.cex = 1, 
                                                       rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                                       show.base.legend = TRUE, debugMode = DEFAULTDEBUGMODE,  singleEndMode = plotter$res@singleEnd, 
                                                       plot = TRUE,...){

  res <- plotter$res
  isPlottable <- !(is.null(plotter$res@qc.data[["NVC.lead.clip.R1"]]) )
  if(! plot) return(isPlottable);

  plot.name <- "Lead Clipping NVC";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    if(rasterize.plotting.area){
      check.rasterize.or.die("rasterize.plotting.area");
    }
    res <- plotter$res;
    #print(clip.amt)
    #print(r2.buffer)
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else if(singleEndMode){
      xlim <- c(1,clip.amt)
      ylim <- c(0,1);
      bases <- names(plotter$nvc.colors);

      get.data.list <- function(dl){
         lapply(dl, function(df){
            out <- data.frame(base = rep(bases,clip.amt), 
                              CT = unlist(lapply(1:clip.amt, function(i){
                                      sapply(bases, function(b){
                                         sum( df$CT[ df$base == b & df$tailClipLen <= clip.amt & ((df$tailClipLen) - (df$readLen - df$readPos - 1) == i)] );
                                      });
                                   })),
                              stringsAsFactors = FALSE);
            #names(out) <- c("base","CT");
            return(out);
         })
      }
      #data.list.r1 <- get.data.list(plotter$res@qc.data[["NVC.tail.clip.R1"]]);
      data.list.r1 <- get.data.list(plotter$res@calc.data[["NVC.tail.clip.R1"]]);
      makePlot.generic.NVC.single(plot.name,    data.list.r1, 
                               plotter = plotter, 
                               x.name = "readPos",
                               y.name = "CT",
                               xlim=xlim,
                               ylim=ylim,
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = "-",label.major.cex = label.majority.bases.cex,
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                               ...);
      internal.plot.main.title(paste0("Clipped Nucleotide Rate, Trailing Clipped bases\n(1-",clip.amt,"), by distance to clip start"), plotter, plot.type = "NVC", ...);
      title(xlab = "Distance from Clip Start");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    } else {
      xlim <- c(1,clip.amt)
      ylim <- c(0,1);
      bases <- names(plotter$nvc.colors);

      get.data.list <- function(dl){
         lapply(dl, function(df){
            out <- data.frame(base = rep(bases,clip.amt), 
                              CT = unlist(lapply(1:clip.amt, function(i){
                                      sapply(bases, function(b){
                                         sum( df$CT[ df$base == b & df$tailClipLen <= clip.amt & ((df$tailClipLen) - (df$readLen - df$readPos - 1) == i)] );
                                      });
                                   })),
                              stringsAsFactors = FALSE);
            #names(out) <- c("base","CT");
            return(out);
         })
      }
      #data.list.r1 <- get.data.list(plotter$res@qc.data[["NVC.tail.clip.R1"]]);
      #data.list.r2 <- get.data.list(plotter$res@qc.data[["NVC.tail.clip.R2"]]);
      data.list.r1 <- get.data.list(plotter$res@calc.data[["NVC.tail.clip.R1"]]);
      data.list.r2 <- get.data.list(plotter$res@calc.data[["NVC.tail.clip.R2"]]);
      
      makePlot.generic.NVC.pair(plot.name,    data.list.r1, 
                               data.list.r2, 
                               plotter = plotter, 
                               x.name = "readPos",
                               y.name = "CT",
                               xlim=xlim,
                               ylim=ylim,
                               r2.buffer = r2.buffer,
                               label.major = label.majority.bases, label.cutoff = label.majority.bases.threshold, blank.label.char = "-",label.major.cex = label.majority.bases.cex,
                               show.base.legend = show.base.legend, nvc.colors = plotter$nvc.colors, nvc.colors.light = plotter$nvc.colors.light,
                               rasterize.plotting.area = rasterize.plotting.area, raster.height = raster.height, raster.width = raster.width, debugMode = debugMode,
                               ...);
      internal.plot.main.title(paste0("Clipped Nucleotide Rate, Trailing Clipped bases\n(1-",clip.amt,"), by distance to clip start"), plotter, plot.type = "NVC", ...);
      internal.plot.read.labels(plotter,"bottom",r2.buffer,xlim[2]);
      title(xlab = "Distance from Clip Start");
      title(ylab = "Nucleotide Rate");
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}


####################################################################################
####################################################################################
####################################################################################

makePlot.gene.assignment.rates <- function(plotter, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, plot = TRUE,...) {
  res <- plotter$res
  isPlottable <- !(is.null(plotter$res@qc.data[["summary"]]) | ! all(c("ReadPairs_UniqueGene","ReadPairs_AmbigGene","ReadPairs_NoGene","ReadPairs_NoGene_Intron","ReadPairs_NoGene_OneKbFromGene","ReadPairs_NoGene_TenKbFromGene", "ReadPairs_NoGene_MiddleOfNowhere","READ_PAIR_OK" ) %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);

  plot.name <- "Gene Assignment Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        #generic.points.tf <- function(data.list, x.names, x.titles, norm.by = NULL, field.name = "FIELD", count.name = "COUNT", offset.horiz = 0, horiz.offsets = NULL){
        #makePlot.generic.points <- function(tf.list, plotter, plot.type = "points", ylim = NULL, offset.horiz = 0, ...){
        pre.plot.func <- function(){
             #abline(v=1.5,col="gray",lty=3);
             abline(v=2.5,col="gray",lty=3);
             abline(v=3.5,col="gray",lty=3);
             abline(h=0.0,col="gray",lty=3);
        }

        tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("ReadPairs_UniqueGene","ReadPairs_UniqueGene_UTR","ReadPairs_AmbigGene","ReadPairs_NoGene","ReadPairs_NoGene_Intron","ReadPairs_NoGene_OneKbFromGene","ReadPairs_NoGene_TenKbFromGene", "ReadPairs_NoGene_MiddleOfNowhere"),
                                     x.titles = c("Unique\nGene","Unique\nGene\nUTR","Ambig\nGene","No Gene","No Gene\nIntronic","No Gene\n1kb\nFrom\nGene","No Gene\n10kb\nFrom\nGene","No Gene\nMiddle\nOf\nNowhere"),
                                     norm.by = "READ_PAIR_OK",
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name, tf.list,plotter, pre.plot.func = pre.plot.func,...)
      #title(xlab="Rate");
      title(ylab="Rate");
      internal.plot.main.title("Read Mapping Location Rates", plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}
####  if(is.null(plotter$res@qc.data[["summary"]]) | ! all(c() %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) ){
makePlot.splice.junction.loci.counts <- function(plotter, calc.rate = FALSE, high.low.cutoff = 4, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["summary"]]) | ! all(c("SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads","SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);

  plot.name <- "Splice Junction Loci Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        norm.by <- NULL;
        if(calc.rate){
          norm.by <- c("SpliceLoci_Covered");
        }

        low.label <- paste0("1-",high.low.cutoff-1,"\nReads");
        high.label <- paste0(high.low.cutoff,"+\nReads");
        hl.label <- c(low.label,high.label);

          pre.plot.func <- function(){
             abline(v=2.5,col="gray",lty=3);
             abline(v=4.5,col="gray",lty=3);
             categoryText <- c("All Loci","Known SJ Loci","Novel SJ Loci");
             categoryText.cex <- fit.character.vector(categoryText, min.width = 1.2, max.width = 1.9, max.width.per.char = 0.3);

             text(c(1.5,3.5,5.5),par("usr")[4],categoryText, cex = categoryText.cex, adj=c(0.5,0.95));
          }

        tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("SpliceLoci_Known_Covered","SpliceLoci_Novel","SpliceLoci_Known_FewReads","SpliceLoci_Known_ManyReads","SpliceLoci_Novel_FewReads","SpliceLoci_Novel_ManyReads"),
                                     x.titles = c("Known\nTotal","Novel\nTotal", paste0("",hl.label),paste0("",hl.label)),
                                     norm.by = norm.by,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list,plotter,pre.plot.func=pre.plot.func,...)

      if(calc.rate){
        title(ylab="Proportion of Splice Loci");
        internal.plot.main.title("Observed Splice Junction Loci, by type", plotter, ...);
      } else {
        title(ylab="# Splice Loci");
        internal.plot.main.title("# Observed Splice Junction Loci, by type", plotter, ...);
      }
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.splice.junction.event.ratesPerRead <- function(plotter, high.low.cutoff = 4, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["summary"]]) | ! all(c("SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads","SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);

  plot.name <- "Splice Junction Event Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    calc.ratePerRead = TRUE;
    calc.proportion = FALSE;

    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(! isPlottable){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        low.label <- paste0("1-",high.low.cutoff-1,"\nReads");
        high.label <- paste0(high.low.cutoff,"+\nReads");
        hl.label <- c(low.label,high.label);
        norm.by <- NULL;
        if(calc.ratePerRead){
          norm.by <- "READ_PAIR_OK";
        } else if(calc.proportion){
          norm.by <- "SpliceEvents";
        }
        pre.plot.func <- function(){
             abline(v=2.5,col="gray",lty=3);
             abline(v=4.5,col="gray",lty=3);
             categoryText <- c("All SJ Loci","Known SJ Loci","Novel SJ Loci");
             categoryText.cex <- fit.character.vector(categoryText, min.width = 1.2, max.width = 1.9, max.width.per.char = 0.3);
             text(c(1.5,3.5,5.5),par("usr")[4],categoryText, cex = categoryText.cex, adj=c(0.5,0.95));
        }
        tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("SpliceEvents_KnownLoci","SpliceEvents_NovelLoci","SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads","SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads"),
                                     x.titles = c("Known","Novel", paste0("",hl.label),paste0("",hl.label) ),
                                     norm.by = norm.by,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list,plotter,pre.plot.func=pre.plot.func,...)

      title(ylab="Events per Read-Pair");
      internal.plot.main.title("Splice Junction Event Rates per Read-Pair", plotter, ...);

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}


makePlot.splice.junction.event.proportions <- function(plotter, high.low.cutoff = 4, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["summary"]]) | ! all(c("SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads","SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);
  

  plot.name <- "Splice Junction Event Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    calc.ratePerRead = FALSE;
    calc.proportion = TRUE;

    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        low.label <- paste0("1-",high.low.cutoff-1,"\nReads");
        high.label <- paste0(high.low.cutoff,"+\nReads");
        hl.label <- c(low.label,high.label);

        norm.by <- NULL;
        if(calc.ratePerRead){
          norm.by <- "READ_PAIR_OK";
        } else if(calc.proportion){
          norm.by <- "SpliceEvents";
        }

        pre.plot.func <- function(){
             abline(v=2.5,col="gray",lty=3);
             abline(v=4.5,col="gray",lty=3);
             categoryText <- c("All SJ Loci","Known SJ Loci","Novel SJ Loci");
             categoryText.cex <- fit.character.vector(categoryText, min.width = 1.2, max.width = 1.9, max.width.per.char = 0.3);

             text(c(1.5,3.5,5.5),par("usr")[4],categoryText, cex = categoryText.cex, adj=c(0.5,0.95));
        }

        tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("SpliceEvents_KnownLoci","SpliceEvents_NovelLoci","SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads","SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads"),
                                     x.titles = c("Known","Novel", paste0("",hl.label),paste0("",hl.label) ),
                                     norm.by = norm.by,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
                                     
      #print("1:");
      #print(tf.list[[1]]);
      #print("2:");
      #print(tf.list[[2]]);
      
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list,plotter,pre.plot.func=pre.plot.func,...)


      title(ylab="Proportion of all splice events");
      internal.plot.main.title("Breakdown of Splice Junction Events", plotter, ...);

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}


makePlot.splice.junction.event.counts <- function(plotter, high.low.cutoff = 4, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["summary"]]) | ! all(c("SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads","SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);
  
  plot.name <- "Splice Junction Event Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    calc.ratePerRead = FALSE;
    calc.proportion = FALSE;

    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        low.label <- paste0("1-",high.low.cutoff-1,"\nReads");
        high.label <- paste0(high.low.cutoff,"+\nReads");
        hl.label <- c(low.label,high.label);

        norm.by <- NULL;
        if(calc.ratePerRead){
          norm.by <- "READ_PAIR_OK";
        } else if(calc.proportion){
          norm.by <- "SpliceEvents";
        }

        pre.plot.func <- function(){
             abline(v=2.5,col="gray",lty=3);
             abline(v=4.5,col="gray",lty=3);
             categoryText <- c("All SJ Loci","Known SJ Loci","Novel SJ Loci");
             categoryText.cex <- fit.character.vector(categoryText, min.width = 1.2, max.width = 1.9, max.width.per.char = 0.3);

             text(c(1.5,3.5,5.5),par("usr")[4],categoryText, cex = categoryText.cex, adj=c(0.5,0.95));
        }

        tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("SpliceEvents_KnownLoci","SpliceEvents_NovelLoci","SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads","SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads"),
                                     x.titles = c("Known","Novel", paste0("",hl.label),paste0("",hl.label) ),
                                     norm.by = norm.by,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list,plotter,pre.plot.func=pre.plot.func,...)

      title(ylab="# Splice Events");
      internal.plot.main.title("# Observed Splice Events, by type", plotter, ...);

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}

makePlot.splice.junction.event.proportionsByType <- function(plotter, high.low.cutoff = 4, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, plot = TRUE,...){
  res <- plotter$res 
  isPlottable <- ! (is.null(plotter$res@qc.data[["summary"]]) | ! all(c("SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads","SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);
    
  plot.name <- "Splice Junction Event Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    calc.rate = TRUE
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        low.label <- paste0("1-",high.low.cutoff-1,"\nReads");
        high.label <- paste0(high.low.cutoff,"+\nReads");
        hl.label <- c(low.label,high.label);

        norm.by <- NULL;
        if(calc.rate) norm.by <- "SpliceEvents";
        tf.list.1 <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("SpliceEvents_KnownLoci","SpliceEvents_NovelLoci"),
                                     x.titles = c("Known","Novel"),
                                     norm.by = norm.by,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
        if(calc.rate) norm.by <- "SpliceEvents_KnownLoci";
        tf.list.2 <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("SpliceEvents_KnownLociWithFewReads","SpliceEvents_KnownLociWithManyReads"),
                                     x.titles = hl.label,
                                     norm.by = norm.by,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets, 
                                     x.start = 2);
        if(calc.rate) norm.by <- "SpliceEvents_NovelLoci";
        tf.list.3 <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("SpliceEvents_NovelLociWithFewReads","SpliceEvents_NovelLociWithManyReads"),
                                     x.titles = hl.label,
                                     norm.by = norm.by,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets,
                                     x.start = 4);
        tf.list <- lapply(1:length(tf.list.1), function(i){
            rbind(tf.list.1[[i]],tf.list.2[[i]], tf.list.3[[i]]);
        })
        names(tf.list) <- names(tf.list.1);

        pre.plot.func <- function(){
             abline(h=0.0,col="gray",lty=3);

             abline(v=2.5,col="gray",lty=3);
             abline(v=4.5,col="gray",lty=3);
             categoryText <- c("All SJ Loci","Known SJ Loci","Novel SJ Loci");
             categoryText.cex <- fit.character.vector(categoryText, min.width = 1.2, max.width = 1.9, max.width.per.char = 0.3);

             text(c(1.5,3.5,5.5),par("usr")[4],categoryText, cex = categoryText.cex, adj=c(0.5,0.95));
        }
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list,plotter,pre.plot.func=pre.plot.func,...)
           #usr <- par("usr");
           #text(1.5,usr[4],"All Events",adj=c(0.5,1.1));
           #text(3.5,usr[4],"Known SJ Loci",adj=c(0.5,1.1));
           #text(5.5,usr[4],"Novel SJ Loci",adj=c(0.5,1.1));


      title(ylab="Proportion of Events, by type");
      internal.plot.main.title("Breakdown of Splice Junction Events, by type", plotter, ...);

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}


makePlot.strandedness.test <- function(plotter, plot.target.boxes = FALSE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["summary"]]) | ! all(c("StrandTest_frFirstStrand","StrandTest_frSecondStrand") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);
    
  plot.name <- "Strandedness Test";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        pre.plot.func <- function(){
          abline(h=0.5,col="gray",lty=3);
          #abline(h=0,col="gray",lty=3);
          #abline(h=1,col="gray",lty=3);
        }
        tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("StrandTest_frFirstStrand","StrandTest_frSecondStrand"),
                                     x.titles = c("fr_firstStrand","fr_secondStrand"),
                                     norm.by = c("StrandTest_frFirstStrand","StrandTest_frSecondStrand"),
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list,plotter,ylim=c(0,1),pre.plot.func=pre.plot.func,...)
      title(ylab="Rate");
      internal.plot.main.title("Strandedness Test", plotter, ...);
      if(plot.target.boxes){
        stranded.rule.code <- sapply(plotter$res@qc.data[["summary"]], function(dl){ dl$COUNT[ dl$FIELD == "Stranded_Rule_Code"] });
        if(all(stranded.rule.code == stranded.rule.code[1])){
           #Draw boxes around the target areas:
           if(stranded.rule.code[1] == 0){
             rect(0.6, 0.475, 1.4, 0.525, lwd = 2, lty = 2, border = "Green");
             rect(1.6, 0.475, 2.4, 0.525, lwd = 2, lty = 2, border = "Green");
           } else if(stranded.rule.code[1] == 1){
             rect(0.6, 0.95, 1.4, 1.035, lwd = 2, lty = 2,  border = "Green");
             rect(1.6, -0.035, 2.4, 0.05, lwd = 2, lty = 2, border = "Green");       
           } else if(stranded.rule.code[1] == 2){
             rect(0.6, -0.035, 1.4, 0.05, lwd = 2, lty = 2, border = "Green");
             rect(1.6, 0.95, 2.4, 1.035, lwd = 2, lty = 2,  border = "Green");
           }
        } else {
           #do nothing, stranded rule is not consistant.
        }
      }
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
}})}

makePlot.dropped.rates <- function(plotter, dropAlwaysZeroRows = FALSE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE,...){
  res <- plotter$res
    x.names <- c(#"DROPPED_UNALIGNED",
                 "DROPPED_NOT_PROPER_PAIR",
                 "DROPPED_READ_FAILS_VENDOR_QC",
                 "DROPPED_MARKED_NOT_VALID", 
                 "DROPPED_CHROMS_MISMATCH", 
                 "DROPPED_PAIR_STRANDS_MISMATCH", 
                 "DROPPED_IGNORED_CHROMOSOME", 
                 "DROPPED_NOT_UNIQUE_ALIGNMENT") ;
    x.titles <- c(#"Unaligned",
                  "Not\nProper\nPair",      
                  "Fails\nVendor\nQC", 
                  "Marked\nNot\nValid",
                  "Chroms\nMismatch", 
                  "Pair\nStrands\nDisagree",
                  "On\nIgnored\nChrom",
                  "Multi\nMapped");

  isPlottable <- ! (is.null(plotter$res@qc.data[["summary"]]) | ! all(c(x.names) %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);
  
  plot.name <- "Dropped Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;

    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c(x.names),
                                     x.titles = c(x.titles),
                                     norm.by = c("TOTAL_READ_PAIRS"),
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);

        keep.rows <- 1:nrow(tf.list[[1]]);
        for(i in 1:nrow(tf.list[[1]])){
          currY <- sapply(tf.list, function(tf){ tf$y[i] });
          if(all(currY == -1)){
            keep.rows <- keep.rows[keep.rows != i];
          } else if(dropAlwaysZeroRows & all(currY == 0)){
            keep.rows <- keep.rows[keep.rows != i];
          } else if(! all(currY != -1)){
            message("WARNING WARNING WARNING: read drop rules for ",tf.list[[1]]$x.names[i]," are not consistent! Samples may not be comparable!");
          }
        }
        tf.list <- lapply(tf.list, function(tf){
          tf[keep.rows, , drop=F];
        });
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      if(nrow(tf.list[[1]]) == 0){
        blank.plot("No Reads Dropped");
      } else {
        makePlot.generic.points(plot.name,tf.list,plotter,...)
        title(ylab="Rate");
      }
      internal.plot.main.title("Read Drop Rate, by Reason", plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}


#unique.ID  total.reads mapped.reads  mapping.rate  multi.mapped.reads  multi.mapping.rate
#c("total.reads","mapped.reads","mapping.rate","mm.reads","mm.rate")
#names(plotter$res@calc.data)
makePlot.mapping.rates <- function(plotter, plot.mm = NULL, y.counts.in.millions = TRUE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
    res <- plotter$res
  isPlottable <- ! is.null(plotter$res@calc.data[["map.rates"]])
  if(! plot) return(isPlottable);
    
  plot.name <- "Mapping Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){
      message("WARNING: Skipping ",plot.name," plotting. Data not found!\n",
              "    This may because input.read.pair.count was not included in the decoder.\n",
              "    Unless you supply the input (pre-alignment) read-pair count,\n",
              "    QoRTs cannot calculate mapping rates.");
      blank.plot(c(plot.name,"Data Not Found","Skipping..."));
    } else {
      if(is.null(plot.mm)){
        if("mm.reads" %in% plotter$res@calc.data[["map.rates"]]$FIELD){
          message("    Note: Multimapping data not found.");
          plot.mm <- FALSE;
        } else {
          plot.mm <- TRUE;
        }
        #plot.mm <- ! is.null(plotter$res@calc.data[["map.rates"]]$mm.reads);
      } else if(plot.mm){
         if(is.null(plotter$res@calc.data[["map.rates"]]$mm.reads)){
           message("WARNING: plot.mm is TRUE but multimapping rates not found!");
           plot.mm <- FALSE;
         }
      }
      oldmar <- par("mar");
      if(plot.mm){
         tryCatch({
           #par(mar = c(5,4,4,4) + 0.1);
           old.mar <- par("mar");
           curr.mar <- par("mar");
           curr.mar[4] <- curr.mar[2];
           par(mar = curr.mar);
           
           tf.list <- generic.points.tf(plotter$res@calc.data[["map.rates"]], 
                                        x.names = c("total.reads","mapped.reads","mm.reads"),
                                        x.titles = c("Input\nRead\nPairs","Uniquely\nMapped\nPairs","Multi\nMapped\nPairs"),
                                        norm.by = NULL,
                                        offset.horiz = 0.5,
                                        horiz.offsets = plotter$lanebam.params$horiz.offsets);

           tf.list.2 <- generic.points.tf(plotter$res@calc.data[["map.rates"]], 
                                          x.names = c("mapped.reads","mm.reads"),
                                          x.titles = c("Uniquely\nMapping\nRate","Multi\nMapping\nRate"),
                                          norm.by = c("total.reads"),
                                          offset.horiz = 0.5,
                                          horiz.offsets = plotter$lanebam.params$horiz.offsets);
         }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
         makePlot.generic.points(plot.name,tf.list,plotter, leave.blank.cols = 4, label.y = F,...);
         makePlot.generic.points.right(plot.name,tf.list.2, plotter, section.offset = 4.5, ylim.data = c(0,1), ...);
         usr <- par("usr");
         yrange <- abs(usr[4] - usr[3]);
         rect(4,usr[3] - yrange*0.02,4.5, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
         rect(4,usr[3] - yrange     ,4.5, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);
         
         text(2,usr[4],"Read Count",adj=c(0.5,1.1), ...);
         text(6,usr[4],"Rates",adj=c(0.5,1.1), ...);
         
         axis(4, at = (yrange / 2) + usr[3], labels=c("Rate"),col="transparent", line=1, ...);
         axis(2, at = (yrange / 2) + usr[3], labels=c("Read Counts"),col="transparent", line=1.5, ...);
         par(mar = old.mar);
      } else {
         tryCatch({
           old.mar <- par("mar");
           curr.mar <- par("mar");
           curr.mar[4] <- curr.mar[2];
           par(mar = curr.mar);
           #par(mar = c(5,4,4,4) + 0.1);
           tf.list <- generic.points.tf(plotter$res@calc.data[["map.rates"]], 
                                        x.names = c("total.reads","mapped.reads"),
                                        x.titles = c("Input\nRead\nPairs","Uniquely\nMapped\nPairs"),
                                        norm.by = NULL,
                                        offset.horiz = 0.5,
                                        horiz.offsets = plotter$lanebam.params$horiz.offsets);

           tf.list.2 <- generic.points.tf(plotter$res@calc.data[["map.rates"]], 
                                          x.names = c("mapped.reads"),
                                          x.titles = c("Uniquely\nMapping\nRate"),
                                          norm.by = c("total.reads"),
                                          offset.horiz = 0.5,
                                          horiz.offsets = plotter$lanebam.params$horiz.offsets);
         }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
         makePlot.generic.points(plot.name,tf.list,plotter, leave.blank.cols = 3, label.y = F,...);
         makePlot.generic.points.right(plot.name,tf.list.2, plotter, section.offset = 3.5, ylim.data = c(0,1));
         usr <- par("usr");
         yrange <- usr[4] - usr[3];
         rect.y <- c(usr[3] - yrange, usr[4] + yrange);
         rect(3,rect.y[1],3.5, rect.y[2],col="white", ...);
         rect(3,usr[3] - yrange*0.02,4.5, usr[4] + yrange*0.02, border = "white",col="white", xpd = TRUE,...);
         rect(3,usr[3] - yrange     ,3.5, usr[4] + yrange     , border = "black",col="white", xpd = FALSE,...);
         
         
         text(1.5,usr[4],"Read Counts",adj=c(0.5,1.1), ...);
         text(4.5,usr[4],"Rates",adj=c(0.5,1.1), ...);
         axis(4, at = (yrange / 2) + usr[3], labels=c("Rate"),col="transparent", line=1, ...);
         axis(2, at = (yrange / 2) + usr[3], labels=c("Read Counts"),col="transparent", line=1.5, ...);
         
         par(mar = old.mar);
      }
      if(y.counts.in.millions){
        usr <- par("usr");
        pretty.y <- pretty(c(usr[3],usr[4]))
        pretty.y.labels <- pretty.y;
        pretty.y.labels <- paste0(pretty.y / 1000000, "M");
        axis(2,at=pretty.y,labels=pretty.y.labels,las=1, lwd = -1, lwd.ticks = 1, ...);
      } else {
         axis(2, lwd = -1, lwd.ticks = 1, ...)
      }
      abline(h=0,col="grey",lty=3, ...);
        
      internal.plot.main.title("Mapping Stats", plotter, ...);
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      par(mar = oldmar);
    }
  })
}


UCSC.chrom.category.def <- function(chromName){
   if(is.null(chromName)){
      return(c("Autosomes","X","Y","XY","MT","Other"));
   } else {
      chromID <- substr(chromName,4,nchar(chromName));
      if(chromID == "M"){
         return("MT");
      } else if(chromID == "X") {
         return("X");
      } else if(chromID == "Y"){
         return("Y");
      } else if(chromID == "XY"){
         return("XY");
      } else if(! is.na(   suppressWarnings(as.numeric(chromID))   ))  {
         return("Autosomes");
      } else {
         return("Other");
      }
   }
}
ENSEMBL.chrom.category.def <- function(chromName){
   if(is.null(chromName)){
      return(c("Autosomes","X","Y","XY","MT","Other"));
   } else {
      if(chromName == "M"){
         return("MT");
      } else if(chromName == "X") {
         return("X");
      } else if(chromName == "Y"){
         return("Y");
      } else if(chromName == "XY"){
         return("XY");
      } else if(! is.na(suppressWarnings(as.numeric(chromName)))){
         return("Autosomes");
      } else {
         return("Other");
      }
   }
}
UCSC_WITH_ERCC.chrom.category.def <- function(chromName){
   if(is.null(chromName)){
      return(c("Autosomes","X","Y","XY","MT","ERCC","Other"));
   } else {
      chromID <- substr(chromName,4,nchar(chromName));
      if(chromID == "M"){
         return("MT");
      } else if(chromID == "X") {
         return("X");
      } else if(chromID == "Y"){
         return("Y");
      } else if(chromID == "XY"){
         return("XY");
      } else if(! is.na(   suppressWarnings(as.numeric(chromID))   ))  {
         return("Autosomes");
      } else if(length(chromName) > 4) {
         if(substr(chromName,1,4) == "ERCC"){
           return("ERCC");
         }
      }
      return("Other");
   }
}
ENSEMBL_WITH_ERCC.chrom.category.def <- function(chromName){
   if(is.null(chromName)){
      return(c("Autosomes","X","Y","XY","MT","ERCC","Other"));
   } else {
      if(chromName == "MT"){
         return("MT");
      } else if(chromName == "X") {
         return("X");
      } else if(chromName == "Y"){
         return("Y");
      } else if(chromName == "XY"){
         return("XY");
      } else if(! is.na(   suppressWarnings(as.numeric(chromName))   ))  {
         return("Autosomes");
      } else if(length(chromName) > 4) {
         if(substr(chromName,1,4) == "ERCC"){
           return("ERCC");
         }
      }
      return("Other");
   }
}

DEFAULT_CHROM_CATEGORY_DEF_FUNCTIONS = list(UCSC = UCSC.chrom.category.def,
                                            UCSC_WITH_ERCC = UCSC_WITH_ERCC.chrom.category.def,
                                            ENSEMBL = ENSEMBL.chrom.category.def,
                                            ENSEMBL_WITH_ERCC = ENSEMBL_WITH_ERCC.chrom.category.def  );


makePlot.chrom.type.rates <- function(plotter, 
                                  plot.rates = TRUE,
                                  chromosome.name.style = "UCSC",
                                  exclude.autosomes = FALSE,
                                  chrom.norm.factors = NULL,
                                  custom.chromosome.style.def.function = NULL,
                                  return.table = FALSE, debugMode = DEFAULTDEBUGMODE,  singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
                                      
  res <- plotter$res
  isPlottable <- ! is.null(plotter$res@qc.data[["chrom.counts"]])
  if(! plot) return(isPlottable);
  
  plot.name <- "Chromosome Rates";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if( !isPlottable ){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(is.null(custom.chromosome.style.def.function)){
          if(chromosome.name.style %in% names(DEFAULT_CHROM_CATEGORY_DEF_FUNCTIONS)){
            chrom.style.def.fcn <- DEFAULT_CHROM_CATEGORY_DEF_FUNCTIONS[[chromosome.name.style]];
          } else {
            stop(paste("Unrecognized chromosome name style! chromosome.name.style must be one of: \n        ",paste0("\"",names(DEFAULT_CHROM_CATEGORY_DEF_FUNCTIONS),"\"",collapse= ","), sep=""));
          }
        } else {
          chrom.style.def.fcn <- custom.chromosome.style.def.function;
        }
        chrom.cats <- chrom.style.def.fcn(NULL);

        #Remove categories not found in the dataset:
        #chrom.list <- plotter$res@qc.data[["chrom.counts"]][[1]]$CHROM;
        #chrom.cats <- chrom.cats[chrom.cats %in% sapply(chrom.list, chrom.style.def.fcn)]


        tf.list <- lapply(1:length(plotter$res@qc.data[["chrom.counts"]]), function(i){
           in.df <- plotter$res@qc.data[["chrom.counts"]][[i]];

           df <- data.frame( x = 1,y = sum(in.df$CT[ sapply(in.df$CHROM,function(x){ chrom.style.def.fcn(x) == chrom.cats[1] } ) ]),x.titles = chrom.cats[1], stringsAsFactors=F); 

           for(j in 2:length(chrom.cats)){
              chrom.cat <- chrom.cats[j];
              df <- rbind(df, c(j,sum(in.df$CT[ sapply(in.df$CHROM,function(x){ chrom.style.def.fcn(x) == chrom.cat } ) ]),chrom.cat)); 
           }

           df$y <- as.numeric(df$y);
           df$x <- as.numeric(df$x) ;
           df$x <- df$x + plotter$lanebam.params$horiz.offsets[i] * 0.5
           if(plot.rates) {
             df$y <- 100 * df$y / sum(df$y);
           }
           df;
        });
        names(tf.list) <- names(plotter$res@qc.data[["chrom.counts"]]);

        for(i in 1:length(chrom.cats)){
          chrom.cat <- chrom.cats[i];
          if(chrom.cat == "Autosomes" & exclude.autosomes){
            tf.list <- lapply(tf.list, function(tf.xy){ tf.xy[tf.xy$x.titles != chrom.cat,,drop=F] });
          } else if(all(sapply(tf.list,function(tf.xy){ tf.xy$y[tf.xy$x.titles == chrom.cat] == 0 }))){
            tf.list <- lapply(tf.list, function(tf.xy){ tf.xy[tf.xy$x.titles != chrom.cat,,drop=F] });
          }
        }
        if(dim(tf.list[[1]])[1] == 0){
           stop("> makePlot.chrom.type.rates: All chromosomes excluded! Did zero reads align? Are you using inconsistent chromosome names? Have you set exclude.autosomes on a dataset with no non-autosomal reads?");
        }

        tf.list <- lapply(1:length(tf.list), function(i){
          tf.list[[i]]$x <- 1:length(tf.list[[i]]$x) + plotter$lanebam.params$horiz.offsets[i] * 0.5;
          tf.list[[i]];
        })
        names(tf.list) <- names(plotter$res@qc.data[["chrom.counts"]]);

        if(! is.null(chrom.norm.factors)){
          nonzero.chrom.cats <- tf.list[[1]]$x.titles;  
          chrom.cat.norm.factors <- as.list( rep(0,length(nonzero.chrom.cats) ));
          names(chrom.cat.norm.factors) <- nonzero.chrom.cats;
          for(i in 1:length(chrom.norm.factors)){
            nf <- chrom.norm.factors[[i]];
            chrom.name <- names(chrom.norm.factors)[i];
            chrom.cat <- chrom.style.def.fcn(chrom.name);
            if(chrom.cat %in% names(chrom.cat.norm.factors)){
              chrom.cat.norm.factors[[chrom.cat]] <- chrom.cat.norm.factors[[chrom.cat]] + nf;
            }
          }
          tf.list <- lapply(tf.list, function(tf.xy){
            tf.xy$y <- tf.xy$y / unlist(chrom.cat.norm.factors);
            tf.xy;
          });
        }

        pre.plot.func <- function(){
           abline(h=0.0,col="gray",lty=3);
        }
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list, plotter, pre.plot.func = pre.plot.func,...);

      if(plot.rates) {
          title(ylab="% Reads by Chromosome Type");
      } else {
          title(ylab="# Reads by Chromosome Type");
      }
      if(exclude.autosomes){
        internal.plot.main.title("Chromosome Distribution (Excluding Autosomes)", plotter, ...);
      } else {
        internal.plot.main.title("Chromosome Distribution", plotter, ...);
      }
      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      if(return.table){
        table.df <- as.data.frame(t( sapply(tf.list, function(tf.xy){
          tf.xy$y;
        }) ));
        names(table.df) <- tf.list[[1]]$x.titles;
        return(table.df);
      }
    }
  })
}

makePlot.norm.factors <- function(plotter, by.sample = TRUE, return.table = FALSE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res
  isPlottable <- ! is.null(plotter$res@calc.data[["norm.factors.bySample"]]) 
  if(! plot) return(isPlottable);

  plot.name <- "Normalization Factors";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if( !isPlottable){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(by.sample){
          norm.factor.data <- plotter$res@calc.data[["norm.factors.bySample.splitToLaneBam"]]
        } else {
          norm.factor.data <- plotter$res@calc.data[["norm.factors.byLaneBam"]]
        }

        pre.plot.func <- function(){
           abline(h=1.0,col="gray",lty=3);
        }

        tf.list <- generic.points.tf.from.df(df = norm.factor.data,
                                  x.names = c("Norm_TC","Norm_Geo","Norm_TMM","Norm_UQ","Norm_RLE"),
                                  x.titles = c("Total\nCount","Geometric\n(DESeq)","TMM\n(edgeR\nDefault)","UQ\n(edgeR)","RLE\n(edgeR)"),
                                  skip.missing = TRUE,
                                  offset.horiz = 0.5,
                                  horiz.offsets = plotter$lanebam.params$horiz.offsets
                                  )
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list,plotter,pre.plot.func=pre.plot.func,...);
      internal.plot.main.title("Normalization Factors", plotter, ...);
      title(ylab="Normalization Factor");

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      if(return.table){
        table.df <- as.data.frame(t( sapply(tf.list, function(tf.xy){
          tf.xy$y;
        }) ));
        names(table.df) <- tf.list[[1]]$x.titles;
        return(table.df);
      }
    }
  })
}

makePlot.norm.factors.vs.TC <- function(plotter, by.sample = TRUE, return.table = FALSE, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,plot = TRUE, ...){
  res <- plotter$res

  isPlottable <- ! is.null(plotter$res@calc.data[["norm.factors.bySample"]])
  if(! plot) return(isPlottable);
  
  plot.name <- "Normalization Factors";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if( !isPlottable ){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        if(by.sample){
          norm.factor.data <- plotter$res@calc.data[["norm.factors.bySample.splitToLaneBam"]]
        } else {
          norm.factor.data <- plotter$res@calc.data[["norm.factors.byLaneBam"]]
        }
        tf.list <- generic.points.tf.from.df(df = norm.factor.data,
                                  x.names = c("Norm_TC","Norm_Geo","Norm_TMM","Norm_UQ","Norm_RLE"),
                                  x.titles = c("Total\nCount","Geometric\n(DESeq)","TMM\n(edgeR\nDefault)","UQ\n(edgeR)","RLE\n(edgeR)"),
                                  skip.missing = TRUE,
                                  offset.horiz = 0.5,
                                  horiz.offsets = plotter$lanebam.params$horiz.offsets
                                  )

        tf.tf.list <- lapply(tf.list, function(tf.xy){
          out <- tf.xy[-1,];
          TC.y <- as.numeric(tf.xy$y[1]);
          out$y <- as.numeric(out$y) / TC.y;
          out$x <- out$x - 1;
          out;
        });
        pre.plot.func <- function(){
           abline(h=1.0,col="gray",lty=3);
        }
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.tf.list,plotter,pre.plot.func=pre.plot.func,...);
      internal.plot.main.title("Normalization Factors vs Total-Count Normalization", plotter, ...);
      title(ylab="Alt vs TC ratio");

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }

      if(return.table){
        table.df <- as.data.frame(t( sapply(tf.tf.list, function(tf.xy){
          tf.xy$y;
        }) ));
        names(table.df) <- tf.tf.list[[1]]$x.titles;
        return(table.df);
      }
    }
  })
}

#INCOMPLETE!!!!!!
#Placeholder for future functionality!
makePlot.cigarMismatch <- function(plotter,debugMode = DEFAULTDEBUGMODE,  singleEndMode = plotter$res@singleEnd,plot = TRUE,...){
    res <- plotter$res
  isPlottable <- ! (is.null(plotter$res@qc.data[["summary"]]) | ! all(c("CigChk_staggered","CigChk_noAlignedOverlap") %in% plotter$res@qc.data[["summary"]][[1]]$FIELD) )
  if(! plot) return(isPlottable);
    
  plot.name <- "Cigar Mismatch";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  if(debugMode){ message("Starting: ",plot.name," plot."); }
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){
      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      #Do stuff:
      tryCatch({
        pre.plot.func <- function(){
           abline(h=1.0,col="gray",lty=3);
        }

        tf.list <- generic.points.tf(plotter$res@qc.data[["summary"]], 
                                     x.names = c("StrandTest_frFirstStrand","StrandTest_frSecondStrand"),
                                     x.titles = c("fr_firstStrand","fr_secondStrand"),
                                     norm.by = c("StrandTest_frFirstStrand","StrandTest_frSecondStrand"),
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
      makePlot.generic.points(plot.name,tf.list,plotter,ylim=c(0,1),pre.plot.func=pre.plot.func,...)
      title(ylab="Rate");
      internal.plot.main.title("Strandedness Test", plotter, ...);

      if(debugMode){ message("Finished: ",plot.name," plot.",getTimeAndDiff(ts)); }
    }
  })
}
#INCOMPLETE! Temporary placeholder for future functionality
makePlot.splicing.mismatch <- function(plotter, debugMode = DEFAULTDEBUGMODE,  singleEndMode = plotter$res@singleEnd,plot = TRUE,... ){
  res <- plotter$res
  isPlottable <- ! is.null(plotter$res@qc.data[["insert.size"]])
  if(! plot) return(isPlottable);
  
  plot.name <- "Splice Mismatch";
  readLabel <- if(singleEndMode){ "Reads" } else {"Read-Pairs"}
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    res <- plotter$res;
    plot.name <- "Splicing Mismatch";
    if(debugMode){ message("Starting: ",plot.name," plot."); }
    if(debugMode){ ts <- timestamp() }
    if(!isPlottable){

      message(paste0("Warning: Skipping ",plot.name," plotting. Data not found!"));
      blank.plot(c(plot.name,"Data Not Found\nSkipping..."));
    } else {
      tryCatch({
        tf.list <- generic.points.tf(plot.name,plotter$res@qc.data[["summary"]], 
                                     x.names = c(""),
                                     x.titles = c(""),
                                     norm.by = NULL,
                                     offset.horiz = 0.5,
                                     horiz.offsets = plotter$lanebam.params$horiz.offsets);
      }, error = function(e){ errorPlot(plot.name,e, code = 0); stop(e);});
    }
  })
}

####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################


makePlot.legend.box <- function(plotter,debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd, cex.legend = NULL, ncol = NULL,plot = TRUE, ...){
  plot.name <- "legend box";
  res <- plotter$res
  if(!plot) return(TRUE);
  
  plotter.error.wrapper(plot.name, plotterFcn = function(){
    plotter.NVC <- plotter;
    plot(0,0,col="transparent",xlim=c(0,1),ylim=c(0,1), axes=F,xlab="",ylab="");

    #internal.plot.legend(plotter,"lines","bottomleft");
    #internal.plot.legend(plotter,"points","bottomright");
    if(is.null(plotter$title.annotations[["main.subtitle"]])){
      mainSub <- ""
    } else {
      mainSub <- paste0("\n",plotter$title.annotations[["main.subtitle"]]);
    }

    if(plotter$plot.type == "highlightSample.colorByLane"){
      title.text <- paste0("Sample Highlight:\n",plotter$title.highlight.name,"\nColored by Lane",mainSub)
      title.cex <- fit.character.vector(title.text);

      #text(0.5,1,title.text,adj=c(0.5,1),cex=title.cex, ...);
      #internal.plot.legend(plotter,"lnpt","bottom", ...);
      #internal.plot.legend(plotter.NVC,"NVC.both","center");
    } else if(plotter$plot.type == "highlightSample") {
      #internal.plot.legend(plotter.NVC,"NVC.both","center");
      title.text <- paste0("Sample Highlight:\n",plotter$title.highlight.name)
      title.cex <- fit.character.vector(title.text);

      #internal.plot.legend(plotter,"lnpt","bottom", ...);
      #text(0.5,1,title.text,adj=c(0.5,1),cex=title.cex, ...);
    } else if(plotter$plot.type == "colorBySample") {
      title.text <- paste0("Summary Plots, By Sample")
      title.cex <- fit.character.vector(title.text);
      #text(0.5,1,title.text,adj=c(0.5,1),cex=title.cex, ...);
      #internal.plot.legend(plotter,"lnpt","bottom", ...);
    } else if(plotter$plot.type == "summary"){
      title.text <- paste0("All\nQuality Control\nSummary Plots")
      title.cex <- fit.character.vector(title.text);
      text(0.5,0.5,title.text,adj=c(0.5,0.5),cex=title.cex, ...);
      return("");
    } else if(plotter$plot.type == "colorByLane"){
      title.text <- paste0("Summary Plots, By Lane")
      title.cex <- fit.character.vector(title.text);
      #text(0.5,1,title.text,adj=c(0.5,1),cex=title.cex, ...);
      #internal.plot.legend(plotter,"lnpt","bottom", ...);
    } else if(plotter$plot.type == "colorByGroup"){
      title.text <- paste0("Summary Plots, By Group")
      title.cex <- fit.character.vector(title.text);
      #text(0.5,1,title.text,adj=c(0.5,1),cex=title.cex, ...);
      #internal.plot.legend(plotter,"lnpt","bottom", ...);
    } else if(plotter$plot.type == "colorByX"){
      title.text <- paste0("Summary Plots, By ",plotter$title.highlight.name)
      title.cex <- fit.character.vector(title.text);
    } else {
      title.text <- paste0("Summary Plots, Custom");
      title.cex <- fit.character.vector(title.text);
    }
    text(0.5,1,title.text,adj=c(0.5,1),cex=title.cex, ...);
    
    if(is.null(cex.legend)){
      autofit.limits <- device.limits();
      autofit.limits[4] <- par("usr")[4] - strheight(paste0(title.text,"\n"),cex = title.cex);
      internal.plot.legend(plotter,"lnpt","bottom",autofit.limits=autofit.limits, ncol=ncol, ...);
    } else {
      internal.plot.legend(plotter,"lnpt","bottom", cex = cex.legend, ncol=ncol, ...);
    }
  })
}

makePlot.legend.over <- function(position, plotter, debugMode = DEFAULTDEBUGMODE, singleEndMode = plotter$res@singleEnd,ncol = NULL,plot = TRUE, ...) {
  res <- plotter$res
  plotter.NVC <- plotter;
  internal.plot.legend(plotter, "lnpt", position,ncol=ncol, ...);
}

get.summary.table <- function(res, outfile = NULL, debugMode = DEFAULTDEBUGMODE){
   out.fields <- lapply( res@qc.data[["summary"]], function(dl){
     out.dl <- dl$FIELD;
     out.dl;
   } );
   all.fields <- unique(unlist(out.fields))
   
   out.table <- sapply( res@qc.data[["summary"]], function(dl){
     dl$COUNT[match(all.fields,dl$FIELD)]
   });
   row.names(out.table) <- all.fields
   
   #out.table <- t(out.table);
   out.table <- data.frame(out.table,stringsAsFactors=F);
   if(! is.null(outfile)){
     write.table.with.first.col(out.table,file = outfile, rownamesTitle = "FIELD");
   }
   
   return(out.table);
}

get.size.factors <- function(res, 
                    sf.method = c("DESeq2","DESeq2_GEO","TC","edgeR","edgeR_TMM","edgeR_UQ","edgeR_RLE"), 
                    outfile=NULL, debugMode = DEFAULTDEBUGMODE){
  if(is.null(res@calc.data[["norm.factors.bySample"]])){
    message("Size factor table not found! (Perhaps DESeq2 or edgeR is not installed)?");
    return(NULL);
  } else {
    sf <- match.arg(sf.method);
    
    if(sf == "DESeq2" || sf == "DESeq2_GEO"){
      sfid <- "Norm_Geo";
    } else if(sf == "TC"){
      sfid <- "Norm_TC";
    } else if(sf == "edgeR" || sf == "edgeR_TMM"){
      sfid <- "Norm_TMM";
    } else if(sf == "edgeR_UQ"){
      sfid <- "Norm_UQ";
    } else if(sf == "edgeR_RLE"){
      sfid <- "Norm_RLE";
    }
    
    if(is.null(res@calc.data$norm.factors.bySample[[sfid]])){
      message("Size factor type ",sf," not found in size factor table! (Perhaps DESeq2 or edgeR is not installed)?");
      return(NULL);
    }
    
    sizeFactors <- data.frame(
      sample.ID = res@calc.data$norm.factors.bySample[["sample.ID"]], 
      size.factor = res@calc.data$norm.factors.bySample[[sfid]], 
      stringsAsFactors=F);
    
    if(! is.null(outfile)){
      write.table(sizeFactors,outfile, row.names=F,col.names=T,quote=F, sep='\t');
    }
    return(sizeFactors);
  }
}









