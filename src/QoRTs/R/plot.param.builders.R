
###########################################################################################################################################################################
###########################################################################################################################################################################
######################################################### EXTERNAL-FACING:
###########################################################################################################################################################################
###########################################################################################################################################################################

#highlight.by, color.highlighted.by, plot.type, offset.by = color.highlighted.by, merge.offset.outgroup = TRUE, pch.matched.with.color = TRUE

build.plotter.highlightSample.colorByLane <- function(curr.sample, res, plotter.params = list(), merge.offset.outgroup = TRUE){
   pch.matched.with.color <- TRUE
   if(! (curr.sample %in% res@decoder$sample.ID)){
      stop("FATAL ERROR! ","Cannot find sample ", curr.sample," in the sample list!");
   }
   base.defaultParams <- QoRTs.default.plotting.params;
   defaultParams <- merge.plotting.params(base.defaultParams,list());
   final.params <- merge.plotting.params(defaultParams,plotter.params);
   compiled.params <- compile.plotting.params(final.params);
   
   return(build.plotter.highlight.and.color(curr.highlight = curr.sample, 
                                            res = res, 
                                            compiled.params = compiled.params,
                                            highlight.by = res@decoder$sample.ID, 
                                            color.highlighted.by = res@decoder$lane.ID, 
                                            plot.type = "highlightSample.colorByLane", 
                                            merge.offset.outgroup = merge.offset.outgroup,
                                            pch.matched.with.color = pch.matched.with.color));
}

#########################################################

build.plotter.highlightSample <- function(curr.sample, res, plotter.params = list(), merge.offset.outgroup = TRUE){
   if(! (curr.sample %in% res@decoder$sample.ID)){
      stop("FATAL ERROR! ","Cannot find sample ", curr.sample," in the sample list!");
   }
   base.defaultParams <- QoRTs.default.plotting.params;
   defaultParams <- merge.plotting.params(base.defaultParams,list());
   final.params <- merge.plotting.params(defaultParams,plotter.params);
   compiled.params <- compile.plotting.params(final.params);
   
   return(build.plotter.highlight(curr.sample, 
                                  res = res, 
                                  compiled.params = compiled.params,
                                  highlight.by = res@decoder$sample.ID, 
                                  plot.type = "highlightSample", 
                                  offset.by = res@decoder$lane.ID, 
                                  merge.offset.outgroup = merge.offset.outgroup));
}

#########################################################

build.plotter.colorByLane <- function(res, plotter.params = list()){
   base.defaultParams <- QoRTs.default.plotting.params;
   defaultParams <- merge.plotting.params(base.defaultParams,list(std.lines.alpha = 125));
   final.params <- merge.plotting.params(defaultParams,plotter.params);
   compiled.params <- compile.plotting.params(final.params);
   
   return(build.plotter.color( 
                                            res = res, 
                                            compiled.params = compiled.params,
                                            color.by = res@decoder$lane.ID, 
                                            plot.type = "colorByLane",
                                            color.by.title.name = "Lane"));
}

#########################################################

build.plotter.colorByGroup <- function(res, plotter.params = list()){
   base.defaultParams <- QoRTs.default.plotting.params;
   defaultParams <- merge.plotting.params(base.defaultParams,list(std.lines.alpha = 125));
   final.params <- merge.plotting.params(defaultParams,plotter.params);
   compiled.params <- compile.plotting.params(final.params);
   
   return(build.plotter.color(
                                            res = res, 
                                            compiled.params = compiled.params,
                                            color.by = res@decoder$group.ID, 
                                            plot.type = "colorByGroup",
                                            color.by.title.name = "Group"));
}

#########################################################

build.plotter.colorBySample <- function(res, plotter.params = list()){
   base.defaultParams <- QoRTs.default.plotting.params;
   defaultParams <- merge.plotting.params(base.defaultParams,list(std.lines.alpha = 125));
   final.params <- merge.plotting.params(defaultParams,plotter.params);
   compiled.params <- compile.plotting.params(final.params);
   
   return(build.plotter.color( 
                                            res = res, 
                                            compiled.params = compiled.params,
                                            color.by = res@decoder$sample.ID, 
                                            plot.type = "colorBySample",
                                            color.by.title.name = "Sample"));
}

#########################################################

build.plotter.basic <- function(res, plotter.params = list()){
   base.defaultParams <- QoRTs.default.plotting.params;
   defaultParams <- merge.plotting.params(base.defaultParams,list(alt.lines.alpha = 150, alt.points.alpha = 150, alt.lines.color = "black", alt.lines.lty = 1));
   final.params <- merge.plotting.params(defaultParams,plotter.params);
   compiled.params <- compile.plotting.params(final.params);
   
   build.plotter.basic.helper(res,res@decoder$lane.ID, compiled.params = compiled.params ,plot.type = "summary");
}

#########################################################
#Advanced use:
#########################################################

build.plotter.colorByX <- function(res, color.by.name, color.by.title.name = color.by.name, plotter.params = list()){
  base.defaultParams <- QoRTs.default.plotting.params;
  defaultParams <- merge.plotting.params(base.defaultParams,list(std.lines.alpha = 125));
  final.params <- merge.plotting.params(defaultParams,plotter.params);
  compiled.params <- compile.plotting.params(final.params);
  
  if(! (color.by.name %in% names(res@decoder))){
    stop(paste0("Fatal error: Cannot find field \"",color.by.name,"\" in the decoder! Check the spelling and capitalization?"));
  }
  
   return(build.plotter.color(
                                            res = res, 
                                            compiled.params = compiled.params,
                                            color.by = res@decoder[[color.by.name]], 
                                            plot.type = "colorByX",
                                            color.by.title.name = color.by.title.name));
}

###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
######################################################### INTERNAL:
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################

build.plotter.basic.helper <- function(res, offset.by, compiled.params, plot.type){
  if(! is.null(check.isValid(res))){
    #THROW AN ERROR!
    stop("Error parsing results: ", check.isValid(res));
  }
  
  #plotParams <- new("QoRT_Plotter");
  #plotParams@res <- res;
  #plotParams@plot.type <- plot.type;
  title.highlight.name <- "Summary";
  showLegend <- FALSE;
  nvc.colors <- compiled.params@nvc.colors;
  nvc.colors.light <- compiled.params@nvc.colors;
  
  lanebam.ct <- length(offset.by);
  
  legend.params <- data.frame(name = c(),
                                         lines.col = c(),
                                         lines.lty = c(),
                                         points.pch = c(),
                                         points.col = c(),
                                         stringsAsFactors=F
  );
  
  offset.by.vals <- sort(unique(offset.by));
  offset.ct <- length(offset.by.vals);
  if(offset.ct %% 2 == 0){
    offsets <- (((1:offset.ct) / (offset.ct)) );
    offsets <- offsets - ( offsets[ offset.ct / 2  + 1] + offsets[ offset.ct / 2 ]  ) / 2
  } else {
    offsets <- (((1:offset.ct) / (offset.ct)) - 0.5  );
    offsets <- offsets - offsets[ (offset.ct + 1) / 2 ];
  }
  
  lanebam.offset.indices <- sapply(1:lanebam.ct, FUN=function(i){
    which(offset.by.vals == offset.by[i])
  });
  lanebam.offsets <- offsets[lanebam.offset.indices];

  lanebam.params <- data.frame( plot.priority = rep(2,lanebam.ct),
                                           unique.ID = res@decoder$unique.ID,
                                           lines.col = rep(compiled.params@highlight.color[1], lanebam.ct),
                                           points.col = rep(compiled.params@highlight.color[1], lanebam.ct),
                                           points.pch = rep(compiled.params@highlight.points.pch[1], lanebam.ct),
                                           lines.lty = rep(compiled.params@highlight.lines.lty[1], lanebam.ct),
                                           lines.lwd = rep(compiled.params@highlight.lines.lwd[1], lanebam.ct),
                                           lines.alpha = rep(compiled.params@highlight.lines.alpha[1], lanebam.ct),
                                           points.alpha = rep(compiled.params@highlight.points.alpha[1], lanebam.ct),
                                           horiz.offsets = lanebam.offsets,
                                           vert.offsets = lanebam.offsets,
                                           stringsAsFactors=F
  );
  
  plotParams <- generate.plotter(res = res, plot.type = plot.type, title.highlight.name = title.highlight.name, legend.params = legend.params, showLegend = showLegend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, lanebam.params = lanebam.params)
  return(plotParams);
}

##################################################################################################################
##################################################################################################################
##################################################################################################################

build.plotter.highlight.and.color <- function(curr.highlight, res, compiled.params, highlight.by, color.highlighted.by, plot.type, offset.by = color.highlighted.by, merge.offset.outgroup = TRUE, pch.matched.with.color = TRUE, highlighted.by.name = "Samples"){
  if(! is.null(check.isValid(res))){
    #THROW AN ERROR!
    stop("Error parsing results: ", check.isValid(res));
  }
  if(length(highlight.by) != length(color.highlighted.by)) stop("error!");
  if(length(highlight.by) != length(res@decoder$unique.ID)) stop("error!");
   if(! (curr.highlight %in% highlight.by)){
      stop("FATAL ERROR! ","Cannot find highlighted element ", curr.highlight," in the highlight list!");
   }
  
  #plotParams <- new("QoRT_Plotter");
  #plotParams@res <- res;
  #plotParams@plot.type <- plot.type;
  title.highlight.name <- curr.highlight;
  showLegend <- compiled.params@showLegend;
  nvc.colors <- compiled.params@nvc.colors;
  nvc.colors.light <- compiled.params@nvc.colors.light;
  
  is.highlighted <- highlight.by == curr.highlight;
  lanebam.ct <- length(highlight.by);
  hl.by.factor <- color.highlighted.by[is.highlighted ];
  hl.by.factor.levels <- sort(unique(hl.by.factor));
  
  if(length(hl.by.factor.levels) > length(compiled.params@by.colors)) {
    message("WARNING WARNING WARNING: Too many categories to color! (Max = ",length(compiled.params@by.colors),", curr = ",length(hl.by.factor.levels),") Add more colors? Falling back to a rainbow palette");
    compiled.params@by.colors <- rainbow(length(hl.by.factor.levels));
  } 
  if(pch.matched.with.color){
    if(length(hl.by.factor.levels) > length(compiled.params@by.pch)){
      message("WARNING WARNING WARNING: Too many categories to mark distinct with pch! (Max = ",length(compiled.params@by.pch),", curr = ",length(hl.by.factor.levels),") Add more colors? Falling back to a repeating pattern.");
      compiled.params@by.pch <- rep(compiled.params@by.pch, ceiling( length(hl.by.factor.levels) / length(compiled.params@by.pch)));
    }
    legend.points.pch <- compiled.params@by.pch[1:length(hl.by.factor.levels)];
  } else {
    legend.points.pch <- rep(compiled.params@highlight.points.pch[1], length(hl.by.factor.levels));
  }

  legend.params <- data.frame(name = c(hl.by.factor.levels, paste0("Other ",highlighted.by.name)),
                                         lines.col = c(compiled.params@by.colors[1:length(hl.by.factor.levels)], compiled.params@highlight.color[2]),
                                         lines.lty = c(rep(compiled.params@highlight.lines.lty[1],length(hl.by.factor.levels)), compiled.params@highlight.lines.lty[2]),
                                         points.pch = c(legend.points.pch, compiled.params@highlight.points.pch[2]),
                                         points.col = c(compiled.params@by.colors[1:length(hl.by.factor.levels)], compiled.params@highlight.color[2]),
                                         stringsAsFactors=F
  );
  
  lines.col <- sapply(1:lanebam.ct, FUN=function(i){
    if(! is.highlighted[i]){ return(compiled.params@highlight.color[2]);
    } else {return( legend.params$lines.col[ legend.params$name == color.highlighted.by[i] ] );}
  });
  points.col <- sapply(1:lanebam.ct, FUN=function(i){
    if(! is.highlighted[i]){ return(compiled.params@highlight.color[2]);
    } else {return( legend.params$points.col[ legend.params$name == color.highlighted.by[i] ] );}
  });

  if(pch.matched.with.color){
    #message("!!!");
    if(length(hl.by.factor.levels) > length(compiled.params@by.pch)) stop("to many categories in by.colors! Add more pch values to compiled.params@by.pch!");
    points.pch <- sapply(1:lanebam.ct, FUN=function(i){
      if(! is.highlighted[i]){ return(compiled.params@highlight.points.pch[2]);
      } else {return( legend.params$points.pch[ legend.params$name == color.highlighted.by[i] ] );}
    });
    #message(paste("compiled.params@highlight.points.pch: ",compiled.params@highlight.points.pch,collapse=","));
    #message(paste("compiled.params@highlight.points.pch: ",compiled.params@highlight.points.pch,collapse=","));
    #message(paste("points.pch: ",points.pch,collapse=","));
  } else {
    points.pch <- ifelse(is.highlighted, compiled.params@highlight.points.pch[1], compiled.params@highlight.points.pch[2])
  }


  
  if(merge.offset.outgroup){
    offset.by.vals <- sort(unique(offset.by[is.highlighted]));
    if( all(offset.by %in% offset.by.vals)){
      offset.ct <- length(offset.by.vals) ;
    } else {
      offset.ct <- length(offset.by.vals) + 1;
    }
  } else {
    offset.by.vals <- sort(unique(offset.by));
    offset.ct <- length(offset.by.vals) ;
  }
  if(offset.ct %% 2 == 0){
    offsets <- (((1:offset.ct) / (offset.ct)) );
    offsets <- offsets - ( offsets[ offset.ct / 2  + 1] + offsets[ offset.ct / 2 ]  ) / 2
  } else {
    offsets <- (((1:offset.ct) / (offset.ct)) - 0.5  );
    offsets <- offsets - offsets[ (offset.ct + 1) / 2 ];
  }
  
  lanebam.offset.indices <- sapply(1:lanebam.ct, FUN=function(i){
    if(any(offset.by.vals == offset.by[i])){
      which(offset.by.vals == offset.by[i]);
    } else {
      offset.ct;
    }
  });
  lanebam.offsets <- offsets[lanebam.offset.indices];

  lanebam.params <- data.frame( plot.priority = ifelse(is.highlighted, 2,1),
                                           unique.ID = res@decoder$unique.ID,
                                           lines.col = lines.col,
                                           points.col = points.col,
                                           points.pch = points.pch,
                                           lines.lty = ifelse(is.highlighted, compiled.params@highlight.lines.lty[1], compiled.params@highlight.lines.lty[2]),
                                           lines.lwd = ifelse(is.highlighted, compiled.params@highlight.lines.lwd[1], compiled.params@highlight.lines.lwd[2]),
                                           lines.alpha = ifelse(is.highlighted, compiled.params@highlight.lines.alpha[1], compiled.params@highlight.lines.alpha[2]),
                                           points.alpha = ifelse(is.highlighted, compiled.params@highlight.points.alpha[1], compiled.params@highlight.points.alpha[2]),
                                           horiz.offsets = lanebam.offsets,
                                           vert.offsets = lanebam.offsets,
                                           stringsAsFactors=F
  );
  
  plotParams <- generate.plotter(res = res, plot.type = plot.type, title.highlight.name = title.highlight.name, legend.params = legend.params, showLegend = showLegend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, lanebam.params = lanebam.params)
  return(plotParams);
}

##################################################################################################################
##################################################################################################################
##################################################################################################################

#curr.highlight
#highlight.by
#color.highlighted.by
build.plotter.color <- function(res, compiled.params, color.by, plot.type, offset.by = color.by, merge.offset.outgroup = TRUE, pch.matched.with.color = TRUE, color.by.title.name = ""){
  if(! is.null(check.isValid(res))){
    #THROW AN ERROR!
    stop("Error parsing results: ", check.isValid(res));
  }
  
  #if(length(highlight.by) != length(color.by)) stop("error!");
  #if(length(highlight.by) != length(res@decoder$unique.ID)) stop("error!");

  
  #plotParams <- new("QoRT_Plotter");
  #plotParams@res <- res;
  #plotParams@plot.type <- plot.type;
  #title.highlight.name <- curr.highlight;
  showLegend <- compiled.params@showLegend;
  nvc.colors <- compiled.params@nvc.colors;
  nvc.colors.light <- compiled.params@nvc.colors.light;
  
  is.highlighted <- rep(TRUE,length(color.by))
  lanebam.ct <- length(color.by);
  hl.by.factor <- color.by[is.highlighted ];
  hl.by.factor.levels <- sort(unique(hl.by.factor));
  
  if(length(hl.by.factor.levels) > length(compiled.params@by.colors)) {
    message("WARNING WARNING WARNING: Too many categories to color! (Max = ",length(compiled.params@by.colors),", curr = ",length(hl.by.factor.levels),") Add more colors? Falling back to a rainbow palette");
    compiled.params@by.colors <- rainbow(length(hl.by.factor.levels));
  } 
  
  if(pch.matched.with.color){
    if(length(hl.by.factor.levels) > length(compiled.params@by.pch)){
      message("WARNING WARNING WARNING: Too many categories to mark distinct with pch! Ran out of R symbols and ASCII characters! (Max = ",length(compiled.params@by.pch),", curr = ",length(hl.by.factor.levels),") Add more characters? Falling back to a repeating pattern.");
      compiled.params@by.pch <- rep(compiled.params@by.pch, ceiling( length(hl.by.factor.levels) / length(compiled.params@by.pch)));
    }
    legend.points.pch <- compiled.params@by.pch[1:length(hl.by.factor.levels)];
  } else {
    legend.points.pch <- rep(compiled.params@highlight.points.pch[1], length(hl.by.factor.levels));
  }

  legend.params <- data.frame(name = hl.by.factor.levels,
                                         lines.col = compiled.params@by.colors[1:length(hl.by.factor.levels)],
                                         lines.lty = rep(compiled.params@highlight.lines.lty[1],length(hl.by.factor.levels)),
                                         points.pch = legend.points.pch,
                                         points.col = compiled.params@by.colors[1:length(hl.by.factor.levels)],
                                         stringsAsFactors=F
  );
  
  lines.col <- sapply(1:lanebam.ct, FUN=function(i){
    if(! is.highlighted[i]){ return(compiled.params@highlight.color[2]);
    } else {return( legend.params$lines.col[ legend.params$name == color.by[i] ] );}
  });
  points.col <- sapply(1:lanebam.ct, FUN=function(i){
    if(! is.highlighted[i]){ return(compiled.params@highlight.color[2]);
    } else {return( legend.params$points.col[ legend.params$name == color.by[i] ] );}
  });

  if(pch.matched.with.color){
    #message("!!!");
    if(length(hl.by.factor.levels) > length(compiled.params@by.pch)) stop("to many categories in by.colors! Add more pch values to compiled.params@by.pch!");
    points.pch <- sapply(1:lanebam.ct, FUN=function(i){
      if(! is.highlighted[i]){ return(compiled.params@highlight.points.pch[2]);
      } else {return( legend.params$points.pch[ legend.params$name == color.by[i] ] );}
    });
    #message(paste("compiled.params@highlight.points.pch: ",compiled.params@highlight.points.pch,collapse=","));
    #message(paste("compiled.params@highlight.points.pch: ",compiled.params@highlight.points.pch,collapse=","));
    #message(paste("points.pch: ",points.pch,collapse=","));
  } else {
    points.pch <- ifelse(is.highlighted, compiled.params@highlight.points.pch[1], compiled.params@highlight.points.pch[2])
  }

  if(merge.offset.outgroup){
    offset.by.vals <- sort(unique(offset.by[is.highlighted]));
    if( all(offset.by %in% offset.by.vals)){
      offset.ct <- length(offset.by.vals) ;
    } else {
      offset.ct <- length(offset.by.vals) + 1;
    }
  } else {
    offset.by.vals <- sort(unique(offset.by));
    offset.ct <- length(offset.by.vals) ;
  }
  if(offset.ct %% 2 == 0){
    offsets <- (((1:offset.ct) / (offset.ct)) );
    offsets <- offsets - ( offsets[ offset.ct / 2  + 1] + offsets[ offset.ct / 2 ]  ) / 2
  } else {
    offsets <- (((1:offset.ct) / (offset.ct)) - 0.5  );
    offsets <- offsets - offsets[ (offset.ct + 1) / 2 ];
  }
  
  lanebam.offset.indices <- sapply(1:lanebam.ct, FUN=function(i){
    if(any(offset.by.vals == offset.by[i])){
      which(offset.by.vals == offset.by[i]);
    } else {
      offset.ct;
    }
  });
  lanebam.offsets <- offsets[lanebam.offset.indices];

  lanebam.params <- data.frame( plot.priority = ifelse(is.highlighted, 2,1),
                                           unique.ID = res@decoder$unique.ID,
                                           lines.col = lines.col,
                                           points.col = points.col,
                                           points.pch = points.pch,
                                           lines.lty = ifelse(is.highlighted, compiled.params@highlight.lines.lty[1], compiled.params@highlight.lines.lty[2]),
                                           lines.lwd = ifelse(is.highlighted, compiled.params@highlight.lines.lwd[1], compiled.params@highlight.lines.lwd[2]),
                                           lines.alpha = ifelse(is.highlighted, compiled.params@highlight.lines.alpha[1], compiled.params@highlight.lines.alpha[2]),
                                           points.alpha = ifelse(is.highlighted, compiled.params@highlight.points.alpha[1], compiled.params@highlight.points.alpha[2]),
                                           horiz.offsets = lanebam.offsets,
                                           vert.offsets = lanebam.offsets,
                                           stringsAsFactors=F
  );
  
  plotParams <- generate.plotter(res = res, plot.type = plot.type, title.highlight.name = color.by.title.name, legend.params = legend.params, showLegend = showLegend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, lanebam.params = lanebam.params)
  return(plotParams);
}

##################################################################################################################
##################################################################################################################
##################################################################################################################

build.plotter.highlight <- function(curr.highlight, res, compiled.params, highlight.by, plot.type, offset.by, merge.offset.outgroup = TRUE){
  if(! is.null(check.isValid(res))){
    #THROW AN ERROR!
    stop("Error parsing results: ", check.isValid(res));
  }
  if(length(highlight.by) != length(res@decoder$unique.ID)) stop("error!");
   if(! (curr.highlight %in% highlight.by)){
      stop("FATAL ERROR! ","Cannot find highlighted element ", curr.highlight," in the highlight list!");
   }
  
  #hl.color <- compiled.params@highlight.color;

  #plotParams <- new("QoRT_Plotter");
  #plotParams@res <- res;
  #plotParams@plot.type <- plot.type;
  title.highlight.name <- curr.highlight;
  showLegend <- compiled.params@showLegend;
  nvc.colors <- compiled.params@nvc.colors;
  nvc.colors.light <- compiled.params@nvc.colors.light;
  
  is.highlighted <- highlight.by == curr.highlight;
  lanebam.ct <- length(highlight.by);
  
  legend.params <- data.frame(name = c(curr.highlight,"Other Samples"),
                                         lines.col = compiled.params@highlight.color,
                                         lines.lty = compiled.params@highlight.lines.lty,
                                         points.pch = compiled.params@highlight.points.pch,
                                         points.col = compiled.params@highlight.color,
                                         stringsAsFactors=F
  );
  if(merge.offset.outgroup){
    offset.by.vals <- sort(unique(offset.by[is.highlighted]));
    if( all(offset.by %in% offset.by.vals)){
      offset.ct <- length(offset.by.vals) ;
    } else {
      offset.ct <- length(offset.by.vals) + 1;
    }
  } else {
    offset.by.vals <- sort(unique(offset.by));
    offset.ct <- length(offset.by.vals) ;
  }
  if(offset.ct %% 2 == 0){
    offsets <- (((1:offset.ct) / (offset.ct)) );
    offsets <- offsets - ( offsets[ offset.ct / 2  + 1] + offsets[ offset.ct / 2 ]  ) / 2
  } else {
    offsets <- (((1:offset.ct) / (offset.ct)) - 0.5  );
    offsets <- offsets - offsets[ (offset.ct + 1) / 2 ];
  }
  
  lanebam.offset.indices <- sapply(1:lanebam.ct, FUN=function(i){
    if(any(offset.by.vals == offset.by[i])){
      which(offset.by.vals == offset.by[i]);
    } else {
      offset.ct;
    }
  });
  lanebam.offsets <- offsets[lanebam.offset.indices];

  lanebam.params <- data.frame( plot.priority = ifelse(is.highlighted, 2,1),
                                           unique.ID = res@decoder$unique.ID,
                                           lines.col = ifelse(is.highlighted, compiled.params@highlight.color[1], compiled.params@highlight.color[2]),
                                           points.col = ifelse(is.highlighted, compiled.params@highlight.color[1], compiled.params@highlight.color[2]),
                                           points.pch = ifelse(is.highlighted, compiled.params@highlight.points.pch[1], compiled.params@highlight.points.pch[2]),
                                           lines.lty = ifelse(is.highlighted, compiled.params@highlight.lines.lty[1], compiled.params@highlight.lines.lty[2]),
                                           lines.lwd = ifelse(is.highlighted, compiled.params@highlight.lines.lwd[1], compiled.params@highlight.lines.lwd[2]),
                                           lines.alpha = ifelse(is.highlighted, compiled.params@highlight.lines.alpha[1], compiled.params@highlight.lines.alpha[2]),
                                           points.alpha = ifelse(is.highlighted, compiled.params@highlight.points.alpha[1], compiled.params@highlight.points.alpha[2]),
                                           horiz.offsets = lanebam.offsets,
                                           vert.offsets = lanebam.offsets,
                                           stringsAsFactors=F
  );
  
  plotParams <- generate.plotter(res = res, plot.type = plot.type, title.highlight.name = title.highlight.name, legend.params = legend.params, showLegend = showLegend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, lanebam.params = lanebam.params)
  return(plotParams);
}


##################################################################################################################
##################################################################################################################
##################################################################################################################
