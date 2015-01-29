#Make sure all requires base packages are installed.
#  This really shouldn't be necessary, but you'd be surprised.
require(stats);
require(graphics);
require(grDevices);
require(utils);
require(datasets);
require(methods);
require(base);


setClass("QoRTs_QC_Results",representation(
                                     lanebam.list="character",
                                     sample.list="character",
                                     lane.list="character",
                                     group.list="character",
                                     decoder="data.frame", #decoder has columns: unique.ID	sample.ID	lane.ID	group.ID	cycle.CT	and then any number of user-defined columns (which are ignored internally)
                                     qc.data="list", #List of Lists. Each element corresponds to one qc test, and is composed of a list, one element for each lanebam.
                                     calc.data="list", #List of Lists. Same as above, except it holds data calculated within R rather than raw scala output. Each element corresponds to one qc test, and is composed of a list, one element for each lanebam.
                                     singleEnd="logical" #Logical. whether any of the samples are single-end. Note that paired-end samples can be plotted in single-end mode, but NOT vice versa.
                                     ));


setMethod("show","QoRTs_QC_Results",
   function(object){
      cat("QoRTs_QC_Results object:\n");
      cat("lanebam.list:            [",paste0(object@lanebam.list,collapse=","),"]\n");
      cat("singleEnd:               [",paste0(object@singleEnd),"]\n");
      cat("sample.list:             [",paste0(object@sample.list,collapse=","),"]\n");
      cat("lane.list:               [",paste0(object@lane.list,collapse=","),"]\n");
      cat("files included:          [", paste0(names(object@qc.data), collapse=","),"]\n" );
      cat("calculated data included:[", paste0(names(object@calc.data), collapse=","),"]\n" );
      #Add more printing stuff here?
   }
);

check.isValid <- function(object){
      #Test the validity of the results set!
      return(NULL);
}

##INCOMPLETE!!!
#get.results.subset <- function(res,
#                   keep.lanebams = NULL, 
#                   keep.samples = NULL, 
#                   keep.lanes = NULL,
#                   keep.groups = NULL,
#                   drop.lanebams = NULL, 
#		   drop.samples = NULL, 
#		   drop.lanes = NULL,
#                   drop.groups = NULL
#                   ){
#  decoder <- res@decoder;
#  final.keep <- rep(TRUE,length(decoder$unique.ID));
#  if(! is.null(keep.lanebams)){
#    final.keep <- final.keep & decoder$unique.ID %in% as.character(keep.lanebams);
#  }
#  if(! is.null(keep.samples)){
#    final.keep <- final.keep & decoder$sample.ID %in% as.character(keep.samples);
#  }
#  if(! is.null(keep.lanes)){
#    final.keep <- final.lanes & decoder$lane.ID %in% as.character(keep.lanes);
#  }
#  if(! is.null(keep.groups)){
#    final.keep <- final.keep & decoder$group.ID %in% as.character(keep.groups);
#  }
#
#
#  if(! is.null(drop.lanebams)){
#    final.keep <- final.keep &  (! (decoder$unique.ID %in% as.character(drop.lanebams)));
#  }
#  if(! is.null(drop.samples)){
#    final.keep <- final.keep &  (! (decoder$sample.ID %in% as.character(drop.samples)));
#  }
#  if(! is.null(drop.lanes)){
#    final.keep <- final.lanes & (! (decoder$lane.ID   %in% as.character(drop.lanes)));
#  }
#  if(! is.null(drop.groups)){
#    final.keep <- final.keep &  (! (decoder$group.ID  %in% as.character(drop.groups)));
#  }
#  
#  #NOW WHAT?
#  #FINISH ME?
#}



#############################################################
# Plotter Parameter Holder:
#############################################################

QoRTs_Plotter <- setRefClass("QoRTs_Plotter",fields = list(
                       res = "QoRTs_QC_Results",
                       plot.type = "character",
                       title.highlight.name = "character",
                       legend.params = "data.frame",   #has columns: name	lines.col	lines.lty	lines.lwd	points.pch	points.col
                       showLegend = "logical",
                       nvc.colors = "list",
                       nvc.colors.light = "list",
                       lanebam.params = "data.frame" #has columns:	plot.priority	unique.ID	lines.col	points.col	points.pch	lines.lty	lines.lwd	lines.alpha	lines.lwd	points.alpha	horiz.offsets	vert.offsets
                       ));

#setMethod("show","QoRT_Plotter",
#   function(object){
#      cat("Plotter for QoRTs:\n");
#      cat("@plot.type: ", object@plot.type,"\n");
#      cat("@title.highlight.name: ", object@title.highlight.name,"\n");
#      cat("@legend.params:\n");
#      show(object@legend.params);
#      cat("@showLegend: ", object@showLegend,"\n");
#      cat("@nvc.colors:\n");
#      show(unlist(object@nvc.colors));
#      cat("@nvc.colors.light:\n");
#      show(unlist(object@nvc.colors.light));
#      cat("@lanebam.params:\n");
#      show(object@lanebam.params);
#      cat("@res:" );
#      show(object@res);
#   }
#);

generate.plotter <- function(res, plot.type, title.highlight.name, legend.params, showLegend, nvc.colors, nvc.colors.light, lanebam.params){
  QoRTs_Plotter$new(res = res, plot.type = plot.type, title.highlight.name = title.highlight.name, legend.params = legend.params, showLegend = showLegend, nvc.colors = nvc.colors, nvc.colors.light = nvc.colors.light, lanebam.params = lanebam.params);
}

#############################################################
#############################################################

setClass("QoRTs_Compiled_Plotting_Params", representation(
                       by.colors = "character",
                       by.pch = "numeric",
                       
                       showLegend  = "logical",
                       nvc.colors  = "list",
                       nvc.colors.light = "list",
                       
                       highlight.lines.lty = "numeric",
                       highlight.lines.lwd = "numeric",
                       highlight.lines.alpha = "numeric",
                       highlight.lines.color = "character",
                       highlight.points.pch = "numeric",
                       highlight.points.alpha = "numeric",
                       highlight.points.color = "character"
                       
                       ));

setMethod("show","QoRTs_Compiled_Plotting_Params",
   function(object){
      cat("Master plotting parameters object for QoRTs:\n");
   }
);




#default.plot.master.params <- new("QoRT_Plotter_Master_Params");
#default.plot.master.params@by.colors <- c("red","blue","magenta","green","orange","cyan","purple","lightgreen","crimson","violet","darkyellow","indigo");
#default.plot.master.params@by.pch <- c(0,1,2,3,4,5,6,11,7,8,9,10,12,13,14);
#default.plot.master.params@showLegend <- TRUE;
#default.plot.master.params@highlight.lines.lty <- c(1,3);
#default.plot.master.params@highlight.lines.lwd <- c(2,1);
#default.plot.master.params@highlight.lines.alpha <- c(255,200);
#default.plot.master.params@highlight.lines.color <- c("red","grey");
#default.plot.master.params@highlight.points.pch <- c(4,1);
#default.plot.master.params@highlight.points.alpha <- c(255,200);
#default.plot.master.params@highlight.points.color <- c("red","grey");

##default.plot.master.params@nvc.colors <- list(A = "darkgreen",T = "red", G = "black", C = "blue");
##default.plot.master.params@nvc.colors.light <- list(A = "palegreen",T = "maroon3", G = "grey", C = "deepskyblue");
#default.plot.master.params@nvc.colors <- list(A = "green3",T = "red", G = "orange", C = "blue");
#default.plot.master.params@nvc.colors.light <- list(A = "olivedrab1",T = "lightpink", G = "lightgoldenrod1", C = "deepskyblue");

################################################
#                Master Param Builder:
################################################

  DEFAULT.PLOT.MASTER.PARAMS_by.colors              <-  c("red","blue","green","orange","magenta","cyan","purple","lightgreen","brown4","yellow3","seagreen4","royalblue","palegreen2", "orangered","lightsalmon","green4","khaki4","darkorchid3","deepskyblue4","darkolivegreen3","cyan4");
  DEFAULT.PLOT.MASTER.PARAMS_by.pch                 <-  c(0,2,5,6,4,3,8,7,9,10,12,13,14,11, 65:90, 97:122);
  DEFAULT.PLOT.MASTER.PARAMS_showLegend             <-  FALSE;
  DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.lty    <-  c(1,3);
  DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.lwd    <-  c(2,1);
  DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.alpha  <-  c(255,200);
  DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.color  <-  c("red","grey");
  DEFAULT.PLOT.MASTER.PARAMS_highlight.points.pch   <-  c(4,1);
  DEFAULT.PLOT.MASTER.PARAMS_highlight.points.alpha <-  c(255,200);
  DEFAULT.PLOT.MASTER.PARAMS_highlight.points.color <-  c("red","grey");
  DEFAULT.PLOT.MASTER.PARAMS_nvc.colors             <-  list(A = "green3",T = "red", G = "orange", C = "blue");
  DEFAULT.PLOT.MASTER.PARAMS_nvc.colors.light       <-  list(A = "olivedrab1",T = "lightpink", G = "lightgoldenrod1", C = "deepskyblue");

################################################

#make.QoRTs.plot.params <- function( contrasting.colors = NULL,
#                                   contrasting.pch    = NULL,
#                                   show.legend      = NULL,
#                                   std.lines.lty    = NULL,
#                                   std.lines.lwd    = NULL,
#                                   std.lines.alpha  = NULL,
#                                   std.lines.color  = NULL,
#                                   std.points.pch   = NULL,
#                                   std.points.alpha = NULL,
#                                   std.points.color = NULL,
#                                   std.NVC.colors   = NULL,
#                                   alt.lines.lty    = NULL,
#                                   alt.lines.lwd    = NULL,
#                                   alt.lines.alpha  = NULL,
#                                   alt.lines.color  = NULL,
#                                   alt.points.pch   = NULL,
#                                   alt.points.alpha = NULL,
#                                   alt.points.color = NULL,
#                                   alt.NVC.colors   = NULL
#                                           ){
#   new.master.params <- list(contrasting.colors = contrasting.colors,
#                                   contrasting.pch    = contrasting.pch ,
#                                   show.legend      = show.legend,
#                                   std.lines.lty    = std.lines.lty,
#                                   std.lines.lwd    = std.lines.lwd,
#                                   std.lines.alpha  = std.lines.alpha,
#                                   std.lines.color  = std.lines.color,
#                                   std.points.pch   = std.points.pch,
#                                   std.points.alpha = std.points.alpha,
#                                   std.points.color = std.points.color,
#                                   std.NVC.colors   = std.NVC.colors,
#                                   alt.lines.lty    = alt.lines.lty,
#                                   alt.lines.lwd    = alt.lines.lwd ,
#                                   alt.lines.alpha  = alt.lines.alpha,
#                                   alt.lines.color  = alt.lines.color,
#                                   alt.points.pch   = alt.points.pch,
#                                   alt.points.alpha = alt.points.alpha,
#                                   alt.points.color = alt.points.color,
#                                   alt.NVC.colors   = alt.NVC.colors);
#   return(new.master.params);   
#}

QoRTs.default.plotting.params <- list(           contrasting.colors = DEFAULT.PLOT.MASTER.PARAMS_by.colors,
                                                 contrasting.pch    = DEFAULT.PLOT.MASTER.PARAMS_by.pch,
                                                 show.legend      = DEFAULT.PLOT.MASTER.PARAMS_showLegend,
                                                 std.lines.lty    = DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.lty[1],
                                                 std.lines.lwd    = DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.lwd[1],
                                                 std.lines.alpha  = DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.alpha[1],
                                                 std.lines.color  = DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.color[1],
                                                 std.points.pch   = DEFAULT.PLOT.MASTER.PARAMS_highlight.points.pch[1],
                                                 std.points.alpha = DEFAULT.PLOT.MASTER.PARAMS_highlight.points.alpha[1],
                                                 std.points.color = DEFAULT.PLOT.MASTER.PARAMS_highlight.points.color[1],
                                                 std.NVC.colors   = DEFAULT.PLOT.MASTER.PARAMS_nvc.colors,
                                                 alt.lines.lty    = DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.lty[2],
                                                 alt.lines.lwd    = DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.lwd[2],
                                                 alt.lines.alpha  = DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.alpha[2],
                                                 alt.lines.color  = DEFAULT.PLOT.MASTER.PARAMS_highlight.lines.color[2],
                                                 alt.points.pch   = DEFAULT.PLOT.MASTER.PARAMS_highlight.points.pch[2],
                                                 alt.points.alpha = DEFAULT.PLOT.MASTER.PARAMS_highlight.points.alpha[2],
                                                 alt.points.color = DEFAULT.PLOT.MASTER.PARAMS_highlight.points.color[2],
                                                 alt.NVC.colors   = DEFAULT.PLOT.MASTER.PARAMS_nvc.colors.light);


merge.plotting.params <- function(old.params, new.params){
   out.params <- old.params;
   if(length(new.params) > 0){
     for(i in 1:length(new.params)){
       p <- names(new.params)[i];
       v <- new.params[[i]];
       out.params[[p]] <- v;
     }
   }
   return(out.params);
}

compile.plotting.params <- function(params){
  master.params <- new("QoRTs_Compiled_Plotting_Params");
  master.params@by.colors  <- params[["contrasting.colors"]];
  master.params@by.pch     <- params[["contrasting.pch"]];
  master.params@showLegend <- params[["show.legend"]];
  master.params@highlight.lines.lty    <- c(params[["std.lines.lty"]],params[["alt.lines.lty"]]);
  master.params@highlight.lines.lwd    <- c(params[["std.lines.lwd"]],params[["alt.lines.lwd"]]);
  master.params@highlight.lines.alpha  <- c(params[["std.lines.alpha"]],params[["alt.lines.alpha"]]);
  master.params@highlight.lines.color  <- c(params[["std.lines.color"]],params[["alt.lines.color"]]);
  master.params@highlight.points.pch   <- c(params[["std.points.pch"]],params[["alt.points.pch"]]);
  master.params@highlight.points.alpha <- c(params[["std.points.alpha"]],params[["alt.points.alpha"]]);
  master.params@highlight.points.color <- c(params[["std.points.color"]],params[["alt.points.color"]]);
  master.params@nvc.colors       <- params[["std.NVC.colors"]];
  master.params@nvc.colors.light <- params[["alt.NVC.colors"]];
  return(master.params);
}

decompile.plotting.params <- function(compiled.params){
  return(list(
     contrasting.colors = compiled.params@by.colors,
     contrasting.pch    = compiled.params@by.pch,
     show.legend      = compiled.params@showLegend,
     std.lines.lty    = compiled.params@highlight.lines.lty[1],
     std.lines.lwd    = compiled.params@highlight.lines.lwd[1],
     std.lines.alpha  = compiled.params@highlight.lines.alpha[1],
     std.lines.color  = compiled.params@highlight.lines.color[1],
     std.points.pch   = compiled.params@highlight.points.pch[1],
     std.points.alpha = compiled.params@highlight.points.alpha[1],
     std.points.color = compiled.params@highlight.points.color[1],
     std.NVC.colors   = compiled.params@nvc.colors,
     alt.lines.lty    = compiled.params@highlight.lines.lty[2],
     alt.lines.lwd    = compiled.params@highlight.lines.lwd[2],
     alt.lines.alpha  = compiled.params@highlight.lines.alpha[2],
     alt.lines.color  = compiled.params@highlight.lines.color[2],
     alt.points.pch   = compiled.params@highlight.points.pch[2],
     alt.points.alpha = compiled.params@highlight.points.alpha[2],
     alt.points.color = compiled.params@highlight.points.color[2],
     alt.NVC.colors   = compiled.params@nvc.colors.light
  ));
}