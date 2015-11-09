
f.na <- function(x){
  ifelse(is.na(x),FALSE,x);
}

wordwrap.string <- function(s, width){
  
  if(nchar(s) > width){
    curr <- s;
    out <- c();
    for(i in 1:floor(nchar(s) / width)){
      temp <- substr(curr,1,width);
      remainder <- substr(curr,width+1,nchar(curr));
      out <- c(out,temp);
      curr <- remainder;
    }
    out <- c(out,curr);
    return(paste0(out,collapse='\n'));
  } else {
    return(s);
  }
}



draw.logyaxis.stdScalePlot <- function(ylim, ylim.truncate, lwd, lwd.mini, tcl = -0.5, side = 2, label.style = c("base10","raw","none"),...){
  label.style <- match.arg(label.style);
  if(missing(ylim)){
    ylim <- c(par("usr")[3], par("usr")[4]);
  }
  if(missing(lwd)){
    lwd <- par("lwd");
  }
  if(missing(lwd.mini)){
    lwd.mini <- lwd / 2;
  }
  decades_log10 <- (floor((ylim[1]))-1):(ceiling((ylim[2]))+1);
  decades_unscaled <- 10 ^ decades_log10;
  
  ticks_unscaled <- c();
  for(i in 1:length(decades_unscaled)){
    ticks_unscaled <- c(ticks_unscaled,
                        decades_unscaled[i]*2,
                        decades_unscaled[i]*3,
                        decades_unscaled[i]*4,
                        decades_unscaled[i]*5,
                        decades_unscaled[i]*6,
                        decades_unscaled[i]*7,
                        decades_unscaled[i]*8,
                        decades_unscaled[i]*9
                        );
  }
  ticks_log10 <- log10(ticks_unscaled);
  
  if(label.style == "base10"){
      decade.labels <- c();
      if(any(decades_log10 == 0)){
        decade.labels <- c(expression(1));
      }
      if(any(decades_log10 == 1)){
        decade.labels <- c(decade.labels,expression(10));
      }
      if(decades_log10[1] < 0){
        decade.labels <- c(sapply(decades_log10[1]:(-1), function(X){
                                    substitute(10 ^ x, list(x = X));
                                 }), 
                          decade.labels);
      }
      if(max(decades_log10) > 1){
        decade.labels <- c(decade.labels,
                          sapply(2:max(decades_log10), function(X){
                                    substitute(10 ^ x, list(x = X));
                                 })
                          );
      }
  } else if(label.style == "raw"){
    decade.labels <- decades_unscaled;
  } else if(label.style == "none"){
    decade.labels <- rep("",length(decades_log10));
  }
  
  if(! missing(ylim.truncate)){
    keep_decades <- decades_log10 >= ylim.truncate[1] & decades_log10 <= ylim.truncate[2];
    decades_unscaled <- decades_unscaled[keep_decades];
    decade.labels <- decade.labels[keep_decades];
    decades_log10 <- decades_log10[keep_decades];
    ticks_unscaled <- ticks_unscaled[ticks_log10 >= ylim.truncate[1] & ticks_log10 <= ylim.truncate[2]];
    ticks_log10 <-  ticks_log10[ticks_log10 >= ylim.truncate[1] & ticks_log10 <= ylim.truncate[2]];
  }
  
  axis(side, at = decades_log10, labels=decade.labels, lwd = -1, lwd.ticks = lwd, tcl = tcl, las = 1, ...);
  axis(side, at = ticks_log10, labels=FALSE, lwd = -1, lwd.ticks = lwd.mini, tcl = tcl / 2, ...);
}

draw.yaxis.logScalePlot <- function(ylim, ylim.truncate, lwd, lwd.mini, tcl = -0.5, side = 2, label.style = c("base10","raw","none"),...){
  label.style <- match.arg(label.style);
  if(missing(ylim)){
    ylim <- c(10^par("usr")[3], 10^par("usr")[4]);
  }
  if(missing(lwd)){
    lwd <- par("lwd");
  }
  if(missing(lwd.mini)){
    lwd.mini <- lwd / 2;
  }
  decades_log10 <- (floor(log10(ylim[1]))-1):(ceiling(log10(ylim[2]))+1);
  decades_unscaled <- 10 ^ decades_log10;
  
  ticks_unscaled <- c();
  for(i in 1:length(decades_unscaled)){
    ticks_unscaled <- c(ticks_unscaled,
                        decades_unscaled[i]*2,
                        decades_unscaled[i]*3,
                        decades_unscaled[i]*4,
                        decades_unscaled[i]*5,
                        decades_unscaled[i]*6,
                        decades_unscaled[i]*7,
                        decades_unscaled[i]*8,
                        decades_unscaled[i]*9
                        );
  }
  
  if(label.style == "base10"){
      decade.labels <- c(expression(1), expression(10));
      if(decades_log10[1] < 0){
        decade.labels <- c(sapply(decades_log10[1]:(-1), function(X){
                                    substitute(10 ^ x, list(x = X));
                                 }), 
                          decade.labels);
      }
      if(max(decades_log10) > 1){
        decade.labels <- c(decade.labels,
                          sapply(2:max(decades_log10), function(X){
                                    substitute(10 ^ x, list(x = X));
                                 })
                          );
      }
  } else if(label.style == "raw"){
    decade.labels <- decades_unscaled;
  } else if(label.style == "none"){
    decade.labels <- rep("",length(decades_log10));
  }
  
  if(! missing(ylim.truncate)){
    keep_decades <- decades_unscaled >= ylim.truncate[1] & decades_unscaled <= ylim.truncate[2];
    decades_unscaled <- decades_unscaled[keep_decades];
    decade.labels <- decade.labels[keep_decades];
    ticks_unscaled <- ticks_unscaled[ticks_unscaled >= ylim.truncate[1] & ticks_unscaled <= ylim.truncate[2]];
  }
  
  axis(side, at = decades_unscaled, labels=decade.labels, lwd = -1, lwd.ticks = lwd, tcl = tcl, las = 1, ...);
  axis(side, at = ticks_unscaled, labels=FALSE, lwd = -1, lwd.ticks = lwd.mini, tcl = tcl / 2, ...);
}

draw.decade.lines <- function(side = 2, lty = 3, col = "gray", ...){
  if(side == 2){
    if(par("ylog")){
      ylim <- c(10^par("usr")[3], 10^par("usr")[4]);
      abline(h= 10 ^ ((floor(log(ylim[1]))-1):(ceiling(log(ylim[2]))+1)), lty = lty, col = col, ...);
    } else {
      ylim <- c(par("usr")[3], par("usr")[4]);
      abline(h=(floor(log(ylim[1]))-1):(ceiling(log(ylim[2]))+1), lty = lty, col = col, ...);
    }
  }
  if(side == 1){
    if(par("xlog")){
      xlim <- c(10^par("usr")[1], 10^par("usr")[2]);
      abline(h=10 ^ ((floor(log(xlim[1]))-1):(ceiling(log(xlim[2]))+1)), lty = lty, col = col, ...);
    } else {
      xlim <- c(par("usr")[1], par("usr")[2]);
      abline(h=(floor(log(xlim[1]))-1):(ceiling(log(xlim[2]))+1), lty = lty, col = col, ...);
    }
  }
}


check.rasterize.or.warn <- function(varname = "rasterize.plotting.area"){
    if(! can.plot.raster()){
       message(paste0("cannot open raster images, because package \"png\" not found. Install package png (with command install.packages(\"png\")) or set ",varname," to FALSE!"));
       warning(paste0("cannot open raster images, because package \"png\" not found. Install package png (with command install.packages(\"png\")) or set ",varname," to FALSE!"));
       return(FALSE);
    }
    if(! can.rasterize()){
       message(paste0("cannot rasterize plots, because no viable raster device found. Install R with png support, install the R package Cairo (with command install.packages(\"Cairo\")), or set ",varname," to FALSE!"));
       warning(paste0("cannot rasterize plots, because no viable raster device found. Install R with png support, install the R package Cairo (with command install.packages(\"Cairo\")), or set ",varname," to FALSE!"));
       return(FALSE);
    }
    return(TRUE);
}
check.rasterize.or.die <- function(varname = "rasterize.plotting.area"){
    if(! can.plot.raster()){
       stop(paste0("ERROR: cannot open raster images, because package png not found. Install package png (with command install.packages(\"png\")) or set ",varname," to FALSE!"));
    }
    if(! can.rasterize()){
       stop(paste0("ERROR: cannot rasterize plots, because no viable raster device found. Install R with png support, install the R package Cairo (with command install.packages(\"Cairo\")), or set ",varname," to FALSE!"));
    }
}

can.plot.raster <- function(){
   if(! suppressWarnings(require(png))){
     #warning("Cannot rasterize plots! Package png not found. Recommend installing package png.");
     return(FALSE);
   }
   return(TRUE);
}
rasterizationDevices <- c("png","CairoPNG");
can.rasterize <- function(){
   for(d in rasterizationDevices){
     if(is.available.raster.device(d)){
        return(TRUE);
     }
   }
   #warning("Cannot rasterize plots! No working raster device found. Install R with png support or install package Cairo.");
   return(FALSE);
}
can.select.raster.device <- function(){
   for(d in rasterizationDevices){
     if(is.available.raster.device(d)){
        return(TRUE);
     }
   }
   return(FALSE);
}
autoselect.raster.device <- function(debugMode){
   if(is.available.raster.device("png")){
     if(debugMode) message(">     (rasterizingPlot) selecting png raster device.");
     rasterDev <- grDevices::png;
     outDev <- function(...){
       rasterDev(units = "px", ...);
     }
     return(outDev);
   } else if(is.available.raster.device("CairoPNG")){
     if(debugMode) message("selecting CairoPNG raster device.");
     rasterDev <- CairoPNG;
     return(rasterDev);
   } else {
     stop("Error: no viable raster (png) device found! Install package CairoPNG or reinstall R with png capabilities activated!");
   }
}
is.available.raster.device <- function(d){
   if(d == "png"){
     capabilities()[["png"]];
   } else if(d == "CairoPNG"){
     if(suppressWarnings(require("Cairo"))){
       return(TRUE);
     } else {
       return(FALSE);
     }
   } else {
     return(FALSE);
   }
}

make.mixed.raster.vector.function <- function(use.raster.device = NULL, raster.height = 1000, raster.width = 1000, debugMode = TRUE, ...){
  if(! require(png)){
    stop("FATAL ERROR: cannot make hybrid raster/vector drawings without package \"png\". Install package png or set rasterize.plotting.area to FALSE");
  }

  if(is.null(use.raster.device)){
    use.raster.device <- autoselect.raster.device(debugMode);
  }
  usetempfile <- tempfile();
  if(debugMode) message(">     (rasterizingPlot) Selecting temp file location: ",usetempfile);
  initRaster <- function(){
    if(debugMode) message(">     (rasterizingPlot) Writing temp file to",usetempfile);
    use.raster.device(filename = usetempfile, height=raster.height, width=raster.width, bg = "transparent", ...);
    if(debugMode) message(">     (rasterizingPlot) Temp file built.");
    par(mar = c(0,0,0,0));
  }
  
  closeRaster <- function(){
    dev.off();
  }
  
  printRaster <- function(){
    #if(debugMode) message(">     (rasterizingPlot) Reading Rasterized image...");
    require(png);
    imgData <- png::readPNG(usetempfile);
    #if(debugMode) message(">     (rasterizingPlot) Rasterized image read.");
    rasterData <- as.raster(imgData);
    #if(debugMode) message(">     (rasterizingPlot) Cast Rasterized image to raster class.");
    usr <- par("usr");
    xlog <- par("xlog");
    ylog <- par("ylog");
    if(xlog){
       usr[1:2] <- 10 ^ usr[1:2];
    }
    if(ylog){
       usr[3:4] <- 10 ^ usr[3:4];
    }
    #if(debugMode) message(">     (rasterizingPlot) Printing raster image...");
    rasterImage(rasterData, usr[1], usr[3], usr[2], usr[4]);
    
    if(debugMode) message(">     (rasterizingPlot) deleting temp raster image: ",usetempfile);
    unlink(usetempfile);
  }
  
  #if(debugMode) message(">     (rasterizingPlot) Built rasterization fcn.");
  
  return(list(initRaster = initRaster, closeRaster = closeRaster, printRaster = printRaster, usetempfile = usetempfile));
}



write.table.with.first.col <- function(t, file, rownamesTitle = "RowNames"){
  out <- cbind.data.frame(row.names(t), t, stringsAsFactors=F);
  names(out)[1] <- rownamesTitle;
  
  write.table(out, file = file, sep="\t", row.names=F,col.names=T, quote=F);
}


plotting.limits <- function(){
   usr <- par("usr");
   x.log <- par("xlog");
   y.log <- par("ylog");
   x.out <- c(usr[1],usr[2]);
   if(x.log) x.out <- 10 ^ x.out;
   y.out <- c(usr[3],usr[4]);
   if(y.log) y.out <- 10 ^ y.out;
   out <- c(x.out, y.out);
   return(out);
}

device.limits <- function(){
    usr <- par("usr");
    plt <- par("plt");
    x.log <- par("xlog");
    y.log <- par("ylog");

    x.plotfrac <- plt[2] - plt[1];
    x.range <- usr[2] - usr[1];
    x.adjust <- x.range / x.plotfrac;
    x.out <- c(usr[1] - (plt[1] * x.adjust), usr[2] + ((1 - plt[2]) * x.adjust));
    if(x.log) x.out <- 10 ^ x.out;

    y.plotfrac <- plt[4] - plt[3];
    y.range <- usr[4] - usr[3];
    y.adjust <- y.range / y.plotfrac;
    y.out <- c(usr[3] - (plt[3] * y.adjust), usr[4] + ((1 - plt[4]) * y.adjust));
    if(y.log) y.out <- 10 ^ y.out;

    out <- c(x.out,y.out);
    return(out);
}



nonsimple.replace <- function(n, replacement){
  ifelse(is.simple.number(n), n, replacement)
}

is.simple.number <- function(n){
   (! is.nan(n)) & (is.finite(n)) & (! is.na(n)); 
}

convertString <- function(s){
   if(suppressWarnings(is.na(as.numeric(s)))){
     return(s);
   } else {
     return(as.numeric(s));
   }
}

timestamp <- function(){
   return(Sys.time());
}

timediff <- function(ts){
   return(as.numeric(Sys.time() - ts, units = "secs"));
}

getTimeAndDiff <- function(ts = NULL){
  if(is.null(ts)){
    return(paste0("[time: ",Sys.time(),"]"));
  } else {
    nts <- Sys.time();
    elapsed <- as.numeric(nts - ts, units = "secs");
    if(elapsed < 1){
      elapsed <- as.character(round(elapsed,digits=2));
    } else {
      elapsed <- floor(elapsed);
    }
    return(paste0("[time: ",nts,"],[elapsed: ",elapsed," secs]"));
  }
}

reportTimeAndDiff <- function(ts = NULL){
  message(getTimeAndDiff(ts));
}


####################################################################################
###   Minor Utility functions:

removeLineBreaks <- function(s){
   gsub("\n"," ",s,fixed=TRUE)
}

overmerge.list <- function(list.old,list.new){
  list.out <- list.old;
  if(length(list.new) > 0){
    for(i in 1:length(list.new)){
      list.out[[names(list.new)[i]]] <- list.new[[i]];
    }
  }
  return(list.out);
}

make.multicolored.phrase <- function(phrase, col){
  #title(expression("hair Color" * phantom(" and Eye color")),col.main="red")
  #title(expression(phantom("hair Color and ") * "Eye color"),col.main="blue")
  #title(expression(phantom("hair Color ") * "and" * phantom("Eye color"),col.main="black"))
}

unCamelCase.words <- function(p){
  #INCOMPLETE!
}

cap.words <- function(p){
   return(paste0(sapply(strsplit(p," "),cap.first), collapse=" "));
}
cap.first <- function(w){
   return( paste0(toupper(substring(w,1,1)), substring(w,2)) );
}

limited.pretty <- function(axisLim){
   pretty.axis <- pretty(axisLim);
   pretty.axis[1] <- max( axisLim[1], pretty.axis[1]);
   pretty.axis[length(pretty.axis)] <- min(pretty.axis[length(pretty.axis)], axisLim[2])
   return(pretty.axis);
}


color2transparent <- function(someColor,alpha){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}





####################################################################################
###   Text fitting functions:


fit.title <- function(title.text){
  plt <- par("plt");
  usr <- par("usr");
  
  x.dist <- abs(usr[2] - usr[1]);
  x.dist.frac <- min(plt[1], 1-plt[2]);
  x.frac <- (1 - x.dist.frac - x.dist.frac);
  extra.frac <- 1 + (1 - x.frac);
  out.dist <- x.dist * extra.frac;  
  
  default.width <- strwidth(title.text, cex = par("cex.main"))
  if(default.width > out.dist){
    return( fit.character.vector.helper(title.text, curr.cex = par("cex.main"), min.width = out.dist * 0.8, max.width = out.dist * 0.975, max.width.per.char = Inf) );
  } else {
    return(par("cex.main"));
  }
}

fit.vert <- function(title.text, default.cex = par("cex.ylab")){
  plt <- par("plt");
  usr <- par("usr");
  
  x.dist <- abs(usr[4] - usr[3]);
  x.frac <- abs(plt[4] - plt[3]);
  extra.frac <- 1 + (1 - x.frac);
  out.dist <- x.dist * extra.frac;  
  
  default.width <- strwidth(title.text, cex = default.cex)
  if(default.width > out.dist){
    return( fit.character.vector.helper(title.text, curr.cex = default.cex, min.width = out.dist * 0.8, max.width = out.dist, max.width.per.char = Inf) );
  } else {
    return(default.cex);
  }
}

fit.character.vector <- function(strs, min.width = 0.6, max.width = 0.95, max.width.per.char = 0.15){
   curr.cex <- 1;
   return(fit.character.vector.helper(strs, curr.cex = curr.cex, min.width = min.width, max.width = max.width, max.width.per.char = max.width.per.char));
   
}

#fit.character.vector.helper <- function(strs, curr.cex, min.width, max.width, max.width.per.char){
#   curr.width <- max(strwidth(strs, cex = curr.cex));
#   curr.nchar <- max(nchar(strs));
#   
#   curr.width.per.char <- curr.width / curr.nchar;
#   
#   if( curr.width > max.width | curr.width.per.char > max.width.per.char){
#     new.cex <- curr.cex * (((max.width-min.width)/2)+min.width) / (max.width);
#     fit.character.vector.helper(strs, new.cex, min.width, max.width, max.width.per.char);
#   } else if(curr.width < min.width){
#     new.cex <- curr.cex * (((max.width-min.width)/2)+min.width) / (min.width);
#     fit.character.vector.helper(strs, new.cex, min.width, max.width, max.width.per.char);
#   } else {
#     return(curr.cex);
#   }
#}

fit.character.vector.helper <- function(strs, curr.cex, min.width, max.width, max.width.per.char){
   curr.width <- max(strwidth(strs, cex = curr.cex));
   strs.nchar <- max(nchar(strs));
   
   curr.width.per.char <- curr.width / strs.nchar;
   
   desired.width <- ((max.width - min.width) / 2) + min.width;
   new.cex <- curr.cex * (desired.width / curr.width);
   
   new.width <- max(strwidth(strs, cex = new.cex));
   new.width.per.char <- new.width / strs.nchar;
   
   if(new.width.per.char > max.width.per.char){
      desired.width.per.char <- max.width.per.char;
      new.cex.perchar <- curr.cex * (desired.width.per.char / curr.width.per.char);
      return(new.cex.perchar);
   } else {
      return(new.cex);
   }
}
  
  
#   if( curr.width > max.width | curr.width.per.char > max.width.per.char){
#     new.cex <- curr.cex * (((max.width-min.width)/2)+min.width) / (max.width);
#     fit.character.vector.helper(strs, new.cex, min.width, max.width, max.width.per.char);
#   } else if(curr.width < min.width){
#     new.cex <- curr.cex * (((max.width-min.width)/2)+min.width) / (min.width);
#     fit.character.vector.helper(strs, new.cex, min.width, max.width, max.width.per.char);
#   } else {
#     return(curr.cex);
#   }
#}


#############################################

#Expansion of axis.break from the plotrix package
# adds the ability to modify additional graphical parameters, and fixes a bug where
# this function globally overrides certain graphical parameters via par().

qorts.axis.break <- function (axis = 1, breakpos = NULL, pos = NA, bgcol = "white", 
    breakcol = "black", style = "slash", cex = 1, xw = NULL, yw = NULL, fill = FALSE,...) 
{
    brw <- strwidth("W",cex = cex, units="figure");
    brh <- strheight("W",cex = cex, units="figure");
    figxy <- par("usr")
    xaxl <- par("xlog")
    yaxl <- par("ylog")
    if(is.null(xw)) xw <- (figxy[2] - figxy[1]) * brw
    if(is.null(yw)) yw <- (figxy[4] - figxy[3]) * brh
    if (!is.na(pos)) 
        figxy <- rep(pos, 4)
    if (is.null(breakpos)) 
        breakpos <- ifelse(axis%%2, figxy[1] + xw * 2, figxy[3] + 
            yw * 2)
    if (xaxl && (axis == 1 || axis == 3)) 
        breakpos <- log10(breakpos)
    if (yaxl && (axis == 2 || axis == 4)) 
        breakpos <- log10(breakpos)
    switch(axis, 
        br <- c(breakpos - xw/2, figxy[3] - yw/2, breakpos + xw/2, figxy[3] + yw/2), 
        br <- c(figxy[1] - xw/2, breakpos - yw/2, figxy[1] + xw/2, breakpos + yw/2), 
        br <- c(breakpos - xw/2, figxy[4] - yw/2, breakpos + xw/2, figxy[4] + yw/2), 
        br <- c(figxy[2] - xw/2, breakpos - yw/2, figxy[2] + xw/2, breakpos + yw/2), 
        stop("Improper axis specification."))
    #old.xpd <- par("xpd")
    #par(xpd = TRUE)
    if (xaxl) 
        br[c(1, 3)] <- 10^br[c(1, 3)]
    if (yaxl) 
        br[c(2, 4)] <- 10^br[c(2, 4)]
    if (style == "gap") {
        if (xaxl) {
            figxy[1] <- 10^figxy[1]
            figxy[2] <- 10^figxy[2]
        }
        if (yaxl) {
            figxy[3] <- 10^figxy[3]
            figxy[4] <- 10^figxy[4]
        }
        if (axis == 1 || axis == 3) {
            rect(breakpos, figxy[3], breakpos + xw, figxy[4], 
                col = bgcol, border = bgcol)
            xbegin <- c(breakpos, breakpos + xw)
            ybegin <- c(figxy[3], figxy[3])
            xend <- c(breakpos, breakpos + xw)
            yend <- c(figxy[4], figxy[4])
            if (xaxl) {
                xbegin <- 10^xbegin
                xend <- 10^xend
            }
        }
        else {
            rect(figxy[1], breakpos, figxy[2], breakpos + yw, 
                col = bgcol, border = bgcol)
            xbegin <- c(figxy[1], figxy[1])
            ybegin <- c(breakpos, breakpos + yw)
            xend <- c(figxy[2], figxy[2])
            yend <- c(breakpos, breakpos + yw)
            if (xaxl) {
                xbegin <- 10^xbegin
                xend <- 10^xend
            }
        }
        #par(xpd = TRUE)
    }
    else {
        if(fill) rect(br[1], br[2], br[3], br[4], col = bgcol, border = bgcol)
        if (style == "slash") {
            if (axis == 1 || axis == 3) {
                xbegin <- c(breakpos - xw, breakpos)
                xend <- c(breakpos, breakpos + xw)
                ybegin <- c(br[2], br[2])
                yend <- c(br[4], br[4])
                if (xaxl) {
                  xbegin <- 10^xbegin
                  xend <- 10^xend
                }
            }
            else {
                xbegin <- c(br[1], br[1])
                xend <- c(br[3], br[3])
                ybegin <- c(breakpos - yw, breakpos)
                yend <- c(breakpos, breakpos + yw)
                if (yaxl) {
                  ybegin <- 10^ybegin
                  yend <- 10^yend
                }
            }
        }
        else {
            if (axis == 1 || axis == 3) {
                xbegin <- c(breakpos - xw/2, breakpos - xw/4, 
                  breakpos + xw/4)
                xend <- c(breakpos - xw/4, breakpos + xw/4, breakpos + 
                  xw/2)
                ybegin <- c(ifelse(yaxl, 10^figxy[3 + (axis == 
                  3)], figxy[3 + (axis == 3)]), br[4], br[2])
                yend <- c(br[4], br[2], ifelse(yaxl, 10^figxy[3 + 
                  (axis == 3)], figxy[3 + (axis == 3)]))
                if (xaxl) {
                  xbegin <- 10^xbegin
                  xend <- 10^xend
                }
            }
            else {
                xbegin <- c(ifelse(xaxl, 10^figxy[1 + (axis == 
                  4)], figxy[1 + (axis == 4)]), br[1], br[3])
                xend <- c(br[1], br[3], ifelse(xaxl, 10^figxy[1 + 
                  (axis == 4)], figxy[1 + (axis == 4)]))
                ybegin <- c(breakpos - yw/2, breakpos - yw/4, 
                  breakpos + yw/4)
                yend <- c(breakpos - yw/4, breakpos + yw/4, breakpos + 
                  yw/2)
                if (yaxl) {
                  ybegin <- 10^ybegin
                  yend <- 10^yend
                }
            }
        }
    }
    segments(xbegin, ybegin, xend, yend, col = breakcol, lty = 1, xpd = NA, ...)
    
    #par(xpd = FALSE)
}
