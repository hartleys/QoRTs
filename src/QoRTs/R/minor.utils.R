
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
      elapsed <- as.character(elapsed);
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

blank.plot <- function(label){
  plot(0,0,col="transparent", xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="");
  text(0.5,0.5,labels=paste0(label,collapse="\n"),cex=3);
  box();
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