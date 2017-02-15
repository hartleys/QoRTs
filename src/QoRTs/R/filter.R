


AUTO.FLAG.FILTERS <- c(
    "min_READ_PAIR_OK"
    );

png.with.layout <- function(filename= "png.png",layout = matrix(1),res=150,height=7,width=7,units="in",...){
    png(filename = filename,res=res,height=height*nrow(layout),width=width*ncol(layout),units=units,...);
    layout(layout);
}

prelim.flag.bad.samples <- function(res,
   thresholds = list(
         indel.rate = 0.012,
         indel.read.rate = 0.01,
         trunc.rate = 0.65,
         refMM = 0.25,
         overMM = 0.125,
         depth = 15
       )
){
    st <- as.data.frame(t(res@summaryTable));
    d <- res@decoder;
    
    if(FALSE){
        d <- decoder;
        i <- 1;
        keepList <- sapply(d$qc.data.dir,function(dir){
              file.exists(paste0("qcdata/",dir,"/QC.QORTS_COMPLETED_OK"))
            })
        
        st <- do.call(cbind.data.frame,lapply(d$qc.data.dir,function(uid){
              if( i %% 100 == 0) cat(".");
              i = i + 1;
              read.table(paste0("qcdata/",uid,"/QC.summary.txt"),header=TRUE,sep='\t',stringsAsFactors=FALSE)[,2];
            }));
        row.names(st) <- read.table(paste0("qcdata/",d$qc.data.dir[1],"/QC.summary.txt"),header=TRUE,sep='\t',stringsAsFactors=FALSE)[,1];
        colnames(st) <- d$qc.data.dir
        st <- as.data.frame(t(st));
    }
    
    ins.rate <- st$PAIR_CONTAINS_INS / st$READ_PAIR_OK;
    del.rate <- st$PAIR_CONTAINS_DEL / st$READ_PAIR_OK;
    indel.rate <- st$PAIR_CONTAINS_INDEL / st$READ_PAIR_OK;
    ins.rate1     <- st$READ_CONTAINS_INS_R1 / st$READ_PAIR_OK;
    del.rate1     <- st$READ_CONTAINS_DEL_R1 / st$READ_PAIR_OK;
    ins.rate2     <- st$READ_CONTAINS_INS_R2 / st$READ_PAIR_OK;
    del.rate2     <- st$READ_CONTAINS_DEL_R2 / st$READ_PAIR_OK;
    ref.mm.1 <- st$HAS_REF_BASE_SWAP_R1 / st$READ_PAIR_OK;
    ref.mm.2 <- st$HAS_REF_BASE_SWAP_R2 / st$READ_PAIR_OK;
    ref.over.mm <- st$OM_overlap_mismatch / st$OM_overlap;
    frac.off <- st$fractionOfPairsOffTarget;
    trunc.1 <- st$NumReadsTruncated_R1 / st$READ_PAIR_OK;
    trunc.2 <- st$NumReadsTruncated_R2 / st$READ_PAIR_OK;
    map.rate <- st$READ_PAIR_OK / st$PREALIGNMENT_READ_CT;
    gb.coverage.median <- st$GeneBodyCoverage_Overall_Median;
    gb.coverage.mean <- st$GeneBodyCoverage_Overall_Mean;
    depth     <- st$TargetReadDepth
    depthPair <- st$TargetReadPairDepth
    
    
    #LANEBAM FILTERS:
    layout(matrix(1:12,byrow=TRUE,ncol=4))
        plot.hist.with.norm(ins.rate,"Insertion Rate")
        plot.hist.with.norm(del.rate,"Deletion Rate")
        plot.hist.with.norm(indel.rate,"INDEL Rate")
        plot.hist.with.norm(ins.rate1,"Insertion Rate (R1)")
        plot.hist.with.norm(del.rate1,"Deletion Rate (R1)")
        plot.hist.with.norm(ins.rate2,"Insertion Rate (R2)")
        plot.hist.with.norm(del.rate2,"Deletion Rate (R2)")
        plot.hist.with.norm(trunc.1,"Truncation Rate (R1)")
        plot.hist.with.norm(trunc.2,"Truncation Rate (R2)")
        plot.hist.with.norm(ref.mm.1,"Ref Mismatch Rate (R1)")
        plot.hist.with.norm(ref.mm.2,"Ref Mismatch Rate (R2)")
        plot.hist.with.norm(ref.over.mm,"Overlap Mismatch Rate",breaks=50)

    f <- NULL;
    f <- addFlag(f=NULL,ins.rate > thresholds$indel.rate,"InsRate.B")
    f <- addFlag(f,del.rate > thresholds$indel.rate,"DelRate.B")
    f <- addFlag(f,ins.rate1 > thresholds$indel.read.rate,"InsRate.R1")
    f <- addFlag(f,del.rate1 > thresholds$indel.read.rate,"DelRate.R1")
    f <- addFlag(f,ins.rate2 > thresholds$indel.read.rate,"InsRate.R2")
    f <- addFlag(f,del.rate2 > thresholds$indel.read.rate,"DelRate.R2")
    f <- addFlag(f,trunc.1 > thresholds$trunc.rate,"trunc.R1")
    f <- addFlag(f,trunc.2 > thresholds$trunc.rate,"trunc.R2")
    f <- addFlag(f,ref.mm.1 > thresholds$refMM,"refMM.R1")
    f <- addFlag(f,ref.mm.2 > thresholds$refMM,"refMM.R2")
    f <- addFlag(f,ref.over.mm > thresholds$overMM,"overMM")
    
    
    FLAG <- f[["FLAG"]];
    FLAG.REASON <- f[["FLAG.REASON"]];
    FLAG.FIRST.REASON <- f[["FLAG.FIRST.REASON"]];
    
    table(FLAG.REASON)

    REASON.TABLE <- data.frame(
               REASON=as.character(names(table(FLAG.FIRST.REASON))),
               CT = as.numeric(table(FLAG.FIRST.REASON))
               );
    REASON.TABLE[1,1] <- "OK";
    REASON.TABLE
    
    ###############################
    #Samplewise filters:
    
    depth     <- st$TargetReadDepth
    depthPair <- st$TargetReadPairDepth
    sample.list <- unique(d$sample.ID);
    
    depth.samp <-     sapply(sample.list,function(samp){ sum( depth[d$sample.ID == samp] ) })
    depthPair.samp <- sapply(sample.list,function(samp){ sum( depthPair[d$sample.ID == samp] ) })
    lanebam.flag <- sapply(sample.list,function(samp){ 
        any(f$FLAG[d$sample.ID == samp])
    })
    lanebam.reason <- sapply(sample.list,function(samp){ 
        x <- f$FLAG.REASON[d$sample.ID == samp];
        x[x == ""] <- "OK";
        paste0(x,collapse="|") 
    })

    layout(matrix(1:4,byrow=TRUE,ncol=2))
        plot.hist.with.norm(depth.samp,"Read Depth",50,xmin = 0)
        plot.hist.with.norm(depthPair.samp,"Pair Depth",50,xmin = 0)
    
    fSamp <- NULL;
    fSamp <- addFlag(f=fSamp,depth.samp < thresholds$depth,"depth")
    fSamp <- addFlag(f=fSamp,lanebam.flag,lanebam.reason)
    
    FLAG <- fSamp[["FLAG"]];
    FLAG.REASON <- fSamp[["FLAG.REASON"]];
    FLAG.FIRST.REASON <- fSamp[["FLAG.FIRST.REASON"]];

    REASON.TABLE <- data.frame(
               REASON=as.character(names(table(FLAG.FIRST.REASON))),
               CT = as.numeric(table(FLAG.FIRST.REASON))
               );
    REASON.TABLE[1,1] <- "OK";
    REASON.TABLE
    
    (sum(REASON.TABLE$CT)-REASON.TABLE$CT[1]) / sum(REASON.TABLE$CT)
    
    ###############################
    #Lane/run filters:
    
    run.list <- unique(d$FLOWCELL);
    lane.list <- unique(d$lane.ID);
    fraction.fail <- 0.33;
    
    layout(matrix(1:3,byrow=TRUE,ncol=1))
    fail <- ins.rate > thresholds$indel.rate
    plot.lane.stat(lane.list,ins.rate,d,0,"Deletion Rate (By Lane)",fail,failCt, thresholds$indel.rate,fraction.fail);
    
    fail <- del.rate > thresholds$indel.rate
    failCt <- flag.lane.calc(lane.list,fail,d)
    plot.lane.stat(lane.list,del.rate,d,0,"Insertion Rate (By Lane)",fail,failCt, thresholds$indel.rate,fraction.fail);
    
    del.rate > thresholds$indel.rate
    ins.rate1 > thresholds$indel.read.rate
    del.rate1 > thresholds$indel.read.rate
    ins.rate2 > thresholds$indel.read.rate
    del.rate2 > thresholds$indel.read.rate
    trunc.1 > thresholds$trunc.rate
    trunc.2 > thresholds$trunc.rate
    ref.mm.1 > thresholds$refMM
    ref.mm.2 > thresholds$refMM
    ref.over.mm > thresholds$overMM
    
    
    
    
    ###############################
    

    
}

flag.run.calc <- function(run.list,failVector,d){
    #run.list <- unique(d$FLOWCELL)
    out <- sapply(seq_along(run.list),function(i){
          lane <- run.list[[i]];
          laneLogic <- d$FLOWCELL == lane;
          laneVector <- failVector[laneLogic];
          sum(laneVector) / length(laneVector);
        });
    #names(out) <- run.list;
    out;
}
flag.lane.calc <- function(lane.list,failVector,d){
    #lane.list <- unique(d$lane.ID)
    out <- sapply(seq_along(lane.list),function(i){
          lane <- lane.list[[i]]
          laneLogic <- d$lane.ID == lane;
          laneVector <- failVector[laneLogic];
          sum(laneVector) / length(laneVector);
        });
    #names(out) <- lane.list;
    out;
}




get.lane.stat <- function(run.list,lane.list,stat, d,ymin, title = "",fail,thresh,laneFractionThresh = 0.33,runFractionThresh = laneFractionThresh){
    runRate <- flag.run.calc(run.list, fail,d)
    laneRate <- flag.lane.calc(lane.list,fail,d)
    laneFail <- f.na(unlist(lapply(run.list,function(runid){
          lanes <- paste0(runid,"_",1:8);
          sapply(lanes,function(laneid){
                if(laneid %in% lane.list){
                    laneRate[[which(lane.list == laneid)]] > laneFractionThresh
                } else {
                    FALSE
                }
              })
        }))) 

    laneFailCt <- f.na(unlist(lapply(run.list,function(runid){
          lanes <- paste0(runid,"_",1:8);
          sapply(lanes,function(laneid){
                if(laneid %in% lane.list){
                    laneRate[[which(lane.list == laneid)]]
                } else {
                    0
                }
              })
        }))) 
        
    #names(laneFail) <- paste0(rep(run.list,each=8),"_",rep(1:8,length(run.list)));
    
    runFail <- f.na(unlist(lapply(run.list,function(runid){
                runRate[[which(run.list == runid)]]  > laneFractionThresh
        })))
    return(list(runFail=runFail,laneFail = laneFail, laneFailCt = laneFailCt, thresh = thresh));
}

plot.lane.stat.withSampleConnect <- function(run.list,stat, d,ymin, title = "",fail,thresh,laneFractionThresh = 0.33,runFractionThresh = laneFractionThresh,
                           label.lane=TRUE, label.runs = TRUE, run.titles = run.list, label.thresh = TRUE, 
                           connect.isolated.samples = TRUE, drop.clean.runs = TRUE,
                           cex.run.label = 0.75, cex.lane.label = 0.75, srt.run.label = 90,badThresh =NULL){
    if(drop.clean.runs){
      keep.run <- sapply(run.list,function(r){ any(fail[ d$FLOWCELL == r ]) });
      if(any(keep.run)){
          run.list <- run.list[ keep.run ];
      } else {
         #else leave it alone.
      }
    }
    lane.list <- unlist(lapply(run.list,function(runid){paste0(runid,"_",1:8)}));
    
    runRate <- flag.run.calc(run.list, fail,d)
    laneRate <- flag.lane.calc(lane.list,fail,d)
    laneFail <- f.na(unlist(lapply(run.list,function(runid){
          lanes <- paste0(runid,"_",1:8);
          sapply(lanes,function(laneid){
                if(laneid %in% lane.list){
                    laneRate[[which(lane.list == laneid)]] > laneFractionThresh
                } else {
                    FALSE
                }
              })
        })))
    #names(laneFail) <- paste0(rep(run.list,each=8),"_",rep(1:8,length(run.list)));
    
    runFail <- f.na(unlist(lapply(run.list,function(runid){
                runRate[[which(run.list == runid)]]  > laneFractionThresh
        })))
    
    lane.offsets <- (1:8 / 8) - 0.5625
    
    x <- sapply(seq_along(stat),function(i){
          which(d$FLOWCELL[i] == run.list) + lane.offsets[d$LANE[i]]
        });
    if(missing(ymin)){
        ymin <- min(stat,na.rm=TRUE);
    }
    par(xaxs = 'i');
    plot.new();
    plot.window(xlim=c(0.5,length(run.list)+0.5), 
                ylim=c(ymin,max(stat,na.rm=TRUE)));

    if(label.thresh){
      failPercentile <- round((sum(fail) / length(fail)) * 100,digits=1)
      title(main=paste0(title,"\nThreshold: ",thresh," (",failPercentile,"%)"));
    } else {
        title(main=title);
    }
    box();
    
    #+ runif(length(stat),-0.1,0.1)
    points(x,stat,col=ifelse(fail,"yellow","green"));
    
    if(!is.null(badThresh)){
      abline(h=badThresh,col="red",lty=3);
    }
    
    if(connect.isolated.samples){
        isolated <- fail & (! d$lane.ID %in% names(laneFail)[laneFail] )
        if(any(isolated)){
            points(x[isolated],stat[isolated],col="red",pch=1);
            #points(x[isolated],stat[isolated],col="red",pch=4);
            isolated.samples <- unique(d$sample.ID[isolated]);
            isolated.samples.multiRun <- isolated.samples[sapply(isolated.samples,function(samp){ sum(d$sample.ID == samp) > 1 })]
            for(samp in isolated.samples.multiRun){
                lines(x[d$sample.ID == samp],stat[d$sample.ID == samp],col="red",lty=3);
                points(x[d$sample.ID == samp],stat[d$sample.ID == samp],col="red",pch=4);
            }
        } 
    }
    
    axis(2);
    
    if(label.runs){
       text(seq_along(run.list),par("usr")[3],run.titles,srt=srt.run.label,
            cex=cex.run.label,adj=c(1.1,0.5),
            xpd=NA,col=ifelse(runFail,"red","black"));
    }
    if(label.lane){
      #Note to self: 
      #  This code SEEMS unnecessarily convoluted. But text() performs REALLY slowly
      #  when handed transparent colors. This hack works around that.
      cex.lane.nums <- 0.75;
      color.lane.good <- "gray88"
      color.lane.bad  <- "red"
      char.ht <- strheight("0",cex=cex.lane.nums)
      lane.label.x <- rep(seq_along(run.list),each=8) + rep(lane.offsets, length(run.list));
      lane.label.y <- rep(rep(seq(par("usr")[3] + char.ht * 3,par("usr")[3],length.out=4),length.out=8),length(run.list));
      lane.label.label <- rep(1:8,length(lane.list))
      color.lane.good <- "gray"
      color.lane.bad  <- "red"
      
      if(sum(! laneFail) > 0){
        text(lane.label.x[!laneFail],
         lane.label.y[!laneFail],
         lane.label.label[!laneFail],
         cex=0.75,
         adj=c(0.5,-0.25),
         xpd=NA,
         col=color.lane.good);
      }
      if(sum(laneFail) > 0){
      text(lane.label.x[laneFail],
         lane.label.y[laneFail],
         lane.label.label[laneFail],
         cex=0.75,
         adj=c(0.5,-0.25),
         xpd=NA,
         col=color.lane.bad);
      }
    }
    
    borders <- c(seq_along(run.list)-0.5);
    borders <- borders[-1];
    axis(1,at=borders,labels=FALSE);
    abline(v=borders,lty=3,col="gray");
    abline(h=thresh,lty=3,col="gray");
}


plot.lane.stat.isolated <- function(run.list, FLAGS,stat, d,ymin = NULL,ymax = NULL, title = "",fail,thresh,
                           label.lane=TRUE, label.runs = TRUE, run.titles = run.list, label.thresh = TRUE, 
                           connect.isolated.samples = TRUE, drop.clean.runs = TRUE,
                           cex.run.label = 0.75, cex.lane.label = 0.75, srt.run.label = 90,badThresh =NULL,
                           ymax.pct = NULL){
    if(drop.clean.runs){
      keep.run <- sapply(run.list,function(r){ any(fail[ d$FLOWCELL == r ]) });
      fail.samples <- unique(d$sample.ID[fail]);
      keep.run <- keep.run | sapply(run.list,function(r){ any((d$sample.ID %in% fail.samples)[ d$FLOWCELL == r ]) });
      if(any(keep.run)){
          run.list <- run.list[ keep.run ];
      } else {
         #else leave it alone.
      }
    }
    lane.list <- unlist(lapply(run.list,function(runid){ paste0(runid,"_",1:8) }));
    laneFail <- FLAGS[names(FLAGS) %in% lane.list];
    laneFailList <- names(FLAGS)[FLAGS];
    
    lane.offsets <- (1:8 / 8) - 0.5625
    
    x <- as.numeric(sapply(seq_along(stat),function(i){
          which(d$FLOWCELL[i] == run.list) + lane.offsets[d$LANE[i]]
        }));
    if(is.null(ymin)){
        ymin <- min(stat,na.rm=TRUE);
    }
    if(! is.null(ymax.pct)){
      ymax <- quantile(stat,prob=ymax.pct,na.rm=TRUE);
    } else if(is.null(ymax)){
      ymax <- max(stat,na.rm=TRUE)
    }
    
    par(xaxs = 'i');
    plot.new();
    plot.window(xlim=c(0.5,length(run.list)+0.5), 
                ylim=c(ymin,ymax));

    if(label.thresh){
      failPercentile <- round((sum(fail) / length(fail)) * 100,digits=1)
      title(main=paste0(title,"\nThreshold: ",thresh," (",failPercentile,"%)"));
    } else {
      title(main=title);
    }
    box();
    
    onBadLane <- d$lane.ID %in% laneFailList
      arrowTop <- mean(c(par("usr")[4],par("usr")[4],par("usr")[4],ymax));
      arrowBottom <- mean(c(par("usr")[4],ymax,ymax,ymax));
      breakPoint <- ymax;
      
    statFull <- stat;
    isOvr <- stat > ymax;
    stat <- ifelse(isOvr,arrowBottom,stat);

    color.laneFail <- "orange";
    color.isoFail <- "red";
    color.ok <- "green";
    
    colVector <- ifelse(onBadLane,color.laneFail,ifelse(fail,color.isoFail,color.ok))

    if(any(isOvr)){
      draw.arrows(x=x[isOvr & (!is.na(x))],y=arrowTop,col=colVector[isOvr & (!is.na(x))],cex=0.5);
      
      segments(x0=as.numeric(x[isOvr & (!is.na(x))]),
               y0=as.numeric(rep(arrowBottom,sum(isOvr & (!is.na(x))))),
               x1=as.numeric(x[isOvr & (!is.na(x))]),
               y1=as.numeric(rep(arrowTop,sum(isOvr & (!is.na(x))))),lty=1,col=colVector[isOvr & (!is.na(x))]);
      text(x[isOvr],rep(arrowBottom,sum(isOvr)),paste0("  (",signif(statFull[isOvr],5),")"),adj=c(0,0.5),col=colVector[isOvr],srt=-45,xpd=NA);
      
      qorts.axis.break(axis=2,breakpos=ymax,fill=TRUE,cex=0.5);

    }
    points(x,stat,col=colVector);

    axis(2);
    if(connect.isolated.samples){
        isolated <- fail & (! onBadLane);
        if(any(isolated)){
            points(x[isolated],stat[isolated],col="red",pch=1);
            #points(x[isolated],stat[isolated],col="red",pch=4);
            isolated.samples <- unique(d$sample.ID[isolated]);
            isolated.samples.multiRun <- isolated.samples[sapply(isolated.samples,function(samp){ sum(d$sample.ID == samp) > 1 })]
            for(samp in isolated.samples.multiRun){
                lines(x[d$sample.ID == samp],stat[d$sample.ID == samp],col="red",lty=3);
                points(x[d$sample.ID == samp],stat[d$sample.ID == samp],col="red",pch=4);
            }
        } 
    }
    

    
    if(label.runs){
       text(seq_along(run.list),par("usr")[3],run.titles,srt=srt.run.label,
            cex=cex.run.label,adj=c(1.1,0.5),
            xpd=NA,col="black");
    }
    if(label.lane){
      #Note to self: 
      #  This code SEEMS unnecessarily convoluted. But text() performs REALLY slowly
      #  when handed transparent colors. This hack works around that.
      cex.lane.nums <- 0.75;
      color.lane.good <- "gray88"
      color.lane.bad  <- "red"
      char.ht <- strheight("0",cex=cex.lane.nums)
      lane.label.x <- rep(seq_along(run.list),each=8) + rep(lane.offsets, length(run.list));
      lane.label.y <- rep(rep(seq(par("usr")[3] + char.ht * 3,par("usr")[3],length.out=4),length.out=8),length(run.list));
      lane.label.label <- rep(1:8,length(lane.list))
      color.lane.good <- "gray"
      color.lane.bad  <- "red"
      
      if(sum(! laneFail) > 0){
        text(lane.label.x[!laneFail],
         lane.label.y[!laneFail],
         lane.label.label[!laneFail],
         cex=0.75,
         adj=c(0.5,-0.25),
         xpd=NA,
         col=color.lane.good);
      }
      if(sum(laneFail) > 0){
      text(lane.label.x[laneFail],
         lane.label.y[laneFail],
         lane.label.label[laneFail],
         cex=0.75,
         adj=c(0.5,-0.25),
         xpd=NA,
         col=color.laneFail);
      }
    }
    
    borders <- c(seq_along(run.list)-0.5);
    borders <- borders[-1];
    axis(1,at=borders,labels=FALSE);
    abline(v=borders,lty=3,col="gray");
    abline(h=thresh,lty=3,col="gray");
}

points.arrowLimited <- function(x,y,ymax,col = "black",pch = 1,
                                label.text = TRUE, label.text.signif = 5,arrow.cex = 0.5,
                                cex = 1,
                                ...){
  x <- as.numeric(x);
  y <- as.numeric(y);
  if(length(x) == 1 && length(y) > 1) x <- rep(x,length(y));
  if(length(y) == 1 && length(x) > 1) y <- rep(y,length(x));
  
  if(length(col) == 1) col <- rep(col,length(x));
  if(length(pch) == 1) pch <- rep(pch,length(x));
  
  keep <- (!is.na(x)) & (! is.na(y));
  if(any(keep)){
  x <- x[keep];
  y <- y[keep];
  col <- col[keep];
  pch <- pch[keep];
  
  isOvr <- y > ymax;
  arrowTop <- mean(c(par("usr")[4],par("usr")[4],par("usr")[4],ymax));
  arrowBottom <- mean(c(par("usr")[4],ymax,ymax,ymax));
  breakPoint <- ymax;
  
    if(any(isOvr)){
      draw.arrows(x=x[isOvr],y=arrowTop,col=col[isOvr],cex=arrow.cex,...);
      
      segments(x0=as.numeric(x[isOvr]),
               y0=as.numeric(rep(arrowBottom,sum(isOvr))),
               x1=as.numeric(x[isOvr]),
               y1=as.numeric(rep(arrowTop,sum(isOvr))),lty=1,col=col[isOvr],...);
      if(label.text){
        text(x[isOvr],rep(arrowBottom,sum(isOvr)),paste0("  (",signif(y[isOvr],label.text.signif),")"),
             adj=c(0,0.5),col=col[isOvr],srt=-45,xpd=NA,...);
      }
      qorts.axis.break(axis=2,breakpos=ymax,fill=TRUE,cex=0.5,...);

    }
    points(x,ifelse(isOvr,arrowBottom,y),col=col,pch=pch,cex=cex,...);
  }
}

plot.lane.stat <- function(run.list,lane.list,stat, d,ymin = NULL, ymax = NULL, title = "",fail,thresh,laneFractionThresh = 0.33,runFractionThresh = laneFractionThresh,
                           label.lane=TRUE, label.runs = TRUE, run.titles = run.list, label.thresh = TRUE,
                           add.vioplot = TRUE,
                           cex.run.label = 0.75, cex.lane.label = 0.75, srt.run.label = 90,
                           preplotfunc = NULL,
                           point.pch = NULL,markLanes=c(),drawLaneBorders = TRUE,
                           laneIsFlagged = c(), colorByFlagOnly = FALSE,laneColors = NULL,lineOnFail=FALSE){
    lane.list <- unlist(lapply(run.list,function(runid){ paste0(runid,"_",1:8) }));
    runRate <- flag.run.calc(run.list, fail,d)
    laneRate <- flag.lane.calc(lane.list,fail,d)
    if(is.null(point.pch)) point.pch <- rep(1,length(stat));
    laneFail <- f.na(unlist(lapply(run.list,function(runid){
          lanes <- paste0(runid,"_",1:8);
          sapply(lanes,function(laneid){
                if(laneid %in% lane.list){
                    laneRate[[which(lane.list == laneid)]] > laneFractionThresh
                } else {
                    FALSE
                }
              })
        })))
    #names(laneFail) <- paste0(rep(run.list,each=8),"_",rep(1:8,length(run.list)));
    laneFailList <- lane.list[laneFail]
    runFail <- f.na(unlist(lapply(run.list,function(runid){
                runRate[[which(run.list == runid)]]  > laneFractionThresh
        })))
    lane.offsets <- (1:8 / 8) - 0.5625

    x <- sapply(seq_along(stat),function(i){
          which(d$FLOWCELL[i] == run.list) + lane.offsets[d$LANE[i]]
        });
    if(is.null(ymin)){
        ymin <- min(stat,na.rm=TRUE);
    }
    if(is.null(ymax)){
        ymax <- max(stat,na.rm=TRUE);
    }
    
    #if(suppressWarnings(requireNamespace("vioplot")) && add.vioplot){
    if(add.vioplot){
       xlim <- c(0.5,length(run.list)+1.5)
    } else {
       xlim <- c(0.5,length(run.list)+0.5)
    }
    
    par(xaxs = 'i');
    plot.new();
    plot.window(xlim=xlim, 
                ylim=c(ymin,ymax + abs(ymax-ymin)*0.05));
    axis(2);
    
    if(lineOnFail){
      lane.label.x <- rep(seq_along(run.list),each=8) + rep(lane.offsets, length(run.list));
      abline(v=lane.label.x[which(laneFail)],col=color2transparent("red",50),lty=3);
    }
    
    if(! is.null(preplotfunc)){
      preplotfunc();
    }
    if(label.thresh){
    failPercentile <- round((sum(fail) / length(fail)) * 100,digits=1)
    title(main=paste0(title,"\nThreshold: ",signif(thresh,6)," (",failPercentile,"%)"));
    } else {
        title(main=title);
    }
    #box();
    rect(par("usr")[1],par("usr")[3],length(run.list)+0.5, par("usr")[4],xpd=NA);
    
    isOvr <- stat > ymax;
    
    #+ runif(length(stat),-0.1,0.1)
    if(!is.null(laneColors)){
      col.pch <- laneColors[d$lane.ID];
    } else if(colorByFlagOnly){
      col.pch <- ifelse(d$lane.ID %in% laneFailList,"red",ifelse(d$lane.ID %in% laneIsFlagged,"orange","green"))
    } else {
      col.pch <- ifelse(fail,"red",ifelse(d$lane.ID %in% laneIsFlagged,"orange","green"))
    }
    
    points.arrowLimited(x=x,y=stat,ymax=ymax,col=col.pch,pch=point.pch); 
    
    #points(x[!isOvr],stat[!isOvr],col=ifelse(fail[!isOvr],"red","green"),pch=point.pch[!isOvr]);
    #if(sum(isOvr) > 0){
    #  points(x[isOvr],rep(ymax,sum(isOvr)),col=ifelse(fail[isOvr],"red","green"),pch = 2);
    #  points(x[isOvr],rep(ymax,sum(isOvr)),col=ifelse(fail[isOvr],"red","green"),pch=point.pch[isOvr]);
    #  text(x[isOvr],rep(ymax,sum(isOvr)),labels=round(stat[isOvr],3),col="red",adj=-0.1);
    #}
    
    if(label.runs){
       text(seq_along(run.list),par("usr")[3],run.titles,srt=srt.run.label,
            cex=cex.run.label,adj=c(1.1,0.5),
            xpd=NA,col=ifelse(runFail,"red","black"));
    }
    

    
    if(label.lane){
      #I have removed the hack. let's see what happens.
      #Note to self: 
      #  This code SEEMS unnecessarily convoluted. But text() performs REALLY slowly
      #  when handed transparent colors. This hack works around that.
      cex.lane.nums <- 0.75;
      color.lane.good <- "gray88"
      color.lane.bad  <- "red"
      char.ht <- strheight("0",cex=cex.lane.nums)
      lane.label.x <- rep(seq_along(run.list),each=8) + rep(lane.offsets, length(run.list));
      lane.label.y <- rep(rep(seq(par("usr")[3] + char.ht,par("usr")[3],length.out=2),length.out=8),length(run.list));
      lane.label.label <- rep(1:8,length(lane.list))
      color.lane.good <- "gray"
      color.lane.bad  <- "red"
      color.lane.flaggedGood <- "orange";
      
      if(is.null(laneColors)){
        laneIDcolors <- ifelse(laneFail,color.lane.bad,ifelse(lane.list %in% laneIsFlagged,color.lane.flaggedGood,color.lane.good))
      } else {
        laneIDcolors <- laneColors[lane.list]
      }
      
      text(lane.label.x,
         lane.label.y,
         lane.label.label,
         cex=0.75,
         adj=c(0.5,-0.25),
         xpd=NA,
         col=laneIDcolors);
      
      #if(sum(! laneFail) > 0){
      #  text(lane.label.x[!laneFail],
      #   lane.label.y[!laneFail],
      #   lane.label.label[!laneFail],
      #   cex=0.75,
      #   adj=c(0.5,-0.25),
      #   xpd=NA,
      #   col=color.lane.good);
      #}
      #if(sum(laneFail) > 0){
      #text(lane.label.x[laneFail],
      #   lane.label.y[laneFail],
      #   lane.label.label[laneFail],
      #   cex=0.75,
      #   adj=c(0.5,-0.25),
      #   xpd=NA,
      #   col=color.lane.bad);
      #   if(length(markLanes) > 0){
      #    star.lane <- lane.list %in% markLanes;
      #    if(sum(star.lane) > 0){
      #      text(lane.label.x[star.lane],rep(max(lane.label.y) + char.ht,sum(star.lane)),rep("*",sum(star.lane)),cex=0.75,adj=c(0.5,-0.25),col="red");
      #    }
      #   }
      # }
    }
    
    #if(suppressWarnings(requireNamespace("vioplot")) && add.vioplot){
    if(add.vioplot){
        #vioplot::vioplot(stat, col="red", add=TRUE, at= length(run.list)+1);
        #plotKD.violin(stat,col="red",add=TRUE,at = length(run.list)+1,type="violin");
        plotKD.violin(stat,col="red",add=TRUE,at = length(run.list)+0.5,type="left",whiskerQuantile=0.05);
    }

    
    borders <- c(seq_along(run.list)-0.5);
    borders <- borders[-1];
    axis(1,at=borders,labels=FALSE);
    if(drawLaneBorders){
      abline(v=borders,lty=3,col="gray");
    } else {
      axis(1,at=borders,labels=FALSE,tcl=1);
      axis(3,at=borders,labels=FALSE,tcl=1);
      
    }
    abline(h=thresh,lty=3,col="gray");
}

plot.lane.stat.isolated.OLD <- function(run.list,lane.list,stat, d, flagged,ymin = NULL, ymax = NULL, title = "",fail,thresh,laneFractionThresh = 0.33,runFractionThresh = laneFractionThresh,
                           label.lane=TRUE, label.runs = TRUE, run.titles = run.list, label.thresh = TRUE,
                           add.vioplot = TRUE, connect.isolated.samples = TRUE, drop.clean.runs = TRUE,
                           cex.run.label = 0.75, cex.lane.label = 0.75, srt.run.label = 90){
    original.run.list <- run.list;
    
    if(drop.clean.runs){
      keep.run <- sapply(run.list,function(r){ any(fail[ d$FLOWCELL == r ]) });
      if(any(keep.run)){
          run.list <- run.list[ keep.run ];
      } else {
         #else leave it alone.
      }
    }
    lane.list <- unlist(lapply(run.list,function(runid){paste0(runid,"_",1:8)}));
    

    runRate <- flag.run.calc(run.list, fail,d)
    laneRate <- flag.lane.calc(lane.list,fail,d)
    laneFail <- f.na(unlist(lapply(run.list,function(runid){
          lanes <- paste0(runid,"_",1:8);
          sapply(lanes,function(laneid){
                if(laneid %in% lane.list){
                    laneRate[[which(lane.list == laneid)]] > laneFractionThresh
                } else {
                    FALSE
                }
              })
        })))
    #names(laneFail) <- paste0(rep(run.list,each=8),"_",rep(1:8,length(run.list)));
    
    runFail <- f.na(unlist(lapply(run.list,function(runid){
                runRate[[which(run.list == runid)]]  > laneFractionThresh
        })))
    lane.offsets <- (1:8 / 8) - 0.5625

    x <- sapply(seq_along(stat),function(i){
          which(d$FLOWCELL[i] == run.list) + lane.offsets[d$LANE[i]]
        });
    if(is.null(ymin)){
        ymin <- min(stat,na.rm=TRUE);
    }
    if(is.null(ymax)){
        ymax <- max(stat,na.rm=TRUE);
    }
    
    #if(suppressWarnings(requireNamespace("vioplot")) && add.vioplot){
    if(add.vioplot){
       xlim <- c(0.5,length(run.list)+1.5)
    } else {
       xlim <- c(0.5,length(run.list)+0.5)
    }
    
    par(xaxs = 'i');
    plot.new();
    plot.window(xlim=xlim, 
                ylim=c(ymin,ymax));

    if(label.thresh){
    failPercentile <- round((sum(fail) / length(fail)) * 100,digits=1)
    title(main=paste0(title,"\nThreshold: ",thresh," (",failPercentile,"%)"));
    } else {
        title(main=title);
    }
    #box();
    rect(par("usr")[1],par("usr")[3],length(run.list)+0.5, par("usr")[4],xpd=NA);
    
    isOvr <- stat > ymax;
    
    #+ runif(length(stat),-0.1,0.1)
    curr <- ! isOvr;
    
    points(x[curr],stat[curr],col=ifelse(flagged[curr],ifelse(fail[curr],"yellow","blue"),ifelse(fail[curr],"red","green")));
    if(sum(isOvr) > 0){
      curr <- isOvr;
      
      points(x[curr],stat[curr],col=ifelse(flagged[curr],ifelse(fail[curr],"yellow","blue"),ifelse(fail[curr],"red","green")));
      text(x[curr],rep(ymax,sum(curr)),labels=round(stat[curr],3),col=ifelse(flagged[curr],"yellow","red"),adj=-0.1);
    }
    
    if(connect.isolated.samples){
        isolated <- fail & (! flagged)
        if(any(isolated)){
            points(x[isolated],stat[isolated],col="red",pch=1);
            #points(x[isolated],stat[isolated],col="red",pch=4);
            isolated.samples <- unique(d$sample.ID[isolated]);
            isolated.samples.multiRun <- isolated.samples[sapply(isolated.samples,function(samp){ sum(d$sample.ID == samp) > 1 })]
            for(samp in isolated.samples.multiRun){
                lines(x[d$sample.ID == samp],stat[d$sample.ID == samp],col="red",lty=3);
                points(x[d$sample.ID == samp],stat[d$sample.ID == samp],col="red",pch=4);
            }
        } 
    }
    
    
    axis(2);
    
    if(label.runs){
       text(seq_along(run.list),par("usr")[3],run.titles,srt=srt.run.label,
            cex=cex.run.label,adj=c(1.1,0.5),
            xpd=NA,col=ifelse(runFail,"red","black"));
    }
    if(label.lane){
      #Note to self: 
      #  This code SEEMS unnecessarily convoluted. But text() performs REALLY slowly
      #  when handed transparent colors. This hack works around that.
      cex.lane.nums <- 0.75;
      color.lane.good <- "gray88"
      color.lane.bad  <- "red"
      char.ht <- strheight("0",cex=cex.lane.nums)
      lane.label.x <- rep(seq_along(run.list),each=8) + rep(lane.offsets, length(run.list));
      lane.label.y <- rep(rep(seq(par("usr")[3] + char.ht,par("usr")[3],length.out=2),length.out=8),length(run.list));
      lane.label.label <- rep(1:8,length(lane.list))
      color.lane.good <- "gray"
      color.lane.bad  <- "red"
      
      if(sum(! laneFail) > 0){
        text(lane.label.x[!laneFail],
         lane.label.y[!laneFail],
         lane.label.label[!laneFail],
         cex=0.75,
         adj=c(0.5,-0.25),
         xpd=NA,
         col=color.lane.good);
      }
      if(sum(laneFail) > 0){
      text(lane.label.x[laneFail],
         lane.label.y[laneFail],
         lane.label.label[laneFail],
         cex=0.75,
         adj=c(0.5,-0.25),
         xpd=NA,
         col=color.lane.bad);
      }
    }
    
    #if(suppressWarnings(requireNamespace("vioplot")) && add.vioplot){
    if(add.vioplot){
        #vioplot::vioplot(stat, col="red", add=TRUE, at= length(run.list)+1);
        #plotKD.violin(stat,col="red",add=TRUE,at = length(run.list)+1,type="violin");
        plotKD.violin(stat,col="red",add=TRUE,at = length(run.list)+0.5,type="left");
    }

    
    borders <- c(seq_along(run.list)-0.5);
    borders <- borders[-1];
    axis(1,at=borders,labels=FALSE);
    abline(v=borders,lty=3,col="gray");
    abline(h=thresh,lty=3,col="gray");
}


addFlag <- function(f=NULL,flagged,title = "UNKFLAG"){
    if(is.null(f)){
       f <- list(FLAG = rep(FALSE,length(flagged)), 
                 FLAG.REASON=rep("",length(flagged)), 
                 FLAG.FIRST.REASON=rep("",length(flagged)));
    }
    FLAG <- f$FLAG
    FLAG.REASON <- f$FLAG.REASON
    FLAG.FIRST.REASON <- f$FLAG.FIRST.REASON;

        
    FLAG.REASON <- ifelse(flagged, ifelse(FLAG,paste0(FLAG.REASON,",",title),title),FLAG.REASON);
    
    FLAG.FIRST.REASON <- ifelse(flagged,ifelse(FLAG,FLAG.FIRST.REASON,title),FLAG.FIRST.REASON);
    FLAG <- FLAG | flagged;
    
    return(list(FLAG = FLAG, FLAG.REASON=FLAG.REASON,FLAG.FIRST.REASON=FLAG.FIRST.REASON))
}

plot.hist.with.norm.withCm <- function(v,title="",xmin = 0,breaks = 100){
    h <- hist(v,breaks=breaks,plot=FALSE);  
    normH <- dnorm(h$mids, mean = mean(v), sd = sd(v))
    plot.new();
    
    if(is.null(xmin)){
        xmin <- min(h$breaks)
    }
    plot.window(xlim=c(xmin,max(h$breaks)), 
                ylim=c(0,1));
    title(main=title);
    box();
    pretty.at <- pretty(par("usr")[1:2],n=10)
    pretty.mini <- seq(from=min(pretty.at),to=max(pretty.at),by=abs(pretty.at[1] - pretty.at[2])/2);
    axis(1,at=pretty.at); 
    axis(1,at=pretty.mini,labels=FALSE,tcl=-0.2); 
    axis(2);
    lines(h$mids,h$density/max(h$density),col="red");
    points(h$mids,h$density/max(h$density),col="red");
    lines(h$mids,normH/max(h$density),col="blue");

    lines(h$mids,cumsum(h$counts)/sum(h$counts),col="red",lty=3);
    points(h$mids,cumsum(h$counts)/sum(h$counts),col="red",cex=0.5);
    lines(h$mids,pnorm(h$mids, mean = mean(v), sd = sd(v)),col="blue",lty=3);

}



plot.hist.with.norm <- function(v,title="",xmin = 0,breaks = 100){
    h <- hist(v,breaks=breaks,plot=FALSE);  
    normH <- dnorm(h$mids, mean = mean(v), sd = sd(v), log = FALSE)
    plot.new();
    if(is.null(xmin)){
        xmin <- min(h$breaks)
    }
    plot.window(xlim=c(xmin,max(h$breaks)), ylim=c(min(h$density),max(h$density)));
    title(main=title);
    box();
    pretty.at <- pretty(par("usr")[1:2],n=10)
    pretty.mini <- seq(from=min(pretty.at),to=max(pretty.at),by=abs(pretty.at[1] - pretty.at[2])/2);
    axis(1,at=pretty.at); 
    axis(1,at=pretty.mini,labels=FALSE,tcl=-0.2); 
    axis(2);
    lines(h$mids,h$density,col="red");
    points(h$mids,h$density,col="red");
    lines(h$mids,normH,col="blue");




}

plotKD.violin <- function(v, type = c("violin","left","right","up","down"),
                             col="magenta",border="black",
                             lty=1,
                             lwd=1,
                             rectCol="black",
                             colMed="white",
                             pchMed=19,
                             at=NULL,
                             add=FALSE,
                             wex=1,
                             drawRect=TRUE,
                             horizontal=FALSE,
                             wd = 0.95,
                             rectWd = wd / 16,
                             truncateToPlot = TRUE,
                             whiskerQuantile = NULL,
                             ...){
  if(type == "up" || type == "down") horizontal <- TRUE;
  if(type == "left" || type == "right") horizontal <- FALSE;
  
  if(!is.null(whiskerQuantile)){
    if(length(whiskerQuantile) == 1){
      whiskerQuantile <- c(whiskerQuantile,1-whiskerQuantile);
      whiskerQuantile <- c(min(whiskerQuantile),max(whiskerQuantile));
    }
  }
  
  #wd <- 0.95;
  #rectWd <- wd / 16;
  h <- density(v);  
  qt <- quantile(v, probs = c(0.25,0.5,0.75));
  if(! is.null(whiskerQuantile)){
    whiskers <- quantile(v, probs = whiskerQuantile);
  }
  
  type <- match.arg(type);
  if(is.null(at)){
    at <- 1
  }
  
  if(! add){
    plot.new();
    if(horizontal){
      plot.window(ylim=c(at-wd,at+wd),xlim=c(min(h$x),max(h$x)));
      axis(1)
    } else {
      plot.window(xlim=c(at-wd,at+wd),ylim=c(min(h$x),max(h$x)));
      axis(2)
    }
    box();
  }
  

  if(truncateToPlot){
    if(horizontal){
      hkeep <- h$x > par("usr")[1] & h$x < par("usr")[2];
      h$x <- c(par("usr")[1],h$x[hkeep],par("usr")[2]);
      h$y <- c(0,h$y[hkeep],0);
      
    } else {
      hkeep <- h$x > par("usr")[3] & h$x < par("usr")[4];
      h$x <- c(par("usr")[3],h$x[hkeep],par("usr")[4]);
      h$y <- c(0,h$y[hkeep],0);
    }
  }
  
  if(! horizontal){

    if(type=="violin"){
      y.left <- h$x;
      y.rght <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x.left <- at - ((h$y / rawXSpan) * (wd/2))
      x.rght <- ((h$y / rawXSpan) * (wd/2)) + at;

      polygon(c(x.left,rev(x.rght)),
              c(y.left,rev(y.rght)),border=border,col=col,...);

      if(drawRect) rect(at-rectWd/2,qt[1],at+rectWd/2,qt[3],border=border,col=rectCol,...);
      lines(c(at-rectWd/2,at+rectWd/2),c(qt[2],qt[2]),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(x0=at-rectWd/4,x1=at+rectWd/4,y0=c(qt[1],qt[2]),y1=c(whiskers[1],whiskers[2]),col=border,...)
    } else if(type == "left"){
      y <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x <- (h$y / rawXSpan) * (wd) + at;
      
      polygon(c(x),
              c(y),border=border,col=col,...);

      if(drawRect) rect(at,qt[1],at+rectWd,qt[3],border=border,col=rectCol,...);
      #points(at+rectWd,qt[2],pch="<",col=colMed,...);
      lines(c(at,at+rectWd),c(qt[2],qt[2]),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(x0=at+rectWd/2,x1=at+rectWd/2,y0=c(qt[1],qt[2]),y1=c(whiskers[1],whiskers[2]),col=border,...)
    } else if(type == "right"){
      y <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x <- at - (h$y / rawXSpan) * (wd)
      
      polygon(c(x),
              c(y),border=border,col=col,...);

      if(drawRect) rect(at,qt[1],at-rectWd,qt[3],border=border,col=rectCol,...);
      #points(at-rectWd,qt[2],pch=">",col=colMed,...);
      lines(c(at,at-rectWd),c(qt[2],qt[2]),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(x0=at-rectWd/2,x1=at-rectWd/2,y0=c(qt[1],qt[2]),y1=c(whiskers[1],whiskers[2]),col=border,...)

    }
  } else {
    if(type=="violin"){
      #NOT SUPPORTED YET!!!
    } else if(type == "up"){
      y <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x <- (h$y / rawXSpan) * (wd) + at;
      
      polygon(c(y),c(x),border=border,col=col,...);

      if(drawRect) rect(qt[1],at,qt[3],at+rectWd,border=border,col=rectCol,...);
      #points(qt[2],at+rectWd,pch="v",col=colMed,...);
      lines(c(qt[2],qt[2]),c(at,at+rectWd),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(y0=at+rectWd/2,y1=at+rectWd/2,x0=c(qt[1],qt[2]),x1=c(whiskers[1],whiskers[2]),col=border,...)

    } else if(type == "down"){
      y <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x <- at - (h$y / rawXSpan) * (wd)
      
      polygon(c(y),c(x),border=border,col=col,...);

      if(drawRect) rect(qt[1],at,qt[3],at-rectWd,border=border,col=rectCol,...);
      #points(qt[2],at-rectWd,pch="^",col=colMed,...);
      lines(c(qt[2],qt[2]),c(at,at-rectWd),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(y0=at-rectWd/2,y1=at-rectWd/2,x0=c(qt[1],qt[2]),x1=c(whiskers[1],whiskers[2]),col=border,...)

    }
  }
}
plot.kd.with.norm <- function(v,title="",xmin = 0,breaks = 100){
    h <- density(v);  
    plot(h);
    
    #title(main=title);
    #box();
    #pretty.at <- pretty(par("usr")[1:2],n=10)
    #pretty.mini <- seq(from=min(pretty.at),to=max(pretty.at),by=abs(pretty.at[1] - pretty.at[2])/2);
    #axis(1,at=pretty.at); 
    #axis(1,at=pretty.mini,labels=FALSE,tcl=-0.2); 
    #axis(2);
    #lines(h$mids,h$density,col="red");
    #points(h$mids,h$density,col="red");
    #lines(h$mids,normH,col="blue");

}








