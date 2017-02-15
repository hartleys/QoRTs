writeHtml.complete <- function(res){

}

writeHtml.multiPlotPage <- function(plotter){
  
}

write.html.table <- function(cath, t, title.class = "base", title.row.class = "base", col.class = rep("base",ncol(t)), class.matrix = NULL){
  cath(paste0("<table>"));
  if(length(col.class == 1)) col.class <- rep(col.class,ncol(t));
  
  if(is.null(class.matrix)){
    class.matrix <- matrix(col.class,nrow=nrow(t),ncol=ncol(t),byrow=TRUE)
  }
  
  cath(paste0("<tr>\n"));
  cath(paste0("<th></th>"));
  for(i in 1:ncol(t)){
    cath(paste0("   <th class=\"",title.class,"\">"));
    cath(paste0("",colnames(t)[i],""));
    cath(paste0("</th>\n"));
  }
  cath(paste0("</tr>\n"));
  
  for(i in 1:nrow(t)){
    cath(paste0("<tr>\n"));
    cath(paste0("   <td class=\"",title.row.class,"\">",row.names(t)[i],"</td>"));
    for(j in 1:ncol(t)){
      cath(paste0("   <td class=\"",class.matrix[i,j],"\">",t[i,j],"</td>"));
    }
    cath(paste0("</tr>\n"));
  }
  cath(paste0("</table>"));
}

write.plot.html <- function(file, title,imgRef,
                            parentRef = NULL, prevRef = NULL, nextRef = NULL, 
                            parentName = "BACK", prevName = "PREV", nextName = "NEXT",
                            indexRefs = NULL, indexNames = 1:length(indexRefs),
                            divClass = "plotPageBody", linkTableClass = "centered", linkTdClass = "navCell",
                            plotContainerClass = "plotContainer", multiPlotClass = "multiPlot",
                            plotWidth = "1850px"){
  
  cat(paste0(
  "<html><head>\n",
  "<title> ",title," </title>\n",
  "<link rel=\"stylesheet\" type=\"text/css\" href=\"../styles.css\"></head>\n",
  "<body>\n",
  "<div class=\"",divClass,"\">\n"
  ),
  file = file, append = FALSE);
  
  cath <- function(...){
    cat(paste0(...),file=file,append=TRUE)
  }
  
  if(! is.null(parentRef)){
    cath("<h4>","<a href=\"",parentRef,"\">",parentName,"</a>","</h4>\n");
  }
  
  if((! is.null(prevRef)) || (! is.null(nextRef))){
    cath("<table class=\"",linkTableClass,"\"><tr>\n");
    if((! is.null(prevRef))){
      cath("  <td class=\"",linkTdClass,"\">","<a href=\"",prevRef,"\">" ,prevName,"</a>","</td>\n");
    }
    if((! is.null(nextRef))){
      cath("  <td class=\"",linkTdClass,"\">","<a href=\"",nextRef,"\">" ,nextName,"</a>","</td>\n");
    }
    cath("</tr></table>\n");
  }
  
  if(! is.null(indexRefs)){
    cath("<table class=\"",linkTableClass,"\"><tr>\n");
    for(i in 1:length(indexRefs)){
      cath("  <td class=\"",linkTdClass,"\">","<a href=\"",indexRefs[[i]],"\">" ,indexNames[[i]],"</a>","</td>\n");
    }
    cath("</tr></table>\n");
  }
  
  cath(
          "<div class=\"",plotContainerClass,"\">\n",
          "  <img src=\"",imgRef,"\"","width=",plotWidth,", class=\"",multiPlotClass,"\">\n",
          "</div>\n"
      );
  
  cath("</div></body></html>\n");

  
}

write.index.html <- function(htmlDir,
                             multiPlotSetList,
                             byLane.titleList, byLane.shortTitleList,byLane.flags,byLane.hasFilter, byLane.failPct, byLane.threshSet,laneFailed.failCtMatrix
                             ){
  file = paste0(htmlDir,"index.html");
  
  cat(paste0(
  "<html><head>",
  "<title> CCSS Quality Control </title>",
  "<link rel=\"stylesheet\" type=\"text/css\" href=\"styles.css\"></head>",
  "<body>",
  "<div class=\"mainTextBody\">"
  ),
  file = file, append = FALSE);
  
  cath <- function(...){
    cat(paste0(...),file=file,append=TRUE)
  }
  
  cath(paste0("<h1> QC SUMMARY </h1>"));
  
  cath("<h2> Summary Multiplots: </h2>");
  multiPlotSetList$plot.title <- paste0("<a href=\"./byGroup/",multiPlotSetList$plot.ID,".html\">",multiPlotSetList$plot.title,"</a>");
  write.html.table(cath=cath, t=multiPlotSetList,
                   title.row.class="numberCellMed",title.class="numberCellMed",
                   col.class = "numberCellMed")
  
  cath(paste0(
     "<h2> Summary Plots, By Lane </h2>",
     "<table>"
  ));
  
  cath(paste0("<tr>\n"));
  cath(paste0("   <th class=\"topLink\"> Flag Filter </th>\n"));
  cath(paste0("   <th class=\"topLink\"> Link </th>\n"));
  cath(paste0("   <th class=\"topLink\"> Flag Threshold </th>\n"));
  cath(paste0("   <th class=\"topLink\"> % BAMs flagged </th>\n"));
  cath(paste0("   <th class=\"topLink\"> # Lanes Flagged </th>\n"));
  cath(paste0("   <th class=\"topLink\"> # New Lanes Flagged </th>\n"));
  cath(paste0("</tr>\n"))
  
  for(j in seq_along(byLane.shortTitleList)){
        cath(paste0("<tr>\n"));
        cath(paste0("   <td class=\"topLink\">"));
        cath(paste0(byLane.titleList[[j]],":"));
        cath(paste0("</td>\n"));
        cath(paste0("   <td class=\"topLink\">"));
        cath(paste0("<a href=\"",paste0("summaryPlots/",byLane.shortTitleList[[j]],".html"),"\">",byLane.shortTitleList[[j]],"</a>"));
        cath(paste0("</td>\n"));
        
        cath(paste0("   <td class=\"numberCell\">"));
        cath(paste0("",byLane.threshSet[[j]]));
        cath(paste0("</td>\n"));
        
        cath(paste0("   <td class=\"numberCell\">"));
        cath(paste0("",round(100*byLane.failPct[[j]],1)));
        cath(paste0("</td>\n"));
        
        if(byLane.hasFilter[[j]]){
          numFailed    <- sum(sapply(byLane.flags$FLAG.REASON,function(r){any(strsplit(r,",",fixed=TRUE)[[1]] == byLane.shortTitleList[[j]])}))
          numNewFailed <- sum(byLane.flags$FLAG.FIRST.REASON == byLane.shortTitleList[[j]]);
        } else {
          numFailed <- "N/A";
          numNewFailed <- "N/A";
        }
        
        cath(paste0("   <td class=\"numberCell\">"));
        cath(paste0("",numFailed));
        cath(paste0("</td>\n"));
        
        cath(paste0("   <td class=\"numberCell\">"));
        cath(paste0("",numNewFailed));
        cath(paste0("</td>\n"));
        
        cath(paste0("</tr>\n"));
  }
  cath(paste0("</table>"));
  
  cath(paste0(
  "<h2> Flagged Lanes: </h2>"
  ));
  #numberCellSmRed
  write.html.table(cath=cath, t=laneFailed.failCtMatrix,
                   title.row.class="numberCellSm",title.class="numberCellSm",
                   class.matrix= apply(laneFailed.failCtMatrix,MARGIN=c(1,2),FUN=function(e){
                       ifelse(f.na(eval(parse(text=e)) > 0.5),"numberCellSmRed","numberCellSm");
                     })
                  )
  

  
  
  
  cath("<a href=\"split/set.html\"> <h2> Outlier ID Plots </h2> </a>");
  
  
  
  cath("</div></body></html>");
}

make.byLane.summary.plot.html <- function(titleList,shortTitleList, plotDir, htmlDir){
    summaryHtmlDir <- paste0(htmlDir,"/summaryPlots/");
    
    indexRefs = sapply(1:length(shortTitleList),function(i){
      paste0("./",shortTitleList[[i]],".html");
    })
    indexNames = sapply(1:length(shortTitleList),function(i){
      paste0(" ",shortTitleList[[i]]," ");
    })
    
    for(i in seq_along(shortTitleList)){
       file <- paste0(summaryHtmlDir,shortTitleList[[i]],".html");
       title <- paste0(titleList[[i]]);
       imgRef <- paste0("../../summaryPlots/byLane/",shortTitleList[[i]],".png");
       if(i > 1){
         prevRef = paste0("./",shortTitleList[[i-1]],".html");
       }  else {
         prevRef = NULL
       }
       if(i < length(shortTitleList)){
         nextRef = paste0("./",shortTitleList[[i+1]],".html");
       } else {
         nextRef = NULL
       }
       write.plot.html(file, title,imgRef,
                            parentRef = "../index.html", prevRef = prevRef, nextRef = nextRef, 
                            parentName = "BACK", prevName = "PREV", nextName = "NEXT",
                            indexRefs = indexRefs, indexNames = indexNames,
                            divClass = "plotPageBody", linkTableClass = "centered", linkTdClass = "navCell",
                            plotContainerClass = "plotContainer", multiPlotClass = "lanePlot",
                            plotWidth = "800px")
    }
}


make.byLane.summary.plot.html.OLD <- function(titleList,shortTitleList, plotDir, htmlDir){
    #summaryPlotDir <- paste0(plotDir,"/summaryPlots/byLane/");
    summaryHtmlDir <- paste0(htmlDir,"/summaryPlots/");
    
    for(i in seq_along(shortTitleList)){
      file = paste0(summaryHtmlDir,shortTitleList[[i]],".html");

      message("writing html: ",file)

      #Start HTML:
      cat(paste0("<html><head>\n",
                 "<title> ", titleList[[i]], " </title>\n",
                 "<link rel=\"stylesheet\" type=\"text/css\" href=\"../styles.css\">",
                 "</head>\n",
                 "<body>\n"
                 ),
          file = file, append = FALSE);
      #start up table:
      
      cat(paste0("<table class=\"centered\">\n"),file = file, append = TRUE);
      
      cat(paste0("<tr>\n"),file = file, append = TRUE);
      if(i > 1){
        cat(paste0("<td class=\"topLink\">"),file = file, append = TRUE);
        cat(paste0("<a href=\"",paste0("./",shortTitleList[[i-1]],".html"),"\">","PREV","</a>\n"),file = file, append = TRUE);
        cat(paste0("</td>"),file = file, append = TRUE);
      } else {
        cat(paste0("<td> </td>"),file = file, append = TRUE);
      }
      
      if(i < length(shortTitleList)){
        cat(paste0("<td class=\"topLink\">"),file = file, append = TRUE);
        cat(paste0("<a href=\"",paste0("./",shortTitleList[[i+1]],".html"),"\">","NEXT","</a>\n"),file = file, append = TRUE);
        cat(paste0("</td>"),file = file, append = TRUE);
      } else {
        cat(paste0("<td> </td>\n"),file = file, append = TRUE);
      }
      cat(paste0("</tr>\n"),file = file, append = TRUE);
      cat(paste0("</table>\n"),file = file, append = TRUE);

      cat(paste0("<table class=\"centered\">\n"),file = file, append = TRUE);
      cat(paste0("<tr>\n"),file = file, append = TRUE);
      for(j in seq_along(shortTitleList)){
        cat(paste0("<td class=\"topLink\">"),file = file, append = TRUE);
        cat(paste0("<a href=\"",paste0("./",shortTitleList[[j]],".html"),"\">",shortTitleList[[j]],"</a>"),file = file, append = TRUE);
        cat(paste0("</td>\n"),file = file, append = TRUE);
      }
      cat(paste0("</tr>\n"),file = file, append = TRUE);
      cat(paste0("</table>\n"),file = file, append = TRUE);
      
      #Add image:
      #if(samp %in% has.sample.plot){
      cat(paste0(
            "<div class=\"plotContainer\">\n",
            "<img src=\"../../summaryPlots/byLane/",shortTitleList[[i]],".png\"","width=1600>\n",
            "</div>"
          ),file = file, append = TRUE);
      #}

      #Finish up:
      cat(paste0("</body></html>"),file = file, append = TRUE);

    }
      
}

generateHTML.samplePage  <- function(res,
                                     samp,
                                     html.dir,
                                     has.sample.plot = c()
                                     ){
    decoder <- res@decoder
    file = paste0(html.dir,"/samp/",samp,".html");
    
    message("writing html: ",file)
    
    #Start HTML:
    cat(paste0("<html><head>\n",
               "<title> Sample Info: ", samp, " </title>\n",
               "<link rel=\"stylesheet\" type=\"text/css\" href=\"../styles.css\">",
               "</head>\n",
               "<body>\n"
               ),
        file = file, append = FALSE);
    
    #start up table:
    cat(paste0(
          "<div class=\"linkTableContainer\">\n",
               "  <table>\n"
        ),file = file, append = TRUE);
    
    unique.IDs <- res@decoder$unique.ID[res@decoder$sample.ID == samp];
    table <- res@summaryTable;
    table <- table[,colnames(table) %in% unique.IDs];
    d <- t(decoder[match(colnames(table),decoder$unique.ID),]);
    colnames(d) <- colnames(table);
    table <- rbind(d,table);
    
    cat(paste0("  <tr>\n"),file = file, append = TRUE); 
    cat(paste0("  <td><b>FIELD</b></td>\n"),file = file, append = TRUE); 
    for(j in 1:ncol(table)){
        cat(paste0("  <td><b>",colnames(table)[j],"</b></td>\n"),file = file, append = TRUE); 
    }
    cat(paste0("  <td><b>FIELD</b></td>\n"),file = file, append = TRUE); 
    for(j in 1:ncol(table)){
        i.col2 <- i + ceiling(nrow(table)/2);
        cat(paste0("  <td><b>",colnames(table)[j],"</b></td>\n"),file = file, append = TRUE); 
    }
    cat(paste0("  </tr>\n"),file = file, append = TRUE); 
    
    for(i in 1:ceiling(nrow(table)/2)){
      i.col2 <- i + ceiling(nrow(table)/2);
      
      cat(paste0("  <tr>\n"),file = file, append = TRUE); 
      cat(paste0("  <td><b>",gsub("_","&#8203;",rownames(table)[[i]],fixed=T),"</b></td>\n"),file = file, append = TRUE); 
      
      for(j in 1:ncol(table)){
        cat(paste0("  <td>",table[i,j],"</td>\n"),file = file, append = TRUE); 
      }
      
      if(i.col2 < nrow(table)){
        cat(paste0("  <td><b>",gsub("_","&#8203;",rownames(table)[[i.col2]],fixed=T),"</b></td>\n"),file = file, append = TRUE); 
        for(j in 1:ncol(table)){
          cat(paste0("  <td>",table[i.col2,j],"</td>\n"),file = file, append = TRUE); 
        }
      } else {
        cat(paste0("  <td><b>"," ","</b></td>\n"),file = file, append = TRUE); 
        for(j in 1:ncol(table)){
          cat(paste0("  <td>"," ","</td>\n"),file = file, append = TRUE); 
        }
      }
      cat(paste0("  </tr>\n"),file = file, append = TRUE); 
    }
    
    #Finish up table:
    cat(paste0("  </table>\n",
               "</div>\n"
        ),file = file, append = TRUE);

    #Add image:
    if(samp %in% has.sample.plot){
    cat(paste0(
          "<div class=\"plotContainer\">\n",
          "<img src=\"./",samp,".png\"","width=1850>\n",
          "</div>"
        ),file = file, append = TRUE);
    }
        
    #Finish up:
    cat(paste0("</body></html>"),file = file, append = TRUE);

}