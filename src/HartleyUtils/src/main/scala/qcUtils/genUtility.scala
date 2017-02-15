package qcUtils

import net.sf.samtools._
import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.JavaConversions._;

import internalUtils.genomicUtils._;
import internalUtils.optionHolder._;

object genUtility {
  

  
   def createOutputFile(outdir : String, fileSuffix : String, fileSummary : String, docWriter : DocWriterUtil, tableLines : (String,String,String)*) : WriterUtil = {
     val writer = openWriterSmart_viaGlobalParam(outdir + "." + fileSuffix);
     
     if(docWriter == null){
       //do nothing
     } else {
       /*docWriter.write("\t\t\n");
       docWriter.write("FILE:\t"+fileSuffix+"\t"+fileSummary+"\n");
       for((a,b,c) <- tableLines){
         docWriter.write(a+"\t"+b+"\t"+c+"\n");
       }*/
       docWriter.registerFile(fileSuffix,fileSummary,tableLines:_*);
     }
     return writer;
   } 
   
   
   
   
}

abstract class genUtility {
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null)

  def getUtilityName : String;
  
   def createOutputFile(outdir : String, fileSuffix : String, fileSummary : String, docWriter : DocWriterUtil, tableLines : (String,String,String)*) : WriterUtil = {
     return(genUtility.createOutputFile(outdir=outdir,fileSuffix=fileSuffix,fileSummary=fileSummary,docWriter=docWriter,tableLines=tableLines:_*));
   }
}











