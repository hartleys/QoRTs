package fileConversionUtils

import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
import java.io.OutputStream
import java.io.FileOutputStream
import java.io.InputStream
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.File

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

import sys.process._; 


object generatePlotsWithR {
  
class genSimplePlots extends CommandLineRunUtil {
     override def priority = 50;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "generateSamplePlots", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This simple function invokes R and generates a simple, single-replicate plots (or a similar pdf report) given a single replicate's QoRTs QC output."+
                        ""+
                        ""+
                        "",
          argList = 
                    new UnaryArgument(    name = "makePdf",
                                         arg = List("--makePdf"), // name of value
                                         argDesc = "Flag to indicate that you want a pdf multi-plot report to be generated." // description
                                       ) ::
                    new UnaryArgument(    name = "noPng",
                                         arg = List("--noPng"), // name of value
                                         argDesc = "Flag to indicate that you do NOT want the primary single-png multi-plot to be generated." // description
                                       ) ::
                    new UnaryArgument(    name = "makeSeparatePngs",
                                         arg = List("--makeSeparatePngs"), // name of value
                                         argDesc = "Flag to indicate that you want a battery of separate pngs to be generated." // description
                                       ) ::
                    new BinaryArgument[String](   name = "uniqueID",
                                                        arg = List("--uniqueID"),  
                                                        valueName = "id", 
                                                        argDesc = "The ID of the replicate. This will be only used for the plot labels.", 
                                                        defaultValue = Some("Untitled")
                                                        ) :: 
                    new FinalArgument[String](
                                         name = "qcdataDir",
                                         valueName = "qcDataDir",
                                         argDesc = "The qc directory in which all the QC files are contained." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS 
                   );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       
       if(out){
         generateSimplePlots(
             parser.get[String]("qcDir"),
             parser.get[String]("uniqueID"),
             ! parser.get[Boolean]("noPng"),
             parser.get[Boolean]("makePdf"),
             parser.get[Boolean]("makeSeparatePngs")
           );
         }
     }
   }
  
  def generateSimplePlots(qcdir : String, uniqueID : String, makePng : Boolean, makePdf : Boolean, makeSeparatePngs : Boolean) {
    //val qcdir = outfileprefix.substring(0,outfileprefix.length - 2)
    
    if(makeSeparatePngs){
       val outDir = new File(qcdir +"/" + "QCplots/");
       if(! outDir.exists()){
         reportln("Creating Directory: "+ outDir,"note");
         outDir.mkdir();
         reportln("Successfully Created Directory: " + outDir, "note");
       }
    }
    
    val rscriptString = makeRscriptString(qcdir, uniqueID, makePng,makePdf, makeSeparatePngs);
    val writer = openWriter(qcdir + "/QC.makeMultiplot.R");
    writer.write(rscriptString);
    writer.close();
    
    val rlogstring = new StringBuilder("");
    
    
    val rlogger = ProcessLogger((line : String) => {
      reportln("   > "+line,"debug");
      rlogstring.append(line + "\n");
    }, (line : String) => {
      reportln("   > " + line,"debug");
      rlogstring.append(line + "\n");
    });
    
    reportln(" > Starting R execution. (generating plots)","note");
    try{
      val exitcode = scala.sys.process.Process("Rscript --no-save --no-restore " + qcdir + "QC.makeMultiplot.R").!(rlogger);
    } catch {
      case e : Exception => {
        reportln("##########################","warn");
        reportln("Rscript execution FAILED.","warn");
        reportln("   Check environment variables, etc?\n   is the correct version of R in the PATH?\n   Is the QoRTs package installed?","warn");
        reportln("##########################","warn");
      }
    }
    reportln(" > Finished R execution.","note");
    
    //val rlogstring = ("Rscript --no-save --no-restore " + qcdir + "QC.makeMultiplot.R").!!
    
    val logwriter = openWriter(qcdir + "/QC.makeMultiplot.Rlog");
    logwriter.write(rlogstring.toString);
    logwriter.close();
  }
  
  def makeRscriptString(qcdir : String, uniqueID : String, makePng : Boolean, makePdf : Boolean, makeSeparatePngs : Boolean) : String = {
    "# This is an automatically-generated R script designed to make a simple multiplot and/or pdf report for a sample.\n"+
    "message(\"STARTING...\");\n"+
    "library(QoRTs);\n"+
    "unique.ID <- c(\"" + uniqueID + "\");\n"+
    "qc.data.dir <- c(\"" + qcdir + "/\");\n"+
    "decoder.raw <- data.frame(unique.ID = as.character(unique.ID), qc.data.dir = as.character(qc.data.dir));\n"+
    "decoder <- completeAndCheckDecoder(decoder = decoder.raw)\n"+
    "message(decoder);\n"+
    "message(lapply(names(decoder), function(n){ class(decoder[[n]]) }));\n"+
    "res <- read.qc.results.data(\"\", decoder = decoder, calc.DESeq2 = FALSE, calc.edgeR = FALSE);\n"+
    (if(makePng)          "makeMultiPlot.basic(res, outfile = \""+qcdir+"/QC.multiPlot.png\", plotter.params = list(std.color = \"blue\", std.lines.lwd = 4), plot.device.name = \"png\");\n" else "\n")+
    (if(makePdf)          "makeMultiPlot.basic(res, outfile = \""+qcdir+"/QC.multiPlot.pdf\", plotter.params = list(std.color = \"blue\", std.lines.lwd = 4), plot.device.name = \"pdf\");\n" else "\n")+
    (if(makeSeparatePngs) "makeMultiPlot.basic(res, outfile.dir = paste0(\""+qcdir+"\",\"QCplots/\"), plotter.params = list(std.color = \"blue\", std.lines.lwd = 4), plot.device.name = \"png\", separatePlots = TRUE);\n" else "\n")+
    "message(\"DONE...\");\n"+
    "\n"+
    "\n"+
    "\n";
  }

  /*

#This is an automatically-generated R script designed to make a simple multiplot and/or pdf report for a sample.
message("STARTING...");
library(QoRTs);
unique.ID <- c("SAMP1_RG1");
qc.data.dir <- c(".//");
decoder.raw <- data.frame(unique.ID = as.character(unique.ID), qc.data.dir = as.character(qc.data.dir));
decoder <- completeAndCheckDecoder(decoder = decoder.raw)
message(decoder);
message(lapply(names(decoder), function(n){ class(decoder[[n]]) }));
res <- read.qc.results.data("", decoder = decoder, calc.DESeq2 = FALSE, calc.edgeR = FALSE);
makeMultiPlot.basic(res, outfile = "QC.multiPlot.png", plotter.params = list(std.color = "blue", alt.color = "blue"), plot.device.name = "png");
message("DONE...")

   */
  
  
}





