package runner


import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

import internalUtils.commonSeqUtils._;
import internalUtils.genomicUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.JavaConversions._

object helpDocs {
  
  val MIT_LICENSE = List[String](
 "The MIT License",
 "Copyright (c) 2009 The Broad Institute",
 "Permission is hereby granted, free of charge, to any person obtaining a copy "+
 "of this software and associated documentation files (the \"Software\"), to deal "+
 "in the Software without restriction, including without limitation the rights "+
 "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell "+
 "copies of the Software, and to permit persons to whom the Software is "+
 "furnished to do so, subject to the following conditions: ",
 
 "The above copyright notice and this permission notice shall be included in "+
 "all copies or substantial portions of the Software. ",
 
 "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "+
 "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, "+
 "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE "+
 "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER "+
 "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, "+
 "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN "+
 "THE SOFTWARE."
  );
  
  val AUTHOR = internalUtils.commandLineUI.DEFAULT_AUTHOR;
  val LEGAL = internalUtils.commandLineUI.DEFAULT_LEGAL;
  
  val DESCRIPTION = List[String](
      "The first module of the software package QoRT, which is intended for use with Paired-End, High-Throughput RNA-Seq data.",
      "This tool can perform a number of different functions to assist in assessing the data quality, detecting errors or biases, "+
      "performing analyses, data cleaning, data visualization, and data formatting."
  );
  
  //val DESCRIPTION = List[String](
  //    "<TODO>: Write description!"
  //);
  
  def runHelp(args : Array[String]){
    //report("Help:","output");
    if(args.length <= 1 || args(1) == "?" || args(1) == "--man" || args(1) == "--help" || args(1) == "help" || args(1) == "-help" || args(1) == "man" || args(1) == "-man"){
      generalHelp;
    } else {
      //Print help for a specific function:
      
      val commandList : Map[String, () => CommandLineRunUtil] = runner.commandList;
      
      val helpCommand = args(1);
      val cmd = commandList.get(helpCommand);
      
      report("HELP: " + helpCommand +"\n","output");
      
      if(helpCommand == "samjdkinfo"){
        val fmt_mit_license = lineseq2string(wrapLinesWithIndent(MIT_LICENSE, internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH, "#    ", false ));
        reportln("The SAM jdk from Picard is licensed under the MIT License:\n" + fmt_mit_license, "output");
      } else {
        cmd match {
          case Some(makerFcn) => {
            val cmdRunner = makerFcn();
            cmdRunner.parser.reportManual();
          }
          case None => {
            val cmdOld = runner.depreciated_commandList.get(helpCommand);
            cmdOld match {
              case Some(c) => reportln("Command " + helpCommand +" is DEPRECIATED. No help info found!","output");
              case None => reportln("Command " + helpCommand + " not found!","output");
            }
          }
        }
      }
    }
  }
  def generalHelp {
    //Print general help.
    reportln("GENERAL HELP:","output");
    
    reportln(getGeneralHelp,"output");
  }
   
  def getGeneralHelp : String = {
      val sb = new StringBuilder("");
      sb.append("QoRTs: Quality Of Rna-seq Tool Set\n")
      sb.append("version: " + runner.QORTS_VERSION + "\n");
      sb.append("\n");
      //sb.append("SYNOPSIS\n");
      //sb.append( lineseq2string(wrapLinesWithIndent(SYNOPSIS, internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH , "    ", false))  + "\n");
      //sb.append("\n");
      sb.append("DESCRIPTION:\n");
      sb.append( lineseq2string(wrapLinesWithIndent(DESCRIPTION, internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH, "    ", false)) + "\n");
      sb.append("    NOTE: if you run into OutOfMemoryExceptions, \n    try adding the java options: \"-Xmx18000M -Xms5000M\""+"\n");
      sb.append("\n");
      
      sb.append("GENERAL SYNTAX:\n");
      sb.append("    java [_java_options_] -jar "+runner.Runner_ThisProgramsExecutableJarFileName +" _COMMAND_ [_options_]"+"\n");
      sb.append("\n");
      
      sb.append("COMMANDS:\n");
      for((arg, cmdMaker) <- runner.utilCommandList){
        
        val parser = cmdMaker().parser;
        //Note: Hack! the line below must be at least as long as COMMAND_MAX_LENGTH!
        val blank = "                                                                                                                                                                           ";
        sb.append("    "+arg+":" + blank.substring(0,runner.COMMAND_MAX_LENGTH - arg.length) + parser.getQuickSynopsis + "\n");
        sb.append(wrapLinesWithIndent(parser.getDescription, internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH, "        ", false) + "\n");
        sb.append(wrapLinesWithIndent(parser.getForMoreHelp, internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH, "        ", false) + "\n");
      }
      sb.append("AUTHORS:\n");
      sb.append(lineseq2string(wrapLinesWithIndent(AUTHOR, internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH, "    ", false)) + "\n");
      
      sb.append("LEGAL:\n");
      sb.append(lineseq2string(wrapLinesWithIndent(LEGAL, internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH, "    ", false)) + "\n");
      
      return sb.toString;
  }
}

class helpDocs extends CommandLineRunUtil {
    val parser : CommandLineArgParser = 
      new CommandLineArgParser(
          command = "", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "",   
          argList =  List()
      );
    
    def run(args : Array[String]){
      helpDocs.runHelp(args);
    }
    
    def generalHelp {
      helpDocs.generalHelp;
    }
 }