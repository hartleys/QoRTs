package runner

//import fileConversionUtils._
//import annotationUtils._
//import internalUtils._
//import miscLittleJobUtils._

import internalUtils.commandLineUI._;

object runner {
  
  final val QORTS_VERSION = "1.0.1"; // REPLACE_THIS_QORTS_VERSION_VARIABLE_WITH_VERSION_NUMBER          (note this exact text is used in a search-and-replace. Do not change it.)
  final val QORTS_COMPILE_DATE = "Tue Nov 10 10:05:13 EST 2015"; // REPLACE_THIS_QORTS_DATE_VARIABLE_WITH_DATE          (note this exact text is used in a search-and-replace. Do not change it.)
  
  //final val FOR_HELP_STRING = "For help, use command: "
  
  final val Runner_ThisProgramsExecutableJarFileName : String = "QoRTs.jar";
  final val allowDepreciated : Boolean = true;
  final val COMMAND_MAX_LENGTH = 30;
  
  //Command name -> (execution call, summary, syntax)
  final val utilCommandList : Map[String, () => CommandLineRunUtil] = 
      Map( //NOTE: All commands MUST be of length < COMMAND_MAX_LENGTH!
           ("QC" -> (() => new qcUtils.runAllQC.allQC_runner)),
           ("makeFlatGff" -> (() => new fileConversionUtils.prepFlatGtfFile.prepFlatGtfFile_runner)),
           ("mergeWig" ->  (() => new fileConversionUtils.SumWigglesFast.SumWigglesFast_runner)),
           ("mergeAllCounts" ->  (() => new fileConversionUtils.mergeQcOutput.multiMerger)),
           ("mergeCounts" ->  (() => new fileConversionUtils.mergeQcOutput.merger)),
           ("bamToWiggle" ->  (() => new fileConversionUtils.bamToWiggle.wiggleMaker)),
           //("makeSpliceBed" ->  (() => new fileConversionUtils.convertSpliceCountsToBed.converter)),
           ("mergeNovelSplices" -> (() => new fileConversionUtils.addNovelSplices.mergeNovelSplices)),
           ("makeJunctionTrack" -> (() => new fileConversionUtils.makeSpliceJunctionBed.converter)),
           ("generateSamplePlots" -> (() => new fileConversionUtils.generatePlotsWithR.genSimplePlots)),
           ("makeOrphanJunctionTrack" -> (() => new fileConversionUtils.makeOrphanJunctionBed.converter))
            
           //(("prepFlatGtfFile",((fileConversionUtils.prepFlatGtfFile.run(_), "", "")))),
           //(("QC", ((qcUtils.runAllQC.run(_)),"",""))),
           //(("convertSoftToHardClipping", ((fileConversionUtils.convertSoftToHardClipping.run(_)),"",""))),
           //(("bamToWiggle", ((fileConversionUtils.bamToWiggle.run(_)),"",""))),
           //(("sumWiggles", ((fileConversionUtils.SumWigglesFast.run(_)),"","")))
         )
  final val sortedCommandList : Seq[(String, () => CommandLineRunUtil)] = utilCommandList.toVector.sortBy(_._2().priority);
         
         
  final val helpCommandList : Map[String, () => CommandLineRunUtil] = 
    (internalUtils.commandLineUI.HELP_COMMAND_LIST ++ internalUtils.commandLineUI.MANUAL_COMMAND_LIST).map( (c : String) => {
      (c, (() => new helpDocs()));
    }).toMap;
  
  final val commandList = utilCommandList ++ helpCommandList;
         
  final val depreciated_commandList : Map[String, () => CommandLineRunUtil] = 
     Map(
         ("makeBedFromGtf" -> (() => new fileConversionUtils.makeBedFromGtf.converter))
         );

  def main(args: Array[String]){
    //println("Initializing...");
    
    if(args.contains("--verbose")){
      internalUtils.Reporter.init_base();
      internalUtils.Reporter.init_stderrOnly(internalUtils.Reporter.VERBOSE_CONSOLE_VERBOSITY);
    } else if(args.contains("--quiet")){
      internalUtils.Reporter.init_base();
      internalUtils.Reporter.init_stderrOnly(internalUtils.Reporter.QUIET_CONSOLE_VERBOSITY);
    } else {
      internalUtils.Reporter.init_base();
      internalUtils.Reporter.init_stderrOnly();
    }
    
    internalUtils.Reporter.reportln("Starting QoRTs v"+QORTS_VERSION+" (Compiled " + QORTS_COMPILE_DATE + ")","note");
    internalUtils.Reporter.reportln("Starting time: ("+(new java.util.Date()).toString+")","note");
    
    if(System.getProperty("sun.arch.data.model") == "32"){
      internalUtils.Reporter.reportln(
               "> Warning: 32-bit Java detected! 32-bit JVMs generally have a built-in hard ceiling on their memory usage. \n"+
               ">          Usually between 1.5 and 3 gigabytes of RAM (The precise ceiling varies depending on the version, the hardware, and the operating system).\n"+
               ">          It is generally recommended that you install a 64-bit version of Java, if available.\n"+
               ">          This can be downloaded from \"www.java.com\".","warn")
    } 
    
    try{
	    if(args.length == 0){
	      internalUtils.Reporter.reportln("No command given!","output");
	      helpDocs.generalHelp;
	    } else {
		    val cmd = commandList.get(args(0));
		    cmd match {
		      case Some(makerFcn) => {
		        val cmdRunner = makerFcn();
		        cmdRunner.run(args);
		      }
		      case None => {
		        if(! allowDepreciated) {
		          internalUtils.Reporter.reportln("[runner.runner Error]: Command " + args(0) + " not found, and depreciated tools are deactivated!","output");
		          helpDocs.generalHelp;
		        } else {
		          val cmdOld = depreciated_commandList.get(args(0));
		          cmdOld match {
		            case Some(makerFcn) => {
		              internalUtils.Reporter.reportln("WARNING: Running Beta tool: " + args(0),"warn");
		              val cmdRunner = makerFcn();
		              cmdRunner.run(args);
		            }
		            case None => {
		              internalUtils.Reporter.reportln("[runner.runner Error]: Command " + args(0) + " not found!","output");
		              helpDocs.generalHelp;
		            }
		          }
		        }
		      }
		    }
	    }
    } catch {
      case e : Exception => {
        internalUtils.Reporter.reportln("============================FATAL_ERROR============================\n"+
                                        "QoRTs encountered a FATAL ERROR. For general help, use command:\n"+
                                        "          java -jar path/to/jar/QoRTs.jar --man\n"+
                                        "============================FATAL_ERROR============================\n"+
                                        "Error info:","note");
        throw e;
      }
    }
    
   // helloWorld.run(args);
    //} catch {
    //  case e : Exception => {
    //    internalUtils.Reporter.reportln("Error Caught. General Help:","note");
    //    helpDocs.generalHelp;
    //    throw e;
    //  }
    //}
    
    internalUtils.Reporter.reportln("Done. (" + (new java.util.Date()).toString + ")","note");
    internalUtils.Reporter.closeLogs;
  }
  
}
