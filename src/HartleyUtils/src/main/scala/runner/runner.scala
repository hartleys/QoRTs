package runner

//import fileConversionUtils._
//import annotationUtils._
//import internalUtils._
//import miscLittleJobUtils._

import internalUtils.commandLineUI._;

object runner {
  
  final val QORTS_VERSION = "0.0.18";
  //final val FOR_HELP_STRING = "For help, use command: "
  
  final val Runner_ThisProgramsExecutableJarFileName : String = "QoRT.jar";
  final val allowDepreciated : Boolean = true;
  final val COMMAND_MAX_LENGTH = 20;
  
  //Command name -> (execution call, summary, syntax)
  final val utilCommandList : Map[String, () => CommandLineRunUtil] = 
      Map( //NOTE: All commands MUST be of length < COMMAND_MAX_LENGTH!
           ("QC" -> (() => new qcUtils.runAllQC.allQC_runner)),
           ("makeFlatGtf" -> (() => new fileConversionUtils.prepFlatGtfFile.prepFlatGtfFile_runner)),
           ("mergeWig" ->  (() => new fileConversionUtils.SumWigglesFast.SumWigglesFast_runner)),
           ("mergeAllCounts" ->  (() => new fileConversionUtils.mergeQcOutput.multiMerger)),
           ("mergeCounts" ->  (() => new fileConversionUtils.mergeQcOutput.merger)),
           ("bamToWiggle" ->  (() => new fileConversionUtils.bamToWiggle.wiggleMaker)),
           ("makeSpliceBed" ->  (() => new fileConversionUtils.convertSpliceCountsToBed.converter)),
           ("mergeNovelSplices" -> (() => new fileConversionUtils.addNovelSplices.mergeNovelSplices)),
           ("makeSpliceJunctionBed" -> (() => new fileConversionUtils.makeSpliceJunctionBed.converter))

          //(("prepFlatGtfFile",((fileConversionUtils.prepFlatGtfFile.run(_), "", "")))),
           //(("QC", ((qcUtils.runAllQC.run(_)),"",""))),
           //(("convertSoftToHardClipping", ((fileConversionUtils.convertSoftToHardClipping.run(_)),"",""))),
           //(("bamToWiggle", ((fileConversionUtils.bamToWiggle.run(_)),"",""))),
           //(("sumWiggles", ((fileConversionUtils.SumWigglesFast.run(_)),"","")))
         )
  final val helpCommandList : Map[String, () => CommandLineRunUtil] = 
    Map(
           ("?" -> (() => new helpDocs())),
           ("help" -> (() => new helpDocs())),
           ("man" -> (() => new helpDocs())),
           ("--help" -> (() => new helpDocs())),
           ("--man" -> (() => new helpDocs())),
           ("-help" -> (() => new helpDocs())),
           ("-man" -> (() => new helpDocs()))
       );
  final val commandList = utilCommandList ++ helpCommandList;
         
  final val depreciated_commandList : Map[String, (Array[String]) => Unit] = 
     Map(
         (("bamToWiggleDepreciated", ((fileConversionUtils.bamToWiggleDepreciated.run(_)))  ))
         );

  def main(args: Array[String]){
    //println("Initializing...");
    
    internalUtils.Reporter.init_stderrOnly;
    internalUtils.Reporter.reportln("Starting QoRTs v"+QORTS_VERSION+" (" + (new java.util.Date()).toString + ")","debug");
    
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
            case Some(c) => c(args);
            case None => {
              internalUtils.Reporter.reportln("[runner.runner Error]: Command " + args(0) + " not found!","output");
              helpDocs.generalHelp;
            }
          }
        }
      }
    }}
   // helloWorld.run(args);
    internalUtils.Reporter.reportln("Done. (" + (new java.util.Date()).toString + ")","note");
    internalUtils.Reporter.closeLogs;
  }
  
}
