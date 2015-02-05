package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;
import internalUtils.GtfTool._;

object makeBedFromFlatGff {

  class converter extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "makeSpliceBedFromGff", 
          quickSynopsis = "", 
          synopsis = "", 
          description = ""+
                        ""+
                        ""+
                        ""+
                        ""+
                        ""+
                        ""+
                        ""+
                        "",   
          argList = 
                    new BinaryOptionArgument[String](
                                         name = "rgb", 
                                         arg = List("--rgb"), 
                                         valueName = "r,g,b",  
                                         argDesc = "The rgb color for all the bed file lines."
                                        ) ::
                    new UnaryArgument(name = "flattenedGff",
                              arg = List("--flattenedGff"), // name of value
                              argDesc = "." // description
                              ) :: 
                    new UnaryArgument(name = "stranded",
                              arg = List("--stranded"), // name of value
                              argDesc = "Flag to indicate that data is stranded." // description
                              ) :: 
                    new FinalArgument[String](
                                         name = "gffFile",
                                         valueName = "annotation.gff.gz",
                                         argDesc = "The gff or gtf file" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output bed file. If the filename ends with \".gz\" or \".zip\" then the file will be compressed using the appropriate method." // description
                                        ) :: List() );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         makeBedFromFlatGff.run(
             parser.get[String]("gffFile"),
             parser.get[String]("outfile"),
             parser.get[Boolean]("stranded"),
             parser.get[Boolean]("flattenedGff"),
             parser.get[Option[String]]("rgb")
           );
         }
     }
   }
  
  def run(gffFile : String, outfile : String, stranded : Boolean, flattened : Boolean, rgb : Option[String]) {
    val gffLines : Iterator[FlatGtfLine] = if(flattened){
      GtfReader.getFlatGtfReader(gffFile, stranded, true, "\\s+", new GtfCodes());
    } else {
      fileConversionUtils.prepFlatGtfFile.getFlatGtfLines(gffFile,stranded).iterator;
    }
    
    //INCOMPLETE!
  }
}

















