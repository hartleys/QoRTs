package fileConversionUtils


import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.File;
import internalUtils.stdUtils._;
import internalUtils.Reporter._;
import internalUtils.commandLineUI._;
import internalUtils.GtfTool._;

object gtfConverter {
  //WORK IN PROGRESS
  class GTF_to_GFF3 extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "GTF_to_GFF3", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This simple utility converts a gtf transcript annotation file into a .bed transcript annotation file."+
                        "Note that this conversion may theoretically be \"lossy\", as it is possible for gtf files to contain information that "+
                        "cannot be perfectly translated into the bed format."+
                        ""+
                        "WARNING: THIS SUB-UTILITY IS BETA! NOT FOR GENERAL USE!"+
                        ""+
                        ""+
                        ""+
                        "",   
          argList = 
                    new FinalArgument[String](
                                         name = "gtffile",
                                         valueName = "annotation.gtf.gz",
                                         argDesc = "The gtf file, or '-' to read from stdin. If the filename ends with \".gz\" or \".zip\" then the file will be decompressed using the appropriate method." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output bed file. If the filename ends with \".gz\" or \".zip\" then the file will be compressed using the appropriate method." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         gtfConverter.GTF_to_GFF3_RUN(
             parser.get[String]("gtffile"),
             parser.get[String]("outfile")
           );
         }
     }
   }
  
 // def test() {
  def GTF_to_GFF3_RUN(infile : String, outfile : String){
    def gtflines : Iterator[StdGtfLine] = GtfReader.getStdGtfReader(infile, true, true, "\\s+", internalUtils.GtfTool.GtfCodes());
    for(gtfline <- gtflines){
      //INCOMPLETE!!!!
      //OutputGtfLine(in_chromName : String, in_featureSource : String, in_featureType : String, in_start : Int, in_end : Int, in_score : String, in_strand : Char, in_attributeMap : Map[String,String], in_gtfFmt_attributeBreak : String, in_stranded : Boolean, in_codes : GtfCodes = new GtfCodes())
    }
  }
  
  
  
  

}
