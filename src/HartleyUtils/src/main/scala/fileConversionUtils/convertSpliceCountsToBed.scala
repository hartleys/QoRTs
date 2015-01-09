package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

object convertSpliceCountsToBed {

  class converter extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "makeSpliceBed", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This utility generates a splice-junction 'bed' file from the QoRTs-generated "+
                        "splice junction counts produced by the QC utility. This splice junction bed file "+
                        "can be used to visualize splice junction counts using the UCSC genome browser "+
                        "and other similar utilities."+
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
                    new BinaryOptionArgument[Double](
                                         name = "sizeFactor", 
                                         arg = List("--sizeFactor"), 
                                         valueName = "val",  
                                         argDesc = "A double-precision floating-point value. If this option is set, all counts will be divided by the given normalization factor."
                                        ) ::
                    new UnaryArgument(name = "includeSpliceNames",
                              arg = List("--includeSpliceNames"), // name of value
                              argDesc = "Flag to indicate that splice names should be used as well as splice counts." // description
                              ) :: 
                    new BinaryArgument[Int](   name = "digits",
                                                        arg = List("--digits"),  
                                                        valueName = "num", 
                                                        argDesc = "The number of digits after the decimal to include in counts. CURRENTLY NOT IMPLEMENTED!", 
                                                        defaultValue = Some(1)
                                                        ) ::
                    new BinaryOptionArgument[Double](
                                         name = "filterMin", 
                                         arg = List("--filterMin"), 
                                         valueName = "val",  
                                         argDesc = "If this option is set, then all bed lines with a count LESS THAN the given value will be dropped."
                                        ) ::
                    new BinaryOptionArgument[Double](
                                         name = "filterMax", 
                                         arg = List("--filterMax"), 
                                         valueName = "val",  
                                         argDesc = "If this option is set, then all bed lines with a count GREATER THAN the given value will be dropped."
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "infile",
                                         argDesc = "The input splice counts file, or '-' to read from stdin.  If the filename ends with \".gz\" or \".zip\" then the file will be read using the appropriate decompression method." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output bed file, or '-' to write to stdout. If the filename ends with \".gz\" or \".zip\" then the file will be compressed using the appropriate method." // description
                                        ) :: List() );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         convertSpliceCountsToBed.run(
             parser.get[String]("infile"),
             parser.get[String]("outfile"),
             parser.get[Option[Double]]("sizeFactor"),
             parser.get[Option[Double]]("filterMin"),
             parser.get[Option[Double]]("filterMax"),
             parser.get[Option[String]]("rgb"),
             parser.get[Int]("digits"),
             parser.get[Boolean]("includeSpliceNames")
           );
         }
     }
   }
  
  def run(infile : String, outfile : String, sizeFactor : Option[Double], filterMin : Option[Double], filterMax : Option[Double], opt_rgb : Option[String], digits : Int, includeSpliceNames : Boolean) {
    val rgb = if(opt_rgb.isEmpty) "255,255,255" else opt_rgb.get;
    def scoreFunction(ct : Int) : Int = ct;
    def countFilter(ct : Int) : Boolean = {
      if((! filterMin.isEmpty) && filterMin.get > ct) false;
      else if((! filterMax.isEmpty) && filterMax.get < ct) false;
      else true;
    }
    
    convert(infile,outfile, rgb, countFilter, scoreFunction, sizeFactor, includeSpliceNames);
  }
  
  def convert(infile : String, outfile : String, rgb : String, countFilter : (Int => Boolean), scoreFunction : (Int => Int), sizeFactor : Option[Double], includeSpliceNames : Boolean, delim : String = "	") {
    val lines : Iterator[String] = getLinesSmartUnzip(infile, true);
    val writer : WriterUtil = openWriterSmart(outfile, true);
    
    val titleCells = lines.next.split(delim);
    
    val chromCol = titleCells.indexOf("chrom");
    val strandCol = titleCells.indexOf("strand");
    val startCol = titleCells.indexOf("start");
    val endCol = titleCells.indexOf("end");
    val ctCol = titleCells.indexOf("CT");
    val nameCol = titleCells.indexOf("spliceName");
    
    if(chromCol == -1) error("ERROR: No \"chrom\" column found!");
    if(strandCol == -1) error("ERROR: No \"strand\" column found!");
    if(startCol == -1) error("ERROR: No \"start\" column found!");
    if(endCol == -1) error("ERROR: No \"end\" column found!");
    if(ctCol == -1) error("ERROR: No \"CT\" column found!");
    if(includeSpliceNames && nameCol == -1) error("ERROR: option --includeSpliceNames is TRUE, but no \"spliceName\" column found!");
    
    for(line <- lines){
      val cells : Array[String] = line.split(delim);
      val chrom : String = cells(chromCol);
      val strand : String = cells(strandCol);
      val start : Int = string2int(cells(startCol));
      val end : Int = string2int(cells(endCol));
      val ct : Int = string2int(cells(ctCol));
      
      if(countFilter(ct)){
        val score : Int = math.max(0,math.min(scoreFunction(ct),1000));
        val len : Int = end - start;
        val name : String = (if(includeSpliceNames) cells(nameCol) else "") +
          ( sizeFactor match {
            case Some(sf) => ct.toDouble * sf;
            case None => ct;
          }).toString;
        //NOTE TO SELF: CHECK FOR OFF-BY-ONE ERRORS!
        writer.write(chrom+"	"+start+"	"+end+"	"+name+"	"+score+"	"+strand+"	"+start+"	"+end+"	"+rgb+"	2	1,1	0,"+(len - 1)+"\n");
      } //else skip this line!
    }
    close(writer);
  }
  
  
  
}