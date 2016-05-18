package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;
import java.io.File;


object convertSpliceCountsToBed {

  final val junctionTypeFileMap : Map[String,String] = 
      Map(
         ("known" -> "QC.spliceJunctionCounts.knownSplices.txt"),
         ("novel" -> "QC.spliceJunctionCounts.novelSplices.txt"),
         ("orphan" -> "QC.spliceJunctionCounts.orphanSplices.txt")
      )
  final val junctionTypeList : List[String] = junctionTypeFileMap.keySet.toList.sorted;
  final val junctionTypeOutFileMap : Map[String,String] = 
      Map(
         ("known" -> "QC.spliceJunctionCounts.knownSplices.bed"),
         ("novel" -> "QC.spliceJunctionCounts.novelSplices.bed"),
         ("orphan" -> "QC.spliceJunctionCounts.orphanSplices.bed")
      )
   /*   
    val chromCol = titleCells.indexOf("chrom");
    val strandCol = titleCells.indexOf("strand");
    val startCol = titleCells.indexOf("start");
    val endCol = titleCells.indexOf("end");
    val ctCol = titleCells.indexOf("CT");
    val nameCol = titleCells.indexOf("spliceName");
      */
      
  final val junctionTypeIdxMap : Map[String,Map[String,Int]] = 
      Map(
         ("known" ->       Map("chrom" -> 1, "strand" -> 2, "start" -> 3, "end" -> 4, "CT" -> 5, "spliceName" -> 0)),
         ("novel" ->       Map("chrom" -> 0, "strand" -> 1, "start" -> 2, "end" -> 3, "CT" -> 4, "spliceName" -> -1)),
         ("orphan" ->      Map("chrom" -> 0, "strand" -> 1, "start" -> 2, "end" -> 3, "CT" -> 4, "spliceName" -> 5))
      );
  final val junctionTypeTLMAP : Map[String,Boolean] = 
      Map(
         ("known" -> false),
         ("novel" -> true),
         ("orphan" -> true)
      );
  
  class converter2 extends CommandLineRunUtil {
     override def priority = 100;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "makeSimpleJunctionTrack", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This utility converts QoRTs splice junction count files into "+
                        "bed format. Unlike makeJunctionTrack and makeOrphanJunctionTrack, "+
                        "this utility is not designed to compile multiple samples or replicates "+
                        "together. It is a simple converter from a QoRTs junction count file "+
                        "to a bed file. The count files come in 3 types: known, novel, and orphan.",
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
                    new UnaryArgument( name = "noGzip",
                                         arg = List("--noGzip"), // name of value
                                         argDesc = "Flag to indicate whether whether input and output data is/will be gzip-compressed." // description
                                       ) :: 
                    new BinaryArgument[Int](   name = "maxIdentifierLength",
                                                        arg = List("--maxIdentifierLength"),  
                                                        valueName = "num", 
                                                        argDesc = "The max number of characters a junction ID can have. If longer than this, the middle will be truncated. This is to prevent browser errors when the length of a bed line name is greater than 255 characters.", 
                                                        defaultValue = Some(45)
                                                        ) ::
                    new BinaryOptionArgument[String]( 
                                         name = "outfileSuffix",
                                         arg = List("--outfileSuffix"), // name of value
                                         valueName = "file.bed.gz",
                                         argDesc = "The name of the output file." // description
                                       ) ::
                    new FinalArgument[String](
                                         name = "indir",
                                         valueName = "indir",
                                         argDesc = "The data directory in which to find the splice junction count files." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "fileType",
                                         valueName = "fileType",
                                         argDesc = "The type of splice junction counts to compile. Must be one of: ["+junctionTypeList.mkString(",")+"]." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outdir",
                                         valueName = "outdir",
                                         argDesc = "The location to output the bed file. Traditionally the same as indir." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         convertSpliceCountsToBed.run2(
             parser.get[String]("indir"),
             parser.get[String]("outdir"),
             parser.get[String]("fileType"),
             parser.get[Option[Double]]("sizeFactor"),
             parser.get[Option[Double]]("filterMin"),
             parser.get[Option[Double]]("filterMax"),
             parser.get[Option[String]]("rgb"),
             parser.get[Int]("digits"),
             parser.get[Boolean]("includeSpliceNames"),
             parser.get[Boolean]("noGzip"),
             parser.get[Option[String]]("outfileSuffix"),
             parser.get[Int]("maxIdentifierLength")
           );
         }
     }
   }
  
  def run2(indir : String, outdir : String, fileType : String, sizeFactor : Option[Double], filterMin : Option[Double], filterMax : Option[Double], opt_rgb : Option[String], digits : Int, includeSpliceNames : Boolean, noGzip : Boolean, outfileSuffix : Option[String], maxIdentifierLength : Int) {
    val rgb = if(opt_rgb.isEmpty) "255,255,255" else opt_rgb.get;
    def scoreFunction(ct : Int) : Int = ct;
    def countFilter(ct : Int) : Boolean = {
      if((! filterMin.isEmpty) && filterMin.get > ct) false;
      else if((! filterMax.isEmpty) && filterMax.get < ct) false;
      else true;
    }
    
    if(! junctionTypeFileMap.contains(fileType)){
      error("fileType \""+fileType+"\" not recognized. Must be one of: ["+junctionTypeList.mkString(",")+"].")
    }
    val inDir = new File(indir);
    if(! inDir.exists()){
      error("Input directory \""+indir+"\" does not exist!");
    }
    
    val infile = indir +"/"+ junctionTypeFileMap(fileType) + {if(! noGzip) ".gz" };
    val inFile = new File(infile);
    if(! inFile.exists()){
      error("Input file \""+infile+"\" does not exist!");
    }
    
    val outDir = new File(outdir); 
    if(! outDir.exists()){
       reportln("Creating Directory: "+ outDir,"note");
       outDir.mkdir();
       reportln("Successfully Created Directory: " + outDir, "note");
     }
    
    val outfile = if(outfileSuffix.isEmpty){
      outdir +"/" + junctionTypeOutFileMap(fileType) + {if(! noGzip) ".gz" };
    } else {
      outdir + "/" + outfileSuffix.get;
    }
    
    val idxmap = junctionTypeIdxMap(fileType);
    val TLmap = junctionTypeTLMAP(fileType);
    
    reportln("Reading from file \""+infile+"\"","debug");
    reportln("Writing to file \""+outfile+"\"","debug");

    convert2(infile=infile,
          outfile=outfile, 
          rgb=rgb, 
          countFilter=countFilter, 
          scoreFunction=scoreFunction, 
          sizeFactor=sizeFactor, 
          includeSpliceNames=includeSpliceNames,
          idxmap = Some(idxmap),
          TL = TLmap,
          maxIdentifierLength = maxIdentifierLength
          );
  }
  
  def convert2(infile : String, outfile : String, rgb : String, countFilter : (Int => Boolean), scoreFunction : (Int => Int), sizeFactor : Option[Double], includeSpliceNames : Boolean, 
              idxmap : Option[Map[String,Int]] = None,
              TL : Boolean = false,
              maxIdentifierLength : Int = 45,
              delim : String = "\t") {
    
    val (titleLine,lines) : (String, Iterator[String]) = peekIterator(getLinesSmartUnzip(infile, true));
    val titleCells = titleLine.split(delim);

    var chromCol = -1;
    var strandCol = -1;
    var startCol = -1;
    var endCol = -1;
    var ctCol = -1;
    var nameCol = -1;
    
    if(! idxmap.isEmpty){
      val idxMap = idxmap.get;
      chromCol = idxMap("chrom");
      strandCol = idxMap("strand");
      startCol = idxMap("start");
      endCol = idxMap("end");
      ctCol = idxMap("CT");
      nameCol = idxMap("spliceName");
      
      if(titleCells(startCol) == "start" && titleCells(endCol) == "end"){
        val buffer = lines.next;
      }
    }
    
    if(idxmap.isEmpty){
        chromCol = titleCells.indexOf("chrom");
        strandCol = titleCells.indexOf("strand");
        startCol = titleCells.indexOf("start");
        endCol = titleCells.indexOf("end");
        ctCol = titleCells.indexOf("CT");
        nameCol = titleCells.indexOf("spliceName");
        val buffer = lines.next;
    }

    if(chromCol == -1) error("ERROR: No \"chrom\" column found!");
    if(strandCol == -1) error("ERROR: No \"strand\" column found!");
    if(startCol == -1) error("ERROR: No \"start\" column found!");
    if(endCol == -1) error("ERROR: No \"end\" column found!");
    if(ctCol == -1) error("ERROR: No \"CT\" column found!");
    if(includeSpliceNames && nameCol == -1) error("ERROR: option --includeSpliceNames is TRUE, but no \"spliceName\" column found!");

    val writer : WriterUtil = openWriterSmart(outfile, true);
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
        val id : String = if(! includeSpliceNames) "" else {
          val rawID = cells(nameCol);
          if(rawID.length() > maxIdentifierLength){
            rawID.substring( 0,(maxIdentifierLength/2)-3 ) + "..." + rawID.substring(rawID.length - (maxIdentifierLength/2), rawID.length())
          } else {
            rawID;
          }
        }
        
        val name : String = (if(includeSpliceNames) cells(nameCol) else "") + "(" +
          ( sizeFactor match {
            case Some(sf) => ct.toDouble / sf;
            case None => ct;
          }).toString + ")";
        //NOTE TO SELF: CHECK FOR OFF-BY-ONE ERRORS!
        writer.write(chrom+"\t"+(start-1)+"\t"+(end+1)+"\t"+name+"\t"+score+"\t"+strand+"\t"+(start-1)+"\t"+(end+1)+"\t"+rgb+"\t2\t1,1\t0,"+(len+1)+"\n");
      } //else skip this line!
    }
    close(writer);
  }
  
  
  /*
  class converter extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "convertSpliceCountsToBed", 
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
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
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
            case Some(sf) => ct.toDouble / sf;
            case None => ct;
          }).toString;
        //NOTE TO SELF: CHECK FOR OFF-BY-ONE ERRORS!
        writer.write(chrom+"	"+start+"	"+end+"	"+name+"	"+score+"	"+strand+"	"+start+"	"+end+"	"+rgb+"	2	1,1	0,"+(len - 1)+"\n");
      } //else skip this line!
    }
    close(writer);
  }
  
  */
  
}