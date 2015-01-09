package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

object mergeSpliceBed {
  class merger extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "mergeSpliceBed", 
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
                    new UnaryArgument(name = "quiet",
                              arg = List("--quiet","-q"), // name of value
                              argDesc = "" // description
                              ) :: 
                    new BinaryOptionArgument[String](
                                         name = "rgb", 
                                         arg = List("--rgb"), 
                                         valueName = "r,g,b",  
                                         argDesc = ""
                                        ) ::
                    new UnaryArgument(name = "calcMean",
                              arg = List("--calcMean","-m"), // name of value
                              argDesc = "" // description
                              ) ::
                    new BinaryOptionArgument[Double](
                                         name = "filterMin", 
                                         arg = List("--filterMin"), 
                                         valueName = "val",  
                                         argDesc = ""
                                        ) ::
                    new BinaryOptionArgument[Double](
                                         name = "filterMax", 
                                         arg = List("--filterMax"), 
                                         valueName = "val",  
                                         argDesc = ""
                                        ) ::
                    new BinaryOptionArgument[List[Double]](
                                         name = "sizeFactors", 
                                         arg = List("--sizeFactors"), 
                                         valueName = "sf,sf,sf,...",  
                                         argDesc = "normalization factors for each bed file."
                                        ) ::
            new UnaryArgument(name = "ignoreSizeFactors",
                              arg = List("--ignoreSizeFactors"), // name of value
                              argDesc = "Flag to indicate that this utility should ignore size factors even if they are found in the input listFile or manually included via the --sizeFactors option." // description
                              ) ::
            new FinalArgument[String](
                                         name = "filelist",
                                         valueName = "<filelist.txt | - | file1.bed,file2.bed,...>",
                                         argDesc = "One of three things:"+
                                                   "(1) A comma-delimited list of bed files. Note: NO WHITESPACE!"+
                                                   "(2) A file (or '-' to read from stdin) containing a list of the bed files to merge (one on each line).\n"+
                                                   "Optionally, the file can have a second (tab-delimited) column, containing the normalization factors to use. If this column is present, "+
                                                   "this utility will automatically calculate the normalized totals (or the normalized mean, if --calcMean is set).\n"+
                                                   "Note: wig filenames cannot contain whitespace or commas.\n"+
                                                   "Also Note: if the bed file names end in \".gz\" or \".zip\", they will be read using the appropriate decompression method."
                                        ) ::
            new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The name of the output wiggle file, or \'-\' to write to stdout. \n"+
                                         "If this ends with \".gz\" or \".zip\", then the file will automatically be compressed using the appropriate method." 
            ) :: List()
        );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         mergeSpliceBed.run2( 
             parser.get[String]("filelist"),
             parser.get[String]("outfile"),
             parser.get[Boolean]("calcMean"), 
             parser.get[Boolean]("quiet"),
             parser.get[Boolean]("ignoreSizeFactors"),
             parser.get[Option[List[Double]]]("sizeFactors"),
             parser.get[Option[String]]("rgb"),
             parser.get[Option[Double]]("filterMin"),
             parser.get[Option[Double]]("filterMax")
           );
         }
     }
   }
  
  //NOTE: Add normalization factoring, and other options?
  
  def run2(filelist : String, outfile : String, calcAverage : Boolean, quiet : Boolean, ignoreSizeFactors : Boolean, sizeFactors : Option[List[Double]], rgb : Option[String], filterMin : Option[Double], filterMax : Option[Double]){
     val initialpairlist : (Seq[(String, Double)]) = if(filelist.contains(",")){
       val files = filelist.split(",").toVector;
       val sf = if(ignoreSizeFactors || sizeFactors.isEmpty){
         repToSeq(1.0, files.length);
       } else {
         sizeFactors.get;
       }
       if(sizeFactors.get.length != files.length){
         error("Syntax error: list of wiggle files and list of size factors have different length!");
       }
       files.zip(sf);
     } else {
       val lines = getLinesSmartUnzip(filelist, true).toVector;
       if(lines.head.contains("	")){
         val files = lines.map(_.split("	")(0));
         if(ignoreSizeFactors){
           files.zip(  repToSeq(1.0, lines.length) );
         } else {
           if(! lines.forall( _.split("	").length > 1 )){
             error("Error reading file list: " + filelist + ". Not every line has a size factor listed!");
           }
           files.zip(  lines.map(l => string2double(l.split("	")(1))) );
         }
       } else {
         lines.zip( repToSeq(1.0, lines.length) );
       }
     }
     
     val makeNegativeVal =  1.toDouble;
     val calcAverageVal = if(calcAverage) initialpairlist.size.toDouble else 1.toDouble;
     
     val pairlist = initialpairlist.map( p => (p._1, (p._2 * makeNegativeVal * calcAverageVal)  ) );
     
     //val denominator : Double = if(! calcAverage){ if(makeNegative) -1.0 else 1.0 } else
     //        if(makeNegative){ - sampleList.length.toDouble} else { sampleList.length.toDouble};
     
     mergeBeds2(pairlist , outfile  , rgb, filterMin, filterMax);
  }
  
  def run(listFile : String, infile_prefix : String, infile_suffix : String, outfile : String, sizeFactor : Option[List[Double]], ignoreSizeFactors : Boolean, rgb : Option[String], filterMin : Option[Double], filterMax : Option[Double], calcMean : Boolean) {
    
    val lines = getLines(listFile).toVector;
    val cells = lines.map(_.split("	"));
    val infile_infixes = cells.map(x => x(0));
    val calcMeanVal : Double = if(calcMean) 1.0 else infile_infixes.length.toDouble;
    
    val adjustmentFactors : Seq[Double] = if(ignoreSizeFactors){
      repToSeq[Double](calcMeanVal, infile_infixes.length);
    } else if(! sizeFactor.isEmpty){
      sizeFactor.get.map( _ * calcMeanVal);
    } else if(cells.forall(c => c.length > 1)){
      cells.map(c => string2double(c(1)) * calcMeanVal );
    } else {
      repToSeq[Double](calcMeanVal, infile_infixes.length);
    }
    
    val infiles : Seq[String] = infile_infixes.map( infile_prefix + _ + infile_suffix );
    
    mergeBeds(infiles, outfile, adjustmentFactors, rgb, filterMin, filterMax);    
  }
  
  def mergeBeds2(pairlist : Seq[(String,Double)], outfile : String, rgb_opt : Option[String], filterMin : Option[Double], filterMax : Option[Double]){
     val rgb : String = if(rgb_opt.isEmpty) "255,255,255" else rgb_opt.get;
     
     val outmap : Map[GenomicInterval,Double] = pairlist.foldLeft(Map[GenomicInterval,Double]().withDefault(x => 0.0))((soFar,curr) =>{
       val lines = getLinesSmartUnzip(curr._1);
       val sf = curr._2;
       lines.foldLeft(soFar)((soFar2, line) => {
         val cells = line.split("	");
         val iv = GenomicInterval(cells(0), cells(5).charAt(0), string2int(cells(1)), string2int(cells(2)));
         val ct = (string2double(cells(3)) / sf) + soFar2(iv);
         soFar2.updated(iv,ct);
       })
     });
     
     val writer = openWriterSmart(outfile);
     for(iv <- outmap.keys.toVector.sorted){
       val ct = outmap(iv);
       if(filterMin.isEmpty || filterMin.get < ct ){
         if(filterMax.isEmpty || filterMax.get > ct){
           writer.write(makeBedLine(iv,ct,rgb));
         }
       }
       
     }
     close(writer);    
  }
  
  def mergeBeds(infiles : Seq[String], outfile : String, adjustmentFactors : Seq[Double], rgb_opt : Option[String], filterMin : Option[Double], filterMax : Option[Double]){
     val rgb : String = if(rgb_opt.isEmpty) "255,255,255" else rgb_opt.get;
     
     val outmap : Map[GenomicInterval,Double] = infiles.zip(adjustmentFactors).foldLeft(Map[GenomicInterval,Double]().withDefault(x => 0.0))((soFar,curr) =>{
       val lines = getLinesSmartUnzip(curr._1);
       val sf = curr._2;
       lines.foldLeft(soFar)((soFar2, line) => {
         val cells = line.split("	");
         val iv = GenomicInterval(cells(0), cells(5).charAt(0), string2int(cells(1)), string2int(cells(2)));
         val ct = (string2double(cells(3)) / sf) + soFar2(iv);
         soFar2.updated(iv,ct);
       })
     });
     
     val writer = openWriterSmart(outfile);
     for(iv <- outmap.keys.toVector.sorted){
       val ct = outmap(iv);
       if(filterMin.isEmpty || filterMin.get < ct ){
         if(filterMax.isEmpty || filterMax.get > ct){
           writer.write(makeBedLine(iv,ct,rgb));
         }
       }
       
     }
     close(writer);    
  }
  
  def makeBedLine(iv : GenomicInterval, ct : Double, rgb : String) : String = {
    val score = math.min(1000,math.max(0,ct.toInt));
    val len = iv.end - iv.start;
    return "" + iv.chromName + "	"+iv.start +"	"+iv.end+"	"+ ("%1.0f" format ct)+"	" + score +"	"+iv.strand+ "	" + iv.start+"	"+iv.end +"	"+rgb +"	2	1,1	0,"+(len-1)+"\n";
  }
  
}










