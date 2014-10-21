package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.fileUtils._;
import internalUtils.stdUtils._;
import internalUtils.commandLineUI._;
import java.io.File;
import scala.collection.parallel.immutable.ParSeq;
 
object SumWigglesFast {
  
  class SumWigglesFast_runner extends CommandLineRunUtil {
    val parser : CommandLineArgParser = 
      new CommandLineArgParser(
          command = "mergeWig",  
          quickSynopsis = "", 
          synopsis = "",  
          description = "",   
          argList = 
            new UnaryArgument(name = "makeNegative",
                              arg = List("--makeNegative","-n"), // name of value
                              argDesc = "Flag to indicate that every counting bin value should be multiplied by -1" // description
                              ) :: 
            new UnaryArgument(name = "calcMean",
                              arg = List("--calcMean","-m"), // name of value
                              argDesc = "Flag to indicate that the mean average should be calculated, rather than the sum." // description
                              ) ::
            new BinaryOptionArgument[List[Double]](
                                         name = "sizeFactors", 
                                         arg = List("--sizeFactors"), 
                                         valueName = "val,val,val,...",  
                                         argDesc = "normalization factors for each wig file."
                                        ) ::
             new BinaryArgument[String](   name = "infilePrefix",
                                                        arg = List("--infilePrefix"),  
                                                        valueName = "infilePrefix", 
                                                        argDesc = "A file prefix for all input junction count files. By default the full file path should be specified by the infile parameter.", 
                                                        defaultValue = Some("")
                                                        ) ::
            new BinaryArgument[String](   name = "infileSuffix",
                                                        arg = List("--infileSuffix"),  
                                                        valueName = "infileSuffix", 
                                                        argDesc = "A file suffix for all input junction count files. By default the full file path should be specified by the infile parameter.", 
                                                        defaultValue = Some("")
                                                        ) ::
            new UnaryArgument(name = "ignoreSizeFactors",
                              arg = List("--ignoreSizeFactors"), // name of value
                              argDesc = "Flag to indicate that this utility should ignore size factors even if they are found in the input listFile. MAKE SURE NOT TO APPLY SIZE FACTORS TWICE!" // description
                              ) ::
            new UnaryArgument(name = "quiet",
                              arg = List("--quiet","-q"), // name of value
                              argDesc = "" // description
                              ) :: 
            new FinalArgument[String](
                                         name = "filelist",
                                         valueName = "<filelist.txt | - | file1.wig,file2.wig,...>",
                                         argDesc = "One of three things:"+
                                                   "(1) A comma-delimited list of wig files. Note: NO WHITESPACE!"+
                                                   "(2) A file (or '-' to read from stdin) containing a list of the wig files to merge (one on each line).\n"+
                                                   "Optionally, the file can have a second (tab-delimited) column, containing the normalization factors to use. If this column is present, "+
                                                   "this utility will automatically calculate the normalized totals (or the normalized mean, if --calcMean is set).\n"+
                                                   "Note: wig filenames cannot contain whitespace or commas.\n"+
                                                   "Also Note: if the wig file names end in \".gz\" or \".zip\", they will be read using the appropriate decompression method."
                                        ) ::
            new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The name of the output wiggle file, or \'-\' to write to stdout. \n"+
                                         "If this ends with \".gz\" or \".zip\", then the file will automatically be compressed using the appropriate method." 
            ) :: List()
      );
    
    def run(args : Array[String]){
      parser.parseArguments(args.toList.tail);
      
      SumWigglesFast.run( parser.get[String]("filelist"),
                          parser.get[String]("outfile"), 
                          parser.get[Boolean]("makeNegative"), 
                          parser.get[Boolean]("calcMean"), 
                          parser.get[Boolean]("quiet"),
                          parser.get[Boolean]("ignoreSizeFactors"),
                          parser.get[Option[List[Double]]]("sizeFactors"),
                          parser.get[String]("infilePrefix"), 
                          parser.get[String]("infileSuffix")
                          
      );
    }
  }
  
  def run(filelist : String, outfile : String, makeNegative : Boolean, calcAverage : Boolean, quiet : Boolean, ignoreSizeFactors : Boolean, sizeFactors : Option[List[Double]], infilePrefix : String = "", infileSuffix : String = ""){
    
     reportln("runing sumWigglesFast. ignoreSizeFactors = " + ignoreSizeFactors, "debug");
    
     val initialpairlist : (Seq[(String, Double)]) = if(filelist.contains(",")){
       val files = filelist.split(",").toVector;
       val sf = if(ignoreSizeFactors || sizeFactors.isEmpty){
         repToSeq(1.0, files.length);
       } else {
         sizeFactors.get;
       }
       if(sf.length != files.length){
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
     val secondPairList : Seq[(String,Double)] = initialpairlist.map{ case (fileInfix ,sf ) => (infilePrefix + fileInfix + infileSuffix, sf)};
     
     
     val makeNegativeVal = if(makeNegative) (-1).toDouble else 1.toDouble;
     val calcAverageVal = if(calcAverage) secondPairList.size.toDouble else 1.toDouble;
     
     val pairlist = secondPairList.map( p => (p._1, (p._2 * makeNegativeVal * calcAverageVal)  ) );
     
     //val denominator : Double = if(! calcAverage){ if(makeNegative) -1.0 else 1.0 } else
     //        if(makeNegative){ - sampleList.length.toDouble} else { sampleList.length.toDouble};
     
     runHelper2(pairlist , outfile  , quiet );
  }
  
  def runHelper2(pairlist : Seq[(String,Double)], outfile : String, quiet : Boolean){
    val wigIteratorList : Seq[Iterator[WigLine]] = pairlist.map(pair => {
        report("opening file: " + pair._1 +"\n        with adj factor " + pair._2+"\n","note");
        if(! internalUtils.fileUtils.fileExists(pair._1)){
          report("Error: File does not exist!: " + pair._1,"error");
        }
        getWigLines(pair._1).map(_.div(pair._2))
    });
    
    reportln("Made iterators.","debug");
    
    val sumIterator : Iterator[WigLine] = new Iterator[WigLine](){
      var lnct = 0;
      def next : WigLine = {
        lnct += 1;
        if(lnct % 1000000 == 0){if(lnct % 10000000 == 0) report(". ","progress"); else report(".","progress");}
        wigIteratorList.tail.foldLeft(wigIteratorList.head.next)((sum, iter) => sum.add(iter.next));
      }
      def hasNext : Boolean = wigIteratorList.head.hasNext;
    }
    reportln("Made iterators 2.","debug");
    
    val lineIterator = sumIterator.map((wigLine) => wigLine.toString() + "\n");
    report("\n","progress");
    
    val writer = openWriterSmart(outfile, true);
    lineIterator.foreach(writer.write(_));
    close(writer);
  }
  
  def runHelper(sampleList : Seq[String], infile_prefix : String, infile_suffix : String, adjustmentFactors : Seq[Double], outfile : String, quiet : Boolean){
     val ziplist : Seq[(String,Double)] = sampleList.zip(adjustmentFactors);
     
     val wigIteratorList : Seq[Iterator[WigLine]] = ziplist.map(pair => getWigLines(infile_prefix + pair._1 + infile_suffix).map(_.div(pair._2)));
     
     val sumIterator : Iterator[WigLine] = new Iterator[WigLine](){
       override def next : WigLine = wigIteratorList.tail.foldLeft(wigIteratorList.head.next)((sum, iter) => sum.add(iter.next));
       override def hasNext : Boolean = wigIteratorList.head.hasNext;
     }
     
     val writer = openWriter(outfile);
     
     val lineIterator = sumIterator.map((wigLine) => wigLine.toString() + "\n")
     
     lineIterator foreach writer.write
     
     close(writer);
  }
  
  /*def runSum(sampleListFile : String, infile_prefix : String, infile_suffix : String, outfile : String, makeNegative : Boolean, quiet : Boolean){
     val sampleList = getLines(sampleListFile).filter(_.trim != "").toSeq;
     val fileCt : Int = sampleList.length;
     
     val wigIteratorList : Seq[Iterator[WigLine]] = sampleList.map(x => getWigLines(infile_prefix + x + infile_suffix));
     
     val sumIterator : Iterator[WigLine] = new Iterator[WigLine](){
       override def next : WigLine = wigIteratorList.tail.foldLeft(wigIteratorList.head.next)((sum, iter) => sum.add(iter.next));
       override def hasNext : Boolean = wigIteratorList.head.hasNext;
     }
     
     //writeWiggleList(out.seq,outfile,makeNegative);
  }*/
  
  
  
  
  /*def runOld(args : Array[String]) {
    val decoderFile : String = args(1);
    val fileSuffix : String = args(2);
    val groupName : String = args(3);
    val makeNegative : Boolean = if(args(4) == "-") true else false;
    val calcAverage : Boolean = if(args.length > 5 && args(5) == "average") true else false;
    
    val decoderLines = getLines(decoderFile);
    
    var curr : List[WigLine] = List();
    var fileCt : Int = 0;
    
    for(decoderLine <- decoderLines){
      val cells = decoderLine.split("\\s+");
      val subj = cells(0);
      if(curr.length == 0){
        curr = readWiggleToList("../"+ subj + fileSuffix);
      } else {
        curr = curr.zip(readWiggleToList("../"+ subj + fileSuffix)).map(t => t._1.add(t._2));
      }
      fileCt += 1;
    }
    
    if(calcAverage){
      report("Averaging counts.","debug");
      curr = curr.map(_.mult( ((1).toDouble) / (fileCt.toDouble)));
    }
    
    writeWiggleList(curr,groupName+fileSuffix,makeNegative);
  }*/
  
  def writeWiggleList(wigList : Seq[WigLine], outfile : String, makeNegative : Boolean){
    val writer = openWriter(outfile);
    if(makeNegative) wigList.map(x => writer.write(x.div(-1.0).toString() + "\n"));
    else wigList.map(x => writer.write(x.toString() + "\n"));
    close(writer);
  }
  
  def readWiggleToList(infile : String) : List[WigLine] = {
    getLines(infile).map((x : String) => makeWigLine(x)).toList;
  }
  def readWiggleToSeq(infile : String) : Seq[WigLine] = {
    reportln("reading: " + infile,"note")
    getLines(infile).map((x : String) => makeWigLine(x)).toSeq;
  }
  
  def makeWigLine(s : String) : WigLine = {
    if(s.startsWith("track") || s.startsWith("fixedStep")) StringLine(s);
    else DoubleLine(string2double(s));
  }
  
  def getWigLines(s : String) : Iterator[WigLine] = {
    new Iterator[WigLine](){
       val internalIterator = getLinesSmartUnzip(s);
       override def next : WigLine = makeWigLine(internalIterator.next);
       override def hasNext : Boolean = internalIterator.hasNext;
     }
  }
  
  abstract class WigLine {
    def add(w : WigLine) : WigLine
    def div(d : Double) : WigLine
    override def toString() : String
  }
  case class DoubleLine(value : Double) extends WigLine {
    def add(w : WigLine) : WigLine = w match {
      case DoubleLine(v2) => DoubleLine(value + v2);
      case StringLine(v2) => StringLine(v2);
    }
    def div(d : Double) : WigLine = DoubleLine(value / d);
    override def toString() : String = value.toString();
  }
  case class StringLine(value : String) extends WigLine {
    def add(w : WigLine) : WigLine = this;
    def div(d : Double) : WigLine = this;
    override def toString() : String = value;
  }
  
  
}