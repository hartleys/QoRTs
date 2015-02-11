package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.fileUtils._;
import internalUtils.stdUtils._;
import internalUtils.commandLineUI._;
import java.io.File;
import scala.collection.parallel.immutable.ParSeq;
 
object SumWigglesFast {
  
  class SumWigglesFast_runner extends CommandLineRunUtil {
    override def priority = 30;
    val parser : CommandLineArgParser = 
      new CommandLineArgParser(
          command = "mergeWig",  
          quickSynopsis = "", 
          synopsis = "",  
          description = "This utility merges multiple '.wig' wiggle files into a single summary '.wig' wiggle file. "+
                        "Optionally it can be used to display the mean read-pair coverage of each window across all input wiggle files rather than the sum. "+
                        "Also optionally, the mean/sum can be weighted by a set of user-supplied normalization factors."+
                        ""+
                        ""+
                        ""+
                        ""+
                        ""+
                        "",   
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
            new BinaryArgument[String](   name = "trackTitle",
                                                        arg = List("--trackTitle"),  
                                                        valueName = "options", 
                                                        argDesc = "The title of the new merged track.", 
                                                        defaultValue = Some("UntitledWig")
                                                        ) ::
            new BinaryArgument[String](   name = "additionalTrackOptions",
                                                        arg = List("--additionalTrackOptions"),  
                                                        valueName = "options", 
                                                        argDesc = "Additional track definition options, added to the track definition line. See the UCSC documentation for more information.", 
                                                        defaultValue = Some("")
                                                        ) ::
            new BinaryArgument[String](   name = "infilePrefix",
                                                        arg = List("--infilePrefix"),  
                                                        valueName = "infilePrefix", 
                                                        argDesc = "A file prefix for all input wiggle files. By default the full file path should be specified by the infile parameter.", 
                                                        defaultValue = Some("")
                                                        ) ::
            new BinaryArgument[String](   name = "infileSuffix",
                                                        arg = List("--infileSuffix"),  
                                                        valueName = "infileSuffix", 
                                                        argDesc = "A file suffix for all input wiggle files. By default the full file path should be specified by the infile parameter.", 
                                                        defaultValue = Some("")
                                                        ) ::
            new UnaryArgument(name = "ignoreSizeFactors",
                              arg = List("--ignoreSizeFactors"), // name of value
                              argDesc = "Flag to indicate that this utility should ignore size factors even if they are found in the input listFile." // description
                              ) ::
            new UnaryArgument(name = "quiet",
                              arg = List("--quiet","-q"), // name of value
                              argDesc = "" // description
                              ) :: 
                    new BinaryOptionArgument[String](
                                         name = "sizeFactorFile", 
                                         arg = List("--sizeFactorFile"), 
                                         valueName = "val",  
                                         argDesc = "A file containing (at least) two columns: a list of sample ID's and their double-precision floating-point size factors. "+
                                                   "The first line must include at least two columns: \"sample.ID\" and \"size.factor\""+
                                                   "If this option is set, all counts will be divided by the given normalization factors. The length must be the same as the length of infiles."+
                                                   "If sample.ID's is not specified by the --sampleList or --sampleListFile parameters, then all listed samples will be merged."
                                        ) ::
                    new BinaryOptionArgument[List[Double]](
                                         name = "sizeFactors", 
                                         arg = List("--sizeFactors"), 
                                         valueName = "val",  
                                         argDesc = "A list of double-precision floating-point values. "+
                                                   "If this or any size factor option is set,"+
                                                   " all counts will be divided by the given normalization factors."+
                                                   " The length must be the same as the number of files to merge."
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "filenames", 
                                         arg = List("--filenames"), 
                                         valueName = "file1.wig,file2.wig,file3.wig.gz,...",  
                                         argDesc = "A comma-delimited list of wiggle files to merge. "+
                                                   "This is optional, and filenames can be inferred from --infilePrefix, --infileSuffix, and the --sampleList, if those options are specified."+
                                                   ""+
                                                   ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "sampleList", 
                                         arg = List("--sampleList"), 
                                         valueName = "[sampleList.txt | - | samp1,samp2,samp3,...]",  
                                         argDesc = "Either a comma-delimited list of sample id's or a file containing a list of sample id's."+
                                                   "The file must either contain no title line, or contain a title line that includes a \"sample.ID\" column."+
                                                   ""+
                                                   ""
                                        ) ::
            new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The name of the output wiggle file, or \'-\' to write to stdout. \n"+
                                         "If this ends with \".gz\" or \".zip\", then the file will automatically be compressed using the appropriate method." 
            ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS
      );
    
    def run(args : Array[String]){
      parser.parseArguments(args.toList.tail);
      
      SumWigglesFast.run( parser.get[Option[String]]("sizeFactorFile"),
                          parser.get[Option[List[Double]]]("sizeFactors"),
                          parser.get[Option[List[String]]]("filenames"),
                          parser.get[Option[String]]("sampleList"),
                          parser.get[String]("outfile"), 
                          parser.get[Boolean]("makeNegative"), 
                          parser.get[Boolean]("calcMean"), 
                          parser.get[Boolean]("ignoreSizeFactors"),
                          parser.get[String]("infilePrefix"), 
                          parser.get[String]("infileSuffix"),
                          parser.get[String]("trackTitle"),
                          parser.get[String]("additionalTrackOptions")
      );
    }
  }

  def run(sizeFactorFile : Option[String], sizeFactors : Option[List[Double]], filenames : Option[List[String]], sampleList : Option[String],
         outfile : String, makeNegative : Boolean, calcAverage : Boolean, 
         ignoreSizeFactors : Boolean, 
         infilePrefix : String = "", infileSuffix : String = "", 
         trackTitle : String, additionalTrackOptions : String){

    
    val samples = fileConversionUtils.makeSpliceJunctionBed.getSampleList(sampleList, filenames);
    val sf = fileConversionUtils.makeSpliceJunctionBed.getSizeFactors(samples, ignoreSizeFactors, sizeFactorFile, sizeFactors);
    val infiles = if(filenames.isEmpty){
      samples.map(infilePrefix + _ + infileSuffix);
    } else {
      filenames.get.map(infilePrefix + _ + infileSuffix).toVector;
    }
    
    if(infiles.length != samples.length){
      error("ERROR: filenames length != # samples");
    }
     
     //getSizeFactors(input : String, ignoreSizeFactors : Boolean, sizeFactors : Option[List[Double]])
     val initialpairlist : Vector[(String,Double)] = if(sf.isEmpty){
       samples.map( (_, 1.0) );
     } else {
       samples.zip(sf.get);
     }
     
     val secondPairList : Seq[(String,Double)] = initialpairlist.map{ case (fileInfix ,sf ) => (infilePrefix + fileInfix + infileSuffix, sf)};
     
     val makeNegativeVal = if(makeNegative) (-1).toDouble else 1.toDouble;
     val calcAverageVal = if(calcAverage) secondPairList.size.toDouble else 1.toDouble;
     
     val pairlist = secondPairList.map( p => (p._1, (p._2 * makeNegativeVal * calcAverageVal)  ) );
     
     //val denominator : Double = if(! calcAverage){ if(makeNegative) -1.0 else 1.0 } else
     //        if(makeNegative){ - sampleList.length.toDouble} else { sampleList.length.toDouble};
     
     runHelper2(pairlist , outfile, trackDefLine = Some("track type=wiggle_0 name="+trackTitle+" "+additionalTrackOptions));
  }
  
  //def run(filelist : String, outfile : String, makeNegative : Boolean, calcAverage : Boolean, quiet : Boolean, ignoreSizeFactors : Boolean, sizeFactors : Option[List[Double]], infilePrefix : String = "", infileSuffix : String = "", trackDefLine : Option[String]){
  /*def run(filelist : String, outfile : String, makeNegative : Boolean, calcAverage : Boolean, quiet : Boolean, 
      ignoreSizeFactors : Boolean, sizeFactors : Option[List[Double]], 
      infilePrefix : String = "", infileSuffix : String = "", trackDefLine : Option[String]){
     reportln("runing sumWigglesFast. ignoreSizeFactors = " + ignoreSizeFactors, "debug");
    
     //getSizeFactors(input : String, ignoreSizeFactors : Boolean, sizeFactors : Option[List[Double]])
     //val initialpairlist : Vector[(String,Double)] = fileConversionUtils.mergeQcOutput.getSizeFactors(filelist, ignoreSizeFactors, sizeFactors);
     
     /*val initialpairlist : (Seq[(String, Double)]) = if((! filelist.endsWith(".txt")) & ( filelist != "-" )){
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
     }*/
     val secondPairList : Seq[(String,Double)] = initialpairlist.map{ case (fileInfix ,sf ) => (infilePrefix + fileInfix + infileSuffix, sf)};
     
     
     val makeNegativeVal = if(makeNegative) (-1).toDouble else 1.toDouble;
     val calcAverageVal = if(calcAverage) secondPairList.size.toDouble else 1.toDouble;
     
     val pairlist = secondPairList.map( p => (p._1, (p._2 * makeNegativeVal * calcAverageVal)  ) );
     
     //val denominator : Double = if(! calcAverage){ if(makeNegative) -1.0 else 1.0 } else
     //        if(makeNegative){ - sampleList.length.toDouble} else { sampleList.length.toDouble};
     
     runHelper2(pairlist , outfile  , quiet , trackDefLine = trackDefLine);
  }*/
  
  def runHelper2(pairlist : Seq[(String,Double)], outfile : String, trackDefLine : Option[String]){
    val wigIteratorList : Seq[Iterator[WigLine]] = pairlist.map(pair => {
        report("opening file: " + pair._1 +"\n        with adj factor " + pair._2+"\n","note");
        if(! internalUtils.fileUtils.fileExists(pair._1)){
          report("Error: File does not exist!: " + pair._1,"error");
        }
        val wigIter = getWigLines(pair._1).map(_.div(pair._2));
        
        if(trackDefLine.isEmpty){
          wigIter;
        } else {
          wigIter.dropWhile(_.toString.startsWith("track"));
        }
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
    if(! trackDefLine.isEmpty){
      writer.write(trackDefLine.get + "\n");
    }
    
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