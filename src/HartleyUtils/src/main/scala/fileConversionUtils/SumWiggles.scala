package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.fileUtils._;
import internalUtils.stdUtils._;
import internalUtils.commandLineUI._;
import java.io.File;
import scala.collection.parallel.immutable.ParSeq;

object SumWiggles {

  def run(args0 : Array[String]){
    val args = args0.tail;
    if(args.length < 5 || args.exists(_ == "?")){
      reportln("Syntax:" ,"note");
      reportln("java -jar HartleyUtils.jar sumWiggles [options] <sampleList.txt> <infile_prefix> <infile_suffix> <out_filename>","note");
      reportln("\n","note");
      reportln("<sampleList.txt>           A list with the names of the samples to be added together. 1 name per line.\n" +
      		   "                           The input files are assumed to be of format:\n." +
      		   "                           <infile_prefix><samplename><infile_suffix>","note");
      reportln("<infile_prefix>            The prefix, the start of the input file names, including the directory.","note");
      reportln("<infile_suffix>            The suffix, the end of the input file names.","note");
      reportln("<out_filename>             The name of the output file","note");

      reportln("\nOptions:","note");
      reportln("--makeNegative             multiply the result by -1.","note");
      reportln("--calcAverage              Average the sum of all supplied wiggle files.","note");
      reportln("--quiet                    Do not report progress.","note");

      reportln("Note: if you get memory problems, add \"-Xmx15000m -Xms2500m\" after the \"java\" but before the \"-jar\"","note");
      reportln("(Written by Stephen Hartley at NISC (NHGRI, NIH). stephen.hartley@nih.gov)","note");
    } else {
      val makeNegative : Boolean = args.exists(_ == "--makeNegative");
      val calcAverage : Boolean = args.exists(_ == "--calcAverage")
      val quiet : Boolean = args.exists(_ == "--quiet")
        
      val sampleListFile : String = args(args.length - 4);
      val infile_prefix : String = args(args.length - 3);
      val infile_suffix : String = args(args.length - 2);
      val outfile : String = args(args.length -1);
      
      run(sampleListFile , infile_prefix , infile_suffix , outfile , makeNegative , calcAverage , quiet );
    }
  }

  def run2(sampleListFile : String, infile_prefix : String, infile_suffix : String, outfile : String, makeNegative : Boolean, calcAverage : Boolean, quiet : Boolean){
     val sampleList = getLines(sampleListFile);
     var curr : List[WigLine] = List();
     var fileCt : Int = 0;
     
     for(sample <- sampleList){
       val subj = sample;
       if(curr.length == 0){
         curr = readWiggleToList(infile_prefix+ subj + infile_suffix);
       } else {
         curr = curr.zip(readWiggleToList(infile_prefix+ subj + infile_suffix)).map(t => t._1.add(t._2));
       }
       fileCt += 1;
       if(! quiet) report("Done with" + subj,"debug");
     }
     if(calcAverage){
       if(! quiet) report("Averaging counts.","debug");
       curr = curr.map(_.mult( ((1).toDouble) / (fileCt.toDouble)));
     }
     
     writeWiggleList(curr,outfile,makeNegative);

  }
  
  def run(sampleListFile : String, infile_prefix : String, infile_suffix : String, outfile : String, makeNegative : Boolean, calcAverage : Boolean, quiet : Boolean){
     val sampleList = getLines(sampleListFile).filter(_.trim != "").toSeq;
     val fileCt : Int = sampleList.length;
     val sums = (sampleList.par.foldLeft(Seq[WigLine]().par)( (currSeq, sample) => {
       currSeq.zipAll(readWiggleToSeq(infile_prefix+ sample + infile_suffix), DoubleLine(0.0), DoubleLine(0.0)).map(t => t._1.add(t._2)).toSeq
     }))

     val out = if(! calcAverage) sums else {
       if(! quiet) report("Averaging counts.","debug");
       sums.map(_.mult( ((1).toDouble) / (fileCt.toDouble)));
     }
     
     writeWiggleList(out.seq,outfile,makeNegative);
  }
  
  
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
    if(makeNegative) wigList.map(x => writer.write(x.mult(-1.0).toString() + "\n"));
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
  
  abstract class WigLine {
    def add(w : WigLine) : WigLine
    def mult(d : Double) : WigLine
    override def toString() : String
  }
  case class DoubleLine(value : Double) extends WigLine {
    def add(w : WigLine) : WigLine = w match {
      case DoubleLine(v2) => DoubleLine(value + v2);
      case StringLine(v2) => StringLine(v2);
    }
    def mult(d : Double) : WigLine = DoubleLine(value * d);
    override def toString() : String = value.toString();
  }
  case class StringLine(value : String) extends WigLine {
    def add(w : WigLine) : WigLine = this;
    def mult(d : Double) : WigLine = this;
    override def toString() : String = value;
  }
}