package fileConversionUtils


import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
import java.io.OutputStream
import java.io.FileOutputStream
import java.io.InputStream
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.File
import scala.collection.JavaConversions._
import scala.collection.mutable.HashMap

import net.sf.samtools._

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;

object bamToWiggleLazy {

  def run(args0 : Array[String]){
    val args = args0.tail;
    if(args.length == 0 || args(0) == "?"){
      reportln("Syntax:" ,"note");
      reportln("java -jar HartleyUtils.jar bamToWiggle [options] <input.bam> <trackNamePrefix> <span> <genomeIndex.mfa.fai> <outputFilePrefix>","note");
      reportln("\n","note");
      reportln("<input.bam>                An aligned bam file.","note");
      reportln("<trackNamePrefix>          A prefix to go at the start of the track name (usually the sample name)","note");
      reportln("<span>                     The size of the summation windows.","note");
      reportln("<genomeIndex.mfa.fai>      The genome index.","note");
      reportln("<outputFilePrefix>         The prefix to the output file (including the directory location, optionally).","note");
      
      reportln("\nOptions:","note");
      reportln("--stranded                 Used for strand-specific sequencing. Two wiggle plots will be produced, " +
      		   "                           trackNamePrefix-plus-span and trackNamePrefix-minus-span. The \'plus\' " +
      		   "                           wiggle plot corresponds to read-pairs for which the second read was from " +
      		   "                           the forward (+) strand.","note");
      reportln("--noFilter                 Do NOT filter for quality metrics (map quality, mate mapped, etc.)","note");
      reportln("--noRecountOverlaps        Don't double-count paired end reads if they overlap. UNIMPLEMENTED! " +
      		   "                           Note: It's not clear whether this is a good idea.","note");
      reportln("--countByRead              By default, Bam2Wiggle counts the total number of read bases that intersect " +
      		   "                           with each window. Alternatively, when this flag is used, it instead counts " +
      		   "                           the total number of reads that intersects at any point with each window. " +
      		   "                           In most cases, the results should be very similar.","note");
      reportln("--singleEnd                By default, Bam2Wiggle is designed for use with paired-end reads. " +
      		   "                           Use this flag for unpaired reads.","note");
      reportln("--truncate                 By default, Bam2Wiggle includes a final step that is smaller than the " +
      		   "                           others. Note that wigToBigWig does not accept this, so this flag needs " +
      		   "                           to be used if you intend to use wigToBigWig.","note");
      reportln("--stranded_firstRead       If this is set, Bam2Wiggle will use the strand of the first read" +
      		   "                           to determine the strandedness of the fragment " +
      		   "                           Note that the normal \"--stranded\" option will" +
      		   "                           use the 2nd read to determine strandedness" +
      		   "                           (which is what works on my datasets)","note");
      reportln("--adjustmentFactor=factor  By default, Bam2Wiggle reports the total number of reads that covers each " +
      		   "                           window of size <span>, divided by the size of the window. Thus, if the "+
      		   "                           span = 100 and a window only has one read that covers the window from start "+
      		   "                           to finish, this tool will report \'1\' " +
      		   "                           for this window. If the adjustmentFactor is assigned, this will then be divided" +
      		   "                           by this result. The adjustment factor can be any double-precision floating-" +
      		   "                           point value. [default = 1]","note");

      reportln("Note: if you get memory problems, add \"-Xmx5000m -Xms2500m\" after the \"java\" but before the \"-jar\"","note");
      reportln("(Written by Stephen Hartley at NISC (NHGRI, NIH). stephen.hartley@nih.gov)","note");
    } else {
      val stranded = args.exists(_ == "--stranded");
      val countMethod : Int = if(args.exists(_ == "--countByRead")) 1 else 0;
      val noFilter = args.exists(_ == "--noFilter");
      val singleEnd = args.exists(_ == "--singleEnd");
      val truncate = args.exists(_ == "--truncate");
      val adjustmentFactor : Option[Double] = args.find(_.startsWith("--adjustmentFactor")) match {
        case Some(afs) =>{
          val afsa = afs.split("=");
          if(afsa.length != 2) error("Syntax error! --adjustmentFactor option must be in format \'--adjustmentFactor=somefloatvalue\' (with no whitespace!)");
          Some(string2double(afsa(1)))
        }
        case None => None
      }
      val negativeStrand = args.exists(_ == "--negativeStrandIsNegative");
        
      val inputSamOrBamFile : File = new File(args(args.length - 5));
      val trackName : String = args(args.length - 4);
      val span : Int = string2int(args(args.length - 3));
      val chromLengthFile : String = args(args.length - 2);
      val outfile : String = args(args.length -1);
      
      val getStrandFromFirstRead = args.exists(_ == "--stranded_firstRead");
    
      createWiggle(inputSamOrBamFile , trackName, countMethod , (stranded || getStrandFromFirstRead), noFilter, truncate, singleEnd, span , chromLengthFile , outfile , adjustmentFactor, negativeStrand, getStrandFromFirstRead);
    }
  }
  
  case class Chrom(chromName :  String, chromStrand : Char, windowCounts : HashMap[Int,Int], span : Int, spanCt : Int){
    def countSamRecord(samRecord : SAMRecord, method : Int, truncate : Boolean) {
        val blocks = samRecord.getAlignmentBlocks();
      
        var i = 0;
        while(i < blocks.size()){
          val block : AlignmentBlock = blocks.get(i);
          val start = block.getReferenceStart();
          val length = block.getLength();
          val end = start + length;
          
          if(! truncate && (end - 1 / span) >=  spanCt){
            error("ERROR: read extends outside chromosome length!");
          }
          
          var j = start;
          while(j < end && (j / span) < spanCt){
            windowCounts += (( j / span , windowCounts(j / span) + 1));
            j = j + 1;
          }
          
          i = i + 1;
        }
    }
  }
  
  // Count Methods:
  //    0    [default]    coverage-based. For a given span, counts the total # of times that each base is covered by a read.
  //    1                 read-count. For a given span, counts the total # of reads that lie on that span.
  
  def createWiggle(inputSamOrBamFile : File, trackName : String, countMethod : Int, stranded : Boolean, noFilter : Boolean, truncate : Boolean, singleEnd : Boolean, span : Int, chromLengthFile : String, outfile : String, adjustmentFactor : Option[Double], negativeStrand : Boolean, getStrandFromFirstRead : Boolean) {
    val inputSam : SAMFileReader = new SAMFileReader(inputSamOrBamFile);

    val chromMap : Map[(String,Char),Chrom] = genChrom(chromLengthFile,span,stranded, truncate);
    var unknownChromSet : Set[(String,Char)] = Set[(String,Char)]();
    
    for(samRecord : SAMRecord <- inputSam){
      if(useRead(samRecord, noFilter, singleEnd)){
        val chromName = samRecord.getReferenceName();
        val strand = if(stranded) getStrand(samRecord,singleEnd,getStrandFromFirstRead) else '*';
        
        chromMap.get((chromName,strand)) match {
          case Some(chrom) => {
            chrom.countSamRecord(samRecord,countMethod, truncate);
          }
          case None => {
           // error("Chromosome " + chromName + " not found! Strand: \'" + strand+ "\'. Something is wrong with the bam file read!");
          //  report("!["+chromName+","+strand+"]","deepDebug");
            if(! unknownChromSet.contains((chromName,strand))){
              reportln("WARNING: CHROMOSOME NOT RECOGNIZED: ["+chromName+","+strand+"]","deepDebug");
              unknownChromSet += ((chromName,strand));
            }
          }
        }        
      }
    }
    
    val totalDropped = causeOfDropArray.foldLeft(0)(_ + _);
    
    startReport("Dropped: " + totalDropped + " reads.\n","note");
    report("| Dropped: " + causeOfDropArray(0) + ", for mate unmapped\n","note");
    report("| Dropped: " + causeOfDropArray(1) + ", for improperly paired\n","note");
    report("| Dropped: " + causeOfDropArray(2) + ", for failing vendor quality check flag\n","note");
    report("| Dropped: " + causeOfDropArray(3) + ", for not being a primary alignment\n","note");
    report("| Dropped: " + causeOfDropArray(4) + ", for not being \'valid\' for whatever reason, according to picard tools.\n","note");
    
    if(stranded){
      val writerP = openWriter(outfile+"plus"+".wig");
      val writerM = openWriter(outfile+"minus"+".wig");
      initWigglePlot(trackName + "-plus",span,writerP);
      initWigglePlot(trackName + "-minus",span,writerM);
      val sortedKeyList : List[(String,Char)] = chromMap.keySet.toList.sorted
      for(chromPairs <- sortedKeyList){
        val chrom = chromMap.getOrElse(chromPairs,null);
        if(chrom.chromStrand == '+'){
          writeChrom(chrom,writerP,adjustmentFactor,false);
        } else {
          writeChrom(chrom,writerM,adjustmentFactor,negativeStrand);
        }
      }
    
      close(writerP);
      close(writerM);
    } else {
      val writer = openWriter(outfile+".wig");
      initWigglePlot(trackName+"-eitherStrand",span,writer);
      val sortedKeyList : List[(String,Char)] = chromMap.keySet.toList.sorted
      for(chromPairs <- sortedKeyList){
        val chrom = chromMap.getOrElse(chromPairs,null);
      
        writeChrom(chrom,writer,adjustmentFactor,false);
      }
      close(writer);
    }
    

  }
  
  def initWigglePlot( name : String, span : Int, writer: WriterUtil){
    
    write("track type=wiggle_0 name="+name+"-"+span+" visibility=full\n",writer);
  }
  
  def writeChrom(chrom : Chrom,writer : WriterUtil, adjustmentFactor : Option[Double], negativeStrand : Boolean){
    write("fixedStep  chrom="+chrom.chromName+"  start=1  step="+chrom.span+"\n",writer);
    
    adjustmentFactor match {
      case Some(af) => {
        val intAF : Double = if(negativeStrand ) (0.toDouble - af.toDouble) else (af.toDouble);
        
        var i = 0;
         while(i < chrom.spanCt){
            write(((chrom.windowCounts(i).toDouble / intAF) / chrom.span.toDouble) +"\n",writer);
            i = i + 1;
         }
      }
      case None => {
        val intAF : Int = if(negativeStrand) (-1) else (1);
        
        var i = 0;
        while(i < chrom.spanCt){
           write(((chrom.windowCounts(i) * (intAF)) / (chrom.span.toDouble) ) +"\n",writer);
           i = i + 1;
        }
      }
    }

    
  }
  
  def genChrom(chromLengthFile : String, span : Int, stranded : Boolean, truncate : Boolean) : Map[(String,Char),Chrom] = {
    var chromMap = Map[(String,Char),Chrom]();
    
    val lines = getLines(chromLengthFile);
    
    for(line <- lines){
      val cells = line.split("\\s+");
      val chromName = cells(0);
      val chromLength = string2int(cells(1));
      val spanCt =  if(chromLength % span == 0 || truncate){
                       chromLength / span;
                    } else {
                       (chromLength / span) + 1;
                    }
      
      if(stranded){
        val chromP = new Chrom(chromName,'+',new HashMap[Int,Int](){ override def default(key:Int) = 0 },span, spanCt);
        val chromM = new Chrom(chromName,'-',new HashMap[Int,Int](){ override def default(key:Int) = 0 },span, spanCt);
        chromMap = chromMap + (((chromName,'+'), chromP));
        chromMap = chromMap + (((chromName,'-'), chromM));        
      } else {
        val chrom = new Chrom(chromName,'*',new HashMap[Int,Int](){ override def default(key:Int) = 0 },span, spanCt);
        chromMap = chromMap + (((chromName,'*'), chrom));
      }
    }
    chromMap;
  }
  
//Checked.
  def getStrand(samRecord : SAMRecord, singleEnd : Boolean, getStrandFromFirstRead : Boolean) : Char = {
    if(getStrandFromFirstRead){
      if(singleEnd || samRecord.getFirstOfPairFlag()){
        if(samRecord.getReadNegativeStrandFlag) '-' else '+';
      } else {
        if(samRecord.getMateNegativeStrandFlag) '-' else '+';
      }
    } else {
      if(singleEnd || ! samRecord.getFirstOfPairFlag()){
        if(samRecord.getReadNegativeStrandFlag) '-' else '+';
      } else {
        if(samRecord.getMateNegativeStrandFlag) '-' else '+';
      }
    }
    
  }
  
  def isFirstRead(samRecord : SAMRecord) : Boolean = {
    samRecord.getFirstOfPairFlag();
  }
  
  val causeOfDropArray : Array[Int] = Array.ofDim[Int](5);
  def resetCauseOfDropArray = {
    for(i <- 0 until causeOfDropArray.length){
      causeOfDropArray(i) = 0;
    }
  }
  
  def useRead(samRecord : SAMRecord, noFilter : Boolean, singleEnd : Boolean) : Boolean = {
    if(noFilter) true;
    else {
      if(! singleEnd){
        if(samRecord.getMateUnmappedFlag()) { causeOfDropArray(0) = causeOfDropArray(0) + 1; return false; }
        if(! samRecord.getProperPairFlag() ) { causeOfDropArray(1) = causeOfDropArray(1) + 1; return false; }
      }
   // if(! samRecord.getFirstOfPairFlag()) return false;
      if(samRecord.getReadFailsVendorQualityCheckFlag()){ causeOfDropArray(2) = causeOfDropArray(2) + 1; return false; }
      if(samRecord.getNotPrimaryAlignmentFlag()) { causeOfDropArray(3) = causeOfDropArray(3) + 1; return false; }
      if( samRecord.isValid != null) { causeOfDropArray(4) = causeOfDropArray(4) + 1; return false; }
//      if( samRecord.getMappingQuality() < mapqThresh){ causeOfDropArray(5) = causeOfDropArray(5) + 1; return false; }
    //add more checks here?
      return true;
    }
  }
  
}