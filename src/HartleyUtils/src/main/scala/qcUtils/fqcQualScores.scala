package qcUtils

import net.sf.samtools._
import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.JavaConversions._;

import internalUtils.genomicUtils._;
import internalUtils.optionHolder._;

object fqcQualScores {
  


  
}

class fqcQualScores(readLength : Int, isSingleEnd : Boolean, maxQualScore : Int) extends fqcUtility[Nothing] {
  init();
  
  //GC tools:
  val qualByPos_r1 : Array[Array[Int]] = Array.ofDim[Int](readLength,maxQualScore + 1);
  val qualByPos_r2 : Array[Array[Int]] = Array.ofDim[Int](readLength,maxQualScore + 1);
  
  def init(){
    reportln("> Init "+getUtilityName+" Utility","debug");
  }
  
  def runOnReadPair(r1 : FastqLine, r2 : FastqLine, readNum: Int) : Option[Nothing] = {
    r1.qualVector.zipWithIndex.foreach{case (q,i) => {
        qualByPos_r1(i)(q) += 1;
    }}
    if(! isSingleEnd){
      r2.qualVector.zipWithIndex.foreach{case (q,i) => {
        qualByPos_r2(i)(q) += 1;
      }}
    }
    
    return None;
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    val pcts = Seq[Double](0.25,0.5,0.75);
    
    writeOutput(pcts, qualByPos_r1, outfile , "FQ.quals.r1.txt", docWriter = docWriter);
    
    if(! isSingleEnd){
      writeOutput(pcts, qualByPos_r2, outfile , "FQ.quals.r2.txt", docWriter = docWriter);
    }
  }
  
  def writeOutput(pcts : Seq[Double], qualByPos : Array[Array[Int]], outfile : String, outfilesuffix : String, docWriter : DocWriterUtil){
    val cumulativeSum = qualByPos.map((qualCts : Array[Int]) => {
      qualCts.scanLeft(0)( (sum,ct) => sum + ct);
    }).toSeq;
    
    val out_data : Seq[Seq[Int]] = cumulativeSum.map((qualCS : Array[Int]) => {
      val pct_thresh = pcts.map((x : Double) => math.floor(x * qualCS.last.toDouble).toInt);
      (Seq(qualCS.indexWhere(_ > 0)) ++ pct_thresh.map((t : Int) => {
        qualCS.indexWhere( _ >= t );
      }).toSeq ++ Seq(qualCS.indexWhere(_ == qualCS.last)))
    }).toSeq
    
    val writer = createOutputFile(outfile,outfilesuffix,"",docWriter,
             ("readLen","int","The POSITION in the read (an unfortunate misnomer, which can't be fixed without breaking backward compatibility)"),
             ("min","int","The minimum observed PHRED quality score at the given position."),
             ("lowerQuartile","int","The lower quartile bound of the PHRED quality scores at the given position"),
             ("median","int","The median of the PHRED quality scores at the given position."),
             ("upperQuartile","int","The upper quartile bound of the PHRED quality scores at the given position"),
             ("max","int","The maximum observed PHRED quality score at the given position.")
    );
    writer.write("readLen	min	lowerQuartile	median	upperQuartile	max\n");
    (0 until out_data.length).map((i : Int) =>{
      writer.write(i + "");
      out_data(i).map((j : Int) => writer.write("	"+j));
      writer.write("\n");
    })
    
    
    
    close(writer);
  }
  
  def getUtilityName : String = "FastqQualityScore";
}













