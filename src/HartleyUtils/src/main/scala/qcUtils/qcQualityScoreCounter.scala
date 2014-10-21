package qcUtils

import net.sf.samtools._

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

object qcQualityScoreCounter {


  
  
  def run(infile : String, outfile : String){
    /*val (readLength, recordIter) = initSamRecordIterator(infile);
    
    val pairIterator : Iterator[(SAMRecord,SAMRecord)] = samRecordPairIterator(recordIter);
    
    val qcQSC : qcQualityScoreCounter = new qcQualityScoreCounter(readLength, MAX_QUALITY_SCORE);
    pairIterator.map((r12) => {
      qcQSC.runOnReadPair(r12._1,r12._2, 0);
    });
    
    qcQSC.writeOutput(outfile);*/
  }
  
  
  val MAX_QUALITY_SCORE = 41;
}

class qcQualityScoreCounter(readLength : Int, maxQualScore : Int) extends QCUtility[Unit] {
  reportln("Init QualityScoreDistribution","progress");
  
  //var minQualByPos_r1 : Array[Byte] = new Array[Byte](readLen);
 // var maxQualByPos_r1 : Array[Byte] = new Array[Byte](readLen);
 // var minQualByPos_r2 : Array[Byte] = new Array[Byte](readLen);
 // var maxQualByPos_r2 : Array[Byte] = new Array[Byte](readLen);
  //val readLen = readLength;
  
  val qualByPos_r1 : Array[Array[Int]] = Array.ofDim[Int](readLength,maxQualScore + 1);
  val qualByPos_r2 : Array[Array[Int]] = Array.ofDim[Int](readLength,maxQualScore + 1);
  //var readPairCt = 0;
  //val max = 41;

  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    runOnRead_helper(r1,r1.getFirstOfPairFlag());
    runOnRead_helper(r2,r2.getFirstOfPairFlag());
  }
  
  def runOnRead_helper(r : SAMRecord, isFirstRead : Boolean){
    if(isFirstRead){
      if(r.getReadNegativeStrandFlag){
        r.getBaseQualities().zip((0 until readLength).reverse).map((qi : (Byte, Int)) => {
           qualByPos_r1(qi._2)(qi._1.toInt) += 1;
        });
      } else {
        r.getBaseQualities().zip(0 until readLength).map((qi : (Byte, Int)) => {
           qualByPos_r1(qi._2)(qi._1.toInt) += 1;
        });
      }

    } else {
      if(r.getReadNegativeStrandFlag){
        r.getBaseQualities().zip((0 until readLength).reverse).map((qi : (Byte, Int)) => {
          qualByPos_r2(qi._2)(qi._1.toInt) += 1;
        });
      } else {
        r.getBaseQualities().zip(0 until readLength).map((qi : (Byte, Int)) => {
          qualByPos_r2(qi._2)(qi._1.toInt) += 1;
        });
      }
    }
    //readPairCt += 1;
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    val pcts = Seq[Double](0.0,0.25,0.5,0.75,1.0);
    
    writeOutput(pcts, qualByPos_r1, outfile + ".quals.r1.txt");
    writeOutput(pcts, qualByPos_r2, outfile + ".quals.r2.txt");
  }
  
  
  def writeOutput(pcts : Seq[Double], qualByPos : Array[Array[Int]], outfile : String){
    val cumulativeSum = qualByPos.map((qualCts : Array[Int]) => {
      qualCts.scanLeft(0)( (sum,ct) => sum + ct);
    }).toSeq;
    
    val out_data : Seq[Seq[Int]] = cumulativeSum.map((qualCS : Array[Int]) => {
      val pct_thresh = pcts.map((x : Double) => math.floor(x * qualCS.last.toDouble).toInt);
      pct_thresh.map((t : Int) => {
        qualCS.indexWhere( _ >= t );
      }).toSeq
    }).toSeq
    
    val writer = openWriterSmart_viaGlobalParam(outfile);
    writer.write("readLen	min	lowerQuartile	median	upperQuartile	max\n");
    (0 until out_data.length).map((i : Int) =>{
      writer.write(i + "");
      out_data(i).map((j : Int) => writer.write("	"+j));
      writer.write("\n");
    })
    
    close(writer);
  }
  def getUtilityName : String = "QualityScoreDistribution";
}
















