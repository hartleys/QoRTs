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

class qcQualityScoreCounter(isSingleEnd : Boolean, readLength : Int, maxQualScore : Int, adjustPhredScore : Int) extends QCUtility[Unit] {
  reportln("> Init QualityScoreDistribution Utility","debug");
  
  //var minQualByPos_r1 : Array[Byte] = new Array[Byte](readLen);
 // var maxQualByPos_r1 : Array[Byte] = new Array[Byte](readLen);
 // var minQualByPos_r2 : Array[Byte] = new Array[Byte](readLen);
 // var maxQualByPos_r2 : Array[Byte] = new Array[Byte](readLen);
  //val readLen = readLength;
  
  val qualByPos_r1 : Array[Array[Int]] = Array.ofDim[Int](readLength,maxQualScore + 1);
  val qualByPos_r2 : Array[Array[Int]] = Array.ofDim[Int](readLength,maxQualScore + 1);
  //var readPairCt = 0;
  //val max = 41;
  
  //reportln("   DEBUG: qualByPos_r1 dim: (" + qualByPos_r1.length + ","+ qualByPos_r1(0).length +")","debug");
  //reportln("   DEBUG: qualByPos_r2 dim: (" + qualByPos_r2.length + ","+ qualByPos_r2(0).length +")","debug");
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    if(isSingleEnd){
      runOnRead_helper(r1,true);
    } else {
      runOnRead_helper(r1,r1.getFirstOfPairFlag());
      runOnRead_helper(r2,r2.getFirstOfPairFlag());
    }
  }
   
  def runOnRead_helper(r : SAMRecord, isFirstRead : Boolean){
    val ((leadHardClip, leadSoftClip),(tailHardClip, tailSoftClip)) = getAllClipLengthsSimple(r);
    val adjustPhred : Byte = adjustPhredScore.toByte;
    val baseQuals : Seq[Byte] = if(adjustPhredScore == 0){
         r.getBaseQualities().toSeq;
       } else {
         r.getBaseQualities().toSeq.map((b : Byte) => (b -  adjustPhred).toByte)
       }
    val qualIndices = if(r.getReadNegativeStrandFlag()){
      (leadHardClip until (leadHardClip + baseQuals.length)).reverse;
    } else {
      (leadHardClip until (leadHardClip + baseQuals.length));
    }
    
    try{
      if(isFirstRead){
          baseQuals.zip(qualIndices).map((qi : (Byte, Int)) => {
             qualByPos_r1(qi._2)(qi._1.toInt) += 1;
          });
      } else {
          baseQuals.zip(qualIndices).map((qi : (Byte, Int)) => {
             qualByPos_r2(qi._2)(qi._1.toInt) += 1;
          });
      }
    }  catch {
      case e : ArrayIndexOutOfBoundsException => {
        internalUtils.Reporter.reportln("ERROR! ArrayIndexOutOfBoundsException caught while attempting to store quality score metrics.\n" + 
                                        "       Most likely cause: quality score found above "+maxQualScore+"!\n"+
                                        "       You must use the --maxPhredScore parameter to set the maximum legal quality score!\n"+
                                        "       Alternatively, maybe you used the --adjustPhredScore parameter and accidently generated\n"+
                                        "       negative Phred scores?","warn");
        internalUtils.Reporter.reportln("The offending Phred quality score string is:\n"+
                                        "    "+r.getBaseQualityString(),"warn");
        throw e;
      }
      case e : Exception => {
        throw e;
      }
    }
    
    //readPairCt += 1;
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    val pcts = Seq[Double](0.25,0.5,0.75);
    val writer = createOutputFile(outfile,"quals.r1.txt","",docWriter,
             ("readLen","int","The POSITION in the read (an unfortunate misnomer, which can't be fixed without breaking backward compatibility)"),
             ("min","int","The minimum observed PHRED quality score at the given position."),
             ("lowerQuartile","int","The lower quartile bound of the PHRED quality scores at the given position"),
             ("median","int","The median of the PHRED quality scores at the given position."),
             ("upperQuartile","int","The upper quartile bound of the PHRED quality scores at the given position"),
             ("max","int","The maximum observed PHRED quality score at the given position.")
    );
    writeOutput(pcts, qualByPos_r1, writer );

    if(! isSingleEnd){
      val writer2 = createOutputFile(outfile,"quals.r2.txt","",docWriter,
             ("readLen","int","The POSITION in the read (an unfortunate misnomer, which can't be fixed without breaking backward compatibility)"),
             ("min","int","The minimum observed PHRED quality score at the given position."),
             ("lowerQuartile","int","The lower quartile bound of the PHRED quality scores at the given position"),
             ("median","int","The median of the PHRED quality scores at the given position."),
             ("upperQuartile","int","The upper quartile bound of the PHRED quality scores at the given position"),
             ("max","int","The maximum observed PHRED quality score at the given position.")
      );
      
      writeOutput(pcts, qualByPos_r2, writer2);
    }
  }
  
  
  def writeOutput(pcts : Seq[Double], qualByPos : Array[Array[Int]], writer : WriterUtil){
    val cumulativeSum = qualByPos.map((qualCts : Array[Int]) => {
      qualCts.scanLeft(0)( (sum,ct) => sum + ct);
    }).toSeq;
    
    val out_data : Seq[Seq[Int]] = cumulativeSum.map((qualCS : Array[Int]) => {
      val pct_thresh = pcts.map((x : Double) => math.floor(x * qualCS.last.toDouble).toInt);
      (Seq(qualCS.indexWhere(_ > 0)) ++ pct_thresh.map((t : Int) => {
        qualCS.indexWhere( _ >= t );
      }).toSeq ++ Seq(qualCS.indexWhere(_ == qualCS.last)))
    }).toSeq
    
    //val writer = openWriterSmart_viaGlobalParam(outfile);
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
















