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

object qcOverlapMatch {


  val ALL_BASE_BYTES : Seq[Byte] = Vector('A'.toByte,'C'.toByte,'G'.toByte,'T'.toByte,'N'.toByte);
  val BYTE_CODE_MAP : Array[Int] = Array.ofDim(ALL_BASE_BYTES.max + 1);
  BYTE_CODE_MAP('A'.toInt) = 0;
  BYTE_CODE_MAP('C'.toInt) = 1;
  BYTE_CODE_MAP('G'.toInt) = 2;
  BYTE_CODE_MAP('T'.toInt) = 3;
  BYTE_CODE_MAP('N'.toInt) = 4;
}


class qcOverlapMatch(isSingleEnd : Boolean, readLen : Int, anno_holder : qcGtfAnnotationBuilder, stranded : Boolean, fr_secondStrand_bool : Boolean) extends QCUtility[Unit] {
  reportln("> Init OverlapMatch Utility","debug");
  
  var ppoLengthMismatch = 0;
  var ppoReferencePositionMismatch = 0;
  var ppoMatch = 0;
  
  val basePairCount : Array[Array[Array[Array[Int]]]] = Array.ofDim(readLen, readLen,qcOverlapMatch.ALL_BASE_BYTES.length, qcOverlapMatch.ALL_BASE_BYTES.length);
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    val (start1, end1) = (r1.getAlignmentStart() - 1, r1.getAlignmentEnd);
    val (start2, end2) = (r2.getAlignmentStart() - 1, r2.getAlignmentEnd);
    
    val pp1 = getAlignedBasePositionsFromRead(r1);
    val pp2 = getAlignedBasePositionsFromRead(r2);
    
    val ppo1 = pp1.filter((i) => i._1 >= start2 && i._1 < end2);
    val ppo2 = pp2.filter((i) => i._1 >= start1 && i._1 < end1);
    
    val seq1 = r1.getReadBases();
    val seq2 = r2.getReadBases();
    
    if(ppo1.length == ppo2.length){
      val ppp = ppo1.zip(ppo2);
      if( ppp.forall(pair => pair._1._1 == pair._2._1)){
        ppoMatch += 1;
        for(pair <- ppp){
          basePairCount(pair._1._2)(pair._2._2)( qcOverlapMatch.BYTE_CODE_MAP( seq1(pair._1._2).toInt) )( qcOverlapMatch.BYTE_CODE_MAP(seq2(pair._2._2).toInt) ) += 1;
        }
      } else {
        ppoReferencePositionMismatch += 1;
        //do stuff!
      }
    } else {
      ppoLengthMismatch += 1;
      //do stuff!
    }
    
    /*
    val cb1 = getCigarBlocksFromRead(r1);
    val cb2 = getCigarBlocksFromRead(r2);
    
    val cbo1 = cb1.filter((b : CigarBlock) => b.refStart < end2 & b.refEnd > start2);
    val cbo2 = cb2.filter((b : CigarBlock) => b.refStart < end1 & b.refEnd > start1);
    
    val M1 = cbo1.filter(_.op.equals(CigarOperator.MATCH_OR_MISMATCH));
    val M2 = cbo2.filter(_.op.equals(CigarOperator.MATCH_OR_MISMATCH));
    
    */
  }
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    val writer = openWriterSmart_viaGlobalParam(outfile + ".clippedSeqNvc.5prime.txt");
    writer.write("READ_POSITION_R1  READ_POSTION_R2 ");
    
    close(writer);
  }
  
  def getUtilityName : String = "StrandCheck";
}












