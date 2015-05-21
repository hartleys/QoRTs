package qcUtils

import net.sf.samtools._

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

object qcGCContentCount {

}

class qcGCContentCount(isSingleEnd : Boolean, readLength : Int) extends QCUtility[(Int,Int)] {
  reportln("> Init GC counts Utility","debug");
  val readLen = readLength;
  val gcContentCounts : Array[Int] = Array.ofDim[Int](2 * readLen + 1);
  val gcContentCountsR1 : Array[Int] =  Array.ofDim[Int](readLen + 1);
  val gcContentCountsR2 : Array[Int] =  Array.ofDim[Int](readLen + 1);
  val gcContentCountsRB : Array[Int] =  Array.ofDim[Int](readLen + 1);
  
  val gcContentCountsR1_byLength : Array[Array[Int]] =  Array.ofDim[Int](readLen + 1,readLen + 1);
  val gcContentCountsR2_byLength : Array[Array[Int]] =  Array.ofDim[Int](readLen + 1,readLen + 1);
  val gcContentCounts_byLength : Array[Array[Int]] =  Array.ofDim[Int](readLen + 1,readLen + 1);
  //val gcContentCountsR2_byLength : Array[Array[Int]] =  Array.ofDim[Int](readLen + 1,readLen + 1);

  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : (Int,Int) = {
    val r1c = countGCs(r1);
    val r2c = countGCs(r2);
    val gcCount = r1c + r2c;
    gcContentCounts(gcCount) += 1;
    gcContentCountsR1(r1c) += 1;
    gcContentCountsR2(r2c) += 1;
    gcContentCountsRB(r1c) += 1;
    gcContentCountsRB(r2c) += 1;
    
    val r1n = countNs(r1);
    val r2n = countNs(r2);
    
    val len1 = r1.getReadBases().length;
    val len2 = r2.getReadBases().length;
    
    gcContentCountsR1_byLength(len1 - r1n)(r1c) += 1;
    gcContentCountsR2_byLength(len2 - r2n)(r2c) += 1;
    
    gcContentCounts_byLength(len1 - r1n)(r1c) += 1;
    gcContentCounts_byLength(len2 - r2n)(r2c) += 1;
    
    return (r1c,r2c);
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    
    if(! isSingleEnd){
      val writer = openWriterSmart_viaGlobalParam(outfile + ".gc.byPair.txt");
      writer.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCounts.length){
        writer.write(i + "	"+gcContentCounts(i)+"\n");
      } 
      close(writer);
      
      val writer2 = openWriterSmart_viaGlobalParam(outfile + ".gc.R1.txt");
      writer2.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR1.length){
        writer2.write(i + "	"+gcContentCountsR1(i)+"\n");
      } 
      close(writer2);

      val writer3 = openWriterSmart_viaGlobalParam(outfile + ".gc.R2.txt");
      writer3.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR2.length){
        writer3.write(i + "	"+gcContentCountsR2(i)+"\n");
      }
      close(writer3);
    
      val writer4 = openWriterSmart_viaGlobalParam(outfile + ".gc.byRead.txt");
      writer4.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsRB.length){
        writer4.write(i + "	"+gcContentCountsRB(i)+"\n");
      }
      close(writer4);
 
      val writer5 = openWriterSmart_viaGlobalParam(outfile + ".gc.byRead.vsBaseCt.txt");
      writer5.write("NUM_BASES_GC	"+(0 until readLen).map("NONMISSING_CT_"+ _).mkString("	")+"\n");
      for(i <- 0 until gcContentCounts_byLength.length){
        writer5.write(i+"	" + (0 until readLen).map(gcContentCounts_byLength(_)(i)).mkString("	")+"\n");
      }
      close(writer5);
    } else {
      val writer2 = openWriterSmart_viaGlobalParam(outfile + ".gc.byRead.txt");
      writer2.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR1.length){
        writer2.write(i + "	"+gcContentCountsR1(i)+"\n");
      } 
      close(writer2);
      
      val writer5 = openWriterSmart_viaGlobalParam(outfile + ".gc.byRead.vsBaseCt.txt");
      writer5.write("NUM_BASES_GC	"+(0 until readLen).map("NONMISSING_CT_"+ _).mkString("	")+"\n");
      for(i <- 0 until gcContentCountsR1_byLength.length){
        writer5.write(i+"	" + (0 until readLen).map(gcContentCountsR1_byLength(_)(i)).mkString("	")+"\n");
      }
      close(writer5);
    }
  }
  
  val G_BYTE = 'G'.toByte;
  val C_BYTE = 'C'.toByte;
  val N_BYTE = 'N'.toByte;
  
  def countGCs(r : SAMRecord) : Int = {
    r.getReadBases().count( (b : Byte) =>{
      (b == G_BYTE || b == C_BYTE);
    })
  }
  def countNs(r : SAMRecord) : Int = {
    r.getReadBases().count( (b : Byte) =>{
      (b == N_BYTE);
    })
  }
  
  def getUtilityName : String = "GCDistribution";

}