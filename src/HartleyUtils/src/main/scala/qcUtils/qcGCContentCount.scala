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

class qcGCContentCount(readLength : Int) extends QCUtility[Unit] {
  reportln("Init GC counts","progress");
  val readLen = readLength;
  val gcContentCounts : Array[Int] = Array.ofDim[Int](2 * readLen + 1);
  val gcContentCountsR1 : Array[Int] =  Array.ofDim[Int](readLen + 1);
  val gcContentCountsR2 : Array[Int] =  Array.ofDim[Int](readLen + 1);
  val gcContentCountsRB : Array[Int] =  Array.ofDim[Int](readLen + 1);
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    val r1c = countGCs(r1);
    val r2c = countGCs(r2);
    val gcCount = r1c + r2c;
    gcContentCounts(gcCount) += 1;
    gcContentCountsR1(r1c) += 1;
    gcContentCountsR2(r2c) += 1;
    gcContentCountsRB(r1c) += 1;
    gcContentCountsRB(r2c) += 1;
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    val writer = openWriterSmart_viaGlobalParam(outfile + ".gc.txt");
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
    
    val writer4 = openWriterSmart_viaGlobalParam(outfile + ".gc.RB.txt");
    writer4.write("NUM_BASES_GC	CT\n");
    for(i <- 0 until gcContentCountsRB.length){
      writer4.write(i + "	"+gcContentCountsRB(i)+"\n");
    }
    close(writer4);
  }
  
  val G_BYTE = 'G'.toByte;
  val C_BYTE = 'C'.toByte;
  
  def countGCs(r : SAMRecord) : Int = {
    r.getReadBases().count( (b : Byte) =>{
      (b == G_BYTE || b == C_BYTE);
    })
  }
  
  def getUtilityName : String = "GCDistribution";

}