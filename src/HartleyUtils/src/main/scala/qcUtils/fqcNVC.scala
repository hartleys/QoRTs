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

object fqcNVC {
  val ALL_BASE_CHAR : Seq[Char] = Vector('A','C','G','T','N');
  val ALL_BASE_INT : Seq[Int] = ALL_BASE_CHAR.map(_.toInt)

  val BYTE_CODE_MAP : Array[Int] = Array.ofDim(ALL_BASE_INT.max + 1);
  BYTE_CODE_MAP('A'.toInt) = 0;
  BYTE_CODE_MAP('C'.toInt) = 1;
  BYTE_CODE_MAP('G'.toInt) = 2;
  BYTE_CODE_MAP('T'.toInt) = 3;
  BYTE_CODE_MAP('N'.toInt) = 4;

  
  
  def writeCounts(writer : WriterUtil, isFirstRead : Boolean, counter : Array[Array[Int]]){
    val readName = if(isFirstRead) "R1" else "R2";
    //val writer =  openWriterSmart_viaGlobalParam(outfile + ".FQ.NVC."+readName+".txt");
    writer.write("readPos\tbase\tCT\n");
    
    for(i <- 0 until counter.length){
      for((baseChar,baseInt) <- ALL_BASE_CHAR.zip(ALL_BASE_INT)){
        val ct = counter(i)(baseInt);
        writer.write("" + i + "	"+baseChar+"	"+ct+"\n");
      }
    }
    close(writer);
  }
  
}


class fqcNVC(readLength : Int, isSingleEnd : Boolean) extends fqcUtility[Nothing] {
  init();
  val readLn = readLength;
  
  val rawNvcCounter1    : Array[Array[Int]] = Array.ofDim[Int](readLn + 1, fqcNVC.ALL_BASE_INT.max + 1)
  val rawNvcCounter2    : Array[Array[Int]] = Array.ofDim[Int](readLn + 1, fqcNVC.ALL_BASE_INT.max + 1)
  
  def init(){
    reportln("> Init "+getUtilityName+" Utility","debug");
  }
  
  def runOnReadPair(r1 : FastqLine, r2 : FastqLine, readNum: Int) : Option[Nothing] = {
    r1.seqVector.zipWithIndex.foreach{case (b,i) => {
        rawNvcCounter1(i)(b.toInt) += 1;
    }}
    if(! isSingleEnd){
      r2.seqVector.zipWithIndex.foreach{case (b,i) => {
        rawNvcCounter2(i)(b.toInt) += 1;
      }}
    }
    return None;
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    //readPos\tbase\tCT\n
    
    if(true){
      val writer = createOutputFile(outfile, "FQ.NVC."+"R1"+".txt","",docWriter,
             ("readPos","int","The read position"),
             ("base","char","The base-pair"),
             ("CT","int","The number of reads with base <base> at read position <readPos>.")
      );
      fqcNVC.writeCounts(writer, isFirstRead =true, counter =rawNvcCounter1);      
    }
    if(!isSingleEnd){
      val writer = createOutputFile(outfile, "FQ.NVC."+"R2"+".txt","",docWriter,
             ("readPos","int","The read position"),
             ("base","char","The base-pair"),
             ("CT","int","The number of reads with base <base> at read position <readPos>.")
      );
      fqcNVC.writeCounts(writer, isFirstRead =false, counter =rawNvcCounter2);
    }
  }

  def getUtilityName : String = "FastqNVC";
}













