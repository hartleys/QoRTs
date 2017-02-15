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

object fqcGC {

}

class fqcGC(readLength : Int, isSingleEnd : Boolean) extends fqcUtility[(Int,Int)] {
  init();
  
  //GC tools:
  val readLen = readLength;
  val gcContentCounts : Array[Int] = Array.ofDim[Int](2 * readLen + 1);
  val gcContentCountsR1 : Array[Int] =  Array.ofDim[Int](readLen + 1);
  val gcContentCountsR2 : Array[Int] =  Array.ofDim[Int](readLen + 1);
  val gcContentCountsRB : Array[Int] =  Array.ofDim[Int](readLen + 1);
  
  val readLenCounts : Array[Array[Int]] = Array.ofDim[Int](2,readLen + 1);
  val nonMissCounts : Array[Array[Int]] = Array.ofDim[Int](2,readLen + 1);
  val gcContentCounts_byLength : Array[Array[Array[Int]]] =  Array.ofDim[Int](2, readLen + 1,readLen + 1);
  
  
  def init(){
    reportln("> Init "+getUtilityName+" Utility","debug");
  }
  
  def runOnReadPair(r1 : FastqLine, r2 : FastqLine, readNum: Int) : Option[(Int,Int)] = {
    val r1c = countGCs(r1);
    gcContentCountsR1(r1c) += 1;
    gcContentCountsRB(r1c) += 1;
    readLenCounts(0)(r1.readLen) += 1;
    
    val r1n = countNs(r1);
    val len1 = r1.readLen;
    gcContentCounts_byLength(0)(len1 - r1n)(r1c) += 1;
    nonMissCounts(0)(len1 - r1n) += 1;
    
    val r2c = if(! isSingleEnd){
      val r2c = countGCs(r2);
      val gcCount = r1c + r2c;
      val len2 = r2.readLen;
      val r2n = countNs(r2);

      gcContentCounts(gcCount) += 1;
      gcContentCountsR2(r2c) += 1;
      gcContentCountsRB(r2c) += 1;
      readLenCounts(1)(r2.readLen) += 1;
      nonMissCounts(1)(len2 - r2n) += 1;
      gcContentCounts_byLength(1)(len2 - r2n)(r2c) += 1;

      r2c
    } else -1;

    return Some((r1c,r2c));
  }
  
  def countNs(r : FastqLine) : Int = {
    r.seqVector.count( (b : Char) =>{
      (b == 'N');
    })
  }
  
  def countGCs(r : FastqLine) : Int = {
    r.seqVector.count( (b : Char) =>{
      (b == 'G' || b == 'C');
    })
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){

    if(true){
      val writer5 = createOutputFile(outfile, "FQ.gc.byRead.vsBaseCt.txt","",docWriter,
             ("NUM_BASES_GC","int","The number of GC bases"),
             ("NONMISSING_CT_<X>","int","The number of reads with <NUM_BASES_GC> G/C bases and <X> non-missing base-pairs.")
      );
      writer5.write("NUM_BASES_GC	"+(0 to readLen).map("NONMISSING_CT_"+ _).mkString("	")+"\n");
      for(gcCt <- 0 to readLen){
        writer5.write(gcCt+"	" + (0 to readLen).map((nonMissCt) => {
          gcContentCounts_byLength(0)(nonMissCt)(gcCt) + gcContentCounts_byLength(1)(nonMissCt)(gcCt)
        }).mkString("	")+"\n");
      }
      close(writer5);
    }
    if(true){
      val writer = createOutputFile(outfile , "FQ.readLenDist.txt","",docWriter,
             ("LEN","int","The read length."),
             ("CT_R1","int","The number of read-1 reads with read length equal to <LEN>"),
             ("CT_R2","int","The number of read-2 reads with read length equal to <LEN>")
      );
      writer.write("LEN\tCT_R1\tCT_R2\n");
      for(i <- 0 until readLenCounts(0).length){
        writer.write(i+"\t"+readLenCounts(0)(i)+"\t"+readLenCounts(1)(i)+"\n");
      }
      close(writer);
    }
    
    if(true){
      val readCt = gcContentCounts_byLength(0).map(_.sum).sum;
      
        val TEST_VALUE_1 = gcContentCounts_byLength(0).zipWithIndex.tail.map{ case (countArray,nonMissCt) => {
          countArray.zipWithIndex.map{ case (ct,gcCt) =>{
            ct
          }}.sum
        }}.sum.toDouble

      
        val avgGC_r1 = gcContentCounts_byLength(0).zipWithIndex.tail.map{ case (countArray,nonMissCt) => {
          countArray.zipWithIndex.map{ case (ct,gcCt) =>{
            (gcCt.toDouble / nonMissCt.toDouble) * ct.toDouble
          }}.sum
        }}.sum.toDouble / readCt.toDouble;
      
        //report("TEST_VALUE_1 = "+TEST_VALUE_1,"debug");
        //report("readCt = "+readCt,"debug");
        //report("avgGC_r1 = "+avgGC_r1,"debug");
        
      if(! isSingleEnd){
        summaryWriter.write("FQ_AVG_GC_R1\t"+avgGC_r1+"\t"+"Average of the GC fraction for read 1"+"\n");
        
        val avgGC_r2 = gcContentCounts_byLength(1).zipWithIndex.tail.map{ case (countArray,nonMissCt) => {
          countArray.zipWithIndex.map{ case (ct,gcCt) =>{
            (gcCt.toDouble / nonMissCt.toDouble) * ct.toDouble
          }}.sum
        }}.sum.toDouble / readCt.toDouble;
        summaryWriter.write("FQ_AVG_GC_R2\t"+avgGC_r2+"\t"+"Average of the GC fraction for read 2"+"\n");
        
        val avgGC = gcContentCounts_byLength.map(_.zipWithIndex.tail.map{ case (countArray,nonMissCt) => {
          countArray.zipWithIndex.map{ case (ct,gcCt) =>{
            (gcCt.toDouble / nonMissCt.toDouble) * ct.toDouble
          }}.sum
        }}.sum).sum.toDouble / (2.toDouble * readCt.toDouble);
        summaryWriter.write("FQ_AVG_GC\t"+avgGC+"\t"+"Average for the GC fraction across all reads"+"\n");
        
      } else {
        summaryWriter.write("FQ_AVG_GC\t"+avgGC_r1+"\t"+"Average for the GC fraction across all reads"+"\n");
      }
    }
    
    if(! isSingleEnd){
      val writer = createOutputFile(outfile, "FQ.gc.byPair.txt","",docWriter,
             ("NUM_BASES_GC","int","The number of GC bases."),
             ("CT","int","The number of read-pairs with <NUM_BASES_GC> bases that are G/C.")
      );
      writer.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCounts.length){
        writer.write(i + "	"+gcContentCounts(i)+"\n");
      } 
      close(writer); 
      
      val writer2 = createOutputFile(outfile , "FQ.gc.R1.txt","",docWriter,
             ("NUM_BASES_GC","int","The number of GC bases."),
             ("CT","int","The number of reads with <NUM_BASES_GC> bases that are G/C.")
      );
      writer2.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR1.length){
        writer2.write(i + "	"+gcContentCountsR1(i)+"\n");
      } 
      close(writer2);

      val writer3 = createOutputFile(outfile, "FQ.gc.R2.txt","",docWriter,
            ("NUM_BASES_GC","int","The number of GC bases."),
            ("CT","int","The number of reads with <NUM_BASES_GC> bases that are G/C.")
      );
      writer3.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR2.length){
        writer3.write(i + "	"+gcContentCountsR2(i)+"\n");
      }
      close(writer3);

      val writer4 = createOutputFile(outfile , "FQ.gc.byRead.txt","",docWriter,
             ("NUM_BASES_GC","int","The number of GC bases."),
             ("CT","int","The number of reads with <NUM_BASES_GC> bases that are G/C.")
      );
      writer4.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsRB.length){
        writer4.write(i + "	"+gcContentCountsRB(i)+"\n");
      }
      close(writer4);

    } else {
      val writer2 = createOutputFile(outfile, "FQ.gc.byRead.txt","",docWriter,
             ("NUM_BASES_GC","int","The number of GC bases."),
             ("CT","int","The number of reads with <NUM_BASES_GC> bases that are G/C.")
      );
      writer2.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR1.length){
        writer2.write(i + "	"+gcContentCountsR1(i)+"\n");
      } 
      close(writer2);
    }
  }

  def getUtilityName : String = "FastqGC";
}













