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
  
  val gcContentCountsRB_byLength : Array[Array[Array[Int]]] =  Array.ofDim[Int](2, readLen + 1,readLen + 1);
  //val gcContentCountsR1_byLength : Array[Array[Int]] =  Array.ofDim[Int](readLen + 1,readLen + 1);
  //val gcContentCountsR2_byLength : Array[Array[Int]] =  Array.ofDim[Int](readLen + 1,readLen + 1);
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
    
    gcContentCountsRB_byLength(0)(len1 - r1n)(r1c) += 1;
    gcContentCountsRB_byLength(1)(len2 - r2n)(r2c) += 1;
    
    gcContentCounts_byLength(len1 - r1n)(r1c) += 1;
    gcContentCounts_byLength(len2 - r2n)(r2c) += 1;
    
    return (r1c,r2c);
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    
    if(true){
      val readCt = gcContentCountsRB_byLength(0).map(_.sum).sum;
      
        val avgGC_r1 = gcContentCountsRB_byLength(0).zipWithIndex.tail.map{ case (countArray,nonMissCt) => {
          countArray.zipWithIndex.map{ case (ct,gcCt) =>{
            (gcCt.toDouble / nonMissCt.toDouble) * ct.toDouble
          }}.sum
        }}.sum / readCt.toDouble;
      
      
      if(! isSingleEnd){
        summaryWriter.write("AVG_GC_R1\t"+avgGC_r1+"\tAverage GC fraction for read 1"+"\n");
        
        val avgGC_r2 = gcContentCountsRB_byLength(1).zipWithIndex.tail.map{ case (countArray,nonMissCt) => {
          countArray.zipWithIndex.map{ case (ct,gcCt) =>{
            (gcCt.toDouble / nonMissCt.toDouble) * ct.toDouble
          }}.sum
        }}.sum.toDouble / readCt.toDouble;
        summaryWriter.write("AVG_GC_R2\t"+avgGC_r2+"\tAverage GC fraction for read 2"+"\n");
        
        val avgGC = gcContentCountsRB_byLength.map(_.zipWithIndex.tail.map{ case (countArray,nonMissCt) => {
          countArray.zipWithIndex.map{ case (ct,gcCt) =>{
            (gcCt.toDouble / nonMissCt.toDouble) * ct.toDouble
          }}.sum
        }}.sum).sum.toDouble / (2.toDouble * readCt.toDouble);
        summaryWriter.write("AVG_GC\t"+avgGC+"\tAverage GC fraction across all reads"+"\n");
      } else {
        summaryWriter.write("AVG_GC\t"+avgGC_r1+"\tAverage GC fraction across all reads"+"\n");
      }
    }
    
    if(! isSingleEnd){
      val writer = createOutputFile(outfile , "gc.byPair.txt","A table that counts the number of G/C bases observed in the read-pairs.",docWriter,
             ("NUM_BASES_GC","int","Number of GC bases."),
             ("CT","int","Number of read-pairs with the given number of GC bases.")
      );
      writer.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCounts.length){
        writer.write(i + "	"+gcContentCounts(i)+"\n");
      } 
      close(writer);

      val writer2 = createOutputFile(outfile , "gc.R1.txt","A table that counts the number of G/C bases observed in read 1",docWriter,
             ("NUM_BASES_GC","int","Number of GC bases in read 1."),
             ("CT","int","Number of first-read reads with the given number of GC bases.")
      );
      writer2.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR1.length){
        writer2.write(i + "	"+gcContentCountsR1(i)+"\n");
      } 
      close(writer2);

      val writer3 = createOutputFile(outfile , "gc.R2.txt","A table that counts the number of G/C bases observed in read 2",docWriter,
             ("NUM_BASES_GC","int","Number of GC bases in read 2."),
             ("CT","int","Number of 2nd-read reads with the given number of GC bases.")
      );
      writer3.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR2.length){
        writer3.write(i + "	"+gcContentCountsR2(i)+"\n");
      }
      close(writer3);
    
      val writer4 = createOutputFile(outfile ,"gc.byRead.txt","A table that counts the number of G/C bases observed in the reads, counting read 1 and read 2 independantly.",docWriter,
             ("NUM_BASES_GC","int","Number of GC bases."),
             ("CT","int","Number of reads with the given number of GC bases.")
      );
      writer4.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsRB.length){
        writer4.write(i + "	"+gcContentCountsRB(i)+"\n");
      }
      close(writer4);
      
 
      val writer5 = createOutputFile(outfile , "gc.byRead.vsBaseCt.txt","A table that counts the number of G/C bases versus the number of bases that are not N's.",docWriter,
             ("NUM_BASES_GC","int","Number of GC bases."),
             ("NONMISSING_CT_<X>","int","Number of reads with the given number of GC bases, restricted to reads with with <X> nonmissing base-pairs.")
      );
      writer5.write("NUM_BASES_GC	"+(1 to readLen).map("NONMISSING_CT_"+ _).mkString("	")+"\n");
      for(gcCt <- 0 to readLen){
        writer5.write(gcCt+"\t" + (1 to readLen).map(gcContentCounts_byLength(_)(gcCt)).mkString("\t")+"\n");
      }
      close(writer5);
      
    } else {
      val writer2 = createOutputFile(outfile , "gc.byRead.txt","A table that counts the number of G/C bases observed in the reads.",docWriter,
             ("NUM_BASES_GC","int","Number of GC bases."),
             ("CT","int","Number of reads with the given number of GC bases.")
      );
      writer2.write("NUM_BASES_GC	CT\n");
      for(i <- 0 until gcContentCountsR1.length){
        writer2.write(i + "	"+gcContentCountsR1(i)+"\n");
      } 
      close(writer2);
      
      val writer5 = createOutputFile(outfile , "gc.byRead.vsBaseCt.txt","A table that counts the number of G/C bases versus the number of bases that are not N's.",docWriter,
             ("NUM_BASES_GC","int","Number of GC bases."),
             ("NONMISSING_CT_<X>","int","Number of reads with the given number of GC bases, restricted to reads with with <X> nonmissing base-pairs.") 
      );
      writer5.write("NUM_BASES_GC	"+(0 to readLen).map("NONMISSING_CT_"+ _).mkString("	")+"\n");
      for(gcCt <- 0 to readLen){
        writer5.write(gcCt+"\t" + (0 to readLen).map(gcContentCountsRB_byLength(0)(_)(gcCt)).mkString("\t")+"\n");
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