package qcUtils

import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
import java.io.OutputStream
import java.io.FileOutputStream
import java.io.InputStream
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.File
import scala.collection.JavaConversions._
import scala.collection.mutable.HashMap;

import net.sf.samtools._

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;


object qcJunctionDistribution {
   
/*
 * THIS WHOLE OBJECT IS DEPRECIATED!!!!!!!!
 */
  def run(infile : String, outfile : String){
    /*
    val (readLength, recordIter) = initSamRecordIterator(infile);
    val pairedIter : Iterator[(SAMRecord,SAMRecord)] = samRecordPairIterator(recordIter);
    
    reportln("Autodetected read length: " + readLength,"note");
    
    val qcJD : qcJunctionDistribution = new qcJunctionDistribution(readLength);
    
    val numberedPairedIter : Iterator[((SAMRecord,SAMRecord), Int)] = zipIteratorWithCount(pairedIter);
    
    numberedPairedIter.map((next : ((SAMRecord,SAMRecord), Int)) => {
      qcJD.runOnReadPair(next._1._1, next._1._2, next._2);
    })
    
    qcJD.writeOutput(outfile);
    */
    
    /*val causeOfDropArray : Array[Int] = Array.ofDim[Int](5);
    
    val positionOpDistribution1 : Array[HashMap[(CigarOperator,Int),Int]] = new Array[HashMap[(CigarOperator,Int),Int]](readLength);
    for(i <- 0 until positionOpDistribution1.length){
      positionOpDistribution1(i) = new HashMap[(CigarOperator,Int),Int](){ override def default(key : (CigarOperator,Int)) = 0 };
    }
    val positionOpDistribution2 : Array[HashMap[(CigarOperator,Int),Int]] = new Array[HashMap[(CigarOperator,Int),Int]](readLength);
    for(i <- 0 until positionOpDistribution2.length){
      positionOpDistribution2(i) = new HashMap[(CigarOperator,Int),Int](){ override def default(key : (CigarOperator,Int)) = 0 };
    }
    val opLengthDistribution1 : HashMap[(CigarOperator,Int),Int] = new HashMap[(CigarOperator,Int),Int](){ override def default(key : (CigarOperator,Int)) = 0 };
    val opLengthDistribution2 : HashMap[(CigarOperator,Int),Int] = new HashMap[(CigarOperator,Int),Int](){ override def default(key : (CigarOperator,Int)) = 0 };
    
    //val 
    //new HashMap[Int,Int](){ override def default(key:Int) = 0 }
    
    var readCt = 1;
    for((r1 : SAMRecord, r2 : SAMRecord) <- pairedInput){
      val chromName = r1.getReferenceName();
      if(r1.getReadName() != r2.getReadName()){
        error("FATAL ERROR: Paired reads do not match, is dataset sorted by read name?");
      }
      if(r1.getReadLength > readLength || r2.getReadLength > readLength){
        error("FATAL ERROR: Failed to autodetect max read length!\nRead Pair #"+readCt + " has lengths: " + r1.getReadLength +" and " + r2.getReadLength + "!");
      }
      
      //val strand =  getStrand(r1);
      
      if(useRead(r1, causeOfDropArray) & useRead(r2,causeOfDropArray)){
        val r1c = r1.getCigar();
        val r2c = r2.getCigar();
        
        readCigar(r1c,positionOpDistribution1, opLengthDistribution1, r1.getReadNegativeStrandFlag());
        readCigar(r2c,positionOpDistribution2, opLengthDistribution2, r2.getReadNegativeStrandFlag());
      }
      readCt += 1;
    }
    
    
    writePositionOpDistribution( positionOpDistribution1 , outfile + ".cigarOpDistribution.byReadCycle.R1.txt");
    writePositionOpDistribution( positionOpDistribution2 , outfile + ".cigarOpDistribution.byReadCycle.R2.txt");

    writeOpLengthDistribution(opLengthDistribution1, outfile + ".cigarOpLengths.byOp.R1.txt");
    writeOpLengthDistribution(opLengthDistribution2, outfile + ".cigarOpLengths.byOp.R2.txt");*/
  }
  
  val cigarOpList : List[CigarOperator] = List(CigarOperator.DELETION, CigarOperator.HARD_CLIP, CigarOperator.INSERTION, CigarOperator.MATCH_OR_MISMATCH, CigarOperator.PADDING, CigarOperator.SKIPPED_REGION, CigarOperator.SOFT_CLIP);
  def writePositionOpDistribution(pod : Array[HashMap[(CigarOperator,Int),Int]], outfile : String){
    val writer = openWriterSmart_viaGlobalParam(outfile);
    
    writer.write("CYCLE");
    for(op <- cigarOpList){
      writer.write("	"+op.toString + "_S	"+op.toString + "_M	"+op.toString +"_E	"+op.toString+"_B");
    }
    writer.write("\n");
    for(i <- 0 until pod.length){
      writer.write(i + "");
      for(op <- cigarOpList){
        writer.write( "	"+pod(i)((op,0)) + "	"+pod(i)((op,1)) + "	"+pod(i)((op,2))+ "	"+pod(i)((op,3))  );
      }
      writer.write("\n");
    }
    
    writer.close();
  }
  def writeOpLengthDistribution(old : HashMap[(CigarOperator,Int),Int], outfile : String){
    val writer = openWriterSmart_viaGlobalParam(outfile);
    
    writer.write("OP	LEN	CT\n");
    for(op <- cigarOpList){
       val keys = old.keys.toSeq.filter((x) => x._1 == op).sorted;
       
       for( (o,len) <- keys ){
         val ct = old( (o,len) );
         writer.write(op.toString() + "	"+len+"	"+ct+"\n");
       }
    }
    
    writer.close();
  }
  
  def readCigar(c : Cigar, cca : Array[HashMap[(CigarOperator,Int),Int]], tabOps : HashMap[(CigarOperator,Int),Int], reverse : Boolean){
    var readPos = 0;
    
    val elementList : List[CigarElement] = if(reverse) {
      c.getCigarElements().toList.reverse;
    } else {
      c.getCigarElements().toList;
    }
    
    for(ce : CigarElement <- elementList){
      val len : Int = ce.getLength();
      val op : CigarOperator = ce.getOperator();
      
      tabOps( (op,len) ) += 1;
      
      if(op.consumesReadBases()){
        if(len > 1 ){
           cca(readPos) ( (op, 0) ) += 1;
           for(rp <- readPos + 1 until readPos + len - 1){
             cca(rp) ( (op, 1) ) += 1;
           }
           cca(readPos + len - 1)( (op,2) ) += 1;
           readPos += len;
        } else {
          cca(readPos) ( (op, 3) ) += 1;
          readPos += 1;
        }
      } else {
        cca(readPos)( (op,0) ) += 1;
      }
    }
  }

}

class qcJunctionDistribution(readLen : Int)  extends QCUtility[Unit]  {
  val readLength : Int = readLen;
  val positionOpDistribution1 : Array[HashMap[(CigarOperator,Int),Int]] = new Array[HashMap[(CigarOperator,Int),Int]](readLength);
  for(i <- 0 until positionOpDistribution1.length){
    positionOpDistribution1(i) = new HashMap[(CigarOperator,Int),Int](){ override def default(key : (CigarOperator,Int)) = 0 };
  }
  val positionOpDistribution2 : Array[HashMap[(CigarOperator,Int),Int]] = new Array[HashMap[(CigarOperator,Int),Int]](readLength);
  for(i <- 0 until positionOpDistribution2.length){
    positionOpDistribution2(i) = new HashMap[(CigarOperator,Int),Int](){ override def default(key : (CigarOperator,Int)) = 0 };
  }
  val opLengthDistribution1 : HashMap[(CigarOperator,Int),Int] = new HashMap[(CigarOperator,Int),Int](){ override def default(key : (CigarOperator,Int)) = 0 };
  val opLengthDistribution2 : HashMap[(CigarOperator,Int),Int] = new HashMap[(CigarOperator,Int),Int](){ override def default(key : (CigarOperator,Int)) = 0 };
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
      //if(useRead(r1) & useRead(r2)){
        val r1c = r1.getCigar();
        val r2c = r2.getCigar();
        
        qcJunctionDistribution.readCigar(r1c,positionOpDistribution1, opLengthDistribution1, r1.getReadNegativeStrandFlag());
        qcJunctionDistribution.readCigar(r2c,positionOpDistribution2, opLengthDistribution2, r2.getReadNegativeStrandFlag());
      //}
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    qcJunctionDistribution.writePositionOpDistribution( positionOpDistribution1 , outfile + ".cigarOpDistribution.byReadCycle.R1.txt");
    qcJunctionDistribution.writePositionOpDistribution( positionOpDistribution2 , outfile + ".cigarOpDistribution.byReadCycle.R2.txt");

    qcJunctionDistribution.writeOpLengthDistribution(opLengthDistribution1, outfile + ".cigarOpLengths.byOp.R1.txt");
    qcJunctionDistribution.writeOpLengthDistribution(opLengthDistribution2, outfile + ".cigarOpLengths.byOp.R2.txt");
  }
}




































