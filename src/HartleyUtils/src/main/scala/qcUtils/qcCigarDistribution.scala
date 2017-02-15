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

import scala.collection.mutable.AnyRefMap;

object qcCigarDistribution {
   

  def run(infile : String, outfile : String){
    
  }
  
  val cigarOpList : Vector[CigarOperator] = Vector(CigarOperator.DELETION, CigarOperator.HARD_CLIP, CigarOperator.INSERTION, CigarOperator.MATCH_OR_MISMATCH, CigarOperator.PADDING, CigarOperator.SKIPPED_REGION, CigarOperator.SOFT_CLIP);

  
  def writePositionOpDistribution(pod : Array[scala.collection.mutable.Map[(CigarOperator,Int),Int]], writer : WriterUtil){
    //val writer = openWriterSmart_viaGlobalParam(outfile);
    
    writer.write("CYCLE");
    for(op <- cigarOpList){
      writer.write("	"+op.toString + "_S	"+op.toString + "_M	"+op.toString +"_E	"+op.toString+"_B");
    }
    writer.write("\n");
    for(i <- 0 until pod.length){
      writer.write((i) + "");
      for(op <- cigarOpList){
        writer.write( "	"+pod(i)((op,0)) + "	"+pod(i)((op,1)) + "	"+pod(i)((op,2))+ "	"+pod(i)((op,3))  );
      }
      writer.write("\n");
    }
    
    writer.close();
  }
  def writeOpLengthDistribution(old : scala.collection.mutable.Map[(CigarOperator,Int),Int], writer : WriterUtil){
    //val writer = openWriterSmart_viaGlobalParam(outfile);
    
    writer.write("OP\tLEN\tCT\n");
    for(op <- cigarOpList){
       val keys = ((old.keys.toSeq.filter((x) => x._1 == op)).toSet ++ Set((op,1))).toSeq.sorted;
       
       for( (o,len) <- keys ){
         val ct = old( (o,len) );
         writer.write(op.toString() + "	"+len+"	"+ct+"\n");
       }
    }
    
    writer.close();
  }
  
  def readCigar(c : Cigar, cca : Array[scala.collection.mutable.Map[(CigarOperator,Int),Int]], tabOps : scala.collection.mutable.Map[(CigarOperator,Int),Int], reverse : Boolean){
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

class qcCigarDistribution(isSingleEnd : Boolean, readLen : Int,variableReadLen : Boolean)  extends QCUtility[Unit]  {
  reportln("> Init CigarOpDistribution Utility","debug");

  val readLength : Int = readLen;
  val positionOpDistribution1 : Array[scala.collection.mutable.Map[(CigarOperator,Int),Int]] = new Array[scala.collection.mutable.Map[(CigarOperator,Int),Int]](readLength);
  for(i <- 0 until positionOpDistribution1.length){
    positionOpDistribution1(i) = AnyRefMap[(CigarOperator,Int),Int]().withDefault( x => 0);
  }
  val positionOpDistribution2 : Array[scala.collection.mutable.Map[(CigarOperator,Int),Int]] = new Array[scala.collection.mutable.Map[(CigarOperator,Int),Int]](readLength);
  for(i <- 0 until positionOpDistribution2.length){
    positionOpDistribution2(i) = AnyRefMap[(CigarOperator,Int),Int]().withDefault( x => 0);
  }
  val opLengthDistribution1 : scala.collection.mutable.Map[(CigarOperator,Int),Int] = new AnyRefMap[(CigarOperator,Int),Int]().withDefault(x => 0);
  val opLengthDistribution2 : scala.collection.mutable.Map[(CigarOperator,Int),Int] = new AnyRefMap[(CigarOperator,Int),Int]().withDefault(x => 0);
  
  val readLenDistribution1 : Array[Int] = Array.ofDim[Int](readLen+1);
  val readLenDistribution2 : Array[Int] = Array.ofDim[Int](readLen+1);
  
  
  
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
      //if(useRead(r1) & useRead(r2)){
        val r1c = r1.getCigar();
        qcCigarDistribution.readCigar(r1c,positionOpDistribution1, opLengthDistribution1, r1.getReadNegativeStrandFlag());
        if(variableReadLen) readLenDistribution1(r1.getReadLength()) += 1;
        
        if(! isSingleEnd){
          val r2c = r2.getCigar();
          qcCigarDistribution.readCigar(r2c,positionOpDistribution2, opLengthDistribution2, r2.getReadNegativeStrandFlag());
          if(variableReadLen) readLenDistribution2(r2.getReadLength()) += 1;
        }
      //}
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    
    //for(op <- cigarOpList){
    //  writer.write("	"+op.toString + "_S	"+op.toString + "_M	"+op.toString +"_E	"+op.toString+"_B");
    //}
    
          /*
          (List(("CYCLE","String","The read cycle (ie, position in the read)")) ++ qcCigarDistribution.cigarOpList.foldLeft(List()){ case (soFar,op) => 
            soFar ++ List[String]((op.toString+"_S","",""),
                                  (op.toString+"_M","",""),
                                  (op.toString+"_E","",""),
                                  (op.toString+"_B","","")
                                 ) 
          }):_**/
    
    if(true){
      val writer = createOutputFile(outfile,"cigarOpDistribution.byReadCycle.R1.txt","",docWriter,
          ("CYCLE","int","The read cycle (ie, position in the read)"),
          ("<X>_<Y>","int","The number of reads with cigar op <X> of type <Y> at position <CYCLE>. <X> indicates the cigar operation, <Y> indicates whether it is the start, end, middle, or a single-base cigar op."),
          ("<X>_S","int","The number of reads that have cigar op <X> START at position <CYCLE>."),
          ("<X>_M","int","The number of reads that have cigar op <X> that span position <CYCLE>."),
          ("<X>_E","int","The number of reads that have cigar op <X> END at position <CYCLE>."),
          ("<X>_B","int","The number of reads with a single-base-pair-long op <X> at position <CYCLE>"),
          ("D_<Y>","int","The 'D' cigar operation, indicating a DELETION"),
          ("H_<Y>","int","The 'H' cigar operation, indicating HARD CLIP"),
          ("I_<Y>","int","The 'I' cigar operation, indicating INSERTION"),
          ("M_<Y>","int","The 'M' cigar operation, indicating ALIGNED REGION (match or mismatch)"),
          ("P_<Y>","int","The 'P' cigar operation, indicating PADDING"),
          ("N_<Y>","int","The 'N' cigar operation, indicating SKIPPED REGION"),
          ("S_<Y>","int","The 'H' cigar operation, indicating SOFT CLIP")
      );
      qcCigarDistribution.writePositionOpDistribution( positionOpDistribution1 , writer);
    }
    if(true){
      val writer = createOutputFile(outfile,"cigarOpLengths.byOp.R1.txt","",docWriter,
          ("OP","char","The cigar operation"),
          ("LEN","int","The length"),
          ("CT","int","The number of occurances of a cigar operation of type <OP> and length <LEN>.")
      );
      qcCigarDistribution.writeOpLengthDistribution(opLengthDistribution1, writer);
    }
    
    if(! isSingleEnd){
      if(true){
        val writer = createOutputFile(outfile,"cigarOpDistribution.byReadCycle.R2.txt","",docWriter,
          ("CYCLE","int","The read cycle (ie, position in the read)"),
          ("<X>_<Y>","int","The number of reads with cigar op <X> of type <Y> at position <CYCLE>. <X> indicates the cigar operation, <Y> indicates whether it is the start, end, middle, or a single-base cigar op."),
          ("<X>_S","int","The number of reads that have cigar op <X> START at position <CYCLE>."),
          ("<X>_M","int","The number of reads that have cigar op <X> that span position <CYCLE>."),
          ("<X>_E","int","The number of reads that have cigar op <X> END at position <CYCLE>."),
          ("<X>_B","int","The number of reads with a single-base-pair-long op <X> at position <CYCLE>"),
          ("D_<Y>","int","The 'D' cigar operation, indicating a DELETION"),
          ("H_<Y>","int","The 'H' cigar operation, indicating HARD CLIP"),
          ("I_<Y>","int","The 'I' cigar operation, indicating INSERTION"),
          ("M_<Y>","int","The 'M' cigar operation, indicating ALIGNED REGION (match or mismatch)"),
          ("P_<Y>","int","The 'P' cigar operation, indicating PADDING"),
          ("N_<Y>","int","The 'N' cigar operation, indicating SKIPPED REGION"),
          ("S_<Y>","int","The 'H' cigar operation, indicating SOFT CLIP")
        );
        qcCigarDistribution.writePositionOpDistribution( positionOpDistribution2 , writer);
      }
      if(true){
        val writer = createOutputFile(outfile,"cigarOpLengths.byOp.R2.txt","",docWriter,
          ("OP","char","The cigar operation"),
          ("LEN","int","The length"),
          ("CT","int","The number of occurances of a cigar operation of type <OP> and length <LEN>.")
        );
        qcCigarDistribution.writeOpLengthDistribution(opLengthDistribution2, writer);
      }
      //qcCigarDistribution.writePositionOpDistribution( positionOpDistribution2 , outfile + ".cigarOpDistribution.byReadCycle.R2.txt");
      //qcCigarDistribution.writeOpLengthDistribution(   opLengthDistribution2,    outfile + ".cigarOpLengths.byOp.R2.txt");
    }
  }
  
  def getUtilityName : String = "CigarOpDistribution";

}




































