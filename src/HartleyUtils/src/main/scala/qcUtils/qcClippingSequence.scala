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

object qcClippingSequence {
/*
 * THIS ENTIRE OBJECT IS DEPRECIATED!
 */
  
  val CLIPPING_SEQ_BYTE = 'S'.toByte;
  val MISSING_SEQ_BYTE = 'H'.toByte;
  val ALL_BASE_BYTES : Seq[Byte] = Seq('A'.toByte,'C'.toByte,'G'.toByte,'T'.toByte,'N'.toByte);
  
  def writeTabulatedSequence(outfile : String, maxLeadClip : Int, maxTailClip : Int, nvcCounter : NvcCounter, readLengthSet : Set[Int], writeClippedNVC : Boolean){
    writeRawClippingData(true  , outfile , nvcCounter );
    writeRawClippingData(false , outfile , nvcCounter );
    /*var writer = openWriterSmart_viaGlobalParam(outfile + ".ignoreClipping.NVC.raw.txt");
    writer.write("readPos	readNum	base	CT_Aligned_to_Fwd	CT_Aligned_to_Rev	CT\n");
    
    for(i <- 0 until nvcCounter.readLength){
      for(b <- ALL_BASE_BYTES){
        writer.write("" + i + "	1	"+b.toChar+"	"+nvcCounter.get_rawNvcCount(true,true,i,b)+"	"+nvcCounter.get_rawNvcCount(true,false,i,b)+"	"+(nvcCounter.get_rawNvcCount(true,false,i,b) + nvcCounter.get_rawNvcCount(true,false,i,b) )+"\n");
        writer.write("" + i + "	2	"+b.toChar+"	"+nvcCounter.get_rawNvcCount(false,true,i,b)+"	"+nvcCounter.get_rawNvcCount(false,false,i,b)+"	"+(nvcCounter.get_rawNvcCount(false,false,i,b) + nvcCounter.get_rawNvcCount(false,false,i,b) )+"\n");
      }
    }
    writer.close();*/
    
    //writer = openWriter(outfile + ".NVC.raw.txt");
    //writer.write("readPos	readNum	leadClipLen	tailClipLen	base	CT_Aligned_to_Fwd	CT_Aligned_to_Rev	CT\n");
    //for(i <- 0 until nvcCounter.readLength){
    //  for(leadClipLen <- 0 to maxLeadClip){
    //    for(tailClipLen <- 0 to maxTailClip){
    //      for(b <- ALL_BASE_BYTES){
    //        
    //        writer.write("" + i + "	1	"+leadClipLen+"	"+tailClipLen+"	"+b.toChar+"	"+
    //             nvcCounter.get_fullNvcCount(true,true,leadClipLen,tailClipLen,i,b)  + "	" + 
    //             nvcCounter.get_fullNvcCount(true,false,leadClipLen,tailClipLen,i,b) + "	" +
    //            (nvcCounter.get_fullNvcCount(true,true,leadClipLen,tailClipLen,i,b)  +  nvcCounter.get_fullNvcCount(true,false,leadClipLen,tailClipLen,i,b))+"\n");
    //        
    //        writer.write("" + i + "	2	"+leadClipLen+"	"+tailClipLen+"	"+b.toChar+"	"+
    //             nvcCounter.get_fullNvcCount(false,true,leadClipLen,tailClipLen,i,b)  + "	" + 
    //             nvcCounter.get_fullNvcCount(false,false,leadClipLen,tailClipLen,i,b) + "	" +
    //            (nvcCounter.get_fullNvcCount(false,true,leadClipLen,tailClipLen,i,b)  +  nvcCounter.get_fullNvcCount(false,false,leadClipLen,tailClipLen,i,b))+"\n");
    //      }
    //    }
    //  }
    //}
    //writer.close();
    
    val readLengthList = readLengthSet.toList.sorted;
    
    if(writeClippedNVC){
      writeLeadClippingData(true,outfile , maxLeadClip , maxTailClip , nvcCounter , readLengthList );
      writeLeadClippingData(false,outfile , maxLeadClip , maxTailClip , nvcCounter , readLengthList );
      
      writeTailClippingData(true,outfile , maxLeadClip , maxTailClip , nvcCounter , readLengthList );
      writeTailClippingData(false,outfile , maxLeadClip , maxTailClip , nvcCounter , readLengthList );
    }
    writeNoClippingData(true,outfile , maxLeadClip , maxTailClip , nvcCounter , readLengthList );
    writeNoClippingData(false,outfile , maxLeadClip , maxTailClip , nvcCounter , readLengthList );

    //get_clippingNvcCount(isLeadingClip : Boolean, isFirstRead : Boolean, isFwdRefStrand: Boolean, clipAmt : Int, readLength : Int, pos : Int, base : Byte)

  }
  
  def writeRawClippingData(isFirstRead : Boolean, outfile : String, nvcCounter : NvcCounter){
    val readName = if(isFirstRead) "R1" else "R2";
    val writer =  openWriterSmart_viaGlobalParam(outfile + ".NVC.raw."+readName+".txt");
    writer.write("readPos	base	CT_Aligned_to_Fwd	CT_Aligned_to_Rev	CT\n");
    
    for(i <- 0 until nvcCounter.readLength){
      for(b <- ALL_BASE_BYTES){
        writer.write("" + i + "	"+b.toChar+"	"+nvcCounter.get_rawNvcCount(isFirstRead,true,i,b)+"	"+nvcCounter.get_rawNvcCount(isFirstRead,false,i,b)+"	"+(nvcCounter.get_rawNvcCount(isFirstRead,false,i,b) + nvcCounter.get_rawNvcCount(isFirstRead,false,i,b) )+"\n");
      }
    }
    close(writer);
  }
  
  def writeNoClippingData(isFirstRead : Boolean, outfile : String, maxLeadClip : Int, maxTailClip : Int, nvcCounter : NvcCounter, readLengthList : List[Int]){
    val readName = if(isFirstRead) "R1" else "R2";
    val writer =  openWriterSmart_viaGlobalParam(outfile + ".NVC.minus.clipping."+readName+".txt");
    
    //get_noClipNvcCount(isFirstRead : Boolean, isFwdRefStrand: Boolean, pos : Int, base : Byte)
    writer.write("readPos	base	CT_Aligned_to_Fwd	CT_Aligned_to_Rev	CT\n");
    for(readPos <- 0 until nvcCounter.readLength){
      for(b <- ALL_BASE_BYTES){
        val fwdCt = nvcCounter.get_noClipNvcCount(isFirstRead,true, readPos, b);
        val revCt = nvcCounter.get_noClipNvcCount(isFirstRead,false, readPos, b);
        writer.write(""+readPos+"	"+b.toChar+"	"+fwdCt+"	"+revCt+"	"+(fwdCt+revCt)+"\n");
      }
      val b = CLIPPING_SEQ_BYTE;
      val fwdCt = nvcCounter.get_noClipNvcCount(isFirstRead,true, readPos, b);
      val revCt = nvcCounter.get_noClipNvcCount(isFirstRead,false, readPos, b);
      writer.write(""+readPos+"	"+b.toChar+"	"+fwdCt+"	"+revCt+"	"+(fwdCt+revCt)+"\n");
    }
    close(writer);
  }
  
  def writeLeadClippingData(isFirstRead : Boolean, outfile : String, maxLeadClip : Int, maxTailClip : Int, nvcCounter : NvcCounter, readLengthList : List[Int]){
    val readName = if(isFirstRead) "R1" else "R2";
    
    val writer =  openWriterSmart_viaGlobalParam(outfile + ".NVC.lead.clip."+readName+".txt");
    
    writer.write("readPos	leadClipLen	readLen	base	CT_Aligned_to_Fwd	CT_Aligned_to_Rev	CT\n");
    for(leadClipLen <- 1 to maxLeadClip){
      for(pos <- 0 until leadClipLen){
        for(readLen <- readLengthList){
          for(base <- ALL_BASE_BYTES){
            val fwdCt = nvcCounter.get_clippingNvcCount(true, isFirstRead, true,  leadClipLen, readLen, pos, base);
            val revCt = nvcCounter.get_clippingNvcCount(true, isFirstRead, false, leadClipLen, readLen, pos, base);
            writer.write(pos+"	"+leadClipLen+"	"+readLen+"	"+base.toChar+"	"+fwdCt+"	"+revCt+"	"+(fwdCt + revCt)+"\n");
          }
        }
      }
    }
    close(writer);
  }
  
  def writeTailClippingData(isFirstRead : Boolean, outfile : String, maxLeadClip : Int, maxTailClip : Int, nvcCounter : NvcCounter, readLengthList : List[Int]){
    val readName = if(isFirstRead) "R1" else "R2";
    
    val writer =  openWriterSmart_viaGlobalParam(outfile + ".NVC.tail.clip."+readName+".txt");
    
    writer.write("readPos	tailClipLen	readLen	base	CT_Aligned_to_Fwd	CT_Aligned_to_Rev	CT\n");
    for(tailClipLen <- 1 to maxTailClip){
      for(readLen <- readLengthList){
        for(pos <- (readLen - tailClipLen) until readLen){
          for(base <- ALL_BASE_BYTES){
            val fwdCt = nvcCounter.get_clippingNvcCount(false, isFirstRead, true,  tailClipLen, readLen, pos, base);
            val revCt = nvcCounter.get_clippingNvcCount(false, isFirstRead, false, tailClipLen, readLen, pos, base);
            writer.write(pos+"	"+tailClipLen+"	"+readLen+"	"+base.toChar+"	"+fwdCt+"	"+revCt+"	"+(fwdCt + revCt)+"\n");
          }
        }
      }
    }
    close(writer);
  }
  
  
  //def calcNvc(cigarSeq : Seq[CigarElement], seqArray : Seq[Byte], reversed : Boolean, )
  
  class NvcCounter(readLn : Int){
    val readLength : Int = readLn;
    
    val NVC_INTERNAL_COUNTER : HashMap[(Boolean, Boolean, Int, Int, Int, Byte), Int] = new HashMap[(Boolean, Boolean, Int, Int, Int, Byte), Int](){
      override def default(key : (Boolean, Boolean, Int, Int, Int, Byte)) = 0;
    }
    
    def get_fullNvcCount(isFirstRead : Boolean, isFwdRefStrand : Boolean, leadClipLen : Int, tailClipLen : Int, pos : Int, base : Byte) : Int = {
      return NVC_INTERNAL_COUNTER((isFirstRead, isFwdRefStrand, leadClipLen, tailClipLen, pos, base));
    }
    def count_fullNvcCount(isFirstRead : Boolean, isFwdRefStrand : Boolean, leadClipLen : Int, tailClipLen : Int, pos : Int, base : Byte) {
      NVC_INTERNAL_COUNTER((isFirstRead, isFwdRefStrand, leadClipLen, tailClipLen, pos, base)) += 1;
    }
    
    val rawNvcCounter : HashMap[(Boolean, Boolean, Int, Byte), Int] = new HashMap[(Boolean, Boolean, Int, Byte), Int](){
      override def default(key : (Boolean, Boolean, Int, Byte)) = 0;
    }
    def get_rawNvcCount(isFirstRead : Boolean, isFwdRefStrand: Boolean, pos : Int, base : Byte) : Int = {
      return rawNvcCounter(isFirstRead, isFwdRefStrand, pos, base);
    }
    def count_rawNvcCount(isFirstRead : Boolean, isFwdRefStrand: Boolean, pos : Int, base : Byte){
      rawNvcCounter(isFirstRead, isFwdRefStrand, pos, base) += 1;
    }
    
    
    val noClipNvcCounter : HashMap[(Boolean, Boolean, Int,  Byte), Int] = new HashMap[(Boolean, Boolean, Int,  Byte), Int](){
      override def default(key : (Boolean, Boolean, Int, Byte)) = 0;
    }
    val clippingNvcCounter : HashMap[(Boolean, Boolean, Boolean, Int, Int, Int, Byte), Int] = new HashMap[(Boolean, Boolean, Boolean, Int, Int, Int, Byte), Int](){
      override def default(key : (Boolean, Boolean, Boolean, Int, Int, Int, Byte)) = 0;
    }
    def get_noClipNvcCount(isFirstRead : Boolean, isFwdRefStrand: Boolean, pos : Int, base : Byte) : Int = {
      return noClipNvcCounter(isFirstRead, isFwdRefStrand, pos, base);
    }
    def count_noClipNvcCount(isFirstRead : Boolean, isFwdRefStrand: Boolean, pos : Int, base : Byte){
      noClipNvcCounter(isFirstRead, isFwdRefStrand, pos, base) += 1;
    }
    
    def get_clippingNvcCount(isLeadingClip : Boolean, isFirstRead : Boolean, isFwdRefStrand: Boolean, clipAmt : Int, readLength : Int, pos : Int, base : Byte) : Int = {
      return clippingNvcCounter(isLeadingClip, isFirstRead, isFwdRefStrand, clipAmt, readLength, pos, base);
    }
    def count_clippingNvcCount(isLeadingClip : Boolean, isFirstRead : Boolean, isFwdRefStrand: Boolean, clipAmt : Int, readLength : Int, pos : Int, base : Byte){
      clippingNvcCounter(isLeadingClip, isFirstRead, isFwdRefStrand, clipAmt, readLength, pos, base) += 1;
    }
  }
  
  def getClipLengths(r : SAMRecord) : (Int, Int) = {
    val cigarSeq = getCigarSeq(r, r.getReadNegativeStrandFlag());
    
    val leadClipLen = if(cigarSeq.head.getOperator() == CigarOperator.SOFT_CLIP) { cigarSeq.head.getLength() } else 0;
    val tailClipLen = if(cigarSeq.last.getOperator() == CigarOperator.SOFT_CLIP) { cigarSeq.last.getLength() } else 0;
    
    return (leadClipLen, tailClipLen);
  }
  
  def getBasePairs(r : SAMRecord) : Seq[Byte] = {
    if(r.getReadNegativeStrandFlag()){
      return r.getReadBases().toSeq.reverse.map(reverseBaseByteMap(_));
    } else {
      return r.getReadBases().toSeq;
    }
  }
  def countAll(r : SAMRecord, isFirstRead : Boolean, nvcCounter : NvcCounter) : (Int, Int, Int) = {
    
    val (leadClipLen, tailClipLen) = getClipLengths(r);
    val basePairs = getBasePairs(r);
    val isFwdRefStrand = ! r.getReadNegativeStrandFlag();
    val readLength = basePairs.length;
    
    (0 until basePairs.length).map( (i : Int) => nvcCounter.count_rawNvcCount(isFirstRead, isFwdRefStrand, i, basePairs(i)) );
    //(0 until basePairs.length).map( (i : Int) => nvcCounter.count_rawNvcCount(isFirstRead, isFwdRefStrand, i, basePairs(i)) );    
    
    //note:
    //count_nvcCount(isFirstRead : Boolean, isFwdRefStrand : Boolean, leadClipLen : Int, tailClipLen : Int, pos : Int, base : Byte)
    //count_noClipNvcCount(isFirstRead : Boolean, isFwdRefStrand: Boolean, pos : Int, base : Byte)
    //count_clippingNvcCount(isLeadingClip : Boolean, isFirstRead : Boolean, isFwdRefStrand: Boolean, clipAmt : Int, readLength : Int, pos : Int, base : Byte)
    for(i : Int <- 0 until readLength){
      //nvcCounter.count_fullNvcCount(isFirstRead, isFwdRefStrand, leadClipLen, tailClipLen, i, basePairs(i));
      
      if(leadClipLen <= i && i < (readLength - tailClipLen)){
        nvcCounter.count_noClipNvcCount(isFirstRead, isFwdRefStrand, i, basePairs(i));
      } else {
        nvcCounter.count_noClipNvcCount(isFirstRead, isFwdRefStrand, i, CLIPPING_SEQ_BYTE);
      }
    }
    
    for(i : Int <- 0 until leadClipLen){
      nvcCounter.count_clippingNvcCount(true, isFirstRead, isFwdRefStrand, leadClipLen, readLength, i, basePairs(i));
    }
    
    for(i : Int <- (readLength - tailClipLen) until readLength){
      nvcCounter.count_clippingNvcCount(false, isFirstRead, isFwdRefStrand, tailClipLen, readLength, i, basePairs(i));
    }
    
    return (leadClipLen, tailClipLen, readLength);
  }

  def writeTabulatedClippedSequence(outfile : String, maxClip : Int, clippingNvcMap : HashMap[(Boolean,Boolean,Int,Int,Char), Int]){
    val writer = openWriterSmart_viaGlobalParam(outfile + ".clippedSeqNvc.5prime.txt");
    writer.write("CLIP_SIZE	POSITION	BASE	CT_F	CT_R	CT\n");
    for(i <- 1 to maxClip){
      for(j <- 1 to i){
        for(B <- List('A','T','C','G','N')){
          writer.write("" + i + "	"+j+"	"+B+"	"+ clippingNvcMap((true,true,i,j,B)) +"	"+ clippingNvcMap((false,true,i,j,B)) +"	"+ ( clippingNvcMap((true,true,i,j,B)) + clippingNvcMap((false,true,i,j,B)) ) + "\n");
        }
      }
    }
    close(writer);
    
    val writer2 = openWriterSmart_viaGlobalParam(outfile + ".clippedSeqNvc.3prime.txt");
    writer2.write("CLIP_SIZE	POSITION	BASE	CT_F	CT_R	CT\n");
    for(i <- 1 to maxClip){
      for(j <- 1 to i){
        for(B <- List('A','T','C','G','N')){
          writer2.write("" + i + "	"+j+"	"+B+"	"+ clippingNvcMap((true,false,i,j,B)) +"	"+ clippingNvcMap((false,false,i,j,B)) +"	"+ ( clippingNvcMap((true,false,i,j,B)) + clippingNvcMap((false,false,i,j,B)) ) + "\n");
        }
      }
    }
    close(writer2);
  }
  

  
  def tabulateLeadClippedSequence(r : SAMRecord, reverse : Boolean, clippingNvcMap : HashMap[(Boolean,Boolean,Int,Int,Char), Int]) : Int = {
    val cigarSeq = getCigarSeq(r, r.getReadNegativeStrandFlag());
    val clipLength = if(cigarSeq.head.getOperator() == CigarOperator.SOFT_CLIP) { cigarSeq.head.getLength() } else 0;
    val seqArray : Array[Byte] = r.getReadBases();
    val clipSeqArray = if(r.getReadNegativeStrandFlag()){
      seqArray.takeRight(clipLength).reverse.map(reverseBaseByteMap(_)).map(_.toChar);
    } else {
      seqArray.take(clipLength).map(_.toChar);
    }
    for(i <- 0 until clipLength){
      clippingNvcMap( (reverse, true, clipLength, i + 1, clipSeqArray(i)) ) += 1;
    }
    
    return clipLength;
  }
  def tabulateTailClippedSequence(r : SAMRecord, reverse : Boolean, clippingNvcMap : HashMap[(Boolean,Boolean,Int,Int,Char), Int]) : Int = {
    val cigarSeq = getCigarSeq(r, r.getReadNegativeStrandFlag());
    val clipLength = if(cigarSeq.last.getOperator() == CigarOperator.SOFT_CLIP) { cigarSeq.last.getLength() } else 0;
    val seqArray : Array[Byte] = r.getReadBases();
    val clipSeqArray = if(r.getReadNegativeStrandFlag()){
      seqArray.take(clipLength).map(_.toChar);
    } else {
      seqArray.takeRight(clipLength).reverse.map(reverseBaseByteMap(_)).map(_.toChar);
    }
    for(i <- 0 until clipLength){
      clippingNvcMap( (reverse, false, clipLength, i + 1, clipSeqArray(i)) ) += 1;
    }
    
    return clipLength;
  }
  
  def getCigarSeq(r : SAMRecord, reverse : Boolean) : Seq[CigarElement] = {
    if(reverse){
      r.getCigar().getCigarElements().toSeq.reverse;
    } else {
      r.getCigar().getCigarElements().toSeq;
    }
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
  /*
  def getStrand(samRecord : SAMRecord) : Char = {
    if(samRecord.getFirstOfPairFlag()){
      if(samRecord.getReadNegativeStrandFlag) '-' else '+';
    } else {
      if(samRecord.getMateNegativeStrandFlag) '-' else '+';
    }
  }
  
  def isFirstRead(samRecord : SAMRecord) : Boolean = {
    samRecord.getFirstOfPairFlag();
  }
  
  def useRead(r1 : SAMRecord, r2 : SAMRecord, causeOfDropArray : Array[Int]) : Boolean = {
    if(r1.getMateUnmappedFlag() || r2.getMateUnmappedFlag()) { causeOfDropArray(0) = causeOfDropArray(0) + 1; return false; }
    if((! r1.getProperPairFlag()) || ( ! r2.getProperPairFlag())) { causeOfDropArray(1) = causeOfDropArray(1) + 1; return false; }
   // if(! samRecord.getFirstOfPairFlag()) return false;
    if(r1.getReadFailsVendorQualityCheckFlag() || r2.getReadFailsVendorQualityCheckFlag()){ causeOfDropArray(2) = causeOfDropArray(2) + 1; return false; }
    if(r1.getNotPrimaryAlignmentFlag() || r2.getNotPrimaryAlignmentFlag()) { causeOfDropArray(3) = causeOfDropArray(3) + 1; return false; }
    if( r1.isValid != null || r2.isValid != null) { causeOfDropArray(4) = causeOfDropArray(4) + 1; return false; }
    
    //add more checks here?
    return true;
  }*/
  
}

class qcClippingSequence(readLen : Int, writeClippedNVC : Boolean)  extends QCUtility[Unit] {
  val readLength : Int = readLen;
  val nvcCounter : qcClippingSequence.NvcCounter = new qcClippingSequence.NvcCounter(readLength);
  var readLengthSet : Set[Int] = Set[Int]();
  var maxLeadClip = -1;
  var maxTailClip = -1;
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    //if(useRead(r1) && useRead(r2) ){
      val (r1_leadClip, r1_tailClip, r1_curr_readLength) = qcClippingSequence.countAll(r1, true, nvcCounter);
      val (r2_leadClip, r2_tailClip, r2_curr_readLength) = qcClippingSequence.countAll(r2, false, nvcCounter);
      maxLeadClip = math.max(maxLeadClip, math.max(r1_leadClip,r2_leadClip));
      maxTailClip = math.max(maxTailClip, math.max(r1_tailClip,r2_tailClip));
        
      readLengthSet += r1_curr_readLength;
      readLengthSet += r2_curr_readLength;
    //}
  }
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
     qcClippingSequence.writeTabulatedSequence(outfile, maxLeadClip, maxTailClip, nvcCounter, readLengthSet, writeClippedNVC);
  }
  
  def getUtilityName : String = "DEPRECIATED";

  
}



















