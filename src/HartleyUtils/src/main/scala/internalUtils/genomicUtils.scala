package internalUtils

import scala.collection.immutable.TreeMap;
import scala.collection.immutable.HashSet;
import scala.collection.immutable.HashMap;
import internalUtils.stdUtils._;
import internalUtils.Reporter._;

import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._; 

import net.sf.samtools._
//import scala.collection.JavaConversions._;
import scala.collection.JavaConverters._;


object genomicUtils {

  def getBlockStart(b : AlignmentBlock) : Int = {
    b.getReferenceStart - 1;
  }
  def getBlockEnd(b : AlignmentBlock) : Int = {
    b.getReferenceStart - 1 + b.getLength();
  } 
  
  def mergeGenomicIntervalSets(ivs : Seq[GenomicInterval]) : Seq[GenomicInterval] = {
    val sortivs = ivs.sorted; 
    
    sortivs.tail.foldLeft(Seq[GenomicInterval](sortivs.head))((soFar,iv) => {
      if(soFar.last.chromName != iv.chromName || soFar.last.strand != iv.strand){
        soFar :+ iv;
      } else if(soFar.last.end >= iv.start){
        soFar.init :+ new GenomicInterval(iv.chromName,iv.strand,soFar.last.start,math.max(iv.end,soFar.last.end));
      } else {
        soFar :+ iv;
      }
    })
  }
  
  def getOverlap(iv1 : GenomicInterval, iv2 : GenomicInterval) : Int = {
    if(iv1.chromName != iv2.chromName) return 0;
    else if(iv1.strand != iv2.strand) return 0;
    else if(iv1.start < iv2.end && iv1.end > iv2.start){
      math.min(iv1.end,iv2.end) - math.max(iv1.start,iv2.start);
    } else {
      return 0;
    }
  }
  
  def getGenomicIntervalsFromRead(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean) : Iterator[GenomicInterval] = {
    val strand = getStrand(r,stranded,fr_secondStrand);
    val blocks : Iterator[AlignmentBlock] = r.getAlignmentBlocks().asScala.iterator;
    val chromName = r.getReferenceName();
     
    return new Iterator[GenomicInterval](){
      def hasNext = blocks.hasNext;
      def next = {
        val nextBlock : AlignmentBlock = blocks.next;
        new GenomicInterval(chromName,strand,getBlockStart(nextBlock), getBlockEnd(nextBlock));
      }
    }
  }
  def getExpandedGenomicIntervalsFromRead(d : Int, r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean, reverseStrand : Boolean = false) : Iterator[GenomicInterval] = {
    val rawStrand = getStrand(r,stranded,fr_secondStrand);
    val strand = if((! stranded) || (! reverseStrand)){ rawStrand } else { 
      if(rawStrand == '.') '.';
      else if(rawStrand == '+') '-';
      else '+';
    } 
    val blocks : Iterator[AlignmentBlock] = r.getAlignmentBlocks().asScala.iterator;
    val chromName = r.getReferenceName();
    
    return new Iterator[GenomicInterval](){
      def hasNext = blocks.hasNext;
      def next = {
        val nextBlock : AlignmentBlock = blocks.next;
        new GenomicInterval(chromName,strand,getBlockStart(nextBlock) - d, getBlockEnd(nextBlock) + d);
      }
    }
  }

  case class CigarBlock(refStart : Int, refEnd : Int, readStart : Int, readEnd : Int, op : CigarOperator, len : Int){
    override def toString() : String = {
      "[CigarBlock("+len+op.toString()+"):ref("+refStart+","+refEnd+"):read("+readStart+","+readEnd+")]"
    }
  }
  

  def getCigarBlocksFromRead(r : SAMRecord) : Iterator[CigarBlock] = {
    val out = new Iterator[CigarBlock](){
      var currReadPos = 0;
      var currRefPos = r.getAlignmentStart() - 1;
      val cigIter = r.getCigar().getCigarElements.iterator;
      
      def next : CigarBlock = {
        val ce : CigarElement = cigIter.next;
        val op = ce.getOperator();
        val len = ce.getLength();
        
        val refEnd = if(op.consumesReferenceBases()) currRefPos + len else currRefPos;
        val readEnd = if(op.consumesReadBases()) currReadPos + len else currReadPos;
        val refStart = currRefPos;
        val readStart = currReadPos;
        currReadPos = readEnd;
        currRefPos = refEnd; 
        
        CigarBlock(refStart, refEnd, readStart, readEnd, op, len);
      }
      def hasNext : Boolean = cigIter.hasNext;
    }
    return out;
  }
  
  def getSeqStringFromBlock(r : SAMRecord, cb : CigarBlock) : String = {
    r.getReadString().substring(cb.readStart,cb.readEnd);
  }
  
  //UNTESTED!
  def truncateCigarBlocks(cb : Iterator[CigarBlock], start : Int, end : Int) : Iterator[CigarBlock] = {
    val filtcb = cb.filter(b => b.refStart < end && b.refEnd > start);
    
    if(filtcb.isEmpty){
      return Iterator[CigarBlock]();
    } else {
      return filtcb.map(b => truncateBlock(b,start,end));
    }
  }
  
  //UNTESTED!
  def truncateBlock(cb : CigarBlock,start : Int, end : Int) : CigarBlock = {
    if(! cb.op.consumesReferenceBases()){
      return cb;
    } else {
      var out = cb;
      if(out.refStart < start && out.refEnd > start){
        out = CigarBlock(start, out.refEnd, out.readStart + (start - out.refStart), out.readEnd, out.op, out.len - (start - out.refStart));
      }
      if(out.refStart < end && out.refEnd > end){
        out = CigarBlock(out.refStart, end, out.readStart, out.readEnd - (out.refEnd - end), out.op, out.len - (out.refEnd - end));
      }
      return out;
    }
  }
  
  def getAlignedBasePositionsFromRead(r : SAMRecord) : Vector[(Int, Int)] = {
    val cb = getCigarBlocksFromRead(r);
    //val cb = cbAll.filter(_.op.consumesReadBases());
    
    cb.filter(_.op.equals(CigarOperator.MATCH_OR_MISMATCH)).foldLeft(Vector[(Int,Int)]())((soFar, b) => {
      soFar ++ (b.refStart until b.refEnd).zip((b.readStart until b.readEnd)).toVector;
    });
  }
  
  /*def getBasePositionsFromRead(r : SAMRecord) : Iterator[(Int, Int, CigarOperator)] = {
    val cb = getCigarBlocksFromRead(r);
    //val cb = cbAll.filter(_.op.consumesReadBases());
    
    cb.foldLeft(Iterator[(Int,CigarOperator)]())((soFar, b) => {
      if(b.op.consumesReadBases){
        soFar ++ (b.refStart until b.refEnd).zip(b.readStart until b.readEnd).iterator.map( (_,b.op) );
      } else {
        soFar;
      }
    });
  }*/
  

  
}


















