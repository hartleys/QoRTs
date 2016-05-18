package internalUtils

import scala.collection.immutable.TreeMap;
import scala.collection.immutable.HashSet;
import scala.collection.immutable.HashMap;
import internalUtils.stdUtils._;
import internalUtils.Reporter._;

import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._; 

import net.sf.samtools._
import scala.collection.JavaConversions._;


object genomicUtils {

  def getBlockStart(b : AlignmentBlock) : Int = {
    b.getReferenceStart - 1;
  }
  def getBlockEnd(b : AlignmentBlock) : Int = {
    b.getReferenceStart - 1 + b.getLength();
  } 
      
  def getGenomicIntervalsFromRead(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean) : Iterator[GenomicInterval] = {
    val strand = getStrand(r,stranded,fr_secondStrand);
    val blocks : Iterator[AlignmentBlock] = r.getAlignmentBlocks().iterator;
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
    val blocks : Iterator[AlignmentBlock] = r.getAlignmentBlocks().iterator;
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
    
  }
  //FIX ME FIX ME FIX ME!
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
        currReadPos = readEnd;
        currRefPos = refEnd; 
        
        CigarBlock(currRefPos, refEnd, currReadPos, readEnd, op, len);
      }
      def hasNext : Boolean = cigIter.hasNext;
    }
    
    if(r.getReadNegativeStrandFlag()){
      return out
      //FIX ME FIX ME FIX ME!
    } else {
      return out;
    }

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
  def getAlignedBasePositionsFromRead(r : SAMRecord) : Vector[(Int, Int)] = {
    val cb = getCigarBlocksFromRead(r);
    //val cb = cbAll.filter(_.op.consumesReadBases());
    
    cb.filter(_.op.equals(CigarOperator.MATCH_OR_MISMATCH)).foldLeft(Vector[(Int,Int)]())((soFar, b) => {
      soFar ++ (b.refStart until b.refEnd).zip((b.readStart until b.readEnd)).toVector;
    });
  }
  
}


















