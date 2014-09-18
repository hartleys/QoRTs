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
  def getExpandedGenomicIntervalsFromRead(d : Int, r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean) : Iterator[GenomicInterval] = {
    val strand = getStrand(r,stranded,fr_secondStrand);
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
  
}