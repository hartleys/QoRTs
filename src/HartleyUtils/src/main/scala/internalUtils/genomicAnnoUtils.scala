package internalUtils

import scala.collection.immutable.TreeMap;
//import scala.collection.immutable.HashSet;
//import scala.collection.immutable.HashMap;
import internalUtils.stdUtils._;
import internalUtils.Reporter._;

import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.JavaConversions._

object genomicAnnoUtils {

  /*
   * Useful classes:
   */
 
   
   
  object GenomicArrayOfSets {
    def apply[B](isStranded : Boolean) : GenomicArrayOfSets[B] = {
      if(isStranded){
        //new GenomicArrayOfSets_Stranded[B]();
        new GenomicArrayOfSets_Generic[B](isStranded);
      } else {
        //new GenomicArrayOfSets_Unstranded[B]();
        new GenomicArrayOfSets_Generic[B](isStranded);
      }
    }  
     
    def printGenomicArrayToFile[B](outfile : String, ga : GenomicArrayOfSets[B]){
      val writer : internalUtils.fileUtils.WriterUtil = internalUtils.fileUtils.openWriter(outfile);
      reportln("Starting print of gtf. Memory usage: " + MemoryUtil.memInfo,"note");
      for((iv, stepSet) <- ga.finalizeStepVectors.getSteps){
        //  For reference: case class OutputGtfLine(in_chromName : String, in_featureSource : String, in_featureType : String, in_start : Int, in_end : Int, in_score : String, in_strand : Char, in_attributeMap : Map[String,String], in_gtfFmt_attributeBreak : String, in_stranded : Boolean) extends GtfLine
        val attrMap : Map[String,String] = Map("gene_id" -> stepSet.mkString("+"), "gene_ct" -> stepSet.size.toString);
        val gtfLine : GtfLine = new FlatOutputGtfLine(iv.chromName, FlatOutputGtfLine.DEF_FEATURESOURCE, "exonic_part", iv.start + 1, iv.end, FlatOutputGtfLine.DEF_SCORE, iv.strand, attrMap, FlatOutputGtfLine.DEF_ATTRBREAK, ga.isStranded );
        
        writeGtfLine(gtfLine, writer);
      }
      reportln("Finished print of gtf. Memory usage: " + MemoryUtil.memInfo,"note");
      internalUtils.fileUtils.close(writer);
    }
  }
  
  
  abstract class GenomicArrayOfSets[B] {
    def hasChrom(chromName : String) : Boolean;
    def isStranded : Boolean;
    def isFinalized : Boolean;
    
    def addSpan(span : GenomicInterval, element : B);
    
    def finalizeStepVectors : GenomicArrayOfSets[B];
    def findWhollyContainedSteps(iv : GenomicInterval) : Iterator[(GenomicInterval, Set[B])];
    def findIntersectingSteps(iv : GenomicInterval) : Iterator[(GenomicInterval, Set[B])];
    
    def getChroms : Iterator[(String,Char)];
    def getSteps(chromName : String, strand : Char) : Iterator[(GenomicInterval, Set[B])];
    def getSteps : Iterator[(GenomicInterval, Set[B])];
    
    def getValueSet : Set[B];
    
  }
  
  object GenomicSpliceJunctionSet {
    def apply[B](isStranded : Boolean) : GenomicSpliceJunctionSet[B] = {
      return new GenomicSpliceJunctionSet[B](isStranded,  Map[GenomicInterval,B]());
    }
  }
    
  case class GenomicSpliceJunctionSet[B](isStranded : Boolean, sjMap : Map[GenomicInterval,B]) {
    def getSpliceJunction(iv : GenomicInterval) : Option[B] = {
      sjMap.get(iv);
    }
    def addSpliceJunction(iv : GenomicInterval, element : B) : GenomicSpliceJunctionSet[B] = {
      if(iv.strand == '.' && isStranded) error("!");
      return new GenomicSpliceJunctionSet[B](isStranded, sjMap + ((iv.usingStrandedness(isStranded), element)));
    }
  }
  
  /*
   * Implementation:
   */ 
  
  class GenomicArrayOfSets_Generic[B](in_isStranded : Boolean) extends GenomicArrayOfSets[B] {
    private val chromSet : scala.collection.mutable.Set[String] = scala.collection.mutable.Set[String]();
    private val chromMap : scala.collection.mutable.Map[(String,Char),InternalStepVector[B]] = scala.collection.mutable.AnyRefMap[(String,Char),InternalStepVector[B]]();
    private var isFinalizedFlag : Boolean = false;
    private val gsvMap : scala.collection.mutable.AnyRefMap[(String,Char),GenomicStepVector[B]] = scala.collection.mutable.AnyRefMap[(String,Char),GenomicStepVector[B]]();
    private var valueSet : Set[B] = Set[B]();
    
    def getValueSet : Set[B] = valueSet;


    val isStranded : Boolean = in_isStranded;
    
    def hasChrom(chromName : String) : Boolean = if(isStranded){
      chromMap.contains((chromName,'+'));
    } else {
      chromMap.contains((chromName,'.'));
    }
    
    def isFinalized : Boolean = isFinalizedFlag;
     
    def finalizeStepVectors : GenomicArrayOfSets[B] = {
      if(isFinalizedFlag) return this;
      //gsvMap = Map[(String,Char),GenomicStepVector[B]]();
      chromMap.foreach((x) => {
        val chromName = x._1._1;
        val strand = x._1._2;
        val isv = x._2;
        gsvMap((chromName, strand)) = (new GenomicStepVector[B](chromName , strand , isv));
      });
      isFinalizedFlag = true;
      gsvMap.repack;
      return this;
      //gsv = new GenomicStepVector[B](chromName : String, strand : Char, isv : InternalStepVector[B])
    }
    
    def addSpan(iv_in : GenomicInterval, element : B){
      val iv = iv_in.usingStrandedness(isStranded);
      chromMap.get((iv.chromName, iv.strand)) match {
        case Some(isv : InternalStepVector[B]) => {
          chromMap((iv.chromName, iv.strand)) = isv.addInterval(iv.start,iv.end,element);
        }
        case None => {
          chromMap((iv.chromName,swapStrand(iv.strand))) = InternalStepVector[B]();
          chromMap((iv.chromName,iv.strand)) = InternalStepVector[B]().addInterval(iv.start,iv.end,element) ;
        }
      }
      valueSet += element;
    }
    
    /*
     * WARNING: the following methods can ONLY be called AFTER the array of sets is finalized!
     */
    def findWhollyContainedSteps(iv : GenomicInterval) : Iterator[(GenomicInterval, Set[B])] = {
      gsvMap.get((iv.chromName, iv.strandStranded(isStranded))) match {
            case Some(gsv : GenomicStepVector[B]) => {
              gsv.findWhollyContainedSteps(iv.start,iv.end);
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    }
    def findIntersectingSteps(iv : GenomicInterval) : Iterator[(GenomicInterval, Set[B])] = {
      gsvMap.get((iv.chromName, iv.strandStranded(isStranded))) match {
            case Some(gsv : GenomicStepVector[B]) => {
              gsv.findIntersectingSteps(iv.start,iv.end);
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    }
    def getChroms : Iterator[(String,Char)] = {
      return gsvMap.keys.toSeq.sorted.iterator;
    }
    def getSteps(chromName : String, strand : Char) : Iterator[(GenomicInterval, Set[B])] = {
      val strandedStrand = if(isStranded) strand else '.';
      gsvMap.get((chromName, strandedStrand)) match {
            case Some(gsv : GenomicStepVector[B]) => {
              gsv.getSteps
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    }
    def getSteps : Iterator[(GenomicInterval, Set[B])] = {
      //reportln("Getting steps from: " + gsvMap.size + " chromosomes.","note");
      
      return gsvMap.keys.toSeq.sorted.foldLeft( (Seq[(GenomicInterval,Set[B])]()).iterator )( (soFar,currKey) =>{
        val gsv = gsvMap(currKey);
        soFar ++ gsv.getSteps;
      });
    }
  }
  
  
  class GenomicArrayOfSets_Stranded[B] extends GenomicArrayOfSets[B] {
    private var chromMap : Map[(String,Char),InternalStepVector[B]] = Map[(String,Char),InternalStepVector[B]]();
    private var isFinalizedFlag : Boolean = false;
    private var gsvMap : Map[(String,Char),GenomicStepVector[B]] = null;
    private var valueSet : Set[B] = Set[B]();
    def getValueSet : Set[B] = valueSet;

    def hasChrom(chromName : String) : Boolean = {
      chromMap.contains((chromName,'+'));
    }
    def isStranded : Boolean = true;
    def isFinalized : Boolean = isFinalizedFlag;
    
    def finalizeStepVectors : GenomicArrayOfSets[B] = {
      if(isFinalizedFlag) return this;
      gsvMap = Map[(String,Char),GenomicStepVector[B]]();
      chromMap.map((x) => {
        val chromName = x._1._1;
        val strand = x._1._2;
        val isv = x._2;
        gsvMap += ((((chromName, strand)) , new GenomicStepVector[B](chromName , strand , isv)));
      });
      isFinalizedFlag = true;
      return this;
      //gsv = new GenomicStepVector[B](chromName : String, strand : Char, isv : InternalStepVector[B])
    }
    
    def addSpan(iv : GenomicInterval, element : B){
      chromMap.get((iv.chromName, iv.strand)) match {
        case Some(isv : InternalStepVector[B]) => {
          chromMap = chromMap + (( (iv.chromName, iv.strand), isv.addInterval(iv.start,iv.end,element)  ));
        }
        case None => {
          chromMap = chromMap + (( (iv.chromName,iv.strand            ),  InternalStepVector[B]().addInterval(iv.start,iv.end,element)  ));
          chromMap = chromMap + (( (iv.chromName,swapStrand(iv.strand)),  InternalStepVector[B]()                                 ));
        }
      }
      valueSet = valueSet + element;
    }
    
    /*
     * WARNING: the following methods can ONLY be called AFTER the array of sets is finalized!
     */
    
    def findWhollyContainedSteps(iv : GenomicInterval) : Iterator[(GenomicInterval, Set[B])] = {
      gsvMap.get((iv.chromName, iv.strand)) match {
            case Some(gsv : GenomicStepVector[B]) => {
              gsv.findWhollyContainedSteps(iv.start,iv.end);
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    }
    def findIntersectingSteps(iv : GenomicInterval) : Iterator[(GenomicInterval, Set[B])] = {
      gsvMap.get((iv.chromName, iv.strand)) match {
            case Some(gsv : GenomicStepVector[B]) => {
              gsv.findIntersectingSteps(iv.start,iv.end);
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    }
    def getChroms : Iterator[(String,Char)] = {
      return gsvMap.keys.toSeq.sorted.iterator;
    }
    def getSteps(chromName : String, strand : Char) : Iterator[(GenomicInterval, Set[B])] = {
      gsvMap.get((chromName, strand)) match {
            case Some(gsv : GenomicStepVector[B]) => {
              gsv.getSteps
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    }
    def getSteps : Iterator[(GenomicInterval, Set[B])] = {
      //reportln("Getting steps from: " + gsvMap.size + " chromosomes.","note");
      
      return gsvMap.keys.toSeq.sorted.foldLeft( (Seq[(GenomicInterval,Set[B])]()).iterator )( (soFar,currKey) =>{
        val gsv = gsvMap(currKey);
        soFar ++ gsv.getSteps;
      });
    }
  }

  class GenomicArrayOfSets_Unstranded[B] extends GenomicArrayOfSets[B] {
    private val internalArray = new GenomicArrayOfSets_Stranded[B];
    
    private var chromMap : Map[(String,Char),InternalStepVector[B]] = Map[(String,Char),InternalStepVector[B]]();
    private var isFinalizedFlag : Boolean = false;
    private var gsvMap : Map[(String,Char),GenomicStepVector[B]] = null;
    private var valueSet : Set[B] = Set[B]();
    
    def getValueSet : Set[B] = internalArray.getValueSet;

    def hasChrom(chromName : String) : Boolean = {
      chromMap.contains((chromName,'.'));
    }
    def isStranded : Boolean = true;
    def isFinalized : Boolean = isFinalizedFlag;
    
    def finalizeStepVectors : GenomicArrayOfSets[B] = {
      if(isFinalizedFlag) return this;
      gsvMap = Map[(String,Char),GenomicStepVector[B]]();
      chromMap.map((x) => {
        val chromName = x._1._1;
        val strand = x._1._2;
        val isv = x._2;
        gsvMap += ((((chromName, strand)) , new GenomicStepVector[B](chromName , strand , isv)));
      });
      isFinalizedFlag = true;
      return this;
      //gsv = new GenomicStepVector[B](chromName : String, strand : Char, isv : InternalStepVector[B])
    }
    
    def addSpan(iv : GenomicInterval, element : B){
      chromMap.get((iv.chromName, '.')) match {
        case Some(isv : InternalStepVector[B]) => {
          chromMap = chromMap + (( (iv.chromName, '.'), isv.addInterval(iv.start,iv.end,element)  ));
        }
        case None => {
          chromMap = chromMap + (( (iv.chromName,'.'            ),  InternalStepVector[B]().addInterval(iv.start,iv.end,element)  ));
        }
      }
      valueSet = valueSet + element;
    }
    
    /*
     * WARNING: the following methods can ONLY be called AFTER the array of sets is finalized!
     */
    
    def findWhollyContainedSteps(iv : GenomicInterval) : Iterator[(GenomicInterval, Set[B])] = {
      gsvMap.get((iv.chromName, '.')) match {
            case Some(gsv : GenomicStepVector[B]) => {
              gsv.findWhollyContainedSteps(iv.start,iv.end);
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    } 
    def findIntersectingSteps(iv : GenomicInterval) : Iterator[(GenomicInterval, Set[B])] = {
      gsvMap.get((iv.chromName, '.')) match {
            case Some(gsv : GenomicStepVector[B]) => {
              gsv.findIntersectingSteps(iv.start,iv.end);
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    }
    def getChroms : Iterator[(String,Char)] = {
      return gsvMap.keys.toSeq.sorted.iterator;
    }
    def getSteps(chromName : String, strand : Char) : Iterator[(GenomicInterval, Set[B])] = {
      gsvMap.get((chromName, strand)) match {
            case Some(gsv : GenomicStepVector[B]) => {
              return gsv.getSteps
            }
            case None => return (Seq[(GenomicInterval,Set[B])]()).iterator;
          }
    }
    def getSteps : Iterator[(GenomicInterval, Set[B])] = {
      //reportln("Getting steps from: " + gsvMap.size + " chromosomes.","note");
      
      return gsvMap.keys.toVector.sorted.foldLeft( (Seq[(GenomicInterval,Set[B])]()).iterator )( (soFar,currKey) =>{
        val gsv = gsvMap(currKey);
        soFar ++ gsv.getSteps;
      });
    }
  }
  
  /*
   * Internally-used classes and methods:
   */
  
  class GenomicStepVector[B](chromName : String, strand : Char, isv : InternalStepVector[B]){
    private val ivMap : Map[Int,GenomicInterval] = {
      var out = Map[Int,GenomicInterval]();
      var prev = 0;
      val iter = isv.steps.iterator;
      
      while(iter.hasNext) {
        val next = iter.next._1;
        out += ((next, GenomicInterval(chromName, strand,prev,next)));
        prev = next;
      }
      out;
    }
    
    private def wrapIterator(iter : Iterator[(Int, Set[B])]) : Iterator[(GenomicInterval, Set[B])] = {
      return new Iterator[(GenomicInterval,Set[B])]{
        def hasNext = iter.hasNext;
        def next : (GenomicInterval, Set[B]) = {
          val (nextEnd, nextElement) = iter.next;
          return (ivMap(nextEnd), nextElement);
        }
      }
    }
    
    def findWhollyContainedSteps(start : Int, end : Int) : Iterator[(GenomicInterval, Set[B])] = {
      wrapIterator(isv.findWhollyContainedSteps(start, end));
    }
    def findIntersectingSteps(start : Int, end : Int) : Iterator[(GenomicInterval, Set[B])] = {
      wrapIterator(isv.findIntersectingSteps(start, end));
    }
    def getSteps : Iterator[(GenomicInterval, Set[B])] = {
      wrapIterator(isv.getSteps);
    }

  }
  
  /*
   * ******************************************************************
   * Internal Step Vector:
   */
  
  object InternalStepVector {
    
    def apply[B]() : InternalStepVector[B] = {
      val tm = (new TreeMap[Int,Set[B]]()) + ((Int.MaxValue, Set[B]()));
      return new InternalStepVector[B](tm);
    }
  } 
   
  /* InternalStepVector: a tool for searching intervals across a genome.
   * 
   * (note: InternalStepVector is immutable)
   */
  
  case class InternalStepVector[B](steps : TreeMap[Int,Set[B]]) {

    /*  
     * Construction and modification methods: 
     */
    def addInterval(start : Int, end : Int, obj : B) : InternalStepVector[B] = {
      val curr : InternalStepVector[B] = this.insertBreakAt(start).insertBreakAt(end);
      
      val newSteps1 : Iterator[(Int,Set[B])] = curr.findIntersectingSteps(start,end);
      
      return newSteps1.foldLeft[InternalStepVector[B]]( curr )( (acc, step) =>{
        acc.addElementToSpanByEndPosition(step._1, obj);
      })
    }
    
    private def addElementToSpanByEndPosition(end : Int, element : B) : InternalStepVector[B] = {
      //Note: This method is unchecked!
      //Add this check if safe behavior is desired:
      if(! steps.contains(end)) error("ERROR! InternalStepVector.addElementToSpanByEndPosition("+end+")! for InternalStepVector. Impossible state?");
      
      return new InternalStepVector(steps + ((end, steps(end) + element )))
    }
    
    def insertBreakAt(pos : Int) : InternalStepVector[B] = {
      if(steps.contains(pos)){
        return this;
      } else {
        val (oldEnd, breakSet) : (Int, Set[B]) = steps.iteratorFrom(pos).next;
        return new InternalStepVector( steps + ((pos, breakSet)));
      }
    }
    
    /*
     * Search methods:
     */
    
    def getSteps : Iterator[(Int, Set[B])] = {
      steps.iterator;
    }
    final val InternalStepVector_NULL_PAIR = (Int.MaxValue, Set[B]());
     
    def findWhollyContainedSteps(start : Int, end : Int) : Iterator[(Int, Set[B])] = {
      val iter = steps.iteratorFrom(start);
      new Iterator[(Int,Set[B])]{
        var next_holder : (Int, Set[B]) = {
          val b = iter.next;
          if(iter.hasNext){
            iter.next;
          } else {
            InternalStepVector_NULL_PAIR
          }
          //If not assigned already, then the iterator must now be empty, and the curr value is irrelevant:
        }
        def hasNext : Boolean = {
          next_holder._1 <= end && iter.hasNext
        }
        def next : (Int,Set[B]) = {
          val buff = next_holder;
          next_holder = iter.next;
          return next_holder;
        }
      }
    }
    
    def findIntersectingSteps(start : Int, end : Int) : Iterator[(Int,Set[B])] = {
      val iter = steps.iteratorFrom(start+1);
      
      new Iterator[(Int,Set[B])]{
        var prev : (Int, Set[B]) = (-1, Set[B]());
        def hasNext : Boolean = {
          return prev._1 < end && iter.hasNext
        }
        def next : (Int,Set[B]) = {
          prev = iter.next;
          return prev;
        }
      }
    }
  }
  
  private def swapStrand(s : Char) : Char = {
    if(s == '+') '-';
    else if(s == '-') '+';
    else '.';
  }
}