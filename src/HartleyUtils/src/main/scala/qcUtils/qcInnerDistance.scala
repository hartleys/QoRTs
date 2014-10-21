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
import scala.collection.immutable.HashMap;

import net.sf.samtools._

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.immutable.TreeSet;
import scala.collection.immutable.TreeMap;

import scala.util.control.Breaks._;
import internalUtils.optionHolder._;

import scala.collection.GenMap;

object qcInnerDistance {

  
  /* 
   * *****************************************************************************************************************************
   */
  
  def readSplicesFromGtfFile(gtffile : String, stranded : Boolean) : HashMap[(String,Char),TreeSet[(Int,Int)]] = {
    return GtfReader.getGtfReader(gtffile, stranded, true, "\\s+").foldLeft( new HashMap[(String,Char),TreeSet[(Int,Int)]]() )(
        (acc, gtfLine) => {
          readGtfLine(gtfLine,acc);
        }
    )
  }
  
  def readGtfLine(gtfLine : GtfLine, acc :  HashMap[(String,Char),TreeSet[(Int,Int)]]) : HashMap[(String,Char),TreeSet[(Int,Int)]] = {
    if(gtfLine.featureType == "splice_site" || gtfLine.featureType == "novel_splice_site"){
      acc.get((gtfLine.chromName, gtfLine.strand)) match {
        case Some(chromTree) => {
          return acc + (((gtfLine.chromName, gtfLine.strand), chromTree + ((gtfLine.start, gtfLine.end)) ));
        }
        case None => {
          return acc + (((gtfLine.chromName, gtfLine.strand), new TreeSet[(Int,Int)]() + ((gtfLine.start, gtfLine.end))));
        }
      }
    } else return acc;
  }
  
  /*
   * *****************************************************************************************************************************
   */
  
  private def findShortestPath(stranded : Boolean, chromName : String, strand : Char, start : Int, end : Int, spliceAnnotation : GenMap[(String,Char),TreeSet[(Int,Int)]]) : (Int, StringBuilder) = {
    val s = if(stranded) strand else '.';
    val sb = new StringBuilder();
    //reportln("###	["+start + "	"+end+") = " + (end - start),"note");
    
    val shortestPath = if(spliceAnnotation.contains((chromName,s))){
      val spliceAnno = spliceAnnotation((chromName,s));
      val spliceJunctionsBetween = findSplicesBetween(start,end,spliceAnno);
      
      //DEBUG:
      spliceJunctionsBetween.foreach((sj) => sb.append("#	[	"+sj._1+",	"+sj._2+"	]\n"))
      
      findShortestPath_helper(start,end,spliceJunctionsBetween);
    } else {
      //warning? chromosome not found?
      //DEBUG:
      sb.append("#	[	NA,	NA	]\n");
      
      end - start;
    }
    
    return (shortestPath,sb);
  }
  
  private def findShortestPath_helper(start : Int, end : Int, spliceJunctionsBetween : TreeSet[(Int,Int)]) : Int = {
    
    val spliceDistance : scala.collection.mutable.HashMap[(Int,Int), Int] = new scala.collection.mutable.HashMap[(Int,Int),Int]()
    {
      override def default(sj : (Int,Int)) = sj._1 + 1 - start;
    }
    
    //val spliceJunctions = spliceJunctionsBetween + ((start,start)) + ((end,end+1));
    
    spliceDistance((start,start)) = 0;
    var minDistSoFar = end - start;
    
    for(sj <- spliceJunctionsBetween){
      val splicesAfter = findSplicesAfter(sj._2,spliceJunctionsBetween);
      val distanceSoFar = spliceDistance(sj);
      for(sjnext <- splicesAfter){
        spliceDistance(sjnext) =  math.min( spliceDistance(sjnext) , distanceSoFar + sjnext._1 + 1 - sj._2);
        //val newDist = math.min( spliceDistance(sjnext) )
          /*spliceDistance.get(sjnext) match {
          case Some(d) =>  math.min( d, distanceSoFar + sjnext._1 - sj._2 );
          case None => distanceSoFar + sjnext._1 - sj._2;
        }*/
        //spliceDistance += ((sjnext,newDist));
      }
      minDistSoFar = math.min(minDistSoFar, distanceSoFar + (end - sj._2));
    }
    
    return minDistSoFar;
  }
  
  private def findSplicesAfter(start : Int,  spliceJunctions : TreeSet[(Int,Int)]) : TreeSet[(Int,Int)] = {
    spliceJunctions.from( (start, -1) );
  }
  private def findSplicesBetween(start : Int, end : Int, spliceJunctions : TreeSet[(Int,Int)]) : TreeSet[(Int,Int)] = {
    spliceJunctions.range((start,-1),(end,-1)).filter( ab => ab._2 <= end  );
  }
  
  /*
   * *****************************************************************************************************************************
   */
  
  val OVERLAP_CIGAR_MISMATCH_INTERNAL = -999999;
  
  val STAGGERED_NO_OVERLAP = -111111;
  val OVERLAP_CIGAR_MISMATCH_PARTIAL_OVERLAP = -222222;
  val STAGGERED_TOO_MUCH_ADAPTOR_ALIGNED = -333333;
  val STAGGERED_FULL_BLOCK_OF_ADAPTOR_ALIGNED = -444444;
  
  val STAGGERED_ADAPTOR_ALIGNMENT_LIMIT = 5;
  
  def getInsertSize(r1 : SAMRecord, r2 : SAMRecord, spliceAnnotation : GenMap[(String,Char),TreeSet[(Int,Int)]], stranded : Boolean, fr_secondStrand : Boolean) : (Int, Int) = {
    val (rf,rr) = if(r1.getReadNegativeStrandFlag) ((r2,r1)) else ((r1,r2));
    
    //CHECK FOR OFF-BY-ONE ERRORS!
    //reportln("rf: ["+getAlignmentStart(rf)+","+getAlignmentEnd(rf)+")\nrr: ["+getAlignmentStart(rr)+","+getAlignmentEnd(rr)+")","note");

    if(getAlignmentStart(rf) > getAlignmentEnd(rr)){
      //reportln("Method Exit: STAGGERED_BAD! ","note");
      return (STAGGERED_NO_OVERLAP,  2);
    } else if(getAlignmentEnd(rf) <= getAlignmentStart(rr)){
      //reportln("Method Entry: getInnerDistance_noOverlap: ","note");
      return (getInsertSize_noOverlap(rf,rr, spliceAnnotation, stranded, fr_secondStrand), 1);
    } else if(getAlignmentStart(rf) >= getAlignmentStart(rr)){
      //reportln("Method Entry: getInnerDistance_staggeredOverlap: ","note");
      val insertSize = getInsertSize_staggeredOverlap(rf,rr);
      return (insertSize, 0);
    } else {
      //reportln("Method Entry: getInnerDistance_partialOverlap: ","note");
      return (getInsertSize_partialOverlap(rf,rr), 3);
    }
  }
  
  var DEBUG_INTERNAL_InnerDistanceCalc_reportct = 0;
  var DEBUG_INTERNAL_InnerDistanceCalc_reportLimit = 10000;
  
  def getInsertSize_noOverlap(rf : SAMRecord,rr : SAMRecord, spliceAnnotation : GenMap[(String,Char),TreeSet[(Int,Int)]], stranded : Boolean, fr_secondStrand : Boolean) : Int = {
    val endF = getAlignmentEnd(rf);
    val startR = getAlignmentStart(rr);

    val clipF = getTailClipping(rf);
    val clipR = getLeadClipping(rr);
    if(endF == startR) return rf.getReadLength + rr.getReadLength - clipF - clipR;

    //report("noOverlap: clipF="+clipF+", clipR="+clipR+", endF="+endF+", startR="+startR +", maxDistance="+ (startR - endF) ,"note");
    
    val (minDistance, debugString) = findShortestPath(stranded, rf.getReferenceName(), getStrand(rf,stranded,fr_secondStrand), endF,startR,spliceAnnotation);
    //findShortestPath(stranded : Boolean, chromName : String, strand : Char, start : Int, end : Int, spliceAnnotation : HashMap[(String,Char),TreeSet[(Int,Int)]])
    
    val insertSize = rf.getReadLength() + rr.getReadLength() + minDistance - clipF - clipR 
    
    //if(DEBUG_INTERNAL_InnerDistanceCalc_reportct < DEBUG_INTERNAL_InnerDistanceCalc_reportLimit && (insertSize == 202 || insertSize == 201 || insertSize == 203)) {
    //  report("noOverlap: clipF="+clipF+", clipR="+clipR+", endF="+endF+", startR="+startR +", maxDistance="+ (startR - endF) + ", minDistance = "+minDistance + ". Insert Size = " + insertSize + "\n"+debugString.toString+"#	"+rf.getSAMString()+"#	"+rr.getSAMString(),"debug");
    //  DEBUG_INTERNAL_InnerDistanceCalc_reportct += 1;
    //}

    return insertSize; //CHECK FOR OFF-BY-ONE ERRORS!
  }
  
  def getInsertSize_staggeredOverlap(rf : SAMRecord, rr : SAMRecord) : Int = {
    val clipF = rf.getAlignmentBlocks.head.getReadStart() - 1;
    //val clipR = rr.getAlignmentBlocks.head.getReadStart() - 1;
    
    val startF = getAlignmentStart(rf);
    val startR = getAlignmentStart(rr);
    
    val lengthToIntersect = getReadLengthToIntersect(rr,rf);

    if(lengthToIntersect - clipF > STAGGERED_ADAPTOR_ALIGNMENT_LIMIT ){
      return STAGGERED_TOO_MUCH_ADAPTOR_ALIGNED
      //return (STAGGERED_TOO_MUCH_ADAPTOR_ALIGNED, 2, lengthToIntersect - clipF);
    } else {
      val insertSize = rf.getReadLength() - lengthToIntersect + clipF;
      return insertSize;
      //return (insertSize,2,-1);
    }
    
    //return insertSize;
    
    //OLD VERSION:
    //val rrFirstBlockEnd = rr.getAlignmentBlocks.head.getReferenceStart() - 1 + rr.getAlignmentBlocks.head.getLength();
    //if(rrFirstBlockEnd >= startF ) return STAGGERED_FULL_BLOCK_OF_ADAPTOR_ALIGNED;
    //if(startF == startR) return rr.getReadLength() + (clipF - clipR);
    //val distToOverlap = startF - startR;
    //val insertSize = if(distToOverlap - clipF  >= STAGGERED_ADAPTOR_ALIGNMENT_LIMIT) STAGGERED_TOO_MUCH_ADAPTOR_ALIGNED;
    //else rf.getReadLength() - distToOverlap + (clipF - clipR);
 
    //if(DEBUG_INTERNAL_InnerDistanceCalc_reportct < DEBUG_INTERNAL_InnerDistanceCalc_reportLimit && (insertSize == 101 | insertSize == 102 | insertSize == 103)) {
    //  reportln("StaggeredOverlap: clipF="+clipF+", clipR="+clipR+", startF="+startF+", startR="+startR + " equals = " + insertSize + "\n#	" + rf.getSAMString() +"#	"+rr.getSAMString(),"debug");
    //  DEBUG_INTERNAL_InnerDistanceCalc_reportct += 1;
    //}
    //return insertSize;
    
    ////EVEN OLDER VERSION:
    //////report("StaggeredOverlap: clipF="+clipF+", clipR="+clipR+", startF="+startF+", startR="+startR + " equals = " + (rr.getReadLength() + clipF - clipR - ( startF - startR )),"note");
    //////val lengthToIntersect = startF - startR;
    ////val insertSize = if(startF == startR){
    ////  rr.getReadLength() + clipF - clipR; //CHECK FOR OFF-BY-ONE ERRORS
    ////} else {
    ////  rr.getReadLength() + clipF - clipR - ( startF - startR );
    ////}
    ////if(insertSize == 101 | insertSize == 102 | insertSize == 103) {
    ////  reportln("StaggeredOverlap: clipF="+clipF+", clipR="+clipR+", startF="+startF+", startR="+startR + " equals = " + insertSize + "\n#	" + rf.getSAMString() +"\n#	"+rr.getSAMString(),"debug");
    ////}
    ////return insertSize;
  }
  def getInsertSize_partialOverlap(rf : SAMRecord,rr : SAMRecord) : Int = {
    //val abF = rf.getAlignmentBlocks.toSeq;
    ////val abR = rr.getAlignmentBlocks.toSeq;
    
    //val startR = getAlignmentStart(rr);
    val clipR = getLeadClipping(rr);//rr.getAlignmentBlocks.head.getReadStart() - 1;
    
    val lengthToIntersect = getReadLengthToIntersect(rf,rr);
    
    val insertSize = if(lengthToIntersect == OVERLAP_CIGAR_MISMATCH_INTERNAL) OVERLAP_CIGAR_MISMATCH_PARTIAL_OVERLAP;
    else rr.getReadLength + lengthToIntersect - clipR;
    
    return insertSize;
    //for(ab <- abF){
    //  val currStart = ab.getReferenceStart - 1;
    //  val currEnd = currStart + ab.getLength();
    //  
    //  if(currStart > startR){
    //    return OVERLAP_CIGAR_MISMATCH;
    //  } else if(currEnd > startR){
    //    val lengthToIntersect = startR - currStart;
    //    return (ab.getReadStart() + lengthToIntersect - clipR ) + rr.getReadLength();  //CHECK FOR OFF-BY-ONE ERRORS. NOTE: FOUND ONE. THIS COMMENTED-OUT VERSION IS UNFIXED.
    //  }
    //}
    //error("ERROR! qcInnerDistance.getInnerDistance_partialOverlap, Never found an intersecting block! Impossible State!");
    //-1
  }
  
  def getReadLengthToIntersect(rFirst : SAMRecord, rLast : SAMRecord) : Int = {
    val abFirst = rFirst.getAlignmentBlocks.toSeq;
    val startLast = getAlignmentStart(rLast);
    //val clipLast = getLeadClipping(rLast);
    for(ab <- abFirst){
      val currStart = ab.getReferenceStart - 1;
      val currEnd = currStart + ab.getLength();
      if(currStart > startLast){
        return OVERLAP_CIGAR_MISMATCH_INTERNAL;
      } else if(currEnd > startLast){
        val blockLengthToIntersect = startLast - currStart;
        return (ab.getReadStart() - 1) + blockLengthToIntersect;
      } //else do nothing.
    }
    error("ERROR! qcInnerDistance.getReadLengthToIntersect, Never found an intersecting block! Impossible State!");
    -1;
  }
  
  

  def getBlockReferenceSpan(b : AlignmentBlock) : (Int,Int) = {
    val start = b.getReferenceStart - 1;
    (start, b.getLength() + start);
  }
  def getBlockReferenceStart(b : AlignmentBlock) : Int = {
    b.getReferenceStart - 1;
  }
  def getBlockReferenceEnd(b : AlignmentBlock) : Int = {
    b.getReferenceStart - 1 + b.getLength();
  }
  def getBlockReadStart(b : AlignmentBlock) : Int = {
    b.getReadStart - 1;
  }
  def getBlockReadEnd(b : AlignmentBlock) : Int = {
    b.getReadStart - 1 + b.getLength();
  }
  
  /*
   * FINISH ME!
   * (lengthToIntersect, overlapAmt, overlapSpliceCt, code)
   */  
  def getReadLengthToIntersectAndOtherInfo(r1 : SAMRecord, r2 : SAMRecord) : (Int, Int, Int, String) = {
    val r1blocks = r1.getAlignmentBlocks.toVector;
    val r2blocks = r2.getAlignmentBlocks.toVector;
    val (r1start, r1end) = ( getAlignmentStart(r1), getAlignmentEnd(r1) );
    val (r2start, r2end) = ( getAlignmentStart(r2), getAlignmentEnd(r2) );
    //val clipLast = getLeadClipping(rLast);
    
    val r1Overlap = r1blocks.filter( (b : AlignmentBlock) => { 
      val (currStart, currEnd) = getBlockReferenceSpan(b);
      currStart < r2end && r2start < currEnd;
    });
    val r2Overlap = r2blocks.filter( (b : AlignmentBlock) => { 
      val (currStart, currEnd) = getBlockReferenceSpan(b);
      currStart < r1end && r1start < currEnd;
    });
    
    if(r1Overlap.isEmpty || r2Overlap.isEmpty){
      
    } else if(r1Overlap.length != r2Overlap.length) {
      
    } else {
      
      val blockLengthToIntersect = r2start - getBlockReferenceStart(r1Overlap.head);
      
      if(r1Overlap.length == 1){
        
        //val overlapAmt =  math.max()
        //val lengthToOverlap = getBlockReadStart(r1Overlap.head) + blockLengthToIntersect;
        //return (lengthToOverlap, overlapAmt, 0, "SingleBlockOverlap" );
      } else {
        //val blockLengthToIntersect = r2start - getBlockReferenceStart(r1Overlap.head);
      
       // return (lengthToOverlap, overlapAmt, 0, "" );
      }
    } 
    

    
    null;
  }
    

}

class qcInnerDistance(annoHolder : qcGtfAnnotationBuilder, stranded : Boolean, fr_secondStrand : Boolean, readLength : Int)  extends QCUtility[Unit] {
  reportln("Init InsertSize","progress");
  val spliceAnnotation : GenMap[(String,Char),TreeSet[(Int,Int)]] = annoHolder.spliceJunctionTreeMap;
  
  //var insertSizeMap : Map[Int,Int] = Map[Int,Int]().withDefault(i => 0);
  //var insertSizeMap_noOverlap : Map[Int,Int] = Map[Int,Int]().withDefault(i => 0);
  //var insertSizeMap_partialOverlap : Map[Int,Int] = Map[Int,Int]().withDefault(i => 0);
  //var insertSizeMap_staggeredOverlap : Map[Int,Int] = Map[Int,Int]().withDefault(i => 0);
  
  val insertSizeMap : scala.collection.mutable.Map[Int,Int] = scala.collection.mutable.Map[Int,Int]().withDefault(i => 0);
  val insertSizeMap_noOverlap : scala.collection.mutable.Map[Int,Int] = scala.collection.mutable.Map[Int,Int]().withDefault(i => 0);
  val insertSizeMap_partialOverlap : scala.collection.mutable.Map[Int,Int] = scala.collection.mutable.Map[Int,Int]().withDefault(i => 0);
  val insertSizeMap_staggeredOverlap : scala.collection.mutable.Map[Int,Int] = scala.collection.mutable.Map[Int,Int]().withDefault(i => 0);
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    val inSize = qcInnerDistance.getInsertSize(r1,r2, spliceAnnotation, stranded, fr_secondStrand);
    //reportln("InsertSize = " + inSize + "\n################################################","note");
    insertSizeMap(inSize._1) += 1;
    //insertSizeMap = insertSizeMap.updated(inSize._1, insertSizeMap(inSize._1) + 1);
    
    if(inSize._2 == 1){
      insertSizeMap_noOverlap(inSize._1) += 1;
      //insertSizeMap_noOverlap = insertSizeMap_noOverlap.updated(inSize._1, insertSizeMap_noOverlap(inSize._1) + 1);
    } else if(inSize._2 == 2){
      insertSizeMap_staggeredOverlap(inSize._1) += 1;
      //insertSizeMap_staggeredOverlap = insertSizeMap_staggeredOverlap.updated(inSize._1, insertSizeMap_staggeredOverlap(inSize._1) + 1);
    } else {
      insertSizeMap_partialOverlap(inSize._1) += 1;
      //insertSizeMap_partialOverlap = insertSizeMap_partialOverlap.updated(inSize._1, insertSizeMap_partialOverlap(inSize._1) + 1);
    }
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    val maxVal = insertSizeMap.keys.max;
    
    val sizeSeqAll = ((0 until (readLength * 3)).toSet ++  insertSizeMap.keySet).toSeq.sorted;
    val (sizeSeqBad, sizeSeq) = sizeSeqAll.span(_ < 0);
    
    
    val writer = openWriterSmart_viaGlobalParam(outfile+".insert.size.txt");
    writer.write("InsertSize	Ct\n");
    //writer.write("STAGGERED_BAD	"+insertSizeMap(qcInnerDistance.STAGGERED_BAD)+"\n")
    //writer.write("OVERLAP_CIGAR_MISMATCH	"+insertSizeMap(qcInnerDistance.OVERLAP_CIGAR_MISMATCH)+"\n");
    
    for(s  <- sizeSeq){
        writer.write(s + "	"+insertSizeMap(s)+"\n");
    }
    
    //for(i <- 0 to maxVal){
    //  writer.write(i + "	"+insertSizeMap(i)+"\n");
    //}
    writer.close();
    
    val writer2 = openWriterSmart_viaGlobalParam(outfile+".insert.size.debug.txt");
    writer2.write("InsertSize	Ct_noOverlap	Ct_staggered	Ct_partial	Ct\n");
    for(i <- sizeSeq){
      writer2.write(i + "	"+insertSizeMap_noOverlap(i) + "	"+insertSizeMap_staggeredOverlap(i) + "	"+insertSizeMap_partialOverlap(i)+"	"+insertSizeMap(i)+"\n");
    }
    writer2.close();
    
    val writer3 = openWriterSmart_viaGlobalParam(outfile+".insert.size.debug.dropped.txt");
    writer3.write("Category	Ct_noOverlap	Ct_staggered	Ct_partial	Ct_total\n");
    
    var i = qcInnerDistance.STAGGERED_NO_OVERLAP;
    writer3.write("DROPPED_FwdRead_Appears_After_RevRead_Ends	" +insertSizeMap_noOverlap(i) + "	"+insertSizeMap_staggeredOverlap(i) + "	"+insertSizeMap_partialOverlap(i)+"	"+insertSizeMap(i)+"\n");
    i = qcInnerDistance.OVERLAP_CIGAR_MISMATCH_PARTIAL_OVERLAP;
    writer3.write("DROPPED_Reads_Overlap_But_Splice_Inconsistantly_On_Overlapping_Sections	" +insertSizeMap_noOverlap(i) + "	"+insertSizeMap_staggeredOverlap(i) + "	"+insertSizeMap_partialOverlap(i)+"	"+insertSizeMap(i)+"\n");
    i = qcInnerDistance.STAGGERED_TOO_MUCH_ADAPTOR_ALIGNED;
    writer3.write("DROPPED_RevRead_Starts_Before_Fwd_And_Too_Many_Adaptor_Bases_Align_Limit_Is_"+qcInnerDistance.STAGGERED_ADAPTOR_ALIGNMENT_LIMIT+"	" +insertSizeMap_noOverlap(i) + "	"+insertSizeMap_staggeredOverlap(i) + "	"+insertSizeMap_partialOverlap(i)+"	"+insertSizeMap(i)+"\n");
    //i = qcInnerDistance.STAGGERED_FULL_BLOCK_OF_ADAPTOR_ALIGNED;
    //writer3.write("RevRead_Starts_Before_Fwd_And_RevRead_Splices_Before_Overlapping	"+"	" +insertSizeMap_noOverlap(i) + "	"+insertSizeMap_staggeredOverlap(i) + "	"+insertSizeMap_partialOverlap(i)+"	"+insertSizeMap(i)+"\n");
    close(writer3);
    
    val numGood = sizeSeq.map(i => insertSizeMap(i)).sum;
    val numBad = sizeSeqBad.map(i => insertSizeMap(i)).sum;
    summaryWriter.write("InsertSizeCalc_Kept	"+numGood+"\n");
    summaryWriter.write("InsertSizeCalc_lt_readLen	"+ Range(0,readLength).foldLeft(0)((soFar,i) => insertSizeMap(i) + soFar)+"\n");
    summaryWriter.write("InsertSizeCalc_eq_readLen	"+ insertSizeMap(readLength) +"\n");
    summaryWriter.write("InsertSizeCalc_readLen_to_2xreadLen	"+ Range(readLength+1,readLength*2).foldLeft(0)((soFar,i) => insertSizeMap(i) + soFar) +"\n");
    summaryWriter.write("InsertSizeCalc_ge_2xreadLen	"+ insertSizeMap.filter{case(s,ct) => s >= readLength*2}.foldLeft(0){case(soFar,(s,ct)) => ct + soFar} +"\n");
    
    summaryWriter.write("InsertSizeCalc_Dropped	"+numBad+"\n");
    summaryWriter.write("InsertSizeCalc_Dropped_FwdReadAppearsAfterRevReadEnds	"+insertSizeMap(qcInnerDistance.STAGGERED_NO_OVERLAP)+"\n");
    summaryWriter.write("InsertSizeCalc_Dropped_ReadsOverlapButSpliceInconsistantlyOnOverlappingSections	"+insertSizeMap(qcInnerDistance.OVERLAP_CIGAR_MISMATCH_PARTIAL_OVERLAP)+"\n");
    summaryWriter.write("InsertSizeCalc_Dropped_RevReadStartsBeforeFwdAndTooManyAdaptorBasesAlign	"+insertSizeMap(qcInnerDistance.STAGGERED_TOO_MUCH_ADAPTOR_ALIGNED)+"\n");
    //summaryWriter.write("InsertSizeCalc_Dropped_RevReadStartsBeforeFwdAndRevReadSplicesBeforeOverlapping	"+insertSizeMap(qcInnerDistance.STAGGERED_FULL_BLOCK_OF_ADAPTOR_ALIGNED)+"\n");    
  }
  /*
   *   STAGGERED_ADAPTOR_ALIGNMENT_LIMIT
   *   val OVERLAP_CIGAR_MISMATCH_INTERNAL = -999999;
  val STAGGERED_NO_OVERLAP = -111111;
  val OVERLAP_CIGAR_MISMATCH_PARTIAL_OVERLAP = -222222;
  val STAGGERED_TOO_MUCH_ADAPTOR_ALIGNED = -333333;
  val STAGGERED_FULL_BLOCK_OF_ADAPTOR_ALIGNED = -444444;
   */
  
  def getUtilityName : String = "InsertSize";
}























/*
OLD SCRAP WORK:
  def findShortestPath(chromName : String, strand : Char, start : Int, end : Int, spliceJunctionSets : HashMap[(String,Char),TreeSet[(Int,Int)]]) : Int = {
    spliceJunctionSets.get((chromName,strand)) match {
      case Some(treeSet) => findShortestPath(start,end,treeSet);
      case None => end - start;
    }
  }
  

  
  def findShortestPath(start : Int, end : Int, spliceJunctions : TreeSet[(Int,Int)]): Int = {
    val internalSplices : TreeSet[(Int,Int)] = spliceJunctions.range((start,-1),(end,-1)).filter( ab => ab._2 <= end  );
    
    if(internalSplices.size == 0) {
      return end - start;
    } else {
     // findShortestPath_helper(start,end,0,internalSplices);
    }
    -1;
  }
  private def findShortestPath_init(start : Int, end : Int, spliceJunctions : TreeSet[(Int,Int)]){
    
    val spliceJunctionsBetween = findSplicesBetween(start,end, spliceJunctions);
    
    var minSoFar = end - start;
    
    for(sj <- spliceJunctionsBetween){
      val currDistToSplice = sj._1 - start;
      if(minSoFar > currDistToSplice){
//        minSoFar = 
      }
    }
    
  }
  
  private def findShortestPath_helper(currPos : Int, end : Int, minSoFar : Int, distanceSoFar : Int, spliceJunctionsBetween : TreeSet[(Int,Int)], spliceJunctionsChecked : scala.collections.mutable.HashMap[(Int,Int),Int] ) : Int = {
    if(spliceJunctionsBetween.size == 0) return math.min(end - currPos, minSoFar);
    if(distanceSoFar > minSoFar) return minSoFar;
    
    val currSplice = spliceJunctionsBetween.head
    if(spliceJunctionsChecked.contains(currSplice)){
      val currDist = spliceJunctionsChecked(currSplice);
      if(currDist < distanceSoFar){
        
      }
      
      
      //return findShortestPath_helper()
      
      //return findShortestPath_helper(currPos, end, minSoFar, distanceSoFar, spliceJunctionsBetween.tail, spliceJunctionsChecked);
    } 
    
    
  }*/







