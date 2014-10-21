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


object qcCigarMatch {

  val NO_OVERLAPPING_CIGAROP = 0;
  val OVERLAP = 1;
  val OVERLAP_MATCH = 2;
  val OVERLAP_MISMATCH = 3;
  val OVERLAP_MISMATCH_DIFFNUM = 4;
  val OVERLAP_MISMATCH_DIFFOP = 5;
  val NUM_BINS = 6;
  
  def mergeComparisons(allCompareList : Seq[Seq[Int]], opCheckArray : Array[Int]) : Boolean = {
      if( allCompareList.forall( comp => comp.contains(NO_OVERLAPPING_CIGAROP) ) ){
        opCheckArray(NO_OVERLAPPING_CIGAROP) += 1;
        return true;
      } else {
        opCheckArray(OVERLAP) += 1;
        if(allCompareList.forall(comp => ! comp.contains(OVERLAP_MISMATCH)  )  ){
          opCheckArray(OVERLAP_MATCH) += 1;
          return true;
        } else {
          opCheckArray(OVERLAP_MISMATCH) += 1;
          if( allCompareList.exists( comp => comp.contains(OVERLAP_MISMATCH_DIFFNUM) ) ){
            opCheckArray(OVERLAP_MISMATCH_DIFFNUM) += 1;
          }
          if( allCompareList.exists( comp => comp.contains(OVERLAP_MISMATCH_DIFFOP) ) ) {
            opCheckArray(OVERLAP_MISMATCH_DIFFOP) += 1;
          }
          return false;
        }
      }
  }
  
  
  def addComparison(comp : Seq[Int], opCheckArray : Array[Int]){
    comp.foreach(x => opCheckArray(x) += 1);
  }
  
  def compareOverlappingOps(cigFoverlapping : Vector[CigOp], cigRoverlapping : Vector[CigOp], op : CigarOperator) : Seq[Int] = {
      
      val spliceF = cigFoverlapping.filter((c)=> c.op == op);
      val spliceR = cigRoverlapping.filter((c)=> c.op == op);
    
      if(spliceF.length != spliceR.length){
        //spliceCheck_overlap += 1;
        //spliceCheck_overlap_mismatch_differentNumberOfSplices += 1;
        
        return List(OVERLAP, OVERLAP_MISMATCH, OVERLAP_MISMATCH_DIFFNUM);
        
      } else if(spliceF.length > 1){
        //spliceCheck_overlap += 1;
        if(spliceF.zip(spliceR).forall{ case(cF, cR)  => { cF.ref_start == cR.ref_start && cF.ref_end == cR.ref_end } }  ){
          //spliceCheck_overlap_match += 1;
          return List(OVERLAP, OVERLAP_MATCH);
        } else {
          //spliceCheck_overlap_mismatch += 1;
          //spliceCheck_overlap_mismatch_differentSplices += 1;
          return List(OVERLAP, OVERLAP_MISMATCH, OVERLAP_MISMATCH_DIFFOP);
        }
      } else {
        //spliceCheck_noOverlappingSplices += 1;
        return List(NO_OVERLAPPING_CIGAROP);
      }
  }
  
  /*
  def getOverlapLengthOLD(cigFoverlapping : Vector[CigOp], cigRoverlapping : Vector[CigOp]) : Int = {
    val miF : Array[(Int,Int)] = cigFoverlapping.filter(c => c.op == CigarOperator.INSERTION || c.op == CigarOperator.MATCH_OR_MISMATCH).map(c => (c.ref_start,c.ref_end) ).toArray;
    val miR : Array[(Int,Int)] = cigRoverlapping.filter(c => c.op == CigarOperator.INSERTION || c.op == CigarOperator.MATCH_OR_MISMATCH).map(c => (c.ref_start,c.ref_end) ).toArray;
    
    miF(0) = (math.max(miF(0)._1, miR(0)._1) , miF(0)._2);
    miR(0) = (math.max(miF(0)._1, miR(0)._1) , miR(0)._2);
    
    miF(miF.length - 1) = ( miF.last._1, math.min(miF.last._2, miR.last._2));
    miR(miR.length - 1) = ( miR.last._1, math.min(miF.last._2, miR.last._2));
    
    val overlapLen = miF.foldLeft(0){ case(soFar, (s,e)) => (soFar + (e - s)) }
    
    if(overlapLen != miR.foldLeft(0){ case(soFar, (s,e)) => (soFar + (e - s)) }){
      error("IMPOSSIBLE STATE: Overlap lengths do not match!");
    }
    
    return overlapLen;
  }
  
  def getOverlapLength(cigFoverlapping : Vector[CigOp], cigRoverlapping : Vector[CigOp]) : Int = {
    val miF                     = cigFoverlapping.filter(c => c.op == CigarOperator.INSERTION || c.op == CigarOperator.MATCH_OR_MISMATCH);
    val miR : Vector[(Int,Int)] = cigRoverlapping.filter(c => c.op == CigarOperator.INSERTION || c.op == CigarOperator.MATCH_OR_MISMATCH).map(c => (c.ref_start,c.ref_end) );
    
    //val miB : Vector[(Int,Int)] = miF.filter{ case(fa,fb) => miR.exists{case(ra,rb) => (fa < rb && ra < fb) } }.foldLeft()
    val overlapLen = miF.foldLeft(0){case(soFar,op) => {
      val fullOverlap = miR.filter{ case(ra,rb) => op.ref_start <= ra && op.ref_end <= rb};
      
      
      0;
    }};
    
    return overlapLen;
  }*/
   
  /*
   * NO! WRONG WRONG WRONG!
   */
  def getOverlapLength(bf : Vector[AlignmentBlock], br : Vector[AlignmentBlock], readLength : Int) : (Int,Int,Int,Int) = {
    

    
    //val (startF, endF) = (getAlignmentStart(rf), getAlignmentEnd(rf));
    //val (startR, endR) = (getAlignmentStart(rr), getAlignmentEnd(rr));
    
    val (posAtStartF, posAtStartR) = getFirstOverlapPositions(bf,br,readLength);
    val (posAtEndF, posAtEndR) = getLastOverlapPositions(bf,br,readLength);
    
    //val insertSize : Int = posAtLocusF + readLength - (readLength - posAtLocusR) + 1;
    val insertSize : Int = posAtStartF + posAtStartR + 1;
    val insertSize_alt : Int = posAtEndF + posAtEndR + 1;
    
    val overlapLength : Int = 1 + posAtEndF - posAtStartF;
    val overlapLength_alt : Int = 1 + posAtStartR - posAtEndR;
    
    return (insertSize,insertSize_alt, overlapLength,overlapLength_alt);
    
    //val firstOverlapF = bf.find( b =>  startR < getBlockReferenceEnd(b) );
    //val firstOverlapR = br.find( b =>  startF < getBlockReferenceEnd(b) );
    //val lastOverlapF = bf.reverseIterator.find( b =>  getBlockReferenceStart(b) < endR );
    //val lastOverlapR = br.reverseIterator.find( b =>  getBlockReferenceStart(b) < endF );
    
    //val distToOverlapF = getAlignmentStart(rr) - getBlockReferenceStart(firstOverlapF);
    
  }
    def doBlocksIntersect(a : AlignmentBlock, b : AlignmentBlock) : Boolean = {
      getBlockReferenceStart(a) < getBlockReferenceEnd(b) && getBlockReferenceStart(b) < getBlockReferenceEnd(a);
    }  
  def getFirstOverlapPositions(b1 :  Vector[AlignmentBlock], b2 :  Vector[AlignmentBlock], readLength : Int) : (Int,Int) = {

    val over1 = b1.find( a => b2.exists( b => doBlocksIntersect(a,b) )  ).get;
    val over2 = b2.find( a => doBlocksIntersect(a,over1) ).get;
    
    val intersectLocus = math.max(getBlockReferenceStart(over1), getBlockReferenceStart(over2));
    
    val distToLocus1 = intersectLocus - getBlockReferenceStart(over1);
    val distToLocus2 = intersectLocus - getBlockReferenceStart(over2);
    
    val posAtLocus1 = getBlockReadStart(over1) + distToLocus1;
    val posAtLocus2 = readLength - (getBlockReadStart(over2) + distToLocus2+1);
    
    return (posAtLocus1, posAtLocus2);
  }
  def getLastOverlapPositions(b1 :  Vector[AlignmentBlock], b2 :  Vector[AlignmentBlock], readLength : Int) : (Int,Int) = {
    def doBlocksIntersect(a : AlignmentBlock, b : AlignmentBlock) : Boolean = {
      getBlockReferenceStart(a) < getBlockReferenceEnd(b) && getBlockReferenceStart(a) < getBlockReferenceEnd(b);
    }
    
    val over1 = b1.reverse.find( a => b2.exists( b => doBlocksIntersect(a,b) )  ).get;
    val over2 = b2.reverse.find( a => doBlocksIntersect(a,over1) ).get;
    
    val intersectLocus = math.min(getBlockReferenceEnd(over1) - 1, getBlockReferenceEnd(over2) - 1);
    
    val distToLocus1 = intersectLocus - getBlockReferenceStart(over1);
    val distToLocus2 = intersectLocus - getBlockReferenceStart(over2);
    
    val posAtLocus1 = getBlockReadStart(over1) + distToLocus1;
    val posAtLocus2 = readLength -  (getBlockReadStart(over2) + distToLocus2+1);
    
    return (posAtLocus1, posAtLocus2);
  }
  
}



class qcCigarMatch(readLen : Int) extends QCUtility[Unit] {
  reportln("Init cigarMatch","progress");
  reportln("WARNING: This function is INCOMPLETE, and not ready for general use!","progress");
  
  //COUNTING BINS:
  var staggeredBad = 0;
  var noOverlap = 0;
  var overlap = 0;
  var misOverlap = 0;
  
  var insertSize_agree = 0;
  var insertSize_disagree = 0;
  var insertSize_noMismatch_disagree = 0;
  var overlap_agree = 0;
  var overlap_disagree = 0;
  var overlap_noMismatch_disagree = 0;
  
  var insertSizeMap = Array.ofDim[Int](readLen * 2 + 1,11);
  var overlapMap = Array.ofDim[Int](readLen + 1,9);
  
  val spliceCheckArray = new Array[Int](qcCigarMatch.NUM_BINS);
  val delCheckArray = new Array[Int](qcCigarMatch.NUM_BINS);
  val insCheckArray = new Array[Int](qcCigarMatch.NUM_BINS);
  ////val clipCheckArray = new Array[Int](qcCigarMatch.NUM_BINS);
  val allCheckArray = new Array[Int](qcCigarMatch.NUM_BINS);
  
  ////val byOverlapCheckArray = Array.ofDim[Int](4,qcCigarMatch.NUM_BINS,readLen + 1);
  ////val overlapBaseCt = new Array[Int](readLen + 1);

  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    
    summaryWriter.write("CigChk_staggered	"+staggeredBad+"\n");
    summaryWriter.write("CigChk_noAlignedOverlap	"+noOverlap+"\n");
    summaryWriter.write("CigChk_AlignedOverlapMiss	"+misOverlap+"\n");
    summaryWriter.write("CigChk_AlignedOverlap	"+overlap+"\n");
    
    /*
    summaryWriter.write("CigChk_SJ_noOverlap	"+spliceCheckArray(qcCigarMatch.NO_OVERLAPPING_CIGAROP)+"\n");
    summaryWriter.write("CigChk_SJ_overlap	"+spliceCheckArray(qcCigarMatch.OVERLAP)+"\n");
    summaryWriter.write("CigChk_SJ_overlap_match	"+spliceCheckArray(qcCigarMatch.OVERLAP_MATCH)+"\n");
    summaryWriter.write("CigChk_SJ_overlap_mismatch	"+spliceCheckArray(qcCigarMatch.OVERLAP_MISMATCH)+"\n");
    summaryWriter.write("CigChk_SJ_overlap_mismatch_diffnum	"+spliceCheckArray(qcCigarMatch.OVERLAP_MISMATCH_DIFFNUM)+"\n");
    summaryWriter.write("CigChk_SJ_overlap_mismatch_diffops	"+spliceCheckArray(qcCigarMatch.OVERLAP_MISMATCH_DIFFOP)+"\n");
    */
    
    summaryWriter.write("CigChk_I_noOverlap	"+                insCheckArray(qcCigarMatch.NO_OVERLAPPING_CIGAROP)+"\n");
    summaryWriter.write("CigChk_I_overlap	"+                insCheckArray(qcCigarMatch.OVERLAP)+"\n");
    summaryWriter.write("CigChk_I_overlap_match	"+            insCheckArray(qcCigarMatch.OVERLAP_MATCH)+"\n");
    summaryWriter.write("CigChk_I_overlap_mismatch	"+        insCheckArray(qcCigarMatch.OVERLAP_MISMATCH)+"\n");
    summaryWriter.write("CigChk_I_overlap_mismatch_diffnum	"+insCheckArray(qcCigarMatch.OVERLAP_MISMATCH_DIFFNUM)+"\n");
    summaryWriter.write("CigChk_I_overlap_mismatch_diffops	"+insCheckArray(qcCigarMatch.OVERLAP_MISMATCH_DIFFOP)+"\n");
    
    summaryWriter.write("CigChk_D_noOverlap	"+                delCheckArray(qcCigarMatch.NO_OVERLAPPING_CIGAROP)+"\n");
    summaryWriter.write("CigChk_D_overlap	"+                delCheckArray(qcCigarMatch.OVERLAP)+"\n");
    summaryWriter.write("CigChk_D_overlap_match	"+            delCheckArray(qcCigarMatch.OVERLAP_MATCH)+"\n");
    summaryWriter.write("CigChk_D_overlap_mismatch	"+        delCheckArray(qcCigarMatch.OVERLAP_MISMATCH)+"\n");
    summaryWriter.write("CigChk_D_overlap_mismatch_diffnum	"+delCheckArray(qcCigarMatch.OVERLAP_MISMATCH_DIFFNUM)+"\n");
    summaryWriter.write("CigChk_D_overlap_mismatch_diffops	"+delCheckArray(qcCigarMatch.OVERLAP_MISMATCH_DIFFOP)+"\n");
    
    summaryWriter.write("CigChk_ALL_noOverlap	"+                allCheckArray(qcCigarMatch.NO_OVERLAPPING_CIGAROP)+"\n");
    summaryWriter.write("CigChk_ALL_overlap	"+                    allCheckArray(qcCigarMatch.OVERLAP)+"\n");
    summaryWriter.write("CigChk_ALL_overlap_match	"+            allCheckArray(qcCigarMatch.OVERLAP_MATCH)+"\n");
    summaryWriter.write("CigChk_ALL_overlap_mismatch	"+        allCheckArray(qcCigarMatch.OVERLAP_MISMATCH)+"\n");
    summaryWriter.write("CigChk_ALL_overlap_mismatch_diffnum	"+allCheckArray(qcCigarMatch.OVERLAP_MISMATCH_DIFFNUM)+"\n");
    summaryWriter.write("CigChk_ALL_overlap_mismatch_diffops	"+allCheckArray(qcCigarMatch.OVERLAP_MISMATCH_DIFFOP)+"\n");
    
    /*
    summaryWriter.write("CigChk_InsertSize_agree	"+insertSize_agree+"\n");
    summaryWriter.write("CigChk_InsertSize_disagree	"+insertSize_disagree+"\n");
    summaryWriter.write("CigChk_InsertSize_noMismatch_disagree	"+insertSize_noMismatch_disagree+"\n");
    summaryWriter.write("CigChk_overlap_agree	"+overlap_agree+"\n");
    summaryWriter.write("CigChk_overlap_disagree	"+overlap_disagree+"\n");
    summaryWriter.write("CigChk_overlap_noMismatch_disagree	"+overlap_noMismatch_disagree+"\n");

    val writer =  openWriterSmart_viaGlobalParam(outfile + ".overlap.insert.size.debug.txt");
    writer.write("InsertSize	COUNT_FwdFirst_Start	COUNT_FwdFirst_End	COUNT_RevFirst_Start	COUNT_RefFirst_End	COUNT_FwdFirst_match	COUNT_FwdFirst_mismatch	COUNT_RevFirst_match	COUNT_RevFirst_mismatch	COUNT_match\n"); 
    for(i <- Range(0,insertSizeMap.length)){
      val curr = insertSizeMap(i);
      writer.write(i+"	"+curr.mkString("	")+"\n");
    }
    close(writer);
    
    val writer2 =  openWriterSmart_viaGlobalParam(outfile + ".overlap.lengths.txt");
    writer2.write("overlap	COUNT_FwdFirst_Start	COUNT_FwdFirst_End	COUNT_RevFirst_Start	COUNT_RefFirst_End	COUNT_FwdFirst_match	COUNT_FwdFirst_mismatch	COUNT_RevFirst_match	COUNT_RevFirst_mismatch	COUNT_match\n"); 
    for(i <- Range(0,overlapMap.length)){
      val curr = overlapMap(i);
      writer2.write(i+"	"+curr.mkString("	")+"\n");
    }
    close(writer2);
    */
  }
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    checkCigarMatch(r1,r2);
  }
  
  def checkCigarMatch(r1 : SAMRecord, r2 : SAMRecord){
    val (rf,rr) = if(r1.getReadNegativeStrandFlag) ((r2,r1)) else ((r1,r2));
    
    val (startF, endF) = (getAlignmentStart(rf), getAlignmentEnd(rf));
    val (startR, endR) = (getAlignmentStart(rr), getAlignmentEnd(rr));

    if(startF > endR){
      //reportln("Method Exit: STAGGERED_BAD! ","note");
      staggeredBad += 1;
      //Count STAGGERED_BAD
    } else if(endF <= startR){
      //reportln("Method Entry: getInnerDistance_noOverlap: ","note");
      noOverlap += 1;
      //Count NO_OVERLAP
    } else {
      //COUNT OVERLAP
      val bf : Vector[AlignmentBlock] = rf.getAlignmentBlocks().toVector;
      val br : Vector[AlignmentBlock] = rr.getAlignmentBlocks().toVector;      
      val blocksOverlap = bf.exists( a => br.exists( b => qcCigarMatch.doBlocksIntersect(a,b) )  );
      if(! blocksOverlap){
        misOverlap += 1;  
      } else {
    
      overlap += 1;
      
      val cigF = CigarHolder(rf,false,false).cigOps.toVector;   
      val cigR = CigarHolder(rr,false,false).cigOps.toVector;   
      //cigop variables: op, length, ref_start, ref_end, read_start, read_end, chromName, strand, ref_iv, read_iv; 
      
      ////Old version, includes all ops that intersect with the shared region.
      //val cigFoverlapping = cigF.filter((c) => (c.ref_start < endR && startR < c.ref_end)); 
      //val cigRoverlapping = cigR.filter((c) => (c.ref_start < endF && startF < c.ref_end)); 
      ////New version, includes only ops that are entirely inside the shared region.
      val cigFoverlapping = cigF.filter((c) => (c.ref_end <= endR && startR <= c.ref_start)); 
      val cigRoverlapping = cigR.filter((c) => (c.ref_end <= endF && startF <= c.ref_start));       
      
      //val spliceCompare = qcCigarMatch.compareOverlappingOps(cigFoverlapping,cigRoverlapping, CigarOperator.SKIPPED_REGION);
      val delCompare =    qcCigarMatch.compareOverlappingOps(cigFoverlapping,cigRoverlapping, CigarOperator.DELETION);
      val insCompare =    qcCigarMatch.compareOverlappingOps(cigFoverlapping,cigRoverlapping, CigarOperator.INSERTION);
      //val clipCompare = qcCigarMatch.compareOverlappingOps(cigFoverlapping,cigRoverlapping, CigarOperator.SOFT_CLIP);

      //qcCigarMatch.addComparison(spliceCompare, spliceCheckArray);
      qcCigarMatch.addComparison(delCompare, delCheckArray);
      qcCigarMatch.addComparison(insCompare, insCheckArray);
      //qcCigarMatch.addComparison(clipCompare, clipCheckArray);
      
      
      //val allCompareList : Seq[Seq[Int]] = Vector[Seq[Int]](spliceCompare, delCompare, insCompare)//, clipCompare);
      val allCompareList : Seq[Seq[Int]] = Vector[Seq[Int]](delCompare, insCompare);
      qcCigarMatch.mergeComparisons(allCompareList, allCheckArray);
      
      /*
      val isMatch = qcCigarMatch.mergeComparisons(allCompareList, allCheckArray);
      
  
      val (insertSize,insertSize_alt, overlapLength,overlapLength_alt) = qcCigarMatch.getOverlapLength(bf , br , readLen);
      
      if(startF < startR){ 
        insertSizeMap(insertSize)(0) += 1;
        insertSizeMap(insertSize_alt)(1) += 1;
        overlapMap(overlapLength)(0) += 1;
        overlapMap(overlapLength_alt)(1) += 1;
        if(insertSize == insertSize_alt){
          insertSizeMap(insertSize)(4) += 1;
          insertSizeMap(insertSize)(8) += 1;
        } else {
          insertSizeMap(insertSize)(5) += 1;
        }
        if(overlapLength == overlapLength_alt){
          overlapMap(overlapLength)(4) += 1;
          overlapMap(overlapLength)(8) += 1;
        } else {
          overlapMap(overlapLength)(5) += 1;
        }
        
        val otherMethodInsertSize = qcInnerDistance.getInsertSize_partialOverlap(rf,rr);
        if(mismatchCt < 10000 && insertSize != otherMethodInsertSize){
          reportln("NONMATCHING INSERT SIZE, PARTIALOVERLAP: " + insertSize + " vs " + otherMethodInsertSize,"debug");
          reportln(rf.getSAMString().trim,"debug");
          reportln(rr.getSAMString().trim,"debug");
          mismatchCt += 1;
        }
        insertSizeMap(insertSize)(9) += 1;
        
      } else {
        insertSizeMap(insertSize)(2) += 1;
        insertSizeMap(insertSize_alt)(3) += 1;
        overlapMap(overlapLength)(2) += 1;
        overlapMap(overlapLength_alt)(3) += 1;
        if(insertSize == insertSize_alt){
          insertSizeMap(insertSize)(6) += 1;
          insertSizeMap(insertSize)(8) += 1;
        } else {
          insertSizeMap(insertSize)(7) += 1;
        }
        if(overlapLength == overlapLength_alt){
          overlapMap(overlapLength)(6) += 1;
          overlapMap(overlapLength)(8) += 1;
        } else {
          overlapMap(overlapLength)(7) += 1;
        }
      }
      
      if(insertSize == insertSize_alt){
        insertSize_agree += 1;
      } else {
        insertSize_disagree += 1;
        if(isMatch){
          insertSize_noMismatch_disagree += 1;
        }
      }
      if(overlapLength == overlapLength_alt){
        overlap_agree += 1;
      } else {
        overlap_disagree += 1;
        if(isMatch){
          overlap_noMismatch_disagree += 1;
        }
      }
      
        val otherMethodInsertSize = qcInnerDistance.getInsertSize_staggeredOverlap(rf,rr);
        if(mismatchCt < 10000 && insertSize != otherMethodInsertSize){
          reportln("NONMATCHING INSERT SIZE, PARTIALOVERLAP: " + insertSize + " vs " + otherMethodInsertSize,"debug");
          reportln(rf.getSAMString().trim,"debug");
          reportln(rr.getSAMString().trim,"debug");
          mismatchCt += 1;
        }
        insertSizeMap(insertSize)(10) += 1;
      
      
      //val overlapCt = qcCigarMatch.getOverlapLength(cigFoverlapping,cigRoverlapping);
      //overlapBaseCt(overlapCt) += 1;
      */
    }}
  }
  
  var mismatchCt = 0;

  
  def getUtilityName : String = "cigarMatch";

  
}


