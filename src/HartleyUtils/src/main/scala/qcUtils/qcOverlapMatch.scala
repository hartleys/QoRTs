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

object qcOverlapMatch {

  val ALL_REAL_BASE_CHAR : Seq[Char] = Vector('A','C','G','T');
  
  val ALL_BASE_BYTES : Array[Byte] = Array('A'.toByte,'C'.toByte,'G'.toByte,'T'.toByte,'N'.toByte);
  val BYTE_CODE_MAP : Array[Int] = Array.ofDim(ALL_BASE_BYTES.max + 1);
  BYTE_CODE_MAP('A'.toInt) = 0;
  BYTE_CODE_MAP('C'.toInt) = 1;
  BYTE_CODE_MAP('G'.toInt) = 2;
  BYTE_CODE_MAP('T'.toInt) = 3;
  BYTE_CODE_MAP('N'.toInt) = 4;
  //qcOverlapMatch.BYTE_CODE_MAP
  final val NO_OVERLAP_CODE = 0;
  final val NO_OVERLAP_STAGGERED_CODE = 1;
  final val OVERLAP_CIGAR_MISMATCH = 2;
  final val OVERLAP_MATCH = 3;
  final val OVERLAP_MISMATCH = 4;
  final val OVERLAP_CODE_LIST = List(NO_OVERLAP_CODE,NO_OVERLAP_STAGGERED_CODE,OVERLAP_CIGAR_MISMATCH,OVERLAP_MATCH,OVERLAP_MISMATCH);
  
  final val OVERLAP_CODE_MAP = Map[String,Int](
      "NO_OVERLAP_CODE" -> NO_OVERLAP_CODE,
      "NO_OVERLAP_STAGGERED_CODE" -> NO_OVERLAP_STAGGERED_CODE,
      "OVERLAP_CIGAR_MISMATCH"-> OVERLAP_CIGAR_MISMATCH,
      "OVERLAP_MATCH"-> OVERLAP_MATCH,
      "OVERLAP_MISMATCH"-> OVERLAP_MISMATCH
      );
  final val OVERLAP_CODE_MAP_INV : Map[Int,String] = OVERLAP_CODE_MAP.map{case (a,b) => (b,a)}.toMap;
  
  final val INDEL_CODE_MAP = Map[String,Int](
        "None" -> 0,
        "r1" -> 1,
        "r2" -> 2,
        "rBoth" -> 3
      );
  
  
  def getTrimmedOverlapRegion(rA : SAMRecord, rB : SAMRecord, offset : Int) : Option[((Int,Int),(Int,Int))] = {
    val (startB,endB) = (0,math.min(rB.getReadLength(),rA.getReadLength()-offset));
    val overlapLen = endB;
    val (startA,endA) = (offset,offset+overlapLen);
    
    val leadTrimA = math.max(startA,getLeadClipping(rA)) - offset;
    val leadTrimB = getLeadClipping(rB);
    val tailTrimA = endA - math.min(endA,rA.getReadLength() - getTailClipping(rA));
    val tailTrimB = endB - math.min(endB,rB.getReadLength() - getTailClipping(rB));
    
    val leadTrim = math.max(leadTrimA,leadTrimB);
    val tailTrim = math.max(tailTrimA,tailTrimB);
    
    if(leadTrim + tailTrim >= overlapLen){
      //overlap eliminated by trimming!
      return None;
    } else {
      val (startBt, endBt) = (startB+leadTrim,endB - tailTrim);
      val (startAt, endAt) = (startA+leadTrim,endA - tailTrim);
      
      return Some((startAt, endAt),(startBt, endBt));
    }
  }
  
  
  
  
}


class qcOverlapMatch(readLen : Int, mismatchSamFile : Option[String], 
                     maxQualScore : Int, adjustPhredScore : Int, 
                     isSingleEnd : Boolean,
                     calcReferenceMatch : Boolean = false, 
                     genomeFa : Option[Seq[String]] = None,
                     genomeBufferSize : Int = 20000,
                     calcOverlapMismatch : Boolean = true,
                     calcOverlapIndel    : Boolean = true,
                     var writeOverlapCoverage : Boolean = true,
                     var writeQualScoreErrorRate : Boolean = true,
                     var writeQualScoreOverlapCoverage : Boolean = true,
                     var writeReadLengthDistro : Boolean = true) extends QCUtility[Unit] {
  reportln("> Init OverlapMatch Utility","debug");
  
  val samwriter : Option[WriterUtil] = if(mismatchSamFile.isEmpty) None else Some(openWriterSmart(mismatchSamFile.get));
   
  val genomeSeq = if(genomeFa.isEmpty) null else buildEfficientGenomeSeqContainer(genomeFa.get);
  val bufferSize = genomeBufferSize;
  
  var ppoLengthMismatch = 0;
  var ppoReferencePositionMismatch = 0;
  var ppoMatch = 0;
  
  
  //basePairCount(strand)(pos1)(pos2)(base1)(base2) = count
  val basePairCount      : Array[Array[Array[Array[Array[Int]]]]] = Array.ofDim(2,readLen, readLen,qcOverlapMatch.ALL_BASE_BYTES.length, qcOverlapMatch.ALL_BASE_BYTES.length);
  val basePairCountClean : Array[Array[Array[Array[Array[Int]]]]] = Array.ofDim(2,readLen, readLen,qcOverlapMatch.ALL_BASE_BYTES.length, qcOverlapMatch.ALL_BASE_BYTES.length);
  
  val overlapMismatchArrayAlned      : Array[Array[Array[Array[Array[Int]]]]] = Array.ofDim(2,readLen, readLen,qcOverlapMatch.ALL_BASE_BYTES.length, qcOverlapMatch.ALL_BASE_BYTES.length);
  val overlapMismatchArrayAlnedClean : Array[Array[Array[Array[Array[Int]]]]] = Array.ofDim(2,readLen, readLen,qcOverlapMatch.ALL_BASE_BYTES.length, qcOverlapMatch.ALL_BASE_BYTES.length);

  
  val qualBins : Array[Int] = Array(0,10,20,40);
  
  //(quality bin)(hasIndel?)(strand?)(pos1)(pos2)(base1)(base2)
  
  //Redo more efficiently:
  //val basePairCountByQualBin : Array[Array[Array[Array[Array[Array[Array[Int]]]]]]] = qualBins.map(q => {
  //  Array(0,1).map(k => Array.ofDim[Int](2,readLen, readLen,qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1, qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1));
  //})
  
  val overlapCoverage : Array[Array[Long]] = Array.ofDim[Long](2,readLen);
  val overlapCoverageClean : Array[Array[Long]] = Array.ofDim[Long](2,readLen);

  val readLenDist : Array[Array[Long]] = Array.ofDim[Long](2,readLen+1);
  var noOverlap_staggered = 0;
  var noOverlap_normal = 0;
  var overlap = 0;
  
  var BAD_OVERLAP = 0;
  
  var PERFECTMATCH = 0;
  var MISMATCH = 0;
  
  val perfectMatchRate : Array[Long] = Array.ofDim[Long](readLen + 1);
  
  var PAIR_CONTAINS_INDEL : Long = 0;
  var PAIR_CONTAINS_NO_INDEL : Long= 0;
  
  var PAIR_CONTAINS_DEL : Long = 0;
  var PAIR_CONTAINS_INS : Long = 0;
  var PAIR_CONTAINS_INSANDDEL : Long = 0;
  
  //misMatchCount(indel true/false,# bases overlap,# bases mismatch) = count
  val misMatchCount : Array[Array[Array[Long]]] = Array.ofDim[Long](2,readLen+1,readLen+1);
  
  //referenceMismatch(hasIndel)(read num)(strand)(read pos)(refBase)(readBase)
  val referenceMismatch        : Array[Array[Array[Array[Array[Int]]]]] = Array.ofDim[Int](2,2,readLen+1,qcOverlapMatch.ALL_BASE_BYTES.length,qcOverlapMatch.ALL_BASE_BYTES.length);
  //val referenceMismatchNoIndel : Array[Array[Array[Array[Array[Int]]]]] = Array.ofDim[Int](2,2,readLen+1,qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1,qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1);

  
  //referenceMismatchCount(read num)(num bases mismatched)
  val referenceMismatchCount : Array[Array[Long]] = Array.ofDim[Long](2,readLen+1);
  
  val PERFECTREFMATCH   : Array[Long] = Array.ofDim[Long](2);
  val IMPERFECTREFMATCH : Array[Long] = Array.ofDim[Long](2);
  
  def convertCbPosition(r : SAMRecord, pos : Int) = if(r.getReadNegativeStrandFlag()) r.getReadLength() - pos - 1 else pos;
  
  var chromSeq : Vector[String] = Vector[String]();
  
  //flagCountArray(hasDel)(hasIns)(hasMismatchR1, f/r 0,1,2)(hasMismatchR2, f/r 0,1,2)(hasOverlapMismatch)
  val flagCountArray : Array[Array[Array[Array[Array[Long]]]]] = Array.ofDim[Long](qcOverlapMatch.INDEL_CODE_MAP.map(_._2).max+1,
                                                                                 qcOverlapMatch.INDEL_CODE_MAP.map(_._2).max+1,
                                                                                 3,3,qcOverlapMatch.OVERLAP_CODE_LIST.max+1);
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    val (rf,rr) = if(r1.getReadNegativeStrandFlag) ((r2,r1)) else ((r1,r2));
    val (unstartF, unendF) = (getAlignmentStart(rf) - getLeadClipping(rf), getAlignmentEnd(rf) + getTailClipping(rf));
    val (unstartR, unendR) = (getAlignmentStart(rr) - getLeadClipping(rr), getAlignmentEnd(rr) + getTailClipping(rr));
    
    var (hasDel,hasIns,mmidx1,mmidx2,omidx) = (0,0,0,0,0);
    
    if(chromSeq.isEmpty) chromSeq = chromSeq :+ r1.getReferenceName()
    if(r1.getReferenceName() != chromSeq.last){
      chromSeq = chromSeq :+ r1.getReferenceName()
    }
    
    readLenDist(0)(r1.getReadLength()) += 1;
    readLenDist(1)(r2.getReadLength()) += 1;
    val cb1 = getCigarBlocksFromRead(r1).toVector;
    val cb2 = getCigarBlocksFromRead(r2).toVector;
    val (cbf,cbr) = if(r1.getReadNegativeStrandFlag) ((cb2,cb1)) else ((cb1,cb2));
    
    hasDel = if(cb1.exists(cb => cb.op == CigarOperator.DELETION)){
      if(cb2.exists(cb => cb.op == CigarOperator.DELETION ))  2 else 1
    } else 0;
    hasIns = if(cb1.exists(cb => cb.op == CigarOperator.INSERTION)){
      if(cb2.exists(cb => cb.op == CigarOperator.INSERTION )) 2 else 1
    } else 0;
    
    //hasDel = if(cb1.exists(cb => cb.op == CigarOperator.DELETION)  || cb2.exists(cb => cb.op == CigarOperator.DELETION )) 1 else 0;
    //hasIns = if(cb1.exists(cb => cb.op == CigarOperator.INSERTION) || cb2.exists(cb => cb.op == CigarOperator.INSERTION)) 1 else 0;
    if(hasDel > 0) PAIR_CONTAINS_DEL += 1;
    if(hasIns > 0) PAIR_CONTAINS_INS += 1;
    if(hasDel > 0 && hasIns > 0) PAIR_CONTAINS_INSANDDEL += 1;
    
    val hasIndel = if(hasDel > 0 || hasIns > 0){
      PAIR_CONTAINS_INDEL += 1;
      true;
    } else {
      PAIR_CONTAINS_NO_INDEL += 1;
      false;
    }

    if(calcReferenceMatch && ! genomeFa.isEmpty){
      try{
        val minStart = math.min(r1.getAlignmentStart()-1,r2.getAlignmentStart()-1);
        val chrom = r1.getReferenceName()
        
        //if(readNum < 1000 || readNum % 100000 == 0) genomeSeq.reportBufferStatus;
        if( readNum % 1000000 == 1) genomeSeq.reportBufferStatus;
        
        //val mismatchArray = if(hasIndel) referenceMismatch else referenceMismatchNoIndel
        val mismatchArray = referenceMismatch;
        genomeSeq.shiftBufferTo(chrom,minStart-bufferSize);
        val cb1f = cb1.filter(c => c.op.consumesReadBases() && c.op.consumesReferenceBases())
        val cb2f = cb2.filter(c => c.op.consumesReadBases() && c.op.consumesReferenceBases())  
        
        val ircb = if(r1.getReadPairedFlag()){
          Vector((0,r1,cb1f),(1,r2,cb2f))
        } else {
          Vector((0,r1,cb1f))
        }
        
        for((i,r,cb) <- ircb){
          val refMM = calcRefMismatches(chrom , i ,r , cb , mismatchArray);
          //calcQualScoreCoverage(i, r, cb);
          if(i == 0) mmidx1 = if(mmidx1 > 0) mmidx1 else if(refMM){ if(r.getReadNegativeStrandFlag()) 2 else 1 } else 0;
          if(i == 1) mmidx2 = if(mmidx2 > 0) mmidx2 else if(refMM){ if(r.getReadNegativeStrandFlag()) 2 else 1 } else 0;
        }
      } catch {
        case e: Exception => {
          reportln("Caught exception in calcReferenceMatch on read "+readNum+"\n"+
                   "   offending read1:\n"+r1.getSAMString()+""+
                   "   offending read2:\n"+r2.getSAMString(),"note");
          throw e;
        }
      }
    }
    if(calcOverlapMismatch){
      omidx = overlapMismatchCalcs(rf=rf , rr=rr , cbf=cbf, cbr=cbr,hasIndel=hasIndel)
    }
    
    flagCountArray(hasDel)(hasIns)(mmidx1)(mmidx2)(omidx) += 1;
  }
  
  def overlapMismatchCalcs(rf : SAMRecord, rr : SAMRecord, cbf : Vector[CigarBlock], cbr : Vector[CigarBlock],hasIndel : Boolean) : Int = {
    if(getAlignmentStart(rf) >= getAlignmentEnd(rr)){
      noOverlap_staggered += 1;
      return qcOverlapMatch.OVERLAP_CODE_MAP("NO_OVERLAP_STAGGERED_CODE");
    } else if(getAlignmentEnd(rf) <= getAlignmentStart(rr)){
      noOverlap_normal += 1;
      return qcOverlapMatch.OVERLAP_CODE_MAP("NO_OVERLAP_CODE");
    } else {
      val (r1,r2,cb1,cb2) = if(rf.getSecondOfPairFlag()) (rr,rf,cbr,cbf) else (rf,rr,cbf,cbr);
      
      val (rA,rB,offset) = getPairOffset(rf,rr,cbf,cbr);
      
      if(offset == OVERLAP_CIGAR_MISMATCH_INTERNAL){
        BAD_OVERLAP+= 1;
        writeToSam(title = "BADOVERLAP",samwriter=samwriter,r1=r1, r2=r2, seqA=rf.getReadString().toVector, seqB=rr.getReadString().toVector,offset=offset)
        return qcOverlapMatch.OVERLAP_CODE_MAP("OVERLAP_CIGAR_MISMATCH");
      } else {
        //val offset = lenToIX - clipB;
        overlap+=1;
        
        val (startB,endB) = (0,math.min(rB.getReadLength(),rA.getReadLength()-offset));
        val overlapLen = endB;
        val (startA,endA) = (offset,offset+overlapLen);
        //val overlapLen = endA - startA;
        
        //REWRITE THIS, THIS WAY:
        val swap = rA.getSecondOfPairFlag();
        val ((start1,end1),(start2,end2)) = if(swap){
          ((startB,endB),(startA,endA))
        } else {
          ((startA,endA),(startB,endB))
        }
        
        //if(calcOverlapIndel && hasIndel){
        //  calculateOverlapIndel(r1,r2,cb1,cb2);
        //}
        
        val seq1 = r1.getReadString().toVector.slice(start1,end1);
        val seq2 = r2.getReadString().toVector.slice(start2,end2);
        
        val mismatches = seq1.zip(seq2).zipWithIndex.filter{case ((x,y),i) => {(x != 'N' && y != 'N' && x != y)}}
        
        //if(! mismatches.isEmpty){

        //}
        if(writeQualScoreOverlapCoverage){
          calcQualOverlapCoverage(r1=r1,r2=r2,start1=start1,end1=end1,start2=start2,end2=end2,hasIndel = hasIndel);
        }
        
        if(writeOverlapCoverage){
          calcOverlapCoverage(r1=r1, r2=r2 , 
                             start1 =start1, end1 =end1, 
                             start2 =start2, end2 =end2,
                             hasIndel =hasIndel);
        }
        
        if(mismatches.isEmpty){
          perfectMatchRate(overlapLen) += 1;
          PERFECTMATCH += 1;
          return qcOverlapMatch.OVERLAP_CODE_MAP("OVERLAP_MATCH");
          //writeToSam(title = "PERFECTMATCH",samwriter=samwriter,r1=r1, r2=r2, seqA=seqA, seqB=seqB,offset=offset)
        } else {
          MISMATCH += 1;
          val indelIdx = if(hasIndel) 1 else 0;
          misMatchCount(indelIdx)(overlapLen)(mismatches.size) += 1;
          
          if(writeQualScoreErrorRate){
            calcMismatchQuality(mismatches=mismatches,r1,r2,
                                start1=start1,end1=end1,
                                start2=start2,end2=end2,hasIndel = hasIndel);
          }
          
          calculateOverlapMismatch(mismatches =mismatches,
                           r1 =r1,r2=r2,
                           start1=start1,end1=end1,
                           start2=start2,end2=end2, 
                           hasIndel =hasIndel);
          
          if(! mismatchSamFile.isEmpty){
            val (seqA,seqB) = if(swap) (seq2,seq1) else (seq1,seq2)
            writeToSam(title = "MISMATCH",samwriter=samwriter,r1=r1, r2=r2, seqA=seqA, seqB=seqB,offset=offset)
          }
          return qcOverlapMatch.OVERLAP_CODE_MAP("OVERLAP_MISMATCH");
        }
      }
    }
  }
  
  
  
  val NO_OVERLAPPING_INDELS = 0;
  
  //overlapIndelCountArray(ins/del)(strand)(r1Pos)(r2Pos)
  val overlapIndelCountArray = Array.ofDim[Long](2,2,readLen+1,readLen+1);
  
  //overlapMismatchIndelCountArray(ins/del)(strand)(r1/r2)(pos);
  val overlapMismatchIndelCountArray = Array.ofDim[Long](2,2,2,readLen+1);
  
  var OVERLAP_NoOverlappingIndels : Long= 0;
  var OVERLAP_MultipleOverlappingIndels : Long= 0;
  var OVERLAP_MatchingIndels : Long= 0;
  var OVERLAP_WrongTypeMismatchIndel : Long= 0;
  var OVERLAP_SlightMismatchIndel : Long= 0;
  var OVERLAP_TotalMismatchIndel : Long= 0;
  
  var OVERLAP_MissingIndelFromR1 : Long= 0;
  var OVERLAP_MissingIndelFromR2 : Long= 0;
  
  
  //missingOverlapIndel(in/del)(r1/r2)(strand)(readPos)
  val missingOverlapIndel = Array.ofDim[Long](2,2,2,readLen+1);
  
  val totalMismatchMargin = 6;
  
  def calculateOverlapIndel(r1 : SAMRecord,r2: SAMRecord,cb1 : Vector[CigarBlock],cb2: Vector[CigarBlock]){
    val (start1,end1) = (r1.getAlignmentStart()-1,r1.getAlignmentEnd());
    val (start2,end2) = (r2.getAlignmentStart()-1,r2.getAlignmentEnd());
    val (start,end) = (math.max(start1,start2), math.min(end1,end2));
    
    val indel1 = cb1.filter(c => c.op == CigarOperator.DELETION || c.op == CigarOperator.INSERTION).filter(c => c.refStart >= start && c.refEnd <= end);
    val indel2 = cb2.filter(c => c.op == CigarOperator.DELETION || c.op == CigarOperator.INSERTION).filter(c => c.refStart >= start && c.refEnd <= end);
    
    if(indel1.isEmpty && indel2.isEmpty){
      //no overlapping indels!
      OVERLAP_NoOverlappingIndels += 1;
    } else if(indel1.size > 1 || indel2.size > 1){
      OVERLAP_MultipleOverlappingIndels += 1;
    } else if(! indel1.isEmpty){
      val r = r1;
      val indel = indel1.head;
      val indelCode = if(indel.op == CigarOperator.DELETION) 1 else 0;
      val strandCode = if(r.getReadNegativeStrandFlag()) 1 else 0;
      val readCode = 0;
      val pos = if(strandCode == 1) r.getReadLength() - indel.readEnd else indel.readStart;
      overlapMismatchIndelCountArray(indelCode)(readCode)(strandCode)(pos) += 1;
      OVERLAP_MissingIndelFromR1 += 1;
    } else if(! indel2.isEmpty){
      val r = r2;
      val indel = indel2.head;
      val indelCode = if(indel.op == CigarOperator.DELETION) 1 else 0;
      val strandCode = if(r.getReadNegativeStrandFlag()) 1 else 0;
      val readCode = 1;
      val pos = if(strandCode == 1) r.getReadLength() - indel.readEnd else indel.readStart;
      overlapMismatchIndelCountArray(indelCode)(readCode)(strandCode)(pos) += 1;
      OVERLAP_MissingIndelFromR2 += 1;
    } else {
      //check if overlaps match
      val i1 = indel1.head;
      val i2 = indel2.head;
      if(i1.op == i2.op && i1.refStart == i2.refStart && i1.refEnd == i2.refEnd){
        OVERLAP_MatchingIndels += 1;
      } else if(i1.op != i2.op){
        OVERLAP_WrongTypeMismatchIndel += 1;
      } else if(i1.refStart - totalMismatchMargin < i2.refEnd && i2.refStart - totalMismatchMargin < i1.refEnd && i1.len - totalMismatchMargin < i2.len && i2.len - totalMismatchMargin < i1.len){
        OVERLAP_SlightMismatchIndel += 1;
      } else {
        OVERLAP_TotalMismatchIndel += 1;
      }
    }
    
  }
  
  
  //overlapMismatchByQuality(hasIndel)(qual1)(qual2)
  val overlapMismatchByQuality : Array[Array[Array[Long]]] = Array.ofDim[Long](2,maxQualScore+1,maxQualScore+1);
  val overlapQualityCoverage   : Array[Array[Array[Long]]] = Array.ofDim[Long](2,maxQualScore+1,maxQualScore+1);
  
  //overlapMismatchByQualityBySeq(hasIndel)(strand)(qual1)(qual2)(base1)(base2)
  val overlapMismatchByQualityBySeq : Array[Array[Array[Array[Array[Array[Long]]]]]] = Array.ofDim[Long](2).map(buf =>{
    Array.ofDim[Long](2,maxQualScore+1,maxQualScore+1,qcOverlapMatch.ALL_BASE_BYTES.length,qcOverlapMatch.ALL_BASE_BYTES.length);
  })
  //Array.ofDim[Long](2,maxQualScore+1,maxQualScore+1,qcOverlapMatch.ALL_BASE_BYTES.length,qcOverlapMatch.ALL_BASE_BYTES.length);
/*
 * 
qualBins.map(q => {
  //  Array(0,1).map(k => Array.ofDim[Int](2,readLen, readLen,qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1, qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1));
  //})
 */
  
  //(quality bin)(hasIndel?)(strand?)(pos1)(pos2)(base1)(base2)
  //val basePairCountByQualBin : Array[Array[Array[Array[Array[Array[Array[Int]]]]]]] = qualBins.map(q => {
  //  Array(0,1).map(k => Array.ofDim[Int](2,readLen, readLen,qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1, qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1));
  //})
  
  
  def calcMismatchQuality(mismatches : Vector[((Char,Char),Int)],
                           r1 : SAMRecord,r2: SAMRecord,
                           start1:Int,end1:Int,
                           start2:Int,end2:Int, 
                           hasIndel : Boolean){
    val (q1,q2) = (r1.getBaseQualities(),r2.getBaseQualities());
    val strandIdx = if(r1.getReadNegativeStrandFlag()) 1 else 0;
    
    val getMismatchQuals : Vector[(Int,Int,Int,Int,Int)] = mismatches.map{case ((a,b),i) =>{
        //delete later! testing:
        assert( a == r1.getReadString().charAt(start1 + i) );
        assert( b == r2.getReadString().charAt(start2 + i) );
        (q1(start1 +i).toInt - adjustPhredScore , q2(start2 +i).toInt - adjustPhredScore, i, qcOverlapMatch.BYTE_CODE_MAP(a.toInt),qcOverlapMatch.BYTE_CODE_MAP(b.toInt));
    }}
    
    getMismatchQuals.foreach{case (x,y,i,a,b) =>{
      overlapMismatchByQuality(0)(x)(y) += 1;
      overlapMismatchByQualityBySeq(0)(strandIdx)(x)(y)(a)(b) += 1;
    }}
    if(!hasIndel){
      getMismatchQuals.foreach{case (x,y,i,a,b) =>{
        overlapMismatchByQuality(1)(x)(y) += 1;
        overlapMismatchByQualityBySeq(1)(strandIdx)(x)(y)(a)(b) += 1;
      }}
    }
    

    
    
    val mm = mismatches;
    val mismatchQualBins = getMismatchQuals.map{case (q1,q2,i,a,b) =>{
      qualBins.length - 1 - math.max(qualBins.reverse.indexWhere(m => m <= q1),qualBins.reverse.indexWhere(m => m <= q2))
    }}
    
    val mismatchInfo = mm.zip(mismatchQualBins);
    
    
    
    //basePairCountByQualBin(quality bin)(hasIndel?)(strand?)(pos1)(pos2)(base1)(base2)
    //TO DO: REDO more efficiently!
    /*if(r1.getReadNegativeStrandFlag()){
      mismatchInfo.foreach{case (((x,y),i),bin) => {
        basePairCountByQualBin(bin)(0)(1)(r1.getReadLength() - (start1 + i+1))(start2+i)(x.toInt)(y.toInt) += 1;
        if(! hasIndel) basePairCountByQualBin(bin)(1)(1)(r1.getReadLength() - (start1 + i+1))(start2+i)(x.toInt)(y.toInt) += 1;
      }}
    } else {
      mismatchInfo.foreach{case (((x,y),i),bin) => {
        basePairCountByQualBin(bin)(0)(0)(start1+i)(r2.getReadLength() - (start2+i+1))(x.toInt)(y.toInt) += 1;
        if(! hasIndel) basePairCountByQualBin(bin)(1)(0)(start1+i)(r2.getReadLength() - (start2+i+1))(x.toInt)(y.toInt) += 1;
      }}
    }*/
  }
  
  def calculateOverlapMismatch(mismatches : Vector[((Char,Char),Int)],
                           r1 : SAMRecord,r2: SAMRecord,
                           start1:Int,end1:Int,
                           start2:Int,end2:Int, 
                           hasIndel : Boolean){
    
    //basePairCount(readIdx)(pos1)(pos2)(base1)(base2) = count
    if(r1.getReadNegativeStrandFlag()){
      mismatches.foreach{case ((x,y),i) => {
        basePairCount(1)(r1.getReadLength() - (start1 + i+1))(start2+i)(qcOverlapMatch.BYTE_CODE_MAP(x.toInt))(qcOverlapMatch.BYTE_CODE_MAP(y.toInt)) += 1;
        if(! hasIndel) basePairCountClean(1)(r1.getReadLength() - (start1 + i+1))(start2+i)(qcOverlapMatch.BYTE_CODE_MAP(x.toInt))(qcOverlapMatch.BYTE_CODE_MAP(y.toInt)) += 1;
      }}
    } else {
      mismatches.foreach{case ((x,y),i) => {
        basePairCount(0)(start1+i)(r2.getReadLength() - (start2+i+1))(qcOverlapMatch.BYTE_CODE_MAP(x.toInt))(qcOverlapMatch.BYTE_CODE_MAP(y.toInt)) += 1;
        if(! hasIndel) basePairCountClean(0)(start1+i)(r2.getReadLength() - (start2+i+1))(qcOverlapMatch.BYTE_CODE_MAP(x.toInt))(qcOverlapMatch.BYTE_CODE_MAP(y.toInt)) += 1;
      }}
    }
  }
  
  def calcQualOverlapCoverage(r1 : SAMRecord,r2: SAMRecord,
                              start1:Int,end1:Int,start2:Int,end2:Int, hasIndel : Boolean){
    
    val (q1,q2) = (r1.getBaseQualities(),r2.getBaseQualities());
    //if(q1.exists(q => q - adjustPhredScore > maxQualScore) || q2.exists(q => q - adjustPhredScore > maxQualScore)){
    /*  reportln("Qual String:","debug");
      reportln("   "+q1.mkString(","),"debug");
      reportln("   "+q2.mkString(","),"debug");
      reportln("   ("+start1+","+end1+") ("+start2+","+end2+")","debug");
      reportln("   swap = "+swap,"debug");
      reportln("   r1.len = "+r1.getReadLength(),"debug");
      reportln("   r2.len = "+r2.getReadLength(),"debug");
      reportln("   qlen1 = " +r1.getBaseQualities().length,"debug");
      reportln("   qlen2 = " +r2.getBaseQualities().length,"debug");*/
    //}
    (start1 until end1).zip(start2 until end2).foreach{case (i,j) =>{
      overlapQualityCoverage(0)(q1(i).toInt - adjustPhredScore)(q2(j).toInt - adjustPhredScore) += 1;
    }}
    if(!hasIndel){
      (start1 until end1).zip(start2 until end2).foreach{case (i,j) =>{
        overlapQualityCoverage(1)(q1(i).toInt - adjustPhredScore)(q2(j).toInt - adjustPhredScore) += 1;
      }}
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////
  
  //(read1/read2)(qualScore)
  val refMismatchByQuality : Array[Array[Long]] = Array.ofDim[Long](2,maxQualScore+1);
  val refQualityCoverage   : Array[Array[Long]] = Array.ofDim[Long](2,maxQualScore+1);
  
  //refMismatchByQualityByBase(read1/read2)(strand)(qualScore)(refBase)(readBase)
  val refMismatchByQualityByBase : Array[Array[Array[Array[Array[Long]]]]] = Array.ofDim[Long](2,2,maxQualScore+1,qcOverlapMatch.ALL_BASE_BYTES.length+1,qcOverlapMatch.ALL_BASE_BYTES.length+1);
  
  val mappedByPositionCoverage : Array[Array[Long]] = Array.ofDim[Long](2,readLen+1);
  
  def calcQualScoreCoverage(readNum : Int, r : SAMRecord, cb : Vector[CigarBlock]){
    val qs = r.getBaseQualityString();
    
    cb.foreach(b => {
      qs.substring(b.readStart,b.readEnd).foreach( c => {
        refQualityCoverage(readNum)(c.toInt - 33 - adjustPhredScore);
      })
    })
    
    //takes all bases:
    //r.getBaseQualityString().substring(getLeadClipping(r),r.getReadLength() - getTailClipping(r)).foreach( c=>{
    //  refQualityCoverage(readNum)(c.toInt - 33 - adjustPhredScore);
    //})
  }
  
  def calcRefMismatches(chrom : String, i : Int,r : SAMRecord, cb : Vector[CigarBlock], mismatchArray : Array[Array[Array[Array[Array[Int]]]]]) : Boolean = {
    val seqIter = cb.flatMap{ c => {
      genomeSeq.getSeqForInterval(chrom,c.refStart,c.refEnd).toVector.zip(
          getSeqStringFromBlock(r,c).toVector
      ).zip(c.readStart until (c.readStart + c.len) );
      //(refSeq + genomeSeq.getSeqForInterval(chrom,c.refStart,c.refEnd),readSeq+getSeqStringFromBlock(r1,c))
    }};
    val mismatches = seqIter.filter{case ((refchar,readchar),pos) => {readchar != 'N' && refchar != 'N' && refchar != readchar}};
    
    //DEBUGGING:
    if((! mismatches.isEmpty) && (! samwriter.isEmpty)){
      val writer = samwriter.get;
      writer.write("Mismatch report ("+r.getReadName()+"/"+(i+1)+")"+"\n");
      writer.write("   "+r.getReferenceName()+":"+r.getAlignmentStart+", "+r.getCigarString()+"\n");
      val mismatchSet = mismatches.map(_._2).toSet;
      val seqIter = cb.foreach{ c => {
        writer.write("    "+c.toString()+"\n");
        writer.write("        rf: "+genomeSeq.getSeqForInterval(chrom,c.refStart,c.refEnd)+"\n");
        writer.write("        rd: "+getSeqStringFromBlock(r,c)+"\n");
        writer.write("            ");
        (c.readStart until c.readEnd).foreach(j => {
          if(mismatchSet.contains(j))writer.write("^") else writer.write(" ");
        });
        writer.write("\n");
      }};
      //if(mismatches.isEmpty){
      //  writer.write("    No Mismatches!\n");
      //} else {
        writer.write("    Mismatches: "+mismatches.map{case ((a,b),i) =>{"("+a+","+b+","+i+")"}}.mkString(",")+"\n");
      //}
    }
    val qs = r.getBaseQualities();
    
    if(writeQualScoreErrorRate){
      cb.foreach(b => {
        qs.slice(b.readStart,b.readEnd).foreach( c => {
          refQualityCoverage(i)(c.toInt - adjustPhredScore) += 1;
        })
      })
    }
    
    cb.foreach(b => {
      (b.readStart until b.readEnd).foreach(x => {
        mappedByPositionCoverage(i)(x) += 1;
      })
    })
    val strandIdx = if(r.getReadNegativeStrandFlag()) 1 else 0;
    
    if(mismatches.isEmpty){
      PERFECTREFMATCH(i) += 1;
      return false;
    } else {
      IMPERFECTREFMATCH(i) += 1;
      mismatches.foreach{case ((refchar,readchar),pos) =>{
        //referenceMismatchCount(read num)(strand)(read pos)(refBase)(readBase)
        mismatchArray(i)(strandIdx)(convertCbPosition(r,pos))(qcOverlapMatch.BYTE_CODE_MAP(refchar.toInt))(qcOverlapMatch.BYTE_CODE_MAP(readchar.toInt)) += 1;
      }}
      
      if(writeQualScoreErrorRate){
        mismatches.foreach{case ((refchar,readchar),pos) =>{
          val q = qs(pos).toInt - adjustPhredScore;
          refMismatchByQuality(i)(q) += 1;
          refMismatchByQualityByBase(i)(strandIdx)(q)(qcOverlapMatch.BYTE_CODE_MAP(refchar.toInt))(qcOverlapMatch.BYTE_CODE_MAP(readchar.toInt)) += 1;
        }}
      }
      return true;
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////

  def calcOverlapCoverage(r1 : SAMRecord, r2 : SAMRecord, 
                          start1 : Int, end1 : Int, 
                          start2 : Int, end2 : Int,
                          hasIndel : Boolean){
    //val ((r1,start1,end1),(r2,start2,end2)) = if(rA.getSecondOfPairFlag()){
    //  ((rB,startB,endB),(rA,startA,endA))
    //} else {
    //  ((rA,startA,endA),(rB,startB,endB))
    //}
    
    val (rango1,rango2) = if(r1.getReadNegativeStrandFlag()){
      ((r1.getReadLength - end1 until r1.getReadLength - start1),
       (start2 until end2))
    } else {
      ((start1 until end1),
       (r2.getReadLength - end2 until r2.getReadLength - start2))
    }
    
    rango1.foreach(i => overlapCoverage(0)(i) += 1);
    rango2.foreach(i => overlapCoverage(1)(i) += 1);
    if(! hasIndel){
      rango1.foreach(i => overlapCoverageClean(0)(i) += 1);
      rango2.foreach(i => overlapCoverageClean(1)(i) += 1);
    }
  }
  /////////////////////////////////////////////////////////////////////////////

  //UNFINISHED!?
  /*
  def calcOverlapMismatch(mismatches : Vector[((Char,Char),Int)],
                          rA : SAMRecord, rB : SAMRecord, 
                          startA : Int, endA : Int, 
                          startB : Int, endB : Int,
                          hasIndel : Boolean){
    
    val swap = rA.getSecondOfPairFlag();
    val ((r1,start1,end1),(r2,start2,end2)) = if(swap){
      ((rB,startB,endB),(rA,startA,endA))
    } else {
      ((rA,startA,endA),(rB,startB,endB))
    }
    
    val mm = if(swap){
      mismatches.map{case ((a,b),i) => { ((b,a),i) }}
    } else mismatches;
    
    if(r1.getReadNegativeStrandFlag()){
      mismatches.foreach{case ((x,y),i) => {
        basePairCount(1)(r1.getReadLength() - (start1 + i+1))(start2+i)(x.toInt)(y.toInt) += 1;
        if(! hasIndel) basePairCountClean(1)(r1.getReadLength() - (start1 + i+1))(start2+i)(x.toInt)(y.toInt) += 1;
      }}
    } else {
      mismatches.foreach{case ((x,y),i) => {
        basePairCount(0)(start1+i)(r2.getReadLength() - (start2+i+1))(x.toInt)(y.toInt) += 1;
        if(! hasIndel) basePairCountClean(0)(start1+i)(r2.getReadLength() - (start2+i+1))(x.toInt)(y.toInt) += 1;
      }}
    }
  }*/
  
  /////////////////////////////////////////////////////////////////////////////

  final val OVERLAP_CIGAR_MISMATCH_INTERNAL = qcInnerDistance.OVERLAP_CIGAR_MISMATCH_INTERNAL;
  
  def getPairOffset(rA : SAMRecord, rB : SAMRecord, cbAraw : Seq[CigarBlock], cbBraw : Seq[CigarBlock]) : (SAMRecord,SAMRecord,Int) = {
    val cbA = cbAraw.filter(c => c.op.consumesReadBases() && c.op.consumesReferenceBases()).iterator;
    val cbB = cbBraw.filter(c => c.op.consumesReadBases() && c.op.consumesReferenceBases()).iterator;
    var A = cbA.next;
    var B = cbB.next;
    
    //reportln("starting with: "+rA.getReadName(),"debug");
    
    while(true){
      //reportln("   A = "+A.toString + " B = "+B.toString,"debug");
      if(A.refEnd < B.refStart){
        if(cbA.hasNext) {
          A = cbA.next;
          //reportln("   new A = "+A.toString,"debug");
        } else {
          //reportln("   no cbA next!","debug");
          return (rA,rB,OVERLAP_CIGAR_MISMATCH_INTERNAL);
        }
      } else if(B.refEnd < A.refStart){
        if(cbB.hasNext) {
          B = cbB.next;
          //reportln("   new B = "+B.toString,"debug");
        } else {
          //reportln("   no cbB next!","debug");
          return (rA,rB,OVERLAP_CIGAR_MISMATCH_INTERNAL);
        }
      } else {
        val (c1,c2,r1,r2) = if(A.refStart <= B.refStart) (A,B,rA,rB) else (B,A,rB,rA);
        val refOffset = c2.refStart - c1.refStart;
        val offset = c1.readStart + refOffset - c2.readStart;
        if(offset < 0) return (r2,r1,-offset);
        else return (r1,r2,offset);
      }
    }

    return (rA,rB,OVERLAP_CIGAR_MISMATCH_INTERNAL);
  }
  
  
  
  def writeToSam(title : String,samwriter : Option[WriterUtil],r1 : SAMRecord, r2 : SAMRecord, seqA : Vector[Char], seqB : Vector[Char],offset : Int){
    if(! samwriter.isEmpty) {
      samwriter.get.write(title+": \n"+
                           "    "+seqA.mkString("")+"\n"+
                           "    vs\n"+
                           "    "+seqB.mkString("")+"\n"+
                           "    ("+offset+")\n");
      samwriter.get.write( "    "+r1.getSAMString()+"");
      samwriter.get.write( "    "+r2.getSAMString()+"");
      samwriter.get.flush();
    } //else do nothing!
  }
  
  val fwdBases = Vector('A','C','G','T');
  val revBases = Vector('T','G','C','A');
  val bases = fwdBases.zip(revBases);
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    val pairs = bases.flatMap{case (af,ar) => bases.map{case (bf,br) => ((af,ar),(bf,br))}}.filter{case ((af,ar),(bf,br)) => af != bf}.map{case ((af,ar),(bf,br)) => (af,qcOverlapMatch.BYTE_CODE_MAP(af.toInt),ar,qcOverlapMatch.BYTE_CODE_MAP(ar.toInt),bf,qcOverlapMatch.BYTE_CODE_MAP(bf.toInt),br,qcOverlapMatch.BYTE_CODE_MAP(br.toInt))};
    //val mismatchedBasePairs = bases.flatMap{case (af,ar) => bases.map{case (bf,br) => { (af,ar,bf,br) }} }
    
    //for((af,afi,ar,ari,bf,bfi,br,bri) <- pairs){
    //  reportln("   ("+af+","+afi+","+ar+","+ari+"),("+bf+","+bfi+","+br+","+bri+")","debug");
    //}
    
      //(read1/read2)(qualScore)
  //val refMismatchByQuality : Array[Array[Int]] = Array.ofDim[Int](2,maxQualScore+1);
  //val refQualityCoverage   : Array[Array[Int]] = Array.ofDim[Int](2,maxQualScore+1);
  
  
    if(writeQualScoreErrorRate && calcReferenceMatch){
        val writer = createOutputFile(outfile , "referenceMismatch.byScore.txt","",docWriter,
             ("Score","int","The PHRED quality score."),
             ("MismatchCt_R1","int","The number of read-1 reference mismatches with the given PHRED score."),
             ("MismatchCt_R2","int","The number of read-2 reference mismatches with the given PHRED score."),
             ("QualCoverage_R1","int","The total number of aligned read-1 bases that have the given PHRED score."),
             ("QualCoverage_R2","int","The total number of aligned read-2 bases that have the given PHRED score.")
        );
        writer.write("Score\tMismatchCt_R1\tMismatchCt_R2\tQualCoverage_R1\tQualCoverage_R2\n");
        for(i <- 0 to maxQualScore){
            writer.write(i+"\t"+
                         refMismatchByQuality(0)(i)+"\t"+
                         refMismatchByQuality(1)(i)+"\t"+
                         refQualityCoverage(0)(i)+"\t"+
                         refQualityCoverage(1)(i)+""+
                         "\n");
        }
        close(writer);
    }
//refMismatchByQualityByBase(read1/read2)(strand)(qualScore)(refBase)(readBase)
    if(writeQualScoreErrorRate && calcReferenceMatch){
        val writer = createOutputFile(outfile , "referenceMismatch.byScoreAndBP.txt","",docWriter,
             ("Score","int","The PHRED quality score."),
             ("REFBASE","int","The reference base-pair (as it should have appeared on the read)."),
             ("READBASE","int","The read base-pair."),
             ("CT_R1_FWD","int","As CT_R1, but only counting FWD strand aligned reads."),
             ("CT_R1_REV","int","As CT_R1, but only counting REV strand aligned reads."),
             ("CT_R2_FWD","int","As CT_R2, but only counting FWD strand aligned reads."),
             ("CT_R2_REV","int","As CT_R2, but only counting REV strand aligned reads."),
             ("CT_R1","int","The number of base-pair swaps of the given type with the given PHRED quality score for read 1."),
             ("CT_R2","int","The number of base-pair swaps of the given type with the given PHRED quality score for read 2.")
        );
        writer.write("Score\tREFBASE\tREADBASE\tCT_R1_FWD\tCT_R1_REV\tCT_R2_FWD\tCT_R2_REV\tCT_R1\tCT_R2\n");
        for(i <- 0 to maxQualScore){
            for((af,afi,ar,ari,bf,bfi,br,bri) <- pairs){
              writer.write(i+"\t"+af+"\t"+bf+"\t"+
                          refMismatchByQualityByBase(0)(0)(i)(afi)(bfi)+"\t"+
                          refMismatchByQualityByBase(0)(1)(i)(ari)(bri)+"\t"+
                          refMismatchByQualityByBase(1)(0)(i)(afi)(bfi)+"\t"+
                          refMismatchByQualityByBase(1)(1)(i)(ari)(bri)+"\t"+
                         (refMismatchByQualityByBase(0)(0)(i)(afi)(bfi)+
                          refMismatchByQualityByBase(0)(1)(i)(ari)(bri))+"\t"+
                         (refMismatchByQualityByBase(1)(0)(i)(afi)(bfi)+
                          refMismatchByQualityByBase(1)(1)(i)(ari)(bri))+
                         "\n");
            }
        }
        close(writer);
    }
    
    /*
        writer.write("POS\tREFBASE\tREADBASE\tCT_R1_FWD\tCT_R1_REV\tCT_R2_FWD\tCT_R2_REV\tCT_R1\tCT_R2\n");

              for((af,afi,ar,ari,bf,bfi,br,bri) <- pairs){
            writer.write(i+"\t"+af+"\t"+bf+"\t"+
                         referenceMismatch(0)(0)(i)(afi)(bfi)+"\t"+
                         referenceMismatch(0)(1)(i)(ari)(bri)+"\t"+
                         referenceMismatch(1)(0)(i)(afi)(bfi)+"\t"+
                         referenceMismatch(1)(1)(i)(ari)(bri)+"\t"+
                        (referenceMismatch(0)(0)(i)(afi)(bfi) + referenceMismatch(0)(1)(i)(ari)(bri)) + "\t"+
                        (referenceMismatch(1)(0)(i)(afi)(bfi) + referenceMismatch(1)(1)(i)(ari)(bri)) +
                         "\n");
          }
     */
    
    
    if(writeQualScoreErrorRate && calcOverlapMismatch){
      //overlapMismatchByQuality(hasIndel)(qual1)(qual2)
        val writer = createOutputFile(outfile , "overlapMismatch.byScore.txt","",docWriter,
             ("ScoreR1","int","The PHRED quality score of Read 1."),
             ("ScoreR2","int","The PHRED quality score of Read 2."),
             ("OverlapCt","int","The number of read1/read2 overlap base pairs with the given PHRED quality scores."),
             ("OverlapCtNoIndel","int","The number of read1/read2 overlap base pairs with the given PHRED quality scores, if you discount all base-pairs that contain INDELs."),
             ("MismatchCt","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores."),
             ("MismatchCtNoIndel","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores, if you discount all base-pairs that contain INDELs.")
        );
        writer.write("ScoreR1\tScoreR2\tOverlapCt\tOverlapCtNoIndel\tMismatchCt\tMismatchCtNoIndel\n");
        for(i <- 0 to maxQualScore){
          for(j <- 0 to maxQualScore){
            writer.write(i+"\t"+j+"\t"+
                         overlapQualityCoverage(0)(i)(j)+"\t"+
                         overlapQualityCoverage(1)(i)(j)+"\t"+
                         overlapMismatchByQuality(0)(i)(j)+"\t"+
                         overlapMismatchByQuality(1)(i)(j)+""+
                         "\n");
          }
        }
        close(writer);
    }
    
    if(writeQualScoreErrorRate && calcOverlapMismatch){
      //overlapMismatchByQuality(hasIndel)(qual1)(qual2)
      //overlapMismatchByQualityBySeq(hasIndel)(strandr1)(qual1)(qual2)(base1)(base2)
        val writer = createOutputFile(outfile , "overlapMismatch.byScoreAndBP.txt","",docWriter,
             ("ScoreR1","int","The PHRED quality score of Read 1."),
             ("ScoreR2","int","The PHRED quality score of Read 2."),
             ("baseR1","char","The base-pair on read 1"),
             ("baseR2","char","The base-pair on read 2"),
             ("CT","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs."),
             ("CT_FWD","int","As CT, except for the fwd strand only."),
             ("CT_REV","int","As CT, except for the rev strand only."),
             ("CT_NOINDEL_FWD","int","As CT_NOINDEL, except for the fwd strand only."),
             ("CT_NOINDEL_REV","int","As CT_NOINDEL, except for the rev strand only."),
             ("CT","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs."),
             ("CT_NOINDEL","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs, if you discount all base-pairs that contain INDELs.") 
        );
        writer.write("ScoreR1\tScoreR2\tbaseR1\tbaseR2\tCT_FWD\tCT_REV\tCT_NOINDEL_FWD\tCT_NOINDEL_REV\tCT\tCT_NOINDEL\n");
        for(i <- 0 to maxQualScore){
          for(j <- 0 to maxQualScore){
            for((af,afi,ar,ari,bf,bfi,br,bri) <- pairs){
                writer.write(i+"\t"+j+"\t"+af+"\t"+bf+"\t"+
                         overlapMismatchByQualityBySeq(0)(0)(i)(j)(afi)(bfi)+"\t"+
                         overlapMismatchByQualityBySeq(0)(1)(i)(j)(ari)(bri)+"\t"+
                         overlapMismatchByQualityBySeq(1)(0)(i)(j)(afi)(bfi)+"\t"+
                         overlapMismatchByQualityBySeq(1)(1)(i)(j)(ari)(bri)+"\t"+
                         (overlapMismatchByQualityBySeq(0)(0)(i)(j)(afi)(bfi)+overlapMismatchByQualityBySeq(0)(1)(i)(j)(ari)(bri))+"\t"+
                         (overlapMismatchByQualityBySeq(1)(0)(i)(j)(afi)(bfi)+overlapMismatchByQualityBySeq(1)(1)(i)(j)(ari)(bri))+
                         "\n");
            }
          }
        }
        close(writer);
    }
    
    if(calcOverlapMismatch){
      val writer = createOutputFile(outfile , "overlapMismatch.txt","",docWriter,
             ("READ_POSITION_R1","int","The read position on read 1."),
             ("READ_POSITION_R2","int","The read position on read 2."),
             ("baseA","char","The base-pair on read 1"),
             ("baseB","char","The base-pair on read 2"),
             ("CT","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs."),
             ("CT_FWD","int","As CT, except for the fwd strand only."),
             ("CT_REV","int","As CT, except for the rev strand only."),
             ("CT_NOINDEL_FWD","int","As CT_NOINDEL, except for the fwd strand only."),
             ("CT_NOINDEL_REV","int","As CT_NOINDEL, except for the rev strand only."),
             ("CT","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs."),
             ("CT_NOINDEL","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs, if you discount all base-pairs that contain INDELs.")
      );
      writer.write("READ_POSITION_R1\tREAD_POSITION_R2\tbaseA\tbaseB\tCT_FWD\tCT_REV\tCT\tCT_NOINDEL_FWD\tCT_NOINDEL_REV\tCT_NOINDEL\n");
      for(i <- 0 until readLen){
        for(j <- 0 until readLen){
          for((af,afi,ar,ari,bf,bfi,br,bri) <- pairs){
              //val (bfi,bri) = (bf.toInt,br.toInt);
              writer.write(i+"\t"+j+"\t"+af+"\t"+bf+"\t"+
                            basePairCount(0)(i)(j)(afi)(bfi)+"\t"+
                            basePairCount(1)(i)(j)(ari)(bri)+"\t"+
                           (basePairCount(0)(i)(j)(afi)(bfi)+
                            basePairCount(1)(i)(j)(ari)(bri))+"\t"+
                            basePairCountClean(0)(i)(j)(afi)(bfi)+"\t"+
                            basePairCountClean(1)(i)(j)(ari)(bri)+"\t"+
                           (basePairCountClean(0)(i)(j)(afi)(bfi)+
                            basePairCountClean(1)(i)(j)(ari)(bri))+"\n"
              );
          //  }
          }
        }
      }
      close(writer);
    }
    
    if(calcOverlapMismatch){
      val writer = createOutputFile(outfile , "overlapMismatch.byBase.txt","",docWriter,
             ("baseA","char","The base-pair on read 1"),
             ("baseB","char","The base-pair on read 2"),
             ("CT","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs."),
             ("CT_FWD","int","As CT, except for the fwd strand only."),
             ("CT_REV","int","As CT, except for the rev strand only."),
             ("CT_NOINDEL_FWD","int","As CT_NOINDEL, except for the fwd strand only."),
             ("CT_NOINDEL_REV","int","As CT_NOINDEL, except for the rev strand only."),
             ("CT","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs."),
             ("CT_NOINDEL","int","The number of read1/read2 overlap mismatches with the given PHRED quality scores and base-pairs, if you discount all base-pairs that contain INDELs.")
      );
      writer.write("baseA\tbaseB\tCT_FWD\tCT_REV\tCT\tCT_NOINDEL_FWD\tCT_NOINDEL_REV\tCT_NOINDEL\n");
      //for(i <- 0 until readLen){
      //  for(j <- 0 until readLen){
          for((af,afi,ar,ari,bf,bfi,br,bri) <- pairs){
              //val (bfi,bri) = (bf.toInt,br.toInt);
              writer.write(af+"\t"+bf+"\t"+
                            basePairCount(0).map(a => { a.map(b => { b(afi)(bfi)}).sum}).sum+"\t"+
                            basePairCount(1).map(a => { a.map(b => { b(ari)(bri)}).sum}).sum+"\t"+
                           (basePairCount(0).map(a => { a.map(b => { b(afi)(bfi)}).sum}).sum+
                            basePairCount(1).map(a => { a.map(b => { b(ari)(bri)}).sum}).sum)+"\t"+
                            basePairCountClean(0).map(a => { a.map(b => { b(afi)(bfi)}).sum}).sum+"\t"+
                            basePairCountClean(1).map(a => { a.map(b => { b(ari)(bri)}).sum}).sum+"\t"+
                           (basePairCountClean(0).map(a => { a.map(b => { b(afi)(bfi)}).sum}).sum+
                            basePairCountClean(1).map(a => { a.map(b => { b(ari)(bri)}).sum}).sum)+"\n"
              );
          }
         // }
        //}
      //}
      close(writer);
    }
    
    if(calcOverlapMismatch){
      val writer = createOutputFile(outfile , "overlapMismatch.byRead.txt","",docWriter,
             ("POS","int","Read position"),
             ("CT_R1","int","Number of read-1 reads with a mismatch at the given position."),
             ("CT_R2","int","Number of read-2 reads with a mismatch at the given position."),
             ("CT_NOINDEL_R1","int","Number of read-1 reads with a mismatch at the given position, not counting read-pairs with indels."),
             ("CT_NOINDEL_R2","int","Number of read-2 reads with a mismatch at the given position, not counting read-pairs with indels.")
      );
      writer.write("POS\tCT_R1\tCT_R2\tCT_NOINDEL_R1\tCT_NOINDEL_R2\n");
      for(i <- 0 until readLen){
        val missCt1 = (0 until readLen).map(j => {
          pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            basePairCount(0)(i)(j)(afi)(bfi)+basePairCount(1)(i)(j)(ari)(bri)
          }}.sum
        }).sum
        val missCt2 = (0 until readLen).map(j => {
          pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            basePairCount(0)(j)(i)(afi)(bfi)+basePairCount(1)(j)(i)(ari)(bri)
          }}.sum
        }).sum
        val missCt1NI = (0 until readLen).map(j => {
          pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            basePairCountClean(0)(i)(j)(afi)(bfi)+basePairCountClean(1)(i)(j)(ari)(bri)
          }}.sum
        }).sum
        val missCt2NI = (0 until readLen).map(j => {
          pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            basePairCountClean(0)(j)(i)(afi)(bfi)+basePairCountClean(1)(j)(i)(ari)(bri)
          }}.sum
        }).sum
        
        writer.write(i+"\t"+missCt1+"\t"+missCt2+"\t"+missCt1NI+"\t"+missCt2NI+"\n");
      }
      close(writer);
    }
    /*if(calcOverlapMismatch){
      val writer = openWriterSmart_viaGlobalParam(outfile + ".overlapMismatch.byRead.noIndel.txt");
      writer.write("POS\tCT_R1\tCT_R2\n");
      for(i <- 0 until readLen){
        val missCt1 = (0 until readLen).map(j => {
          pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            basePairCountClean(0)(i)(j)(afi)(bfi)+basePairCountClean(1)(i)(j)(ari)(bri)
          }}.sum
        }).sum
        val missCt2 = (0 until readLen).map(j => {
          pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            basePairCountClean(0)(j)(i)(afi)(bfi)+basePairCountClean(1)(j)(i)(ari)(bri)
          }}.sum
        }).sum
        writer.write(i+"\t"+missCt1+"\t"+missCt2+"\n");
      }
      close(writer);
    }*/
    
    
    if(calcOverlapMismatch){
      //misMatchCount(indelIdx)(overlapLen)(mismatches.size) += 1;
      
      val writer = createOutputFile(outfile , "mismatchSizeRates.txt","",docWriter,
             ("OVERLAP_SIZE","int","Number of overlapping bases"),
             ("BASES_MISMATCHED","int","Number of overlapping bases that mismatch one another."),
             ("CT_NOINDEL","int","Number of read-pairs with the given amount of overlap and the given number of mismatching bases, not counting read-pairs that contain INDELs"),
             ("CT","int","Number of read-pairs with the given amount of overlap and the given number of mismatching bases")
      );
      writer.write("OVERLAP_SIZE\tBASES_MISMATCHED\tCT_NOINDEL\tCT\n");
      for(i <- 0 until readLen+1){
        for(j <- 1 to i){
          writer.write(i+"\t"+j+"\t"+misMatchCount(0)(i)(j)+"\t"+(misMatchCount(0)(i)(j) + misMatchCount(1)(i)(j))+"\n");
        }
      } 
      close(writer);
    }
    
    if(writeOverlapCoverage){
      val writer = createOutputFile(outfile , "overlapCoverage.txt","",docWriter,
             ("POS","int","Read position"),
             ("CT_R1","int","Number of read-1 reads that overlap with read 2 at the given position."),
             ("CT_R2","int","Number of read-2 reads that overlap with read 1 at the given position."),
             ("CT_NOINDEL_R1","int","As CT_R1, but discounting read-pairs that contain INDELs"),
             ("CT_NOINDEL_R2","int","As CT_R2, but discounting read-pairs that contain INDELs")
      );
      writer.write("POS\tCT_R1\tCT_R2\tCT_NOINDEL_R1\tCT_NOINDEL_R2\n");
      for(i <- 0 until readLen){
        writer.write(i+"\t"+overlapCoverage(0)(i)+"\t"+overlapCoverage(1)(i)+"\t"+overlapCoverageClean(0)(i)+"\t"+overlapCoverageClean(1)(i)+"\n");
      }
      close(writer);
    }
    
    summaryWriter.write("PAIR_CONTAINS_DEL\t"+PAIR_CONTAINS_DEL+"\tNumber of read-pairs containing one or more deletions\n");
    summaryWriter.write("PAIR_CONTAINS_INS\t"+PAIR_CONTAINS_INS+"\tNumber of read-pairs containing one or more insertions\n");
    summaryWriter.write("PAIR_CONTAINS_INS_AND_DEL\t"+PAIR_CONTAINS_INSANDDEL+"\tNumber of read-pairs containing one or more insertions and deletions\n");
    summaryWriter.write("PAIR_CONTAINS_INDEL\t"+PAIR_CONTAINS_INDEL+"\tNumber of read-pairs containing one or more insertions or deletions\n");
    summaryWriter.write("PAIR_CONTAINS_NO_INDEL\t"+PAIR_CONTAINS_NO_INDEL+"\tNumber of read-pairs containing no indels\n");
    
    //hasDel: Option[Int],hasIns: Option[Int],mmidx1: Option[Int],mmidx2: Option[Int],omidx: Option[Int]
    summaryWriter.write("READ_CONTAINS_DEL_R1\t"+(getFlagSum(hasDel = Some(1)) + getFlagSum(hasDel = Some(3))) +"\tNumber of read pairs containing a del in read 1\n");
    summaryWriter.write("READ_CONTAINS_INS_R1\t"+(getFlagSum(hasIns = Some(1)) + getFlagSum(hasIns = Some(3)))+"\tNumber of read pairs containing a ins in read 2\n");
    summaryWriter.write("READ_CONTAINS_NO_INDEL_R1\t"+(getFlagSum(hasDel = Some(0),hasIns = Some(0))+
                                                  getFlagSum(hasDel = Some(0),hasIns = Some(2))+
                                                  getFlagSum(hasDel = Some(2),hasIns = Some(0))+
                                                  getFlagSum(hasDel = Some(2),hasIns = Some(2)))+
                                                 "\tNumber of read pairs containing no indels in read 1\n");
    summaryWriter.write("READ_CONTAINS_DEL_R2\t"+(getFlagSum(hasDel = Some(2)) + getFlagSum(hasDel = Some(3))) +"\tNumber of read pairs containing a del in read 2\n");
    summaryWriter.write("READ_CONTAINS_INS_R2\t"+(getFlagSum(hasIns = Some(2)) + getFlagSum(hasIns = Some(3)))+"\tNumber of read pairs containing a ins in read 2\n");
    summaryWriter.write("READ_CONTAINS_NO_INDEL_R2\t"+(getFlagSum(hasDel = Some(0),hasIns = Some(0))+
                                                  getFlagSum(hasDel = Some(0),hasIns = Some(1))+
                                                  getFlagSum(hasDel = Some(1),hasIns = Some(0))+
                                                  getFlagSum(hasDel = Some(1),hasIns = Some(1)))+
                                                 "\tNumber of read pairs containing no indels in read 2\n");
    
    if(calcReferenceMatch){
      summaryWriter.write("NO_REF_BASE_SWAP_R1\t"+PERFECTREFMATCH(0)+"\tNumber of read pairs with no reference mismatches in read 1\n");
      summaryWriter.write("HAS_REF_BASE_SWAP_R1\t"+IMPERFECTREFMATCH(0)+"\tNumber of read pairs with one or more reference mismatches in read 1\n");
      
      if(! isSingleEnd){
        summaryWriter.write("NO_REF_BASE_SWAP_R2\t"+PERFECTREFMATCH(1)+"\tNumber of read pairs with no reference mismatches in read 2\n");
        summaryWriter.write("HAS_REF_BASE_SWAP_R2\t"+IMPERFECTREFMATCH(1)+"\tNumber of read pairs with one or more reference mismatches in read 2\n");
      }
      
      val PerfectMatch = flagCountArray(0)(0)(0)(0).sum
      summaryWriter.write("PERFECT_REF_MATCH_PAIR\t"+PerfectMatch+"\tNumber of reads or read pairs with no reference mismatches\n");
      
      //summaryWriter.write("PERFECT_REF_MATCH_R1\t"+getFlagSum(hasIns = Some(1))+"\n");
      
      //referenceMismatch(read num)(strand)(read pos)(refBase)(readBase)
      if(calcReferenceMatch){
        val writer = createOutputFile(outfile , "referenceMismatchRaw.byReadStrand.txt","",docWriter,
             ("POS","int","The read position."),
             ("REFBASE","int","The reference base-pair (as it should have appeared on the read)."),
             ("READBASE","int","The read base-pair."),
             ("CT_R1_FWD","int","As CT_R1, but only counting FWD strand aligned reads."),
             ("CT_R1_REV","int","As CT_R1, but only counting REV strand aligned reads."),
             ("CT_R2_FWD","int","As CT_R2, but only counting FWD strand aligned reads."),
             ("CT_R2_REV","int","As CT_R2, but only counting REV strand aligned reads."),
             ("CT_R1","int","The number of base-pair swaps of the given type at the given position for read 1."),
             ("CT_R2","int","The number of base-pair swaps of the given type at the given position for read 2.")
        );
        writer.write("POS\tREFBASE\tREADBASE\tCT_R1_FWD\tCT_R1_REV\tCT_R2_FWD\tCT_R2_REV\tCT_R1\tCT_R2\n");
        for(i <- 0 until readLen+1){
          for((af,afi,ar,ari,bf,bfi,br,bri) <- pairs){
            writer.write(i+"\t"+af+"\t"+bf+"\t"+
                         referenceMismatch(0)(0)(i)(afi)(bfi)+"\t"+
                         referenceMismatch(0)(1)(i)(ari)(bri)+"\t"+
                         referenceMismatch(1)(0)(i)(afi)(bfi)+"\t"+
                         referenceMismatch(1)(1)(i)(ari)(bri)+"\t"+
                        (referenceMismatch(0)(0)(i)(afi)(bfi) + referenceMismatch(0)(1)(i)(ari)(bri)) + "\t"+
                        (referenceMismatch(1)(0)(i)(afi)(bfi) + referenceMismatch(1)(1)(i)(ari)(bri)) +
                         "\n");
          }
        }
        close(writer);
      }
      
      /*
      if(calcReferenceMatch){
        val writer = createOutputFile(outfile, "referenceMismatchRaw.byRefStrand.txt");
        writer.write("POS\tREFBASE\tREADBASE\tCT_R1_FWD\tCT_R1_REV\tCT_R2_FWD\tCT_R2_REV\tCT_R1\tCT_R2\n");
        for(i <- 0 until readLen+1){
          for((af,afi,ar,ari,bf,bfi,br,bri) <- pairs){
            writer.write(i+"\t"+af+"\t"+bf+"\t"+
                         referenceMismatch(0)(0)(i)(afi)(bfi)+"\t"+
                         referenceMismatch(0)(1)(i)(afi)(bfi)+"\t"+
                         referenceMismatch(1)(0)(i)(afi)(bfi)+"\t"+
                         referenceMismatch(1)(1)(i)(afi)(bfi)+"\t"+
                        (referenceMismatch(0)(0)(i)(afi)(bfi) + referenceMismatch(0)(1)(i)(afi)(bfi)) + "\t"+
                        (referenceMismatch(1)(0)(i)(afi)(bfi) + referenceMismatch(1)(1)(i)(afi)(bfi)) +
                         "\n");
          }
        }
        close(writer);
      }*/
      
      if(calcReferenceMatch){
        val writer = createOutputFile(outfile, "referenceMismatchCounts.txt","",docWriter,
             ("POS","int","The read position."),
             ("CT_R1","int","The number of base-pair swaps at the given position for read 1."),
             ("CT_R2","int","The number of base-pair swaps at the given position for read 2."),
             ("MAPPEDCT_R1","int","The number of mapped bases at the given positon for read 1."),
             ("MAPPEDCT_R2","int","The number of mapped bases at the given positon for read 2.")
        );
        writer.write("POS\tCT_R1\tCT_R2\tMAPPEDCT_R1\tMAPPEDCT_R2\n");
        for(i <- 0 until readLen+1){
          val mm1 = pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            referenceMismatch(0)(0)(i)(afi)(bfi)+referenceMismatch(0)(1)(i)(ari)(bri)
          }}.sum
          val mm2 = pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            referenceMismatch(1)(0)(i)(afi)(bfi)+referenceMismatch(1)(1)(i)(ari)(bri)
          }}.sum
          writer.write(i+"\t"+mm1+"\t"+mm2+"\t"+mappedByPositionCoverage(0)(i)+"\t"+mappedByPositionCoverage(1)(i)+"\n");
        }
        close(writer);
      }
      
      /*if(true){
        val writer = openWriterSmart_viaGlobalParam(outfile + ".referenceMismatchCounts.noIndel.txt");
        writer.write("POS\tREFBASE\tREADBASE\tCT_R1_FWD\tCT_R1_REF\tCT_R2_FWD\tCT_R2_REV\tCT_R1\tCT_R2\n");
        for(i <- 0 until readLen+1){
        
          val mm1 = pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            referenceMismatchNoIndel(0)(0)(i)(afi)(bfi)+referenceMismatchNoIndel(0)(1)(i)(ari)(bri)
          }}.sum
          val mm2 = pairs.map{case (af,afi,ar,ari,bf,bfi,br,bri) =>{
            referenceMismatchNoIndel(1)(0)(i)(afi)(bfi)+referenceMismatchNoIndel(1)(1)(i)(ari)(bri)
          }}.sum
          writer.write(i+"\t"+mm1+"\t"+mm2+"\n");
        }
        close(writer);
      }*/
      
    }
    
    
    
    if(calcOverlapMismatch){
      summaryWriter.write("OM_noOverlap_staggered\t"+noOverlap_staggered+"\tNumber of read-pairs with no overlap, mis-staggered such that the fwd read occurs second.\n");
      summaryWriter.write("OM_noOverlap_normal\t"+noOverlap_normal+"\tNumber of read-pairs with no overlap, arranged normally\n");
      summaryWriter.write("OM_overlap\t"+overlap+"\tNumber of read-pairs with good overlap.\n");
      summaryWriter.write("OM_BAD_OVERLAP\t"+BAD_OVERLAP+"\tNumber of read-pairs with a bad overlap\n");
      summaryWriter.write("OM_overlap_Match\t"+PERFECTMATCH+"\tNumber of read-pairs with a good overlap, and that match perfectly.\n");
      summaryWriter.write("OM_overlap_mismatch\t"+MISMATCH+"\tNumber of read-pairs with a good overlap, but that contain one or more base mismatches\n");
    }

    //flagCountArray(hasDel)(hasIns)(mmidx1)(mmidx2)(omidx)
    
    if(calcReferenceMatch || calcOverlapMismatch){
        val writer = createOutputFile(outfile, "mismatchSummary.txt","",docWriter,
             ("hasDel","boolean","Whether the read has deletions"),
             ("hasIns","boolean","Whether the read has insertions"),
             ("refMismatchR1","int","Number of read-1 reference-mismatches"),
             ("refMismatchR2","int","Number of read-2 reference-mismatches"),
             ("overlapMismatch","int","Number of overlap-mismatches")
        );
        writer.write("hasDel\thasIns\trefMismatchR1\trefMismatchR2\toverlapMismatch\n");
        for((delCode,hasDel) <- qcOverlapMatch.INDEL_CODE_MAP){
          for((insCode,hasIns) <- qcOverlapMatch.INDEL_CODE_MAP){
            for( (mmCode1,mmidx1) <- List(("MATCH",0),("MM_FWD",1),("MM_REV",2))){
              for( (mmCode2,mmidx2) <- List(("MATCH",0),("MM_FWD",1),("MM_REV",2))){
                for( (omCode,omidx) <- qcOverlapMatch.OVERLAP_CODE_MAP){
                  writer.write(delCode+"\t"+
                               insCode+"\t"+
                               mmCode1+"\t"+
                               mmCode2+"\t"+
                               omCode+"\t"+
                               flagCountArray(hasDel)(hasIns)(mmidx1)(mmidx2)(omidx)+"\n");
                }
              }
            }
          }
        }
        close(writer);
    }
    
    if(! samwriter.isEmpty) samwriter.get.close();
  }

    def getFlagSum(hasDel: Option[Int] = None,
                   hasIns: Option[Int] = None,
                   mmidx1: Option[Int] = None,
                   mmidx2: Option[Int] = None,
                   omidx: Option[Int]  = None) : Long = {
      //flagCountArray(hasDel)(hasIns)(mmidx1)(mmidx2)(omidx)
      val s1 = omidx match {
        case Some(i) => flagCountArray.map(a => a.map(b => b.map(c => c.map(d => d(i)))));
        case None    => flagCountArray.map(a => a.map(b => b.map(c => c.map(d => d.sum))));
      }
      val s2 = mmidx2 match {
        case Some(i) => s1.map(a => a.map(b => b.map(c => c(i))));
        case None    => s1.map(a => a.map(b => b.map(c => c.sum)));
      }
      val s3 = mmidx1 match {
        case Some(i) => s2.map(a => a.map(b => b(i)));
        case None    => s2.map(a => a.map(b => b.sum));
      }
      val s4 = hasIns match {
        case Some(i) => s3.map(a => a(i));
        case None    => s3.map(a => a.sum);
      }
      val out = hasDel match {
        case Some(i) => s4(i)
        case None    => s4.sum
      }
      return out;
    }
  
         /*
  
  val referenceMismatch : Array[Array[Array[Array[Array[Int]]]]] = Array.ofDim[Int](2,2,readLen+1,qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1,qcOverlapMatch.ALL_BASE_BYTES.max.toInt+1);
  //referenceMismatchCount(read num)(num bases mismatched)
  val referenceMismatchCount : Array[Array[Int]] = Array.ofDim[Int](2,readLen+1);
  
  val PERFECTREFMATCH   : Array[Int] = Array.ofDim[Int](2);
  val IMPERFECTREFMATCH : Array[Int] = Array.ofDim[Int](2);
  
      val ixRefStart = math.max(getAlignmentStart(rf),getAlignmentStart(rr));
      val ixRefEnd   = math.min(getAlignmentEnd(rf),getAlignmentEnd(rr));
      //val ixRefStart = math.max(unstartF,unstartR);
      //val ixRefEnd = math.min(unendF,unendR);
      
      val cb1 = getCigarBlocksFromRead(r1).toVector;
      val cb2 = getCigarBlocksFromRead(r2).toVector;
      val cbix1 = truncateCigarBlocks(cb1.iterator,ixRefStart,ixRefEnd).toVector;
      val cbix2 = truncateCigarBlocks(cb2.iterator,ixRefStart,ixRefEnd).toVector;
      
      if(cbix1.isEmpty || cbix2.isEmpty){
        BAD_OVERLAP += 1;
      } else {
        val rawOffset = cbix1.head.readStart - cbix2.head.readStart;
        if(cbix1.head.refStart != cbix2.head.refStart){
          BAD_OVERLAP_A += 1;
        } else if(getAlignmentStart(r1) == ixRefStart && cbix2.head.refStart != ixRefStart){
          BAD_OVERLAP_B += 1;
        } else if(getAlignmentStart(r2) == ixRefStart && cbix1.head.refStart != ixRefStart){
          BAD_OVERLAP_C += 1;
        } else {
          val (rA,rB,offset) = if(rawOffset >= 0) (r1,r2,rawOffset) else (r2,r1,-rawOffset);
          val seqA = rA.getReadString().toVector.drop(offset);
          val seqB = rB.getReadString().toVector.dropRight();
          
        }
      }*/ 
    /*val cb1 = getCigarBlocksFromRead(r1);
    val cb2 = getCigarBlocksFromRead(r2);
    
    //get overlapping blocks:
    val cb1f = cb1.filter(block1 => cb2.exists(block2 => {
      block1.refStart < block2.refEnd && block1.refEnd > block2.refStart
    }));
    val cb2f = cb2.filter(block2 => cb1.exists(block1 => {
      block1.refStart < block2.refEnd && block1.refEnd > block2.refStart
    }));*/

  
  /*
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    val (start1, end1) = (r1.getAlignmentStart() - 1, r1.getAlignmentEnd);
    val (start2, end2) = (r2.getAlignmentStart() - 1, r2.getAlignmentEnd);
    
    val pp1 = getAlignedBasePositionsFromRead(r1);
    val pp2 = getAlignedBasePositionsFromRead(r2);
    
    val ppo1 = pp1.filter((i) => i._1 >= start2 && i._1 < end2);
    val ppo2 = pp2.filter((i) => i._1 >= start1 && i._1 < end1);
    
    val seq1 = r1.getReadBases();
    val seq2 = r2.getReadBases();
    
    if(ppo1.length == ppo2.length){
      val ppp = ppo1.zip(ppo2);
      if( ppp.forall(pair => pair._1._1 == pair._2._1)){
        ppoMatch += 1;
        for(pair <- ppp){
          basePairCount(pair._1._2)(pair._2._2)( qcOverlapMatch.BYTE_CODE_MAP( seq1(pair._1._2).toInt) )( qcOverlapMatch.BYTE_CODE_MAP(seq2(pair._2._2).toInt) ) += 1;
        }
      } else {
        ppoReferencePositionMismatch += 1;
        //do stuff!
      }
    } else {
      ppoLengthMismatch += 1;
      //do stuff!
    }
    
    /*
    val cb1 = getCigarBlocksFromRead(r1);
    val cb2 = getCigarBlocksFromRead(r2);
    
    val cbo1 = cb1.filter((b : CigarBlock) => b.refStart < end2 & b.refEnd > start2);
    val cbo2 = cb2.filter((b : CigarBlock) => b.refStart < end1 & b.refEnd > start1);
    
    val M1 = cbo1.filter(_.op.equals(CigarOperator.MATCH_OR_MISMATCH));
    val M2 = cbo2.filter(_.op.equals(CigarOperator.MATCH_OR_MISMATCH));
    
    */
    
    //REDO:
  }*/

  def getUtilityName : String = "StrandCheck";
}












