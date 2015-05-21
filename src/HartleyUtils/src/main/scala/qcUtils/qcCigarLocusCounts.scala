package qcUtils

import scala.collection.immutable.TreeMap;
import scala.collection.immutable.TreeSet;

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
import scala.collection.GenMap;


object qcCigarLocusCounts {

  def getCigarLocusSetFromRead(r : SAMRecord, op : CigarOperator, stranded : Boolean, fr_secondStrand : Boolean) : Set[GenomicInterval] = {
    val cigOps : Stream[CigOp] = CigarHolder(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean).cigOps;
    
    cigOps.filter( _.op == op).map(_.ref_iv).toSet;
  }
  
}

class qcCigarLocusCounts(stranded : Boolean, fr_secondStrand : Boolean, 
                         writeAllCountVariants : Boolean, writeHighCountVariants : Boolean,
                         inclusionThreshold : Int)  extends QCUtility[(Int,Int)] {
  
  reportln("> Init qcCigarLocusCounts Utility","debug");
  
  val deletionCountMap : scala.collection.mutable.Map[GenomicInterval,Int] = scala.collection.mutable.AnyRefMap[GenomicInterval,Int]().withDefault(k => 0);
  val insertionCountMap : scala.collection.mutable.Map[GenomicInterval,Int] = scala.collection.mutable.AnyRefMap[GenomicInterval,Int]().withDefault(k => 0);
  
  var insertionCount = 0;
  var deletionCount = 0;
   //(Map[GenomicInterval,Int]()).withDefault(k => 0);
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : (Int,Int) = {
    
    val r1Del = qcCigarLocusCounts.getCigarLocusSetFromRead(r1,CigarOperator.DELETION, stranded, fr_secondStrand);
    val r2Del = qcCigarLocusCounts.getCigarLocusSetFromRead(r1,CigarOperator.DELETION, stranded, fr_secondStrand);
    val r1Ins = qcCigarLocusCounts.getCigarLocusSetFromRead(r2,CigarOperator.INSERTION, stranded, fr_secondStrand);
    val r2Ins = qcCigarLocusCounts.getCigarLocusSetFromRead(r2,CigarOperator.INSERTION, stranded, fr_secondStrand);
    
    val del = r1Del ++ r2Del;
    val ins = r1Ins ++ r2Ins;
    
    del.foreach( (d) => {
      deletionCountMap(d) += 1;
    })
    ins.foreach( (i) => {
      insertionCountMap(i) += 1;
    });
    
    insertionCount += ins.size;
    deletionCount += del.size;
    
    //if(del.size > 0 | ins.size > 0) reportln("Read "+r1.getReadName +" del="+del.size+", ins="+ins.size,"deepDebug");
    
    return(del.size, ins.size);
  }
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    
    summaryWriter.write("deletionLoci  "+deletionCountMap.size+"\n");
    summaryWriter.write("insertionLoci  "+insertionCountMap.size+"\n");
    summaryWriter.write("deletionEventCt "+deletionCount+"\n");
    summaryWriter.write("insertionEventCt  "+insertionCount+"\n");
    
    val delKeyList = deletionCountMap.keys.toVector.sorted;
    val insKeyList = insertionCountMap.keys.toVector.sorted;
    
    if(writeAllCountVariants){
      val writer1 = openWriterSmart_viaGlobalParam(outfile + ".cigarLoci.deletionCounts.all.txt");
      writer1.write("chrom  strand  start end CT\n");
      for(iv <- delKeyList){
        val ct = deletionCountMap(iv);
        writer1.write(iv.chromName+"  "+iv.strand+" "+iv.start+"  "+iv.end+"  "+ct+"\n");
      }
      writer1.close;
      
      val writer2 = openWriterSmart_viaGlobalParam(outfile + ".cigarLoci.insertionCounts.all.txt");
      writer2.write("chrom  strand  start end CT\n");
      for(iv <- insKeyList){
        val ct = insertionCountMap(iv);
        writer2.write(iv.chromName+"  "+iv.strand+" "+iv.start+"  "+iv.end+"  "+ct+"\n");
      }
      writer2.close;
    }
    
    
    
    if(writeHighCountVariants){
      var highCountDel = 0;
      var highCountIns = 0;
      
      val writer3 = openWriterSmart_viaGlobalParam(outfile + ".cigarLoci.deletionCounts.highCoverage.txt");
      writer3.write("chrom  strand  start end CT\n");
      for(iv <- delKeyList){
        val ct = deletionCountMap(iv);
        if(ct > inclusionThreshold){
          writer3.write(iv.chromName+"  "+iv.strand+" "+iv.start+"  "+iv.end+"  "+ct+"\n");
          highCountDel += 1;
        }
      }
      writer3.close;
      
      val writer4 = openWriterSmart_viaGlobalParam(outfile + ".cigarLoci.insertionCounts.highCoverage.txt");
      writer4.write("chrom  strand  start end CT\n");
      for(iv <- insKeyList){
        val ct = insertionCountMap(iv);
        if(ct > inclusionThreshold){
          writer4.write(iv.chromName+"  "+iv.strand+" "+iv.start+"  "+iv.end+"  "+ct+"\n");
          highCountIns += 1;
        }
      }
      writer4.close;
      
      summaryWriter.write("highCoverageDeletionLoci "+highCountDel+"\n");
      summaryWriter.write("highCoverageInsertionLoci  "+highCountIns+"\n");
    }
  }
  def getUtilityName : String = "qcCigarLocusCounts";
  
}

        



















