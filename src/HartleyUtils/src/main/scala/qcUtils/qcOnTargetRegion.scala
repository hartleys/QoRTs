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
import internalUtils.commonSeqUtils._;

object qcOnTargetRegion {
  
  def isOnTarget(r1 : SAMRecord, r2 : SAMRecord, anno_holder : qcGtfAnnotationBuilder) : Boolean = {
    val iv1 = getGenomicIntervalsFromRead(r1, false,false).toVector;
    val iv2 = getGenomicIntervalsFromRead(r2, false,false).toVector;
    val ivB = mergeGenomicIntervalSets(iv1 ++ iv2);
    
    ivB.exists(iv => {
      anno_holder.targetArray.get.findIntersectingSteps(iv).exists{case (iv,tset) => tset.nonEmpty}
    })
  }
  

}


class qcOnTargetRegion(anno_holder : qcGtfAnnotationBuilder, isSingleEnd : Boolean, copyBed : Boolean = true, extended : Boolean = false) extends QCUtility[Boolean] {
  reportln("> Init OnTargetRegion Utility","debug");
  
  val targetArray : GenomicArrayOfSets[Int] = anno_holder.targetArray.get;
  val targetInfo = anno_holder.targetInfo.get;

  val targetCt = Array.ofDim[Long](targetInfo.keys.max+1);
  val targetCovR = Array.ofDim[Long](targetInfo.keys.max+1);
  val targetCovRP = Array.ofDim[Long](targetInfo.keys.max+1);

  val targetIvs = Array.ofDim[GenomicInterval](targetInfo.keys.max+1);
  for(i <- 0 until targetInfo.keys.max+1){
    val cells = targetInfo(i).split("\t");
    targetIvs(i) = new GenomicInterval(cells(0),'.',string2int(cells(1)),string2int(cells(2)));
  }
  
  //geneArray.addSpan(gtfLine.getGenomicInterval.usingStrandedness(stranded), element);
  var onTargetCt : Long = 0;
  var offTargetCt : Long = 0;
  
  var onTargetReadBases : Long = 0;
  var onTargetReadPairBases : Long = 0;
  

  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : Boolean = {
    
    val iv1 = getGenomicIntervalsFromRead(r1, false,false).toVector;
    val iv2 = getGenomicIntervalsFromRead(r2, false,false).toVector;
    val ivB = mergeGenomicIntervalSets(iv1 ++ iv2);
    
    val targets = getTargetSet(ivB);
    
    if(isSingleEnd){
      targets.foreach(t => {
        targetCt(t) += 1;
        val (readCov,readPairCov) = getTargetCoverage(iv1,iv2,ivB,t);
        targetCovR(t) += readPairCov;
        targetCovRP(t) += readPairCov;
        onTargetReadBases += readPairCov;
        onTargetReadPairBases += readPairCov;
      });
    } else {
      targets.foreach(t => {
        targetCt(t) += 1;
        val (readCov,readPairCov) = getTargetCoverage(iv1,iv2,ivB,t);
        targetCovR(t) += readCov;
        targetCovRP(t) += readPairCov;
        onTargetReadBases += readCov;
        onTargetReadPairBases += readPairCov;
      });
    }
    
    if(targets.isEmpty){
      offTargetCt += 1;
      return false;
    } else {
      onTargetCt += 1;
      return true;
    }
  }
  
  def getTargetCoverage(iv1: Seq[GenomicInterval], iv2 : Seq[GenomicInterval], ivB : Seq[GenomicInterval], t : Int) : (Int,Int) = {
      val tarIV = targetIvs(t);
      val ix = iv1.map(iv => {
        getOverlap(iv,tarIV);
      }).sum + iv2.map(iv => {
        getOverlap(iv,tarIV);
      }).sum
      val ixB = ivB.map(iv => {
        getOverlap(iv,tarIV);
      }).sum
      (ix,ixB);
  }
  
  def getTargetSet(readIntervals : Seq[GenomicInterval]) : Set[Int] = {
    return readIntervals.foldLeft(Set[Int]())((setSoFar, iv) => {
      targetArray.findIntersectingSteps(iv).foldLeft(setSoFar)((ssf,featureSet) => {
        ssf ++ featureSet._2;
      })
    })
  }

    //val targetCovR = Array.ofDim[Int](targetInfo.keys.max);
  //val targetCovRP = Array.ofDim[Int](targetInfo.keys.max);

  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    val readDepth = Array.ofDim[Float](targetCt.length);
    val pairDepth = Array.ofDim[Float](targetCt.length);
    val spans = Array.ofDim[Int](targetCt.length);
    for(i <- 0 until targetCt.length){
      val info = targetInfo(i);
      val cells = info.split("\t");
      spans(i) = string2int(cells(3))
      readDepth(i) = targetCovR(i).toFloat  / spans(i).toFloat;
      pairDepth(i) = targetCovRP(i).toFloat / spans(i).toFloat;
    }
    
    
    if(extended){
      val writer = createOutputFile(outfile , "onTarget.extended.txt","",docWriter,
             ("CHROM","String","The chromosome."),
             ("START","int","The start genomic position (0-based."),
             ("END","int","The end position"),
             ("SPAN","int","The length of the target contig"),
             ("CT","int","The number of reads or read-pairs that cover the target contig"),
             ("READ_COVERAGE","long","The number of read-base-pairs that overlap with the target contig. This counts each read independantly, so with paired-end data, if read-pairs overlap with one another, then the overlapping region is counted twice."),
             ("READPAIR_COVERAGE","long","The number of read-pair-base-pairs that overlap with the target contig. If read-pairs overlap with one another, then the overlapping region is counted ONCE."),
             ("READ_DEPTH","double","The READ_COVERAGE divided by the SPAN."),
             ("READPAIR_DEPTH","double","The READPAIR_COVERAGE divided by the SPAN.")
      );
      writer.write("CHROM\tSTART\tEND\tSPAN\tCT\tREAD_COVERAGE\tREADPAIR_COVERAGE\tREAD_DEPTH\tREADPAIR_DEPTH\n");
      for(i <- 0 until targetCt.length){
        val info = targetInfo(i);
        val cells = info.split("\t");
        writer.write(cells(0) +"\t"+ cells(1)+ "\t"+  cells(2) +"\t"+cells(3) +"\t"+
                     targetCt(i)+"\t"+
                     targetCovR(i)+"\t"+
                     targetCovRP(i)+"\t"+
                     readDepth(i) +"\t"+
                     pairDepth(i) +""+
                     "\n");
      }
      close(writer);
      
    } else {
      val writer = createOutputFile(outfile , "onTarget.txt","",docWriter,
             ("CHROM","String","The chromosome."),
             ("START","int","The start genomic position (0-based."),
             ("END","int","The end position"),
             ("CT","int","The number of reads or read-pairs that cover the target contig"),
             ("READ_COVERAGE","long","The number of read-base-pairs that overlap with the target contig. This counts each read independantly, so with paired-end data, if read-pairs overlap with one another, then the overlapping region is counted twice."),
             ("READPAIR_COVERAGE","long","The number of read-pair-base-pairs that overlap with the target contig. If read-pairs overlap with one another, then the overlapping region is counted ONCE.")
      );
      writer.write("CHROM\tSTART\tEND\tCT\tREAD_COVERAGE\tREADPAIR_COVERAGE\n");
      for(i <- 0 until targetCt.length){
        val info = targetInfo(i);
        val cells = info.split("\t");
        writer.write(cells(0) +"\t"+ cells(1)+ "\t"+  cells(2)+"\t"+
                     targetCt(i)+"\t"+
                     targetCovR(i)+"\t"+
                     targetCovRP(i)+"\n");
      }
      close(writer);
    }
     /* if(copyBed){
        val bedWriter = openWriterSmart_viaGlobalParam(outfile + ".target.bed");
        for(i <- 0 until targetCt.length){
          val info = targetInfo(i);
          val cells = info.split("\t");
          bedWriter.write(cells(0) +"\t"+ cells(1)+ "\t"+  cells(2) +"\t"+i+"\n");
        }
        bedWriter.close();
      }
      val writer = openWriterSmart_viaGlobalParam(outfile + ".onTarget.txt");
      writer.write("TARGET\tCT\tREAD_COVERAGE\tREADPAIR_COVERAGE\n");
      for(i <- 0 until targetCt.length){
        val info = targetInfo(i);
        val cells = info.split("\t");
        writer.write(i +"\t"+
                     targetCt(i)+"\t"+
                     targetCovR(i)+"\t"+
                     targetCovRP(i)+"\n");
      }
      close(writer);*/
      
    summaryWriter.write("OnTargetCount\t"+onTargetCt+"\tNumber of reads or read-pairs that intersect with the target region\n");
    summaryWriter.write("OffTargetCount\t"+offTargetCt+"\tNumber of reads or read-pairs that do not intersect with the target region\n");
    summaryWriter.write("OnTargetReadBases\t"+onTargetReadBases+"\tNumber of read-bases that appear over the target regions\n");
    summaryWriter.write("OnTargetReadPairBases\t"+onTargetReadPairBases+"\tNumber of read-pair-bases that appear over the target regions. If read-pairs overlap then the overlapping region is only counted once.\n");
    
    summaryWriter.write("AvgTargetReadDepth\t"+readDepth.sum / readDepth.length+"\tAverage of the read depth across all target contigs\n");
    summaryWriter.write("AvgTargetReadPairDepth\t"+pairDepth.sum / pairDepth.length+"\tAverage of the read-pair depth across all target contigs\n");
    
    summaryWriter.write("TargetReadDepth\t"+    onTargetReadBases.toDouble / spans.sum.toDouble+"\tTotal number of read bases found on-target. If read-pairs overlap then the two reads are counted separately.\n");
    summaryWriter.write("TargetReadPairDepth\t"+onTargetReadPairBases.toDouble / spans.sum.toDouble+"\tTotal number of read-pair bases found on-target. If read-pairs overlap then the overlapping region is only counted once.\n");
    
    summaryWriter.write("OnTargetFraction\t"+onTargetCt.toDouble / (onTargetCt.toDouble + offTargetCt.toDouble)+"\tThe fraction of reads that intersect with the target region\n");
        
    summaryWriter.write("TotalTargetSpan\t"+spans.sum+"\tThe total number of base pairs in the target regions.\n");
  }
  
  
  
  
  def getUtilityName : String = "onTargetRegion";

}