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

import scala.collection.immutable.TreeSet; 
import internalUtils.optionHolder._;
import scala.collection.GenMap;

object qcGetGeneCounts {

  /*
  def readGtf(stranded : Boolean, gtffile : String , geneCounts : scala.collection.mutable.HashMap[String,Int]) : GenomicArrayOfSets[String] = {
    report("reading Gtf: " + gtffile + "\n","note");
    val geneArray : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](stranded);
    val gtfReader = GtfReader.getGtfReader(gtffile, stranded, true, "\\s+");
    var lineCt = 0;
    
    for(gtfLine <- gtfReader){
      lineCt = lineCt + 1;
      readGtfLine(geneArray, gtfLine, geneCounts);
    }
    report("done.\nRead " + lineCt + " gtf lines. Found " + geneCounts.size + " genes.\n","note");
    reportln("Memory usage: " + MemoryUtil.memInfo,"note");    
    return geneArray;
  }
  
  def readGtfLine(geneArray : GenomicArrayOfSets[String], gtfLine : GtfLine, geneCounts : scala.collection.mutable.HashMap[String,Int]) {
    if(gtfLine.featureType == GtfCodes.EXON_TYPE_CODE){
      val geneID = gtfLine.getAttributeOrDie(GtfCodes.GENE_ID_ATTRIBUTE_KEY)
      //reportln("found exon for: " + geneID, "note");
      geneCounts(geneID) = 0;
      geneArray.addSpan(new GenomicInterval(gtfLine.chromName, gtfLine.strand, gtfLine.start - 1, gtfLine.end), geneID);
    }
  }*/
  
  
  
  def calculateGeneBodyCoverage(geneAssignment : String, r1 : SAMRecord, r2 : SAMRecord, geneBodyCoverageMap : GenMap[String, Array[Int]], intervalMap : Map[String, Vector[TreeSet[GenomicInterval]]]){
    if(geneAssignment != "_no_feature" && geneAssignment != "_ambiguous"){
      intervalMap.get(geneAssignment) match {
        case Some(intervalVector) => {
          //val intervalVector = intervalMap(geneAssignment);
          val intervalSet = helper_findReadIntervalCoverage(r1,intervalVector) ++ helper_findReadIntervalCoverage(r2,intervalVector);
          val intervalCountArray = geneBodyCoverageMap(geneAssignment);
      
          for(i <- intervalSet){
            intervalCountArray(i) += 1;
          }
        }
        case None => {
          //Do nothing!
        }
      }
    } // Else do nothing!
  } 
  
  private def helper_findReadIntervalCoverage(r : SAMRecord, intervalVector : Vector[TreeSet[GenomicInterval]]) : Set[Int] = {
    return r.getAlignmentBlocks.iterator.foldLeft( Set[Int]() )( (sofar, currBlock) => {
      val (start,end) = (currBlock.getReferenceStart() - 1, currBlock.getReferenceStart - 1 + currBlock.getLength());
      sofar ++ helper_findSpanIntervalCoverage(start,end,intervalVector);
    });
  }
  
  private def helper_findSpanIntervalCoverage(start : Int, end : Int, intervalVector : Vector[TreeSet[GenomicInterval]]) : Set[Int] = {
    intervalVector.zip(0 until intervalVector.size).filter( (z) => {
      val (ts,i) = z;
      ts.exists( (iv) => iv.overlaps(start,end) );
    }).map((z) => z._2).toSet ;
  }
  
  
  /*
   * Note that intervalBreaks must begin with 0.0 and end with 1.0!
   */

   /*
   * FIX ME FOR STRANDEDNESS!!!!
   */
  def makeGeneIntervalMap(intervalBreaks : Seq[Double], geneArray : GenomicArrayOfSets[String]) : Map[String, Vector[TreeSet[GenomicInterval]]] = {
    //val initialMap = geneSet.foldLeft( new scala.collection.immutable.HashMap[String,  TreeSet[GenomicInterval] ]() )((soFar,curr) =>{
    //  soFar + ((curr, new TreeSet[GenomicInterval]() ));
    //})
    val geneMap = helper_calculateGeneAssignmentMap_strict(geneArray);
    reportln("making makeGeneIntervalMap for geneBody calculations. Found: " + geneMap.size + " acceptable genes for gene-body analysis.","debug");
    
    val geneLengths = geneMap.map((cg) => {
      val currGene : String = cg._1;
      val currIvSet : TreeSet[GenomicInterval] = cg._2;
      (currGene, currIvSet.foldLeft(0)((sum,curr) => sum + (curr.end - curr.start) ));
    });
    val geneStrands = geneMap.map((cg) => {
      val (currGene, currIvSet) = cg;
      ((currGene,currIvSet.head.strand));
    })
    
    val out = geneMap.foldLeft(Map[String, Vector[TreeSet[GenomicInterval]]]())((sofar,cg) => {
      val currGene : String = cg._1;
      val currGeneLen : Int = geneLengths(currGene);
      if(geneStrands(currGene) == '+') sofar + ((currGene, helper_calculateBreakMap(cg, currGeneLen, intervalBreaks)));
      else sofar + ((currGene, helper_calculateBreakMap(cg, currGeneLen, intervalBreaks).reverse));
    });
    return out;
  }
  
  def makeGeneBodyCoverageCountArrays(intervalCount : Int, geneSet : Iterable[String] ) : scala.collection.mutable.Map[String,Array[Int]] = {
    val out : scala.collection.mutable.AnyRefMap[String,Array[Int]] = scala.collection.mutable.AnyRefMap[String,Array[Int]]();
    
    geneSet.foreach( geneID => {
      out += (geneID, new Array[Int](intervalCount));
    })
    out.repack();
    
    return out;
    
    //return geneSet.map((geneID : String) => {
    //  (geneID, new Array[Int](intervalCount));
    //}).toMap;
  }
  
  private def helper_calculateGeneAssignmentMap_strict(geneArray : GenomicArrayOfSets[String]) : Map[String, TreeSet[GenomicInterval]] = {

    val badGeneSet = geneArray.getSteps.foldLeft(Set[String]())( (soFar, curr) =>{
      val (iv, geneSet) = curr;
      if(geneSet.size > 1){
        soFar ++ geneSet.toSet;
      } else {
        soFar;
      }
    });
    //val allGeneSet =  geneArray.getSteps.foldLeft(Set[String]())( (soFar, curr) =>{
    //  soFar ++ curr._2.toSet;
    //});
    
    reportln("helper_calculateGeneAssignmentMap_strict. Found: " + geneArray.getValueSet.size + " genes in the supplied annotation.","debug");
    reportln("helper_calculateGeneAssignmentMap_strict. Found: " + badGeneSet.size + " genes with ambiguous segments.","debug");
    
    val buf = geneArray.getSteps.foldLeft(  Map[String, TreeSet[GenomicInterval] ]() )( (soFar, curr) => {
      val currIv = curr._1;
      val currGeneSet = curr._2;
      if(currGeneSet.size == 1){
        val currGene = currGeneSet.head;
        if(! badGeneSet.contains(currGene)){
          soFar.get(currGene) match {
            case Some(ts : TreeSet[GenomicInterval]) => soFar + ((currGene, ts + currIv ));
            case None => soFar + ((currGene, TreeSet[GenomicInterval](currIv) ));
          }
        } else {
          soFar;
        }
      } else {
        soFar;
      }
    })
    reportln("helper_calculateGeneAssignmentMap_strict. Found: " + buf.size + " genes after first-pass filtering","debug");

    return buf.filter( (curr) => {
      val (geneID, geneTree) = curr;
      if(geneTree.size == 0) false;
      else {
        geneTree.forall( (c) => geneTree.head.chromName == c.chromName && geneTree.head.strand == c.strand );
      }
    });
  }
  
  /*
   * FIX ME FOR STRANDEDNESS!!!!
   * Done?
   */
  private def helper_calculateBreakMap( cg : (String, TreeSet[GenomicInterval]), currGeneLen : Int,  intervalBreaks : Seq[Double]) : Vector[TreeSet[GenomicInterval]] = {
    val currGene : String = cg._1;
    val currIvSet : TreeSet[GenomicInterval] = cg._2;
    val breakPoints : Seq[Int] = intervalBreaks.map((b) => ( b * currGeneLen.toDouble ).floor.toInt);
    val chromName : String = currIvSet.head.chromName;
    val strand : Char = currIvSet.head.strand;
    
    val breakSpans : Seq[(Int,Int)] = breakPoints.zip(breakPoints.tail);
    val geneCoordMap : Seq[(GenomicInterval,(Int,Int))] = currIvSet.scanLeft( (GenomicInterval(chromName,strand,0,0), (0,0)) )( (sofar, iv) =>{
      val currSpan : (Int,Int) = ((iv.start,iv.end));
      val currSpanGeneCoord : (Int,Int) = (sofar._2._2, sofar._2._2 + (currSpan._2 - currSpan._1)  )
      ((iv,currSpanGeneCoord));
    }).toSeq;
    
    return breakSpans.map((currSpan : (Int,Int)) => {
      helper_calculateBreaks(currSpan, geneCoordMap, chromName, strand);
    }).toVector;
  }
  private def helper_calculateBreaks( currSpan : (Int,Int), geneCoordMap : Seq[(GenomicInterval,(Int,Int))], chromName : String, strand : Char ) : TreeSet[GenomicInterval] = {
    val (currSpanStart, currSpanEnd) = currSpan;
    val intersectingSpans : Array[(GenomicInterval, (Int,Int))] = geneCoordMap.filter( pp => {
        val (iv, (gs,ge)) = pp;
        ( gs < currSpanEnd ) && (currSpanStart < ge);
    }).toArray
    if(intersectingSpans.length > 0){
      val headStart = intersectingSpans.head._2._1;
      val lastEnd = intersectingSpans.last._2._2;
      val addToStart = currSpanStart - headStart;
      val subtractFromEnd = lastEnd - currSpanEnd;
      intersectingSpans(0) = (( new GenomicInterval(chromName, strand, intersectingSpans.head._1.start + addToStart, intersectingSpans.head._1.end ), (intersectingSpans.head._2._1 + addToStart, intersectingSpans.head._2._2) ));
      intersectingSpans(intersectingSpans.length - 1) = (( new GenomicInterval(chromName, strand, intersectingSpans.last._1.start, intersectingSpans.last._1.end - subtractFromEnd ), (intersectingSpans.last._2._1, intersectingSpans.last._2._2 - subtractFromEnd) ));      
      
      return intersectingSpans.foldLeft(TreeSet[GenomicInterval]())( (sofar, c) =>{
        sofar + c._1;
      })
    } else {
      return TreeSet[GenomicInterval]();
    }
  }
  
  def writeGeneBodyCoverage_genewise(outfile : String, geneBody_intervalBreaks : Seq[Double], geneCounts : scala.collection.mutable.Map[String,Int], geneBody_CoverageCountArrays : scala.collection.mutable.Map[String,Array[Int]]){
    val writer = openWriterSmart_viaGlobalParam(outfile+".geneBodyCoverage.genewise.txt");
    
    writer.write("GENE_ID	" + geneBody_intervalBreaks.tail.mkString("	")+"\n");
    for(gene <- geneBody_CoverageCountArrays.keys.toVector.sorted){
      val countArray = geneBody_CoverageCountArrays(gene)
      writer.write(gene +"	"+ countArray.mkString("	") + "\n");
    }
    close(writer);
  }
  def debug_writeGeneBodySpans(outfile : String, geneBody_intervalBreaks : Seq[Double], geneBody_CoverageIntervalMap : Map[String, Vector[TreeSet[GenomicInterval]]]){
    val writer = openGzipWriter(outfile+".geneBodyCoverage.DEBUG.intervals.txt.gz");
    writer.write("GENE_ID	" + geneBody_intervalBreaks.tail.mkString("	")+"\n");
    for((gene,treeVector) <- geneBody_CoverageIntervalMap){
      writer.write(gene +"	"+ treeVector.map( (pp) => {
        "["+pp.map(iv => iv.start +"," + iv.end).mkString("][")+"]";
      }).mkString("	") + "\n");
    }
    close(writer);
  }
  
  
  
  val default_coverageLevelThresholds = Seq(("1.bottomHalf",0.5),("2.upperMidQuartile",0.75),("3.75to90",0.9),("4.high",1.0));
  
  //UNFINISHED?
  def geneBody_calculateGeneBodyCoverage_summaries(outfile : String, geneBody_intervalBreaks : Seq[Double], coverageLevelThresholds : Seq[(String,Double)], geneBody_CoverageCountArrays : GenMap[String,Array[Int]], geneCounts : scala.collection.mutable.Map[String,Int]){
    val geneBody_IntervalCount = geneBody_intervalBreaks.length - 1;
    val totalGeneBodyCoverage_simple = geneBody_CoverageCountArrays.foldLeft(repToSeq(0,geneBody_IntervalCount))((sofar, currPair) =>{
      (0 until sofar.length).map(i => {
        sofar(i) + currPair._2(i);
      }).toSeq;
    });
    
    val includeGenesSet = geneBody_CoverageCountArrays.keySet;
    val sortedReadCountSeq = geneCounts.toVector.filter( (pair) =>  includeGenesSet.contains(pair._1) && pair._2 > 0).sortBy( (pair) => (pair._2,pair._1) );
    val coverageThresholds = coverageLevelThresholds.map(cl_thresh => (sortedReadCountSeq.size.toDouble * cl_thresh._2).toInt );
    
    reportln("DEBUG NOTE: IncludeGenesSet.size: " + includeGenesSet.size,"debug");
    reportln("DEBUG NOTE: sortedReadCountSeq.size: " + sortedReadCountSeq.size,"debug");
    reportln("DEBUG NOTE: coverageThresholds: "+coverageThresholds.mkString(";")+"","debug");
    val coverageSpans = Seq[(Int,Int)]((0, coverageThresholds.head)) ++ coverageThresholds.zip(coverageThresholds.tail);
    reportln("DEBUG NOTE: coverageSpans: ["+coverageSpans.mkString(";")+"]","debug");
    
    val geneBodyByCoverageLevel =  coverageLevelThresholds.zip(coverageSpans).map(cltpair => {
      val ((coverName, coverLevel), (start,end)) = cltpair;
      
      reportln("DEBUG NOTE:	["+coverName+"]["+coverLevel+"] = ["+start+","+end+"]","debug");
      
      sortedReadCountSeq.slice(start,end).foldLeft(repToSeq(0,geneBody_IntervalCount))((sofar, currPair) => {
        val (geneID, geneCounts) = currPair;
        val cca = geneBody_CoverageCountArrays(geneID);
        sofar.zip(cca).map(p => p._1 + p._2);
      })
    });
    
    val writer = openWriterSmart_viaGlobalParam(outfile+".geneBodyCoverage.by.expression.level.txt");
    writer.write("QUANTILE	" + coverageLevelThresholds.map(_._1).mkString("	") + "\n");
    for(i <- 0 until geneBody_IntervalCount){
      writer.write(geneBody_intervalBreaks.tail(i) + "	" + (0 until coverageLevelThresholds.length).map(j => {
        geneBodyByCoverageLevel(j)(i);
      }).mkString("	") + "\n");
    }
    close(writer);
  }
  
  def pairIsNearFeatures(distance : Int, r1 : SAMRecord, r2 : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean, featureArray : GenomicArrayOfSets[String]) : Boolean = {
    readIsNearFeatures(distance,r1,stranded,fr_secondStrand, featureArray) || readIsNearFeatures(distance,r1,stranded,fr_secondStrand, featureArray);
  }
  def readIsNearFeatures(distance : Int, r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean, featureArray : GenomicArrayOfSets[String]) : Boolean = {
    val readIntervals : Iterator[GenomicInterval] = getExpandedGenomicIntervalsFromRead(distance, r , stranded , fr_secondStrand);
    return readIntervals.exists( (iv) => {
      featureArray.findIntersectingSteps(iv).exists(! _._2.isEmpty);
    });
  }
  def getPairFeaturesNear(distance : Int, r1 : SAMRecord, r2 : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean, featureArray : GenomicArrayOfSets[String]): Set[String] = {
    return getPairFeaturesNear(distance, r1,stranded,fr_secondStrand, featureArray) ++ getReadFeatures(r2,stranded,fr_secondStrand, featureArray)
  }
  def getPairFeaturesNear(distance : Int, r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean, featureArray : GenomicArrayOfSets[String]) : Set[String] = {
    val readIntervals : Iterator[GenomicInterval] = getExpandedGenomicIntervalsFromRead(distance, r , stranded , fr_secondStrand);
    return readIntervals.foldLeft(Set[String]())((setSoFar, iv) => {
      featureArray.findIntersectingSteps(iv).foldLeft(setSoFar)((ssf,featureSet) => {
        ssf ++ featureSet._2;
      })
    })
  }
  
  def getPairFeatures( r1 : SAMRecord, r2 : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean, featureArray : GenomicArrayOfSets[String]): Set[String] = {
    return getReadFeatures(r1,stranded,fr_secondStrand, featureArray) ++ getReadFeatures(r2,stranded,fr_secondStrand, featureArray)
  }
  def getReadFeatures(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean, featureArray : GenomicArrayOfSets[String]) : Set[String] = {
    val readIntervals : Iterator[GenomicInterval] = getGenomicIntervalsFromRead(r , stranded , fr_secondStrand);
    return readIntervals.foldLeft(Set[String]())((setSoFar, iv) => {
      featureArray.findIntersectingSteps(iv).foldLeft(setSoFar)((ssf,featureSet) => {
        ssf ++ featureSet._2;
      })
    })
  }
  
}

/************************************************************************************************************************************************************
 * Class:
 */

class qcGetGeneCounts( stranded : Boolean, 
                       fr_secondStrand : Boolean, 
                       anno_holder : qcGtfAnnotationBuilder, 
                       coda : Array[Int], 
                       coda_options : Array[Boolean] , 
                       geneBodyIntervalCount : Int, 
                       calcRPKM : Boolean, writeGenewiseGeneBody : Boolean, writeDESeq : Boolean, writeGeneCounts : Boolean) extends QCUtility[String] {
  
  reportln("Init GeneCalcs","progress");
  
  val geneArray : GenomicArrayOfSets[String] = anno_holder.geneArray;
  val strandedGeneArray : GenomicArrayOfSets[String] = anno_holder.strandedGeneArray;
  
  val geneBody_IntervalCount : Int = geneBodyIntervalCount;
  val geneBody_intervalBreaks = (0 to geneBody_IntervalCount).map(_.toDouble / geneBody_IntervalCount.toDouble).toSeq
  val geneBody_CoverageIntervalMap = qcGetGeneCounts.makeGeneIntervalMap(geneBody_intervalBreaks, strandedGeneArray);
  //reportln("making geneBody_CoverageIntervalMap for geneBody calculations. Found: " + geneBody_CoverageIntervalMap.size + " acceptable genes.","debug");
  val geneBody_CoverageCountArrays : scala.collection.mutable.Map[String,Array[Int]] = qcGetGeneCounts.makeGeneBodyCoverageCountArrays(geneBody_IntervalCount, geneBody_CoverageIntervalMap.keys);
  
  //val mapLocation_CDS : GenomicArrayOfSets[String]
  val geneArea_cdsArray = anno_holder.qcGetGeneCounts_cdsArray;
  val geneArea_intronsArray = anno_holder.qcGetGeneCounts_intronArray;

  //INITIALIZE COUNTERS:
  val geneCounts : scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](geneArray.getValueSet);
  
  //val utilCounts : scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](Set("_no_feature","_ambiguous"));
  val geneArea_cdsCounts : scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](geneArray.getValueSet);
  val geneArea_intronsCounts : scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](geneArea_intronsArray.getValueSet);
  
  val geneCounts_ambig : scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](geneArray.getValueSet);
  val geneCounts_utr : scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](geneArray.getValueSet);
  
  var readNoFeature : Int = 0;
  var readAmbiguous : Int = 0;
  
  var readExonCount : Int = 0;
  var readUtrCount : Int = 0;
  var readCdsCount : Int = 0;
  var readIntronCount : Int = 0;
  var readOneKb : Int = 0;
  var readTenKb : Int = 0;
  var readMiddleOfNowhere : Int = 0;
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : String = {
    val readGenes = qcGetGeneCounts.getPairFeatures(r1,r2,stranded,fr_secondStrand,geneArray);
    
    if(readGenes.size == 0){
      //utilCounts("_no_feature") += 1;
      readNoFeature += 1;
      
      //val readIntrons = qcGetGeneCounts.getPairFeatures(r1,r2,stranded,fr_secondStrand,geneArea_intronsArray);
      //for(intron <- readIntrons){
      //  geneArea_intronsCounts(intron) += 1;
      //}
      //if(readIntrons.size > 0) readIntronCount += 1;
      if(qcGetGeneCounts.pairIsNearFeatures(0,r1,r2,stranded,fr_secondStrand,geneArea_intronsArray)) readIntronCount += 1;
      else if(qcGetGeneCounts.pairIsNearFeatures(1000,r1,r2,stranded,fr_secondStrand,geneArray)) readOneKb += 1;
      else if(qcGetGeneCounts.pairIsNearFeatures(10000,r1,r2,stranded,fr_secondStrand,geneArray)) readTenKb += 1;
      else readMiddleOfNowhere += 1;
      
      return "_no_feature";
    } else if(readGenes.size > 1){
      //utilCounts("_ambiguous") += 1;
      readAmbiguous += 1;
      for(g <- readGenes){
        geneCounts_ambig(g) += 1;
      }
      return "_ambiguous";
    } else {
      val geneAssignment = readGenes.head
      geneCounts(geneAssignment) += 1;
      readExonCount += 1;
      val cdsGenes = qcGetGeneCounts.getPairFeatures(r1,r2,stranded,fr_secondStrand,geneArea_cdsArray);
      if(cdsGenes.size == 1){
        geneArea_cdsCounts(geneAssignment) += 1;
        readCdsCount += 1;
      } else {
        geneCounts_utr(geneAssignment) += 1;
        readUtrCount += 1;
      }
      
      //calculateGeneBodyCoverage(geneAssignment : String, r1 : SAMRecord, r2 : SAMRecord, geneBodyCoverageMap : Map[String, Array[Int]], intervalMap : Map[String, Vector[TreeSet[GenomicInterval]]])
      qcGetGeneCounts.calculateGeneBodyCoverage(geneAssignment, r1,r2, geneBody_CoverageCountArrays, geneBody_CoverageIntervalMap);
      
      return geneAssignment;
    }
  }
  
  //def getReadGenes(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean) : Set[String] = {
  //  return qcGetGeneCounts.getReadFeatures(r,stranded,fr_secondStrand,geneArray);
  //}
  
  //var readUtrCount : Int = 0;
  //var readCdsCount : Int = 0;
 // var readIntronCount : Int = 0;
  //var readOneKb : Int = 0;
  //var readTenKb : Int = 0;
 // var readMiddleOfNowhere : Int = 0;
  //def geneArea_CDS_getReadGenes(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean) : 
  
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    
    //Start with the summary:
    
    //val writerSummary = openWriter(outfile + ".geneMappingSummary.txt");

    summaryWriter.write(internalUtils.commonSeqUtils.causeOfDropArrayToStringTabbed(coda, coda_options));
    //summaryWriter.write("ReadPairs_AmbigGene	"+utilCounts("_ambiguous")+"\n");
    summaryWriter.write("ReadPairs_AmbigGene	"+readAmbiguous+"\n");
    summaryWriter.write("ReadPairs_UniqueGene	"+readExonCount+"\n");
    summaryWriter.write("ReadPairs_UniqueGene_CDS	"+readCdsCount+"\n");
    summaryWriter.write("ReadPairs_UniqueGene_UTR	"+readUtrCount+"\n");
    //summaryWriter.write("ReadPairs_NoGene	"+utilCounts("_no_feature")+"\n");
    summaryWriter.write("ReadPairs_NoGene	"+readNoFeature+"\n");
    summaryWriter.write("ReadPairs_NoGene_Intron	"+readIntronCount+"\n");
    summaryWriter.write("ReadPairs_NoGene_OneKbFromGene	"+readOneKb+"\n");
    summaryWriter.write("ReadPairs_NoGene_TenKbFromGene	"+readTenKb+"\n");
    summaryWriter.write("ReadPairs_NoGene_MiddleOfNowhere	"+readMiddleOfNowhere+"\n");
    
    val nonzero_genes = geneCounts.count( pair => { pair._2 > 0 });
    val zero_genes = geneCounts.size - nonzero_genes;
    summaryWriter.write("Genes_Total	" + geneCounts.size +"\n");
    summaryWriter.write("Genes_WithZeroCounts	" + zero_genes +"\n");
    summaryWriter.write("Genes_WithNonzeroCounts	" + nonzero_genes +"\n");

    
    //close(writerSummary);
    
    //SUMMARY DONE.
    if(writeDESeq){
      val writer2 = openWriterSmart_viaGlobalParam(outfile + ".geneCounts.formatted.for.DESeq.txt");
      for(key <- geneCounts.keys.toVector.sorted){
        writer2.write(key + "	"+geneCounts(key) +"\n");
      }
      //writer2.write("no_feature	"+utilCounts("_no_feature")+"\n")
      //writer2.write("ambiguous	"+utilCounts("_ambiguous")+"\n")
      writer2.write("no_feature	"+readNoFeature+"\n")
      writer2.write("ambiguous	"+readAmbiguous+"\n")
      writer2.write("too_low_aQual	0\n")
      writer2.write("not_aligned	0\n")
      writer2.write("alignment_not_unique	0\n")
      close(writer2);
    }
    if(writeGeneCounts){
      val writer = openWriterSmart_viaGlobalParam(outfile + ".geneCounts.txt");
      //writer.write("");
      writer.write("GENEID	COUNT	COUNT_CDS	COUNT_UTR	COUNT_AMBIG_GENE\n");
      for(key <- geneCounts.keys.toVector.sorted){
        writer.write(key + "	"+geneCounts(key) +"	"+geneCounts_utr(key)+"	"+geneArea_cdsCounts(key)+"	"+geneCounts_ambig(key)+"\n");
      }
      //writer.write("_total_no_feature	"+ utilCounts("_no_feature") + "	0	0	0\n");
      //writer.write("_total_ambiguous	"+ utilCounts("_ambiguous") + "	0	0	0\n");
      writer.write("_total_no_feature	"+ readNoFeature + "	0	0	0\n");
      writer.write("_total_ambiguous	"+ readAmbiguous + "	0	0	0\n");
      close(writer);
    }
    
    //delete this later:
    debug_writeGeneBodySpans(outfile);
    if(writeGenewiseGeneBody){
      writeGeneBodyCoverage_genewise(outfile);
    }
    writeGeneBodyCoverage_summaryByExpressionLevel(outfile);
    
    if(calcRPKM){
      val writerRPKM = openWriterSmart_viaGlobalParam(outfile + ".FPKM.txt");
      
      val M = readExonCount / 1000000;
      
      writerRPKM.write("GENEID	FPKM\n");
      for(key <- geneCounts.keys.toSeq.sorted){
        val K = anno_holder.geneLengthMap(key) / 1000;
        val F = geneCounts(key);
        val FPKM = (F / K) / M;
        writerRPKM.write(key + "	"+FPKM+"\n");
      }
      
      close(writerRPKM);
    }
  }
  
  def writeGeneBodyCoverage_genewise(outfile :String){
//  def writeGeneBodyCoverage_genewise(outfile : String, geneBody_intervalBreaks : Seq[Double], geneCounts : scala.collection.mutable.Map[String,Int], geneBody_CoverageCountArrays : Map[String,Array[Int]]){
    qcGetGeneCounts.writeGeneBodyCoverage_genewise(outfile, geneBody_intervalBreaks, geneCounts, geneBody_CoverageCountArrays);
  }
  def debug_writeGeneBodySpans(outfile : String){
    qcGetGeneCounts.debug_writeGeneBodySpans(outfile, geneBody_intervalBreaks, geneBody_CoverageIntervalMap);
  }
  
  def writeGeneBodyCoverage_summaryByExpressionLevel(outfile : String){
    //geneBody_calculateGeneBodyCoverage_summaries(outfile : String, geneBody_intervalBreaks : Seq[Double], coverageLevelThresholds : Seq[(String,Double)], geneBody_CoverageCountArrays : Map[String,Array[Int]], geneCounts : scala.collection.mutable.Map[String,Int])
    qcGetGeneCounts.geneBody_calculateGeneBodyCoverage_summaries(outfile, geneBody_intervalBreaks, qcGetGeneCounts.default_coverageLevelThresholds, geneBody_CoverageCountArrays, geneCounts);
  }
}
















