package qcUtils

import internalUtils.commonSeqUtils._;
import internalUtils.genomicUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.JavaConversions._
import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;

import scala.collection.immutable.TreeSet;
import internalUtils.optionHolder._;

import scala.collection.GenMap;


class qcGtfAnnotationBuilder(gtffile : String, flatgtffile : Option[String], stranded : Boolean, stdCodes : GtfCodes, flatCodes : GtfCodes, targetRegionBed : Option[String] = None){
  
  def makeStdReader : Iterator[StdGtfLine] = GtfReader.getStdGtfReader(gtffile, stranded, true, "\\s+", stdCodes);
  
  lazy val makeFlatReader : (() => Iterator[FlatGtfLine]) = if(! flatgtffile.isEmpty){
    () => GtfReader.getFlatGtfReader(flatgtffile.get, stranded, true, "\\s+", flatCodes);
  } else {
    reportln("Compiling flat feature annotation, internally in memory...","progress");
    val flatlines = fileConversionUtils.prepFlatGtfFile.getFlatGtfLines(gtffile,stranded);
    reportln("Internal flat feature annotation compiled!","progress");
    () => flatlines.iterator;
  }
   
  //def makeFlatReader : Iterator[FlatGtfLine] = if(! flatgtffile.isEmpty){
  //   GtfReader.getFlatGtfReader(gtffile, stranded, true, "\\s+", flatCodes);
  //} else {
     
  //}

  lazy val spliceJunctionTreeMap : GenMap[(String,Char),TreeSet[(Int,Int)]] = qcGtfAnnotationBuilder.qcInnerDistance_readSplicesFromGtfFile(makeFlatReader,stranded, flatCodes);
  lazy val knownSpliceJunctionNameMap : GenMap[GenomicInterval,String] = qcGtfAnnotationBuilder.qcJunctionCounts_makeJunctionMap(makeFlatReader, stranded, flatCodes);
  lazy val geneArray : GenomicArrayOfSets[String] = qcGtfAnnotationBuilder.qcGetGeneCounts_readGtf(stranded,gtffile, stdCodes).finalizeStepVectors;
  lazy val qcGetGeneCounts_cdsArray : GenomicArrayOfSets[String]  =  qcGtfAnnotationBuilder.qcGetGeneCounts_geneArea_CDS_readGtf(stranded, gtffile, stdCodes).finalizeStepVectors;
  lazy val qcGetGeneCounts_intronArray : GenomicArrayOfSets[String]  =  qcGtfAnnotationBuilder.qcGetGeneCounts_geneArea_INTRONS_readFlatGtf(stranded, makeFlatReader, flatCodes).finalizeStepVectors;
  lazy val qcGetGeneCounts_spanArray : GenomicArrayOfSets[String] = qcGtfAnnotationBuilder.qcGetGeneCounts_geneArea_regions(stranded,gtffile, stdCodes).finalizeStepVectors;
  lazy val flatFeatureList : IndexedSeq[String] = qcGtfAnnotationBuilder.getFlatFeatureList(makeFlatReader, stranded, flatCodes);
  
  lazy val geneMapStrict  : Map[String,TreeSet[GenomicInterval]] = qcGtfAnnotationBuilder.helper_calculateGeneAssignmentMap_strict(geneArray, strandedGeneArray);
  
  lazy val geneLengthMap : GenMap[String,Int] = qcGtfAnnotationBuilder.getGeneLengthMap(geneArray);
  lazy val flatExonArray : GenomicArrayOfSets[String] = qcGtfAnnotationBuilder.makeFlatExonMap(stranded,makeFlatReader, flatCodes).finalizeStepVectors;
  lazy val flatGeneSet : Set[String] = qcGtfAnnotationBuilder.makeFlatGeneSet(flatExonArray);
  lazy val geneBiotypeMap : Map[String,String] = qcGtfAnnotationBuilder.getGeneBiotypeMap(gtffile, stdCodes);
  lazy val (targetArray,targetInfo) :  (Option[GenomicArrayOfSets[Int]],Option[scala.collection.mutable.Map[Int,String]]) = qcGtfAnnotationBuilder.generateTargetArray(targetRegionBed);
  
  lazy val strandedGeneArray = if(stranded) {
    geneArray; 
  } else {
    qcGtfAnnotationBuilder.qcGetGeneCounts_readGtf(true,gtffile, stdCodes).finalizeStepVectors;
  }
  
  
  def initializationReport() {
    /*reportln("initializationReport: qcInnerDistance_spliceMapTree.size: " + spliceJunctionTreeMap.size + " chromosome/strands with annotated splice Junction sites.", "debug");
    reportln("initializationReport: qcJunctionCounts_knownSpliceMap.size: " + knownSpliceJunctionNameMap.size + " annotated splice junction sites.", "debug");
    reportln("initializationReport: qcGetGeneCounts_geneArray.getValueSet.size: " + geneArray.getValueSet.size + " annotated genes.", "debug");
    reportln("initializationReport: qcGetGeneCounts_cdsArray.getValueSet.size: " + qcGetGeneCounts_cdsArray.getValueSet.size + " annotated genes with annotated CDS.", "debug");
    reportln("initializationReport: qcGetGeneCounts_intronArray.getValueSet.size: " + qcGetGeneCounts_intronArray.getValueSet.size + " genes with annotated splice junctions.", "debug");
    reportln("initializationReport: flatExonArray.getValueSet.size: " + flatExonArray.getValueSet.size + " flattened exons bins.", "debug");
    reportln("initializationReport: geneLengthMap.size: " + geneLengthMap.size + " genes.", "debug");
    reportln("initializationReport: flatGeneSet.size: " + flatGeneSet.size + " aggregate-genes.", "debug");
    reportln("initializationReport: flatFeatureList.size: " + flatFeatureList.size + " gene features.", "debug");*/
  }
  
  initializationReport();
}

object qcGtfAnnotationBuilder {
  
  def generateTargetArray(targetBedFile : Option[String]) : (Option[GenomicArrayOfSets[Int]],Option[scala.collection.mutable.Map[Int,String]]) = {
    targetBedFile match {
      case Some(bedfile) => {
        val targetArray : GenomicArrayOfSets[Int] = GenomicArrayOfSets[Int](false);
        val targetInfo = scala.collection.mutable.Map[Int,String]();
         
        var targetRegionSpan = 0;
        getLinesSmartUnzip(bedfile).zipWithIndex.foreach{case (line,lnct) => {
          val cells = line.split("\t");
          val iv = new GenomicInterval(cells(0),'.',string2int(cells(1)), string2int(cells(2)));
          targetArray.addSpan(iv,lnct);
          val span = string2int(cells(2)) - string2int(cells(1));
          targetRegionSpan = targetRegionSpan + span;
          targetInfo.put(lnct,cells(0)+"\t"+cells(1)+"\t"+cells(2)+"\t"+span);
        }}
        targetArray.finalizeStepVectors;
        
        ((Some(targetArray),Some(targetInfo)));
      }
      case None => ((None,None));
    }
  }
  
  
  //UNIMPLEMENTED: For potential future work:
  final val INDEX_REQUIRE_ANNO_SPLICEJUNCTIONTREEMAP = 0;
  final val INDEX_REQUIRE_ANNO_KNOWNSPLICEJUNCTIONNAMEMAP = 1;
  final val INDEX_REQUIRE_ANNO_GENEARRAY = 2;
  final val INDEX_REQUIRE_ANNO_CDSARRAY = 3;
  final val INDEX_REQUIRE_ANNO_INTRONARRAY = 4;
  final val INDEX_REQUIRE_ANNO_FLATFEATURELIST = 5;
  final val INDEX_REQUIRE_ANNO_GENELENGTHMAP = 6;
  final val INDEX_REQUIRE_ANNO_FLATEXONARRAY = 7;
  final val INDEX_REQUIRE_ANNO_FLATGENESET = 8;
  
  final val LENGTH_REQUIRE_ANNO_ARRAY = 9;
  
  def initializeCounter[B <: AnyRef](bset : Set[B]) : scala.collection.mutable.Map[B,Int] = {
    val out = scala.collection.mutable.AnyRefMap[B,Int]();
    for(b <- bset){
      out(b) = 0;
    }
    return out;
  }
  
  /*
   * Misc other utils:
   */
   
  private def helper_calculateGeneAssignmentMap_strict(stdGeneArray : GenomicArrayOfSets[String], strandedGeneArray : GenomicArrayOfSets[String]) : Map[String, TreeSet[GenomicInterval]] = {

    val badGeneSet = stdGeneArray.getSteps.foldLeft(Set[String]())( (soFar, curr) =>{
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
    
    reportln("helper_calculateGeneAssignmentMap_strict. Found: " + strandedGeneArray.getValueSet.size + " genes in the supplied annotation.","debug");
    reportln("helper_calculateGeneAssignmentMap_strict. Found: " + badGeneSet.size + " genes with ambiguous segments.","debug");
    
    val buf = strandedGeneArray.getSteps.foldLeft(  Map[String, TreeSet[GenomicInterval] ]() )( (soFar, curr) => {
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
  
  private def getFlatFeatureList(makeFlatReader : (() => Iterator[FlatGtfLine]), stranded  : Boolean, codes : GtfCodes) : IndexedSeq[String] = {
    makeFlatReader().foldLeft( IndexedSeq[String]() )( (soFar,gtfLine) =>{
      soFar :+ gtfLine.getFeatureName;
    });
  }

  private def getGeneLengthMap(geneArray : GenomicArrayOfSets[String]) : Map[String,Int] = {
    return geneArray.getSteps.foldLeft( Map[String,Int]() )( (soFar,currPair) =>{
      val len = currPair._1.end - currPair._1.start;
      currPair._2.foldLeft(soFar)((soFar2,currGene) => {
        soFar.get(currGene) match {
          case Some(lenSoFar) => (soFar2 + ((currGene, len + lenSoFar)) );
          case None => (soFar2 + ((currGene, len)) );
        }
      });
    });
  }
  

  /*
   * Inner Distance:
   */
  private def qcInnerDistance_readSplicesFromGtfFile(makeFlatReader : (() => Iterator[FlatGtfLine]), stranded : Boolean, codes : GtfCodes) : Map[(String,Char),TreeSet[(Int,Int)]] = {
    return makeFlatReader().foldLeft( Map[(String,Char),TreeSet[(Int,Int)]]() )(
        (acc, gtfLine) => {
          qcInnerDistance_readGtfLine(gtfLine,acc, codes);
        }
    )
  }
  
  private def getGeneBiotypeMap(gtffile : String, codes : GtfCodes) : Map[String,String] = {
    reportln("      (DEBUG) Generating Biotype Map ["+getDateAndTimeString+"]","debug");
    
    val reader = GtfReader.getStdGtfReader(gtffile, true, true, "\\s+", codes);
    var out = Map[String,String]();
    var bioset = Set[String]();
    
    for(gtfLine <- reader){
      val geneOpt = gtfLine.getAttribute(codes.GENE_ID_ATTRIBUTE_KEY);
      //if(gtfLine.featureType == codes.STD_EXON_TYPE_CODE){
      if(! geneOpt.isEmpty){
        val geneID = geneOpt.get;
        if(! out.containsKey(geneID)){
          gtfLine.getAttribute(codes.BIOTYPE_ATTRIBUTE_KEY) match {
            case Some(b) => {
              bioset += b;
              out += ((geneID,b));
            }
            case None => {
              //do nothing.
            }
          }
        }
      }
    }
    
    reportln("      (DEBUG) Extracted gene BioType using key \""+codes.BIOTYPE_ATTRIBUTE_KEY+"\".","debug");
    reportln("              Found "+bioset.size+" types: ["+bioset.toList.mkString(",")+"]","debug");
    
    reportln("      (DEBUG) Finished Biotype Map ["+getDateAndTimeString+"]","debug");
    return out.withDefault(x => "UNK");
  }
  
  
  private def qcInnerDistance_readGtfLine(gtfLine : FlatGtfLine, acc :  Map[(String,Char),TreeSet[(Int,Int)]], codes : GtfCodes) : Map[(String,Char),TreeSet[(Int,Int)]] = {
    if(gtfLine.isSpliceJunction){
      acc.get((gtfLine.chromName, gtfLine.strand)) match {
        case Some(chromTree) => {
          return acc + (((gtfLine.chromName, gtfLine.strand), chromTree + ((gtfLine.start-1, gtfLine.end)) ));
        }
        case None => {
          return acc + (((gtfLine.chromName, gtfLine.strand), new TreeSet[(Int,Int)]() + ((gtfLine.start-1, gtfLine.end))));
        }
      }
    } else return acc;
  }
  /*
   * Junction counts:
   */
  private def qcJunctionCounts_makeJunctionMap(makeFlatReader : (() => Iterator[FlatGtfLine]), stranded : Boolean, codes : GtfCodes) : Map[GenomicInterval,String] = {
    return makeFlatReader().foldLeft( Map[GenomicInterval, String]() )(
        (acc, gtfLine) => {
          if(gtfLine.isSpliceJunction){
            acc + ((gtfLine.getGenomicInterval, gtfLine.getFeatureName ));
          } else acc;
        }
    )
  }
  
  /*
   * 
   */
  
  private def qcGetGeneCounts_readGtf(stranded : Boolean, gtffile : String, codes : GtfCodes) : GenomicArrayOfSets[String] = {
    return buildGenomicArrayOfSets_fromGtf(stranded, gtffile, (gtfLine : GtfLine) => gtfLine.featureType == codes.STD_EXON_TYPE_CODE, (gtfLine : GtfLine) => extractGeneId(gtfLine, codes) );
  }

  def extractGeneId(gtfLine : GtfLine, codes : GtfCodes) : String = {
    return gtfLine.getAttributeOrDie(codes.GENE_ID_ATTRIBUTE_KEY);
  }
  
  
  def buildGenomicArrayOfSets_fromGtf(stranded : Boolean, gtffile : String, lineFilter : ( GtfLine => Boolean ) , elementExtractor : (GtfLine => String)) : GenomicArrayOfSets[String] = {
    //report("reading Gtf: " + gtffile + "\n","note");
    val geneArray : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](stranded);
    val gtfReader = GtfReader.getGtfReader(gtffile, stranded, true, "\\s+");
    
    for(gtfLine <- gtfReader){
      if(lineFilter(gtfLine)){
        buildGenomicArrayOfSets_fromGtfLine(gtfLine, stranded, geneArray, elementExtractor);
      }
    }
    return geneArray;
  }
  def buildGenomicArrayOfSets_fromGtf(stranded : Boolean, makeReader : (() => Iterator[GtfLine]), lineFilter : ( GtfLine => Boolean ) , elementExtractor : (GtfLine => String)) : GenomicArrayOfSets[String] = {
    //report("reading Gtf: " + gtffile + "\n","note");
    val geneArray : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](stranded);
    val gtfReader = makeReader();
    
    for(gtfLine <- gtfReader){
      if(lineFilter(gtfLine)){
        buildGenomicArrayOfSets_fromGtfLine(gtfLine, stranded, geneArray, elementExtractor);
      }
    }
    return geneArray;
  }
  
  private def buildGenomicArrayOfSets_fromGtfLine( gtfLine : GtfLine, stranded : Boolean, geneArray : GenomicArrayOfSets[String], elementExtractor : (GtfLine => String)){
      val element = elementExtractor(gtfLine);
      geneArray.addSpan(gtfLine.getGenomicInterval.usingStrandedness(stranded), element);
  }
  
  /*
   * JunctionSeq data extract:
   */
  private def makeFlatExonMap(stranded : Boolean, makeFlatReader : (() => Iterator[FlatGtfLine]), codes : GtfCodes) :  GenomicArrayOfSets[String] = {
    return buildGenomicArrayOfSets_fromGtf(stranded, makeFlatReader, (g : GtfLine) => g.featureType == codes.JS_FEATURETYPE_EXON, (g : GtfLine) => getJSFeatureName(g, codes));
  }
  
  private def getJSFeatureName(gtfLine : GtfLine, codes : GtfCodes) : String = {
       val code = codes.JS_FEATURETYPE_CODEMAP(gtfLine.featureType);
       return gtfLine.getAttributeOrDie(codes.GENE_ID_ATTRIBUTE_KEY) + ":" + code + gtfLine.getAttributeOrDie(codes.JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY);
  }
  
  private def makeFlatGeneSet(exonArray : GenomicArrayOfSets[String]) : Set[String] = {
    return exonArray.getSteps.foldLeft(Set[String]())((soFar,step) =>{
      val (iv, exonSet) = step;
      exonSet.foldLeft(soFar)((soFar2,exonID) => {
        soFar + exonID.split(":")(0);
      });
    });
  }
  
  /************************************
   * geneArea:
   */
  
  private def qcGetGeneCounts_geneArea_CDS_readGtf(stranded : Boolean, gtffile : String, codes : GtfCodes) : GenomicArrayOfSets[String] = {
    return buildGenomicArrayOfSets_fromGtf(stranded, gtffile, (gtfLine : GtfLine) => gtfLine.featureType == codes.STD_CDS_TYPE_CODE, (gtfLine : GtfLine) => extractGeneId(gtfLine, codes));
  }
  private def qcGetGeneCounts_geneArea_regions(stranded : Boolean, gtffile : String, codes : GtfCodes) : GenomicArrayOfSets[String] = {
   // return buildGenomicArrayOfSets_fromGtf(stranded, gtffile, (gtfLine : GtfLine) => gtfLine.featureType == codes.STD_CDS_TYPE_CODE, (gtfLine : GtfLine) => extractGeneId(gtfLine, codes));
    val gtfReader = GtfReader.getGtfReader(gtffile, stranded, true, "\\s+");
    
    reportln("      (Loading gene regions)","debug");
    
    val spanArray = scala.collection.mutable.AnyRefMap[String,GenomicInterval]();
    for(gtfLine <- gtfReader){
      if(gtfLine.featureType == codes.STD_EXON_TYPE_CODE || gtfLine.featureType == codes.STD_CDS_TYPE_CODE){
        val geneID = extractGeneId(gtfLine, codes);
        val curriv = gtfLine.getGenomicInterval.usingStrandedness(stranded);

        spanArray.get(geneID) match {
          case Some(newiv) => {
            spanArray(geneID) = GenomicInterval(curriv.chromName,curriv.strand, Math.min(curriv.start,newiv.start), Math.max(curriv.end,newiv.end));
          }
          case None => {
            spanArray(geneID) = GenomicInterval(curriv.chromName,curriv.strand,curriv.start, curriv.end);
          }
        }
      }
    }
    reportln("      (Generating region array)","debug");
    
    val geneArray : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](stranded);
    for(geneID <- spanArray.keySet){
      val iv = spanArray(geneID);
      geneArray.addSpan(iv, geneID);
    }
    
    reportln("      (Region array complete)","debug");
    return geneArray;
  }
  
  
  private def qcGetGeneCounts_geneArea_INTRONS_readFlatGtf(stranded : Boolean, makeFlatReader : (() => Iterator[FlatGtfLine]), codes : GtfCodes) : GenomicArrayOfSets[String] = {
    return buildGenomicArrayOfSets_fromGtf(stranded, makeFlatReader, (gtfLine : GtfLine) => gtfLine.featureType == codes.JS_FEATURETYPE_KNOWNSPLICE || gtfLine.featureType == codes.JS_FEATURETYPE_NOVELSPLICE, (gtfLine : GtfLine) => extractGeneId(gtfLine, codes));
  }
  
}




















