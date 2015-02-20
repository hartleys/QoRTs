package fileConversionUtils

import java.util.zip.GZIPInputStream

import java.util.zip.GZIPOutputStream
import java.io.OutputStream
import java.io.FileOutputStream
import java.io.InputStream
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.File
import scala.collection.JavaConversions._
import scala.collection.mutable.HashMap;

import net.sf.samtools._

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

import internalUtils.commonSeqUtils._;
import internalUtils.genomicUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.JavaConversions._

object prepFlatGtfFile {

  class prepFlatGtfFile_runner extends CommandLineRunUtil {
    override def priority = 50;
    val parser : CommandLineArgParser = 
      new CommandLineArgParser(
          command = "makeFlatGff", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "When running the QC command, QoRT first generates a set of non-overlapping exonic "+
                        "fragments out of all the exons in the genome annotation gtf file. It then assigns each exonic"+
                        " fragment a unique identifier. Similarly, it assigns every splice junction its own unique identifier. "+
                        "This command can be used to write that data to file.\n"+
                        "It can also be used to produce a flattened gff file that adheres to the specifications used by DEXSeq.",   
          argList = 
            new UnaryArgument(name = "stranded",
                              arg = List("--stranded","-r"), // name of value
                              argDesc = "DEPRECIATED" // description
                              ) :: 
            new UnaryArgument(name = "dexseqFmt",
                              arg = List("--DEXSeqFmt"), // name of value
                              argDesc = "Flag to indicate that the output gff file should be formatted for use with DEXSeq." // description
                              ) :: 
            new FinalArgument[String](
                                         name = "infile",
                                         valueName = "gtffile",
                                         argDesc = "The gtf annotation file. This tool was designed to use the standard gtf annotations provided by Ensembl, but other annotations can be used as well. Note: if the file ends in .zip or .gz the compression method will be auto-detected and read accordingly." // description
                                        ) ::
            new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "flatgfffile",
                                         argDesc = "The output destination for the \"flattened\" gff annotation file to be created, or '-' to write to stdout. Note: if the filename ends in \".zip\" or \".gz\" the corresponding compression method will be applied." // description
            ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS
      );
    
    def run(args : Array[String]){
      val out = parser.parseArguments(args.toList.tail);
      
      if(out){      
      prepFlatGtfFile.run(parser.get[String]("infile"),
                          parser.get[String]("outfile"), 
                          parser.get[Boolean]("stranded"),
                          parser.get[Boolean]("dexseqFmt")
                          );
      }
    }
  }
  
  def run(infile : String, outfile : String, stranded : Boolean, dexseqFmt : Boolean){
    buildFlatGtf(infile , stranded , outfile , dexseqFmt);
  }
  
  /*
  private def buildGenomicArrayOfSets_fromGtf(stranded : Boolean, gtffile : String, lineFilter : ( GtfLine => Boolean ) , elementExtractor : (GtfLine => String)) : GenomicArrayOfSets[String] = {
    //report("reading Gtf: " + gtffile + "\n","note");
    val geneArray : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](stranded);
    val gtfReader = GtfReader.getGtfReader(gtffile, stranded, true, "\\s+");
    
    for(gtfLine <- gtfReader){
      if(lineFilter(gtfLine)){
        buildGenomicArrayOfSets_fromGtfLine(gtfLine, geneArray, elementExtractor);
      }
    }
    return geneArray;
  }
  
  private def buildGenomicArrayOfSets_fromGtfLine( gtfLine : GtfLine, geneArray : GenomicArrayOfSets[String], elementExtractor : (GtfLine => String)){
      val element = elementExtractor(gtfLine);
      geneArray.addSpan(gtfLine.getGenomicInterval, element);
  }
  
  private def readGtf(stranded : Boolean, gtffile : String) : GenomicArrayOfSets[String] = {
    return buildGenomicArrayOfSets_fromGtf(stranded, gtffile, (gtfLine) => gtfLine.featureType == GtfCodes.EXON_TYPE_CODE, extractTxId(_) );
  }*/
  //private def extractGeneId(gtfLine : GtfLine) : String = {
  //  return gtfLine.getAttributeOrDie(GtfCodes.GENE_ID_ATTRIBUTE_KEY);
  //}
  //private def extractTxId(gtfLine : GtfLine) : String = {
  //  return gtfLine.getAttributeOrDie(GtfCodes.STD_TX_ID_ATTRIBUTE_KEY);
  //}
  
  private def buildGenomicArrayOfSets_Tx_and_Map(stranded : Boolean, gtffile : String, codes : GtfCodes = new GtfCodes()) : (GenomicArrayOfSets[String], Map[String,String]) = {
    val txArray : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](stranded);
    val gtfReader = GtfReader.getStdGtfReader(gtffile, stranded, true, "\\s+", codes);
    var txMap : Map[String,String] = Map[String,String]();
    
    for(gtfLine <- gtfReader){
      if(gtfLine.isExon){
        val tx = gtfLine.getTxID;
        val gene = gtfLine.getGeneID;
        txMap += (tx -> gene);
        txArray.addSpan(gtfLine.getGenomicInterval, tx);
      }
    }
    
    return ((txArray.finalizeStepVectors, txMap));
  }
  
  private def buildSpliceJunctionMap(txArray : GenomicArrayOfSets[String]) : (Map[GenomicInterval, Set[String]]) = {
    var spliceJunctionMap : Map[GenomicInterval, Set[String]] = Map[GenomicInterval, Set[String]]().withDefault(k => {Set[String]()})
    
    for((chromName, strand) <- txArray.getChroms){
      var lastExonTxMap : Map[String,Int] = Map[String,Int]();
       
      for((iv, txSet) <- txArray.getSteps(chromName,strand)){
        for(tx <- txSet){
          if(lastExonTxMap.containsKey(tx)){
            val prevExonEnd = lastExonTxMap(tx);
            if(iv.start != prevExonEnd){
              val spliceIV = new GenomicInterval(chromName,strand,prevExonEnd,iv.start);
              spliceJunctionMap = spliceJunctionMap.updated(spliceIV , (spliceJunctionMap(spliceIV) + tx) );
            }
          }
          lastExonTxMap = lastExonTxMap.updated(tx , iv.end);
        }
      }
    }
    return spliceJunctionMap;
  }
  private def buildGeneSets(txArray : GenomicArrayOfSets[String], txMap : Map[String,String]) : (Map[String,Set[String]]) = {
    var geneSets : Map[String,Set[String]] = Map[String,Set[String]]().withDefault(k => Set[String]());
    
    for((iv,txSet) <- txArray.getSteps){
      var currGeneSet = Set[String]();
      for(tx <- txSet){
        val g = txMap(tx);
        currGeneSet += g;
        currGeneSet = currGeneSet ++ geneSets(g);
      }
      for(g <- currGeneSet){
        assert(geneSets(g).subsetOf(currGeneSet));
        geneSets = geneSets.updated(g , currGeneSet);
      }
    }
    
    return geneSets;
  }
  private def buildAggregateGeneMap(geneSets : Map[String,Set[String]]) : (Map[String,String],Set[String]) = {
    val aggregateGeneMap : Map[String,String] = geneSets.map( pair => {
      (pair._1 -> pair._2.toSeq.sorted.mkString("+"));
    })
    val aggregateSet = aggregateGeneMap.values.toSet;
    return ((aggregateGeneMap, aggregateSet));
  }
  
  private def buildFeatures2(stranded : Boolean, 
                             txArray : GenomicArrayOfSets[String], 
                             txMap : Map[String,String], 
                             spliceJunctionMap : Map[GenomicInterval,Set[String]], 
                             aggregateGeneMap : Map[String,String], 
                             aggregateSet : Set[String], gtfCodes : GtfCodes = new GtfCodes()) : Iterator[FlatGtfLine] = {
    //var featureCountMap : Map[String,Int] = Map[String,Int]().withDefault(k => 0);
    //var featureList : Seq[GtfLine] = Seq[GtfLine]();
    
    //var featureListMap : Map[String,Seq[GtfLine]] = Map[String,Seq[GtfLine]]().withDefault((k) => Seq[GtfLine]());
    
    //var aggregateIntervalMap : Map[String,GenomicInterval] = Map[String,GenomicInterval]();

    reportln("FlatteningGtf: Iterating through the step-vector...("+getDateAndTimeString+")","debug");
    val featureListMap = txArray.getSteps.foldLeft( Map[String,IndexedSeq[FlatGtfLine]]().withDefault( (k) => IndexedSeq[FlatGtfLine]() ) )( (soFar, pair) =>{
      val (iv,txSet) = pair;
      if(! txSet.isEmpty){
        val geneSet : Set[String] = txSet.map(tx => txMap(tx));
        val aggregateGene : String = aggregateGeneMap(geneSet.head);
        val features : IndexedSeq[FlatGtfLine] = soFar(aggregateGene);
        val featureCt : Int = features.size;
        val gtfLine : FlatGtfLine = FlatOutputGtfLine.makeFlatGtfLine_feature(iv, gtfCodes.JS_FEATURETYPE_EXON, stranded, aggregateGene, txSet, geneSet, featureCt + 1, gtfCodes);
      
        soFar.updated(aggregateGene, features :+ gtfLine);
      } else {
        soFar;
      }
    });
    
    reportln("FlatteningGtf: Adding the aggregate genes themselves...("+getDateAndTimeString+")","debug");
    val featureListMap2 = aggregateSet.foldLeft(featureListMap)((soFar, aggregateGene) =>{
      val features = soFar(aggregateGene);
      val iv = new GenomicInterval(features.head.chromName, features.head.strand, features.head.start, features.last.end);
      val gtfLine = FlatOutputGtfLine.makeFlatGtfLine_aggregateGene(iv, stranded , aggregateGene , gtfCodes);
      soFar.updated(aggregateGene, gtfLine +: features);
    });
    
    reportln("FlatteningGtf: Iterating through the splice junctions...("+getDateAndTimeString+")","debug");
    val featureListMapFinal = spliceJunctionMap.keys.toSeq.sorted.foldLeft(featureListMap2)( (soFar, iv) =>{
      val txSet : Set[String] = spliceJunctionMap(iv);
      val geneSet : Set[String] = txSet.map(tx => txMap(tx));
      val aggregateGene : String = aggregateGeneMap(geneSet.head);
      val features : IndexedSeq[FlatGtfLine] = soFar(aggregateGene);
      val featureCt : Int = features.size;
      val gtfLine : FlatGtfLine = FlatOutputGtfLine.makeFlatGtfLine_feature(iv, gtfCodes.JS_FEATURETYPE_KNOWNSPLICE, stranded, aggregateGene, txSet, geneSet, featureCt, gtfCodes);
      soFar.updated(aggregateGene, features :+ gtfLine);
    })
    
    reportln("FlatteningGtf: Sorting the aggregate genes...("+getDateAndTimeString+")","debug");
    val aggregateGeneList = featureListMapFinal.keys.toSeq.sortBy((ag) => featureListMapFinal(ag).head.getGenomicInterval );
    
    //return new Iterator[GtfLine]{
    //  val linemap = featureListMapFinal;
    //  val geneIter = aggregateGeneList.iterator;
    //  var currIter = linemap(geneIter.next).iterator;
    // 
    //  def hasNext : Boolean = geneIter.hasNext || currIter.hasNext;
    //  def next : GtfLine = {
    //    if(currIter.hasNext){
    //      currIter.next;
    //    } else {
    //      currIter = linemap(geneIter.next).iterator;
    //      currIter.next;
    //    }
    //  }
    //}
    reportln("FlatteningGtf: Folding the FlatGtfLine iterator...("+getDateAndTimeString+")","debug");
    return aggregateGeneList.foldLeft( Iterator[FlatGtfLine]() )((soFar, curr) =>{
      soFar ++ featureListMapFinal(curr).iterator;
    });
  }
  
  /*private def buildFeatures(stranded : Boolean, txArray : GenomicArrayOfSets[String], txMap : Map[String,String], spliceJunctionMap : Map[GenomicInterval,Set[String]], aggregateGeneMap : Map[String,String], aggregateSet : Set[String]) : Seq[GtfLine] = {
    var featureCountMap : Map[String,Int] = Map[String,Int]().withDefault(k => 0);
    var featureList : Seq[GtfLine] = Seq[GtfLine]();
    
    //var featureListMap : Map[String,Seq[GtfLine]] = Map[String,Seq[GtfLine]]().withDefault((k) => Seq[GtfLine]());
    
    var aggregateIntervalMap : Map[String,GenomicInterval] = Map[String,GenomicInterval]();
    
    reportln("Iterating through the step-vector...","progress");
    //Iterate through the transcript step-vector.
    for((iv,txSet) <- txArray.getSteps){
      //If there is at least 1 transcript on the step:
      if(txSet.size > 0){
        val geneSet : Set[String] = txSet.map(tx => txMap(tx));
        val aggregateGene : String = aggregateGeneMap(geneSet.head);
        featureCountMap = featureCountMap.updated(aggregateGene,featureCountMap(aggregateGene) + 1);
        val featureNumber = featureCountMap(aggregateGene);
        
        featureList = featureList :+ OutputGtfLine.makeJSGtfLine(iv, GtfCodes.JS_FEATURETYPE_EXON, stranded, aggregateGene, txSet, geneSet, featureNumber);

        aggregateIntervalMap.get(aggregateGene) match {
          case Some(old_iv) => {
            val new_iv = new GenomicInterval(iv.chromName,iv.strand, math.min(iv.start, old_iv.start), math.max(iv.end,old_iv.end));
            aggregateIntervalMap = aggregateIntervalMap.updated(aggregateGene, new_iv);
          }
          case None => {
            aggregateIntervalMap = aggregateIntervalMap.updated(aggregateGene, iv);
          }
        }
      }
    }
    
    reportln("Iterating through the splice junctions...","progress");
    for(iv <- spliceJunctionMap.keys.toSeq.sorted){
      val txSet : Set[String] = spliceJunctionMap(iv);
      val geneSet : Set[String] = txSet.map(tx => txMap(tx));
      val aggregateGene : String = aggregateGeneMap(geneSet.head);
      featureCountMap = featureCountMap.updated(aggregateGene,featureCountMap(aggregateGene) + 1);
      val featureNumber = featureCountMap(aggregateGene);
      featureList = featureList :+ OutputGtfLine.makeJSGtfLine(iv, GtfCodes.JS_FEATURETYPE_KNOWNSPLICE, stranded, aggregateGene, txSet, geneSet, featureNumber);
    }
    
    reportln("Adding the aggregate genes...","progress");
    for(aggregateGene <- aggregateSet){
      val iv = aggregateIntervalMap(aggregateGene);
      featureList = featureList :+ OutputGtfLine.makeJSGtfLine_aggregateGene(iv, stranded , aggregateGene );
    }
    
    reportln("Sorting the output:","progress")
    val sortedFeatureList = featureList.sortBy(gtfLine => ((gtfLine.getGenomicInterval.chromName, gtfLine.getGenomicInterval.start, gtfLine.getAttributeOrDie(GtfCodes.JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY))) );
    
    return sortedFeatureList;
  }*/
  
        //OutputGtfLine.makeJSGtfLine_aggregateGene(iv : GenomicInterval, stranded : Boolean, aggregateGene : String)
        //OutputGtfLine.makeJSGtfLine(iv : GenomicInterval, featureType : String, stranded : Boolean, aggregateGene : String, txSet : Set[String], geneSet : Set[String], featureNumber : Int)
        //OutputGtfLine(iv : .GenomicInterval, featureType : String, attributeMap : Map[String,String], stranded : Boolean)
  

  
  private def getFlatGtfIterator(infile : String, stranded : Boolean) : Iterator[FlatGtfLine] = {
    reportln("FlatteningGtf: starting...("+getDateAndTimeString+")","debug");
    val (txArray, txMap) : (GenomicArrayOfSets[String], Map[String,String]) = buildGenomicArrayOfSets_Tx_and_Map(stranded,infile);
    reportln("FlatteningGtf: gtf file read complete.("+getDateAndTimeString+")","debug");
    val spliceJunctionMap : Map[GenomicInterval, Set[String]] = buildSpliceJunctionMap(txArray);
    reportln("FlatteningGtf: Splice Junction Map read.("+getDateAndTimeString+")","debug");
    val geneSets : Map[String,Set[String]] = buildGeneSets(txArray, txMap);
    reportln("FlatteningGtf: gene Sets generated.("+getDateAndTimeString+")","debug");
    val (aggregateGeneMap, aggregateSet) : (Map[String,String],Set[String]) = buildAggregateGeneMap(geneSets);
    reportln("FlatteningGtf: Aggregate Sets built.","debug");

    //val sortedFeatureList = buildFeatures(stranded, txArray, txMap , spliceJunctionMap , aggregateGeneMap , aggregateSet );
    val sortedFeatureIterator = buildFeatures2(stranded, txArray, txMap , spliceJunctionMap , aggregateGeneMap , aggregateSet);
    reportln("FlatteningGtf: Features Built.("+getDateAndTimeString+")","debug");
    return(sortedFeatureIterator);
  }
  
  private def writeFlatGtf_stdFormat(sortedFeatureIterator : Iterator[FlatGtfLine], outfile : String){
    val writer = openWriterSmart(outfile, true);
    for(gtfLine <- sortedFeatureIterator){
      writeGtfLine(gtfLine , writer );
    }
    close(writer);
  }
  private def writeFlatGtf_DEXSeqFormat(sortedFeatureIterator : Iterator[FlatGtfLine], outfile : String){
    val writer = openWriterSmart(outfile, true);
    for(gtfLine <- sortedFeatureIterator){
      val dss = gtfLine.getDexSeqStr;
      dss match {
        case Some(dssl : String) => writer.write(dssl + "\n");
        case None => {} //do nothing
      }
    }
    close(writer);
  }
  
  private def buildFlatGtf(infile : String, stranded : Boolean, outfile : String, dexseqFmt : Boolean){
    val sortedFeatureIterator = getFlatGtfIterator(infile , stranded );
    
    if(dexseqFmt){
      writeFlatGtf_DEXSeqFormat(sortedFeatureIterator, outfile);
    } else {
      writeFlatGtf_stdFormat(sortedFeatureIterator, outfile);
    }
    
    reportln("FlatteningGtf: done.("+getDateAndTimeString+")","progress");

    //buildSpliceJunctionMap(txArray : GenomicArrayOfSets[String]) : (Map[GenomicInterval, Set[String]])
  }
  
  def getFlatGtfLines(infile : String, stranded : Boolean) : IndexedSeq[FlatGtfLine] = {
    val sortedFeatureIterator = getFlatGtfIterator(infile,stranded);
    return sortedFeatureIterator.toVector;
  } 
  
}


























