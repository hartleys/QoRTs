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


object qcJunctionCounts {
  
  
  
  def getSpliceJunctionSet(r1 : SAMRecord, r2 : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean) : Set[GenomicInterval] = {
    getSpliceJunctionSetFromRead(r1,stranded,fr_secondStrand) ++ getSpliceJunctionSetFromRead(r2,stranded,fr_secondStrand)
  }
  
  def getSpliceJunctionSetFromRead(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean) : Set[GenomicInterval] = {
    val cigOps : Stream[CigOp] = CigarHolder(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean).cigOps;
    
    cigOps.filter( _.op == CigarOperator.SKIPPED_REGION).map(_.ref_iv).toSet;
  }
  /*
  def makeAnnotatedCountMap(spliceAnnotation : Map[(String,Char),TreeSet[(Int,Int)]]) : scala.collection.mutable.Map[GenomicInterval,Int] = {
    spliceAnnotation.foldLeft( new scala.collection.mutable.HashMap[GenomicInterval,Int]() )( (soFar, curr) =>{
      val ((chromName, strand), ts) = curr;
      soFar ++ ts.map(iv =>{
        (new GenomicInterval(chromName, strand, iv._1, iv._2), 0);
      });
    })
  }
  
  def makeAnnotatedCountMap(gtffile : String, stranded : Boolean) : Map[GenomicInterval,(Int, String)] = {
    return GtfReader.getGtfReader(gtffile, stranded, true, "\\s+").foldLeft( Map[GenomicInterval, (Int, String)]() )(
        (acc, gtfLine) => {
          readGtfLine(gtfLine,acc);
        }
    )
  }
  
  def readGtfLine(gtfLine : GtfLine, acc :  Map[GenomicInterval,(Int, String)]) : Map[GenomicInterval,(Int, String)] = {
    if(gtfLine.featureType == "splice_site" || gtfLine.featureType == "novel_splice_site"){
      val iv = new GenomicInterval(gtfLine.chromName, gtfLine.strand, gtfLine.start - 1, gtfLine.end);
      val name : String = GtfCodes.getJSFeatureName(gtfLine);
      return acc + ((iv, ((0,name))));
    } else return acc;
  }*/
  
  def getExonsAndGenesFromPair(r1 : SAMRecord, r2 : SAMRecord, flatExonArray : GenomicArrayOfSets[String], stranded : Boolean, fr_secondStrand : Boolean) : (Set[String],Set[String]) = {
    val r1e = getExonsFromRead(r1,flatExonArray,stranded, fr_secondStrand);
    val r2e = getExonsFromRead(r2,flatExonArray,stranded, fr_secondStrand);
    
    val exonSet = (r1e ++ r2e);
    val geneSet = exonSet.map( getAggregateGeneID(_) );
    
    return ((exonSet, geneSet));
  }
  
  def getExonsFromRead(r : SAMRecord, flatExonArray : GenomicArrayOfSets[String], stranded : Boolean, fr_secondStrand : Boolean) : Set[String] = {
    val readIntervals : Iterator[GenomicInterval] = getGenomicIntervalsFromRead(r , stranded , fr_secondStrand);
    return readIntervals.foldLeft(Set[String]())((soFar, iv) => {
      //val (exonsSoFar,genesSoFar) = pairSoFar;
      flatExonArray.findIntersectingSteps(iv).foldLeft( soFar )((psf,currPair) => {
        val exonsSoFar = psf;
        val (stepiv,featureSet) = currPair;
        exonsSoFar ++ featureSet;
      })
    })
  }
  
}
//runFunc.contains("writeDEXSeq"), runFunc.contains("writeSpliceExon"), runFunc.contains("writeKnownSplices"), runFunc.contains("writeNovelSplices")
class qcJunctionCounts(anno_holder : qcGtfAnnotationBuilder, stranded : Boolean, fr_secondStrand : Boolean, 
                       writeDEXSeq : Boolean, writeSpliceExon : Boolean, writeKnownSplices : Boolean, writeNovelSplices : Boolean, annotatedSpliceExonCounts : Boolean,
                       limit_high_expressed_sj : Int = 3)  extends QCUtility[Unit] {
  reportln("> Init JunctionCalcs utility","debug");
  
  val knownSpliceMap : GenMap[GenomicInterval, String] = anno_holder.knownSpliceJunctionNameMap;
  val flatFeatureList : IndexedSeq[String] = anno_holder.flatFeatureList;
  var knownCountMap : GenMap[String,(GenomicInterval,Int)] = knownSpliceMap.map( pair => {
      val (iv, spliceID) = pair;
      ((spliceID,((iv,0))));
    }) 
  reportln("length of knownSpliceMap after instantiation: " + knownSpliceMap.size,"debug");  
  reportln("length of knownCountMap after instantiation: "  + knownCountMap.size,"debug");  
  
  var novelCountMap : GenMap[GenomicInterval,Int] = (Map[GenomicInterval,Int]()).withDefault(k => 0);
  
  val flatExonArray : GenomicArrayOfSets[String] = anno_holder.flatExonArray;
  val flatGeneSet   : Set[String] = anno_holder.flatGeneSet;
  
  val exonCountMap :     scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](flatExonArray.getValueSet);
  val flatGeneCountMap : scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](flatGeneSet);
  
  //val flatGeneCountMap : scala.collection.mutable.Map[String,Int] = flatExonArray.getSteps.map(pair => { pair._2.map() })
  //TO DO: junction seq file generation!
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    countSplices( qcJunctionCounts.getSpliceJunctionSet(r1,r2,stranded,fr_secondStrand) );
    countExons(r1,r2);
  }
  
  def countSplices(splices : Set[GenomicInterval]){
    for(iv <- splices){
      if(knownSpliceMap.contains(iv)){
        val spliceID : String = knownSpliceMap(iv);
        val (iv2, ct) = knownCountMap(spliceID);
        knownCountMap = knownCountMap + ((spliceID, (( iv, ct + 1))));
      } else {
        val ct = novelCountMap(iv);
        novelCountMap = novelCountMap + ((iv, ct + 1));
      }
    }
  }
  def countExons(r1 : SAMRecord, r2 : SAMRecord){
    val (readExonSet, readGeneSet) = qcJunctionCounts.getExonsAndGenesFromPair(r1,r2,flatExonArray, stranded, fr_secondStrand);
    if(readGeneSet.size == 1){
      val readGene = readGeneSet.head;
      flatGeneCountMap(readGene) += 1;
      for(e <- readExonSet){
        exonCountMap(e) += 1;
      }
    }
  }
  //writeDEXSeq : Boolean, writeSpliceExon : Boolean, writeKnownSplices : Boolean, writeNovelSplices : Boolean
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    report("length of knownCountMap after run: " + knownCountMap.size,"debug");
    
    if(writeKnownSplices){
      val writer = openWriterSmart_viaGlobalParam(outfile + ".spliceJunctionCounts.knownSplices.txt");
      writer.write("spliceName	chrom	strand	start	end	CT\n");
      for(spliceID <- knownCountMap.keys.toVector.sorted){
        val (iv,ct) = knownCountMap(spliceID);
        writer.write(spliceID + "	"+iv.chromName+"	"+iv.strand+"	"+iv.start+"	"+iv.end+"	"+ct+"\n");
      }
      writer.close();
    }
    
    if(writeNovelSplices){
      val writer2 = openWriterSmart_viaGlobalParam(outfile + ".spliceJunctionCounts.novelSplices.txt");
      writer2.write("chrom	strand	start	end	CT\n");
      for(iv <- novelCountMap.keys.toVector.sorted){
        val ct = novelCountMap(iv);
        writer2.write(iv.chromName+"	"+iv.strand+"	"+iv.start+"	"+iv.end+"	"+ct+"\n");
      }
      writer2.close();
    }
    val featureList = flatFeatureList;

    if(writeSpliceExon){
      val writer3 = openWriterSmart_viaGlobalParam(outfile + ".spliceJunctionAndExonCounts.forJunctionSeq.txt");
      //val featureList = (flatGeneCountMap.keySet.map(_ + ":G000") ++ knownCountMap.keySet ++ exonCountMap.keySet).toSeq.sortBy( k => {
      //  val ks = k.split(":");
      //  val mainID = ks(0);
      //  val subID  = ks(1).takeRight(3);
      //  mainID + subID;
      //});
      for(f <- featureList){
        val featureCode = f.split(":")(1).charAt(0);
        if(featureCode == 'A'){
          writer3.write( f +"	"+ flatGeneCountMap(f.split(":")(0)) +"\n");
        } else if(featureCode == 'J' || featureCode == 'N'){
          writer3.write( f +"	"+ knownCountMap(f)._2 +"\n");
        } else if(featureCode == 'E'){
          writer3.write( f +"	"+ exonCountMap(f) +"\n");
        } else {
          error("IMMPOSSIBLE STATE! FATAL ERROR! qcJunctionCounts.writeOutput, writing forSpliceSeq");
        }
      }
      writer3.close();
    }
    

    if(annotatedSpliceExonCounts){
      val writer4 = openWriterSmart_viaGlobalParam(outfile + ".annoSpliceJunctionAndExonCounts.txt");
      val flatgff = anno_holder.makeFlatReader();
      
      writer4.write("featureID	featurType	chrom	start	end	strand	geneID	binID	readCount\n");
      flatgff.foreach( (gffline : FlatGtfLine) => {
        val f = gffline.getFeatureName;
        val featureCode = f.split(":")(1).charAt(0);
        
        if(featureCode == 'A'){
          writer4.write( f +"	"+gffline.featureType+"	"+gffline.chromName+"	"+gffline.start+"	"+gffline.end+"	"+gffline.strand+"	"+gffline.getFeatureAggregateGene+"	"+gffline.getFeaturePartNumber+"	"+ flatGeneCountMap(f.split(":")(0)) +"\n");
        } else if(featureCode == 'J' || featureCode == 'N'){
          writer4.write( f +"	"+gffline.featureType+"	"+gffline.chromName+"	"+gffline.start+"	"+gffline.end+"	"+gffline.strand+"	"+gffline.getFeatureAggregateGene+"	"+gffline.getFeaturePartNumber+"	"+  knownCountMap(f)._2 +"\n");
        } else if(featureCode == 'E'){
          writer4.write( f +"	"+gffline.featureType+"	"+gffline.chromName+"	"+gffline.start+"	"+gffline.end+"	"+gffline.strand+"	"+gffline.getFeatureAggregateGene+"	"+gffline.getFeaturePartNumber+"	"+  exonCountMap(f) +"\n");
        } else {
          error("IMMPOSSIBLE STATE! FATAL ERROR! qcJunctionCounts.writeOutput, writing forSpliceSeq");
        }
      });
      writer4.close();
    }
    
    if(writeDEXSeq){
      val writer4 = openWriterSmart_viaGlobalParam(outfile + ".exonCounts.formatted.for.DEXSeq.txt");
      for(f <- featureList){
        val featureCode = f.split(":")(1).charAt(0);
        if(featureCode == 'E'){
          val split = f.split(":");
          val ag = split(0);
          val num = split(1).takeRight(3);
          writer4.write( ag +":"+num +"	"+ exonCountMap(f) +"\n");
        }
      }
      writer4.write("_ambiguous	NA\n");
      writer4.write("_empty	NA\n");
      writer4.write("_lowaqual	NA\n");
      writer4.write("_notaligned	NA\n");
      close(writer4);
    }
    
    //SUMMARY DATA:
    val num_aggregate_genes = flatGeneSet.size;
    val num_aggregate_genes_with_no_reads = flatGeneCountMap.count( _._2 == 0 );
    val num_aggregate_genes_with_reads = flatGeneCountMap.count( _._2 != 0 );
    val num_known_sj = knownCountMap.size;
    val num_known_sj_with_no_reads = knownCountMap.count( _._2._2 == 0 );
    val num_known_sj_with_few_reads = knownCountMap.count(x => x._2._2 <= limit_high_expressed_sj && x._2._2 != 0 );
    val num_known_sj_with_many_reads = knownCountMap.count(x => x._2._2 > limit_high_expressed_sj );
    val num_novel_sj = novelCountMap.size;
    val num_novel_sj_with_few_reads = novelCountMap.count(_._2 <= limit_high_expressed_sj);
    val num_novel_sj_with_many_reads = novelCountMap.count(_._2 > limit_high_expressed_sj);
    
    val num_events_known_sj = knownCountMap.map(_._2._2).sum
    val num_events_known_sj_with_few_reads = knownCountMap.map(x => {
        if(x._2._2 <= limit_high_expressed_sj) x._2._2;
        else 0;
    }).sum
    val num_events_known_sj_with_many_reads = knownCountMap.map(x => {
        if(x._2._2 <= limit_high_expressed_sj) 0;
        else x._2._2;
    }).sum;
    
    val num_events_novel_sj = novelCountMap.map(_._2).sum;
    val num_events_novel_sj_with_few_reads = novelCountMap.map(x => {
        if(x._2 <= limit_high_expressed_sj) x._2;
        else 0;
    }).sum
    val num_events_novel_sj_with_many_reads = novelCountMap.map(x => {
        if(x._2 <= limit_high_expressed_sj) 0;
        else x._2;
    }).sum
    
    summaryWriter.write("AggregateGenes	"                  + num_aggregate_genes + "\n");
    summaryWriter.write("AggregateGenes_NoReads	"          + num_aggregate_genes_with_no_reads + "\n");
    summaryWriter.write("AggregateGenes_WithReads	"      + num_aggregate_genes_with_reads + "\n");
    summaryWriter.write("SpliceLoci	"                      + (num_known_sj + num_novel_sj) + "\n");    
    summaryWriter.write("SpliceLoci_Known	"              + num_known_sj + "\n");
    summaryWriter.write("SpliceLoci_Known_NoReads	"      + num_known_sj_with_no_reads + "\n");
    summaryWriter.write("SpliceLoci_Known_FewReads	"      + num_known_sj_with_few_reads + "\n");
    summaryWriter.write("SpliceLoci_Known_ManyReads	"      + num_known_sj_with_many_reads + "\n");
    summaryWriter.write("SpliceLoci_Novel	"              + num_novel_sj + "\n");
    summaryWriter.write("SpliceLoci_Novel_FewReads	"      + num_novel_sj_with_few_reads + "\n");
    summaryWriter.write("SpliceLoci_Novel_ManyReads	"      + num_novel_sj_with_many_reads + "\n");
    
    summaryWriter.write("SpliceEvents	"                         + (num_events_novel_sj + num_events_known_sj) + "\n");
    summaryWriter.write("SpliceEvents_KnownLoci	"                 + num_events_known_sj + "\n");
    summaryWriter.write("SpliceEvents_KnownLociWithFewReads	"     + num_events_known_sj_with_few_reads + "\n");
    summaryWriter.write("SpliceEvents_KnownLociWithManyReads	" + num_events_known_sj_with_many_reads + "\n");
    summaryWriter.write("SpliceEvents_NovelLoci	"                 + num_events_novel_sj + "\n");
    summaryWriter.write("SpliceEvents_NovelLociWithFewReads	"     + num_events_novel_sj_with_few_reads + "\n");
    summaryWriter.write("SpliceEvents_NovelLociWithManyReads	" + num_events_novel_sj_with_many_reads + "\n");
    
  }
  def getUtilityName : String = "JunctionCalcs";
}
















