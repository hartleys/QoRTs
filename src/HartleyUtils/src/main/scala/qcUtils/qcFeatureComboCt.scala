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
//import scala.collection.JavaConversions._;
import scala.collection.JavaConverters._;

import internalUtils.genomicUtils._;
import internalUtils.optionHolder._;
import scala.collection.GenMap;

object qcFeatureComboCt {
  
  
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


class qcFeatureComboCt(anno_holder : qcGtfAnnotationBuilder, stranded : Boolean, fr_secondStrand : Boolean, 
                       
                       limit_high_expressed_sj : Int = 3, outfile : String) extends QCUtility[Unit] {
  reportln("> Init qcFeatureComboCt Utility","note");
  reportln("> WARNING: THIS UTILITY IS BETA! IT IS NOT FOR PRODUCTION USE!","warn");
  
  val goodDistance = List(0,1,2,3,4,5,6,7,8,9,10);
  
  val knownSpliceMap : GenMap[GenomicInterval, String] = anno_holder.knownSpliceJunctionNameMap;
  val flatFeatureList : IndexedSeq[String] = anno_holder.flatFeatureList;
  //var knownCountMap : GenMap[String,(GenomicInterval,Int)] = knownSpliceMap.map( pair => {
  //    val (iv, spliceID) = pair;
  //    ((spliceID,((iv,0))));
  //  }) 
  reportln("length of knownSpliceMap after instantiation: " + knownSpliceMap.size,"debug");
  //reportln("length of knownCountMap after instantiation: "  + knownCountMap.size,"debug");
  
  val flatGtf = anno_holder.makeFlatReader().toVector;
  val txList : Set[String] = flatGtf.foldLeft(Set[String]())((soFar, line) =>  soFar ++ line.getTxSet.split("\\+").toSet  );
  
  //reportln("txList = ["+txList.mkString(",")+"]","debug");
  
  val txMap : Vector[(String,Set[String])] = flatGtf.filter(! _.isAggregateGene).map((line) => { (( line.getFeatureName, line.getTxSet.split("\\+").toSet)) }).toVector;
  val txFeatureMap : Map[String,Set[String]] = txList.map((tx : String) => {
    val featureSet = txMap.filter( _._2.contains(tx) ).map(_._1).toSet;
    ((tx, featureSet));
  }).toMap;
  
  //txList.map((tx : String) => {
  //  reportln("TX:\t"+ tx + "\t[" +  flatGtf.filter(_.getTxSet.split("\\+").toSet.contains(tx)).mkString("   ") +"]","debug")
  //});
  //reportln("-------------","debug")
  //txList.map((tx : String) => {
  //  reportln("TX:\t"+ tx + "\t[" +  txFeatureMap(tx) +"]","debug")
  //});
  //reportln("-------------","debug")
  //txMap.map{ case (feature,txSet) => {
  //  reportln("feature:\t"+ feature + "\t[" +  txSet.toList.mkString("   ") +"]","debug")
  //}}
  //reportln("-------------","debug")
  //flatGtf.map{ (line) => {
  //  reportln("feature:\t"+ line.getFeatureName + "\t[" +  line.getTxSet +"]","debug")
  //}}
  
  val endpointMap : Map[String,(Int,Int)] = txList.map((tx : String) => {
    val txLines = flatGtf.filter(! _.isAggregateGene).filter(_.getTxSet.split("\\+").toSet.contains(tx))
    (tx, (txLines.map(_.start).min, txLines.map(_.end).max))
  }).toMap;
  
  var novelCountMap : GenMap[GenomicInterval,Int] = (Map[GenomicInterval,Int]()).withDefault(k => 0);
  
  val flatExonArray : GenomicArrayOfSets[String] = anno_holder.flatExonArray;
  val flatGeneSet   : Set[String] = anno_holder.flatGeneSet;
  
  //val exonCountMap :     scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](flatExonArray.getValueSet);
  //val flatGeneCountMap : scala.collection.mutable.Map[String,Int] = qcGtfAnnotationBuilder.initializeCounter[String](flatGeneSet);
  
  var readInfoList = Vector[(String,Vector[String],Vector[String],Vector[String])]();
  
  val comboCounts = scala.collection.mutable.AnyRefMap[(Set[String],Set[String],Set[String]),Int]().withDefault( x => 0 );
  val txCounts = scala.collection.mutable.AnyRefMap[String,Int]().withDefault( x => 0 );
  val txGoodCounts = scala.collection.mutable.AnyRefMap[(String,Int),Int]().withDefault( x => 0 );
  
  val writer = openWriterSmart_viaGlobalParam(outfile + ".readFeatureCombo.txt");
  writer.write("readID\tstart\tend\texonList\tknownJctList\tnovelJctList\tinsertions\tdeletions\n");
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    val (readExonSet, readGeneSet) = qcJunctionCounts.getExonsAndGenesFromPair(r1,r2,flatExonArray, stranded, fr_secondStrand);
    val splices : Set[GenomicInterval] = qcJunctionCounts.getSpliceJunctionSet(r1,r2,stranded,fr_secondStrand);
    
    val knownSplices = splices.filter(knownSpliceMap.contains(_)).map(knownSpliceMap(_));
    val novelSpliceIv = splices.filter(! knownSpliceMap.contains(_));
    val novelSplices = novelSpliceIv.map( (iv :GenomicInterval) => iv.chromName + ":"+iv.start+"-"+iv.end+iv.strand );
    
    //readInfoList = readInfoList :+ ((r1.getReadName(),
    //                                 readExonSet.toVector.sorted, 
    //                                 knownSplices.toVector.sorted,
    //                                 novelSplices.toVector.sorted
    //                               ));
    
    val insertions = r1.getCigar().getCigarElements().asScala.filter(_.getOperator() == CigarOperator.INSERTION).map(_.getLength());
    val deletions  = r1.getCigar().getCigarElements().asScala.filter(_.getOperator() == CigarOperator.DELETION).map(_.getLength());
    
    val indelTotal = insertions.sum + deletions.sum;
    val featureSet = readExonSet ++ knownSplices;
    
    val matchesTX : Option[String] = 
      if(novelSplices.size == 0){
          txList.find{ (tx : String) => {
            val currFeatures = txFeatureMap(tx);
            ((currFeatures &~ featureSet) ++ (featureSet &~ currFeatures)).isEmpty
          }}
      } else {
        None;
      }
    val goodMatchesTx = matchesTX match {
      case Some(tx) =>{
        val (txStart,txEnd) = endpointMap(tx);
        val (start,end) = (r1.getAlignmentStart,r1.getAlignmentEnd);
        goodDistance.map((D : Int) => {
          (start - D <= txStart && start + D >= txStart && end - D <= txEnd && end + D >= txEnd);
        })
      }
      case None => {
        List(false,false,false,false,false,false,false,false,false,false)
      }
    }
    
    if(! matchesTX.isEmpty){
      txCounts(matchesTX.get) += 1;
      for(i <- Range(0,goodMatchesTx.length)){
        if(goodMatchesTx(i)){
          txGoodCounts((matchesTX.get,i)) += 1;
        }
      }
    }
    
    
    writer.write(r1.getReadName() + "\t"+
                 r1.getAlignmentStart() +"\t"+
                 r1.getAlignmentEnd() +"\t" +
                 indelTotal+"\t"+
                 novelSplices.size+"\t"+
                 matchesTX +"\t"+
                 readExonSet.toList.sorted.mkString(",") + "\t"+
                 knownSplices.toList.sorted.mkString(",") + "\t"+
                 novelSplices.toList.sorted.mkString(",") + "\t"+
                 insertions.mkString(",")+"\t"+
                 deletions.mkString(",")+
                 "\t"+"\n"
                 );
    
    comboCounts((readExonSet, knownSplices, novelSplices)) += 1;
    
  }

  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    close(writer);

    val writer2 = createOutputFile(outfile , "featureComboCt.txt","",docWriter,
             ("exonList","List[String]","[[TODO]]"),
             ("knownJctList","List[String]","[[TODO]]"),
             ("novelJctList","List[String]","[[TODO]]")
    );
    writer2.write("exonList\tknownJctList\tnovelJctList\tCT\n");
    for((exons,known,novel) <- comboCounts.keysIterator){
      val ct = comboCounts((exons,known,novel));
      writer2.write(exons.toList.sorted.mkString(",") + "\t"+
                   known.toList.sorted.mkString(",") + "\t"+
                   novel.toList.sorted.mkString(",") + "\t"+
                   ct + "\n"
      );
    }
    close(writer2);

    val writer3 = createOutputFile(outfile, "TXcount.txt","",docWriter,
             ("[[TODO]]","[[TODO]]","[[TODO]]")
    );
    writer3.write("TX.id\tCT\t"+goodDistance.map("CT_"+_).mkString("\t")+"\n");
    for(tx <- txList){
      writer3.write(
          tx +"\t"+
          txCounts(tx)+"\t"+
          goodDistance.map(i => txGoodCounts((tx,i))).mkString("\t")+"\n"
      );
    }
    close(writer3);
    
    val writer4 = createOutputFile(outfile, "TXinfo.txt","",docWriter,
             ("[[TODO]]","[[TODO]]","[[TODO]]")
    );
    writer4.write("TX.id\tstart\tend\tfeatureList\n");
    for(tx <- txList){
      writer4.write(
          tx +"\t"+
          endpointMap(tx)._1 +"\t" + endpointMap(tx)._2 +"\t"+
          txFeatureMap(tx).toVector.sorted.mkString(",")+
          "\n"
      );
    }
    close(writer4);
    
    //writer.write("CHROM	FWD_CT	REV_CT	CT\n");
    //summaryWriter.write("NumberOfChromosomesCovered	"+chroms.size+"\n");
  }
  
  def getUtilityName : String = "featureComboCt";

}