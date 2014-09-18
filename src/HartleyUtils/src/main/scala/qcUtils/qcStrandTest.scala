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

object qcStrandTest {

  def hasReadFeatures(r : SAMRecord, stranded : Boolean, fr_secondStrand : Boolean, featureArray : GenomicArrayOfSets[String]) : Boolean = {
    val readIntervals : Iterator[GenomicInterval] = getGenomicIntervalsFromRead(r , stranded , fr_secondStrand);
    return readIntervals.exists( (iv) => {
      featureArray.findIntersectingSteps(iv).exists( ! _._2.isEmpty );
    });
  }
  
  def strandTestPair(r1 : SAMRecord, r2 : SAMRecord, featureArray : GenomicArrayOfSets[String]) : (Int, Int) = {
    // 0 = no feature.
    // 1 = feature on read strand
    // 2 = feature on opposite strand
    // 3 = features on both
    
    val r1A = hasReadFeatures(r1,true,true,  featureArray);
    val r1B = hasReadFeatures(r1,true,false, featureArray);
    val r2A = hasReadFeatures(r2,true,true,  featureArray);
    val r2B = hasReadFeatures(r2,true,false, featureArray);
    
    val r1t = if(r1A && r1B) 3 else if(r1A) 1 else if(r1B) 2 else 0;
    val r2t = if(r2A && r2B) 3 else if(r2A) 2 else if(r2B) 1 else 0;
    
    return ((r1t, r2t));
  }
  
}

class qcStrandTest(anno_holder : qcGtfAnnotationBuilder, stranded : Boolean, fr_secondStrand_bool : Boolean) extends QCUtility[Unit] {
  reportln("Init StrandCheck","progress");

  
  
  val strandedGeneArray = anno_holder.strandedGeneArray;
  
  var fr_firstStrand = 0;
  var fr_secondStrand = 0;
  var ambig_featuresFoundOnBothStrands = 0;
  var ambig_noFeature = 0;
  var ambig_other = 0;
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    val (r1t, r2t) = qcStrandTest.strandTestPair(r1,r2,strandedGeneArray);
    
    if( r1t == 1 && r2t == 2){
      fr_secondStrand += 1;
    } else if(r1t == 2 && r2t == 1){
      fr_firstStrand += 1;
    } else if(r1t == 3 || r2t == 3){
      ambig_featuresFoundOnBothStrands += 1;
    } else if(r1t == 0 || r2t == 0){
      ambig_noFeature += 1;
    } else {
      ambig_other += 1;
    }
  }
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    val total = fr_firstStrand + fr_secondStrand;
    val ratio = fr_firstStrand.toDouble / total.toDouble;

    val inferredStrandType = if( ratio < 0.1 ) "fr_secondStrand";
    else if(ratio > 0.9) "fr_firstStrand";
    else if(ratio < 0.6 && ratio > 0.4) "fr_unstranded";
    else "UNKNOWN_STRANDEDNESS";
    
    val strandedness_ok = if(inferredStrandType == "fr_secondStrand" && stranded && fr_secondStrand_bool) 1;
                          else if(inferredStrandType == "fr_firstStrand" && stranded && ! fr_secondStrand_bool) 1;
                          else if(inferredStrandType == "fr_unstranded" && (! stranded)) 1;
                          else 0;
    
    summaryWriter.write("StrandTest_frFirstStrand	"+fr_firstStrand+"\n");
    summaryWriter.write("StrandTest_frSecondStrand	"+fr_secondStrand+"\n");
    summaryWriter.write("StrandTest_ambig_genesFountOnBothStrands	"+ambig_featuresFoundOnBothStrands+"\n");
    summaryWriter.write("StrandTest_ambig_noGenes	"+ambig_noFeature+"\n");
    summaryWriter.write("StrandTest_ambig_other	"+ambig_other+"\n");
    summaryWriter.write("StrandTest_STRANDEDNESS_MATCHES_INFERRED	" + strandedness_ok +"\n")
  }
}













