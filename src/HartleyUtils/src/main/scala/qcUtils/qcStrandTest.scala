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
    
    return (r1t, r2t);
    
    //(strandTestSingleRead(r1,featureArray),strandTestSingleRead(r2,featureArray));
  }
  
  def strandTestSingleRead(r1 : SAMRecord, featureArray : GenomicArrayOfSets[String]) : Int = {
    // 0 = no feature.
    // 1 = feature on read strand
    // 2 = feature on opposite strand
    // 3 = features on both
    
    val r1A = hasReadFeatures(r1,true,true,  featureArray);
    val r1B = hasReadFeatures(r1,true,false, featureArray);
    
    val r1t = if(r1A && r1B) 3 else if(r1A) 1 else if(r1B) 2 else 0;
    
    return (r1t);
  }
  
}

class qcStrandTest(isSingleEnd : Boolean, anno_holder : qcGtfAnnotationBuilder, stranded : Boolean, fr_secondStrand_bool : Boolean) extends QCUtility[String] {
  reportln("> Init StrandCheck Utility","debug");

  val strandedGeneArray = anno_holder.strandedGeneArray;
  
  var fr_firstStrand = 0;
  var fr_secondStrand = 0;
  var ambig_featuresFoundOnBothStrands = 0;
  var ambig_noFeature = 0;
  var ambig_other = 0;
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : String = {
    if(! isSingleEnd){
      val (r1t, r2t) = qcStrandTest.strandTestPair(r1,r2,strandedGeneArray);
    
      if( r1t == 1 && r2t == 2){
        fr_secondStrand += 1;
        return("fr_secondStrand");
      } else if(r1t == 2 && r2t == 1){
        fr_firstStrand += 1;
        return("fr_firstStrand");
      } else if(r1t == 3 || r2t == 3){
        ambig_featuresFoundOnBothStrands += 1;
        return("ambig_featuresFoundOnBothStrands");
      } else if(r1t == 0 || r2t == 0){
        ambig_noFeature += 1;
        return("ambig_noFeature");
      } else {
        ambig_other += 1;
        return("ambig_other");
      }
    } else {
      val r1t = qcStrandTest.strandTestSingleRead(r1, strandedGeneArray);
      if(r1t == 1){
        fr_secondStrand += 1;
        return("fr_secondStrand");
      } else if(r1t == 2){
        fr_firstStrand += 1;
        return("fr_firstStrand");
      } else if(r1t == 3){
        ambig_featuresFoundOnBothStrands += 1;
        return("ambig_featuresFoundOnBothStrands");
      } else if(r1t == 0){
        ambig_noFeature += 1;
        return("ambig_noFeature");
      } else {
        ambig_other += 1;
        return("ambig_other");
      }
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
    
    
    //Check the strandedness. If it doesn't match the data, send a warning:
    if(stranded){
      if((! fr_secondStrand_bool) && inferredStrandType == "fr_secondStrand"){
        reportln("WARNING: The data appears to be follow the fr_secondStrand rule, but\n"+
                 "         QoRTs ran in fr_firstStrand mode!\n"+
                 "      Are you sure that the data follows the fr_firstStrand rule?\n"+
                 "      Under the fr_firstStrand rule, the first read aligns to the strand\n"+
                 "      opposite the original RNA transcript, and the second read aligns to\n"+
                 "      the same strand as the original RNA transcript.\n"+
                 "      If the data follows the fr_secondStrand rule, you should probably\n"+
                 "      re-run QoRTs with the \"--stranded\" and \"--stranded_fr_secondstrand\" options.","warn");
      } else if(fr_secondStrand_bool && inferredStrandType == "fr_firstStrand") {
        reportln("WARNING: The data appears to be follow the fr_firstStrand rule, but\n"+
                 "         QoRTs ran in fr_secondStrand mode!\n"+
                 "      Are you sure that the data follows the fr_secondStrand rule?\n"+
                 "      Under the fr_firstStrand rule, the first read aligns to the strand\n"+
                 "      opposite the original RNA transcript, and the second read aligns to\n"+
                 "      the same strand as the original RNA transcript.\n"+
                 "      If the data follows the fr_firstStrand rule, you should probably\n"+
                 "      re-run QoRTs WITHOUT the \"--stranded_fr_secondstrand\" option.","warn");
      } else if(inferredStrandType == "fr_unstranded"){
        reportln("WARNING: The data appears to be unstranded, but\n"+
                 "         QoRTs ran in stranded mode!\n"+
                 "      Are you sure that the data is stranded?\n"+
                 "      If not, you should probably\n"+
                 "      re-run QoRTs with the \"--stranded\" option.","warn");
      }
    } else {
      if(inferredStrandType == "fr_firstStrand"){
        reportln("WARNING: The data appears to be STRANDED, following the fr_firstStrand rule.\n"+
                 "         Are you sure this isn't stranded data? If it is stranded, then you should probably\n"+
                 "         re-run QoRTs with the \"--stranded\" option!","warn");
      } else if(inferredStrandType == "fr_secondStrand") {
        reportln("WARNING: The data appears to be STRANDED, following the fr_secondStrand rule.\n"+
                 "         Are you sure this isn't stranded data? If it is stranded, then you should probably\n"+
                 "         re-run QoRTs with the \"--stranded\" and \"--stranded_fr_secondstrand\" options!","warn");
      }
    }
    
    if(inferredStrandType == "UNKNOWN_STRANDEDNESS"){
        reportln("WARNING: QoRTs is unable to infer the strandedness from the data!\n"+
                 "         This isn't a problem per-se, since QoRTs requires that strandedness\n"+
                 "         mode be set manually. However, it might be indicative that something\n"+
                 "         is very wrong with your dataset and/or transcript annotation.","warn");
    }
    
    summaryWriter.write("StrandTest_frFirstStrand	"+fr_firstStrand+"\n");
    summaryWriter.write("StrandTest_frSecondStrand	"+fr_secondStrand+"\n");
    summaryWriter.write("StrandTest_ambig_genesFountOnBothStrands	"+ambig_featuresFoundOnBothStrands+"\n");
    summaryWriter.write("StrandTest_ambig_noGenes	"+ambig_noFeature+"\n");
    summaryWriter.write("StrandTest_ambig_other	"+ambig_other+"\n");
    summaryWriter.write("StrandTest_STRANDEDNESS_MATCHES_INFERRED	" + strandedness_ok +"\n")
  }
  
  def getUtilityName : String = "StrandCheck";
}













