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

object qcChromCount {

}


class qcChromCount(fr_secondStrand : Boolean) extends QCUtility[Unit] {
  reportln("Init chromCount","progress");

  
  var chromMap : Map[(String,Char),Int] = Map[(String,Char),Int]().withDefault((x) => 0);
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
    
    val strand = if(fr_secondStrand != r1.getFirstOfPairFlag) {
      if(r1.getReadNegativeStrandFlag()) '+' else '-';
    } else {
      if(r1.getReadNegativeStrandFlag()) '-' else '+';
    }
    
    val chrom = r1.getReferenceName();
    val cs = (chrom,strand);
    
    chromMap = chromMap.updated(cs, chromMap(cs) + 1);
  }
  def writeOutput(outfile : String, summaryWriter : WriterUtil){
    val writer = openWriterSmart_viaGlobalParam(outfile + ".chromCount.txt");
    
    val chroms = chromMap.keySet.map(_._1).toVector.sorted;
    
    writer.write("CHROM	FWD_CT	REV_CT	CT\n");
    
    for(chrom <- chroms){
      val ctp = chromMap((chrom,'+'));
      val ctm = chromMap((chrom,'-'));
      writer.write(chrom +"	"+ ctp+ "	"+  ctm +"	"+(ctp+ctm)+"\n");
    }
    
    close(writer);
    summaryWriter.write("NumberOfChromosomesCovered	"+chroms.size+"\n");
  }
}