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

object qcChromCount {

}


class qcChromCount(isSingleEnd : Boolean, fr_secondStrand : Boolean) extends QCUtility[String] {
  reportln("> Init chromCount Utility","debug");

  
  var chromMap : Map[(String,Char),Int] = Map[(String,Char),Int]().withDefault((x) => 0);
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : String = {
    val strand = getStrand(r1, true,fr_secondStrand);
    
    val chrom = r1.getReferenceName();
    val cs = (chrom,strand);
    
    chromMap = chromMap.updated(cs, chromMap(cs) + 1);
    return chrom + "(" + strand.toString + ")";
  }
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    val writer = createOutputFile(outfile, "chromCount.txt","",docWriter,
             ("CHROM","String","Chromosome name"),
             ("FWD_CT","int","Number of reads/read-pairs found on FWD genomic strand"),
             ("REV_CT","int","Number of reads/read-pairs found on FWD genomic strand"),
             ("CT","int","Number of reads/read-pairs found on EITHER genomic strand")
    );
    val chroms = chromMap.keySet.map(_._1).toVector.sorted;
    writer.write("CHROM	FWD_CT	REV_CT	CT\n");
    for(chrom <- chroms){
      val ctp = chromMap((chrom,'+'));
      val ctm = chromMap((chrom,'-'));
      writer.write(chrom +"	"+ ctp+ "	"+  ctm +"	"+(ctp+ctm)+"\n");
    }
    close(writer);
    
    summaryWriter.write("NumberOfChromosomesCovered\t"+chroms.size+"\tNumber of chromosomes with 1 or more aligned reads.\n");
    
  }
  
  def getUtilityName : String = "chromCounts";

}