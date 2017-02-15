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

object qcMinorUtils {

}

class qcMinorUtils(isSingleEnd : Boolean, stranded : Boolean, fr_secondStrand_bool : Boolean,readLen : Int,
                   writeChromList : Boolean = true,
                   writeReadLengthDistro : Boolean = true,
                   outputFastqFile : Option[String] = None,
                   outputBAMFile : Option[String] = None) extends QCUtility[Unit] {
  reportln("> Init "+getUtilityName+" Utility","debug");
  
  var chromSeq : Vector[String] = Vector[String]();
  val readLenDist : Array[Array[Int]] = Array.ofDim(2,readLen+1);

  val writeOutputFastq = ! outputFastqFile.isEmpty;
  val (fq1,fq2) : (WriterUtil,WriterUtil) = if(writeOutputFastq){
    if(isSingleEnd){
      ((openWriterSmart_viaGlobalParam(outputFastqFile.get + ".sequence.fq")),null)
    } else {
      ((openWriterSmart_viaGlobalParam(outputFastqFile.get + ".sequence.1.fq")),
       (openWriterSmart_viaGlobalParam(outputFastqFile.get + ".sequence.2.fq")))
    }
  } else {
    (null,null);
  }
  
  val writeOutputBAM = ! outputBAMFile.isEmpty;
  val bamWriter : WriterUtil = if(writeOutputBAM){
    openWriterSmart_viaGlobalParam(outputFastqFile.get + ".PF.bam")
  } else {
    null;
  }
  
  def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) {
    if(chromSeq.isEmpty) chromSeq = chromSeq :+ r1.getReferenceName()
    if(r1.getReferenceName() != chromSeq.last){
      chromSeq = chromSeq :+ r1.getReferenceName()
    }
    readLenDist(0)(r1.getReadLength()) += 1;
    
    if(! isSingleEnd){
      readLenDist(1)(r2.getReadLength()) += 1;
    }
    
    if(writeOutputFastq){
      fq1.write("@"+r1.getReadName()+"\n"+
                  r1.getReadString()+"\n"+
                  "+\n"+
                  r1.getBaseQualityString()+"\n");
      if(!isSingleEnd){
        fq2.write("@"+r2.getReadName()+"\n"+
                  r2.getReadString()+"\n"+
                  "+\n"+
                  r2.getBaseQualityString()+"\n");
      }
    }
    if(writeOutputBAM){
      if(r1.getAlignmentStart() < r2.getAlignmentStart()){
        bamWriter.write(r1.getSAMString()+"\n");
        bamWriter.write(r2.getSAMString()+"\n");
      } else {
        bamWriter.write(r2.getSAMString()+"\n");
        bamWriter.write(r1.getSAMString()+"\n");
      }
      bamWriter.flush();
    }
  }
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    if(writeOutputFastq){
      fq1.close();
      if(!isSingleEnd){
        fq2.close();
      }
    }
    
    if(writeChromList){
      val writer = createOutputFile(outfile , "orderedChromList.txt","",docWriter,
             ("COLUMN1","String","The names of the chromosomes found in the BAM file, in the order that they first appear.")
      );
      chromSeq.foreach(chr => writer.write(chr + "\n"));
      close(writer);
    }
    if(writeReadLengthDistro){
      val writer = createOutputFile(outfile , "readLenDist.txt","",docWriter,
             ("LEN","int","The read-length."),
             ("CT_R1","int","The number of read-1 reads that have the given read-length."),
             ("CT_R2","int","The number of read-2 reads that have the given read-length.")
      );
      writer.write("LEN\tCT_R1\tCT_R2\n");
      for(i <- 0 until readLen+1){
        writer.write(i+"\t"+readLenDist(0)(i)+"\t"+readLenDist(1)(i)+"\n");
      }
      close(writer);
    }
    
    summaryWriter.write("NumReadsAtMaxReadLength_R1\t" + readLenDist(0)(readLen) +"\tNumber of read-1 reads at the max length.\n");
    summaryWriter.write("NumReadsAtMaxReadLength_R2\t" + readLenDist(1)(readLen) +"\tNumber of read-2 reads at the max length.\n");
    summaryWriter.write("NumReadsTruncated_R1\t" + readLenDist(0).init.sum  +"\tNumber of read-1 reads that are not at the max length.\n");
    summaryWriter.write("NumReadsTruncated_R2\t" + readLenDist(1).init.sum +"\tNumber of read-2 reads that are not at the max length.\n");
    
    summaryWriter.write("NumReadsTruncated_25pct_R1\t" + readLenDist(0).slice(0,readLen/4).sum  +"\tNumber of read-1 reads that are truncated to less than 25% of the max length\n");
    summaryWriter.write("NumReadsTruncated_25pct_R2\t" + readLenDist(1).slice(0,readLen/4).sum +"\tNumber of read-2 reads that are truncated to less than 25% of the max length\n");
    summaryWriter.write("NumReadsTruncated_50pct_R1\t" + readLenDist(0).slice(0,readLen/2).sum  +"\tNumber of read-1 reads that are truncated to less than 50% of the max length\n");
    summaryWriter.write("NumReadsTruncated_50pct_R2\t" + readLenDist(1).slice(0,readLen/2).sum +"\tNumber of read-2 reads that are truncated to less than 50% of the max length\n");
    summaryWriter.write("NumReadsTruncated_75pct_R1\t" + readLenDist(0).slice(0,readLen * 3 / 4).sum  +"\tNumber of read-1 reads that are truncated to less than 75% of the max length\n");
    summaryWriter.write("NumReadsTruncated_75pct_R2\t" + readLenDist(1).slice(0,readLen * 3 / 4).sum +"\tNumber of read-2 reads that are truncated to less than 75% of the max length\n");
    
    //summaryWriter.write("StrandTest_STRANDEDNESS_MATCHES_INFERRED	" + strandedness_ok +"\n")
  }
  
  def getUtilityName : String = "MinorUtils";
}













