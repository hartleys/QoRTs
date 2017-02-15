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

object fqcStdUtils {

}

class fqcStdUtils(readLength : Int, isSingleEnd : Boolean) extends fqcUtility[Int] {
  init();
  
  def init(){
    reportln("> Init "+getUtilityName+" Utility","debug");
  }
  
  def runOnReadPair(r1 : FastqLine, r2 : FastqLine, readNum: Int) : Option[Int] = {
    
    return Some(-1);
  }
  def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
    
  }

  def getUtilityName : String = "FastqStdUtils";
}













