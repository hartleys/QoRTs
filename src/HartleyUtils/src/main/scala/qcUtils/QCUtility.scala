package qcUtils

import net.sf.samtools._
import internalUtils.fileUtils._;

abstract class QCUtility[+B <: Any] {
   def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : B;
   def writeOutput(outfile : String, summaryWriter : WriterUtil);
}

object QCUtility {
  def getBlankStringUtil : QCUtility[String] = {
    new blankStringQCUtility;
  }
  def getBlankUnitUtil : QCUtility[Unit] = {
    new blankUnitQCUtility;
  }
  
  //def getBlankUtil[T] : QCUtility[T] = {
    
  //}
  
  class blankStringQCUtility extends QCUtility[String] {
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : String = {
      //do nothing!
      ""
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil){
      //do nothing!
    }
  }
  
  class blankUnitQCUtility extends QCUtility[Unit] {
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
      //do nothing!
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil){
      //do nothing!
    }
  }
}