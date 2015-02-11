package qcUtils

import net.sf.samtools._
import internalUtils.fileUtils._;

abstract class QCUtility[+B <: Any] {
   def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : B;
   def writeOutput(outfile : String, summaryWriter : WriterUtil);
   def getUtilityName : String;
   def priority : Int = 255;
}

object QCUtility {
  trait QCU_Ordering extends Ordering[QCUtility[Any]] {
    def compare(a : QCUtility[Any], b : QCUtility[Any]) = implicitly[Ordering[Int]].compare(a.priority,b.priority);
  }
  implicit object QCUtility extends QCU_Ordering;
  
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
    def getUtilityName : String = "blankStringQCUtility";
  }
  
  class blankUnitQCUtility extends QCUtility[Unit] {
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
      //do nothing!
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil){
      //do nothing!
    }
    def getUtilityName : String = "blankUnitQCUtility";
  }
}