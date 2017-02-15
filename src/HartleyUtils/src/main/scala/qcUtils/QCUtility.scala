package qcUtils

import net.sf.samtools._
import internalUtils.fileUtils._;

abstract class QCUtility[+B <: Any] extends genUtility {
   def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : B;
   def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null);
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
  def getBlankIntUtil : QCUtility[Int] = {
    new blankIntQCUtility;
  }
  
  def getBlankUnitUtil : QCUtility[Unit] = {
    new blankUnitQCUtility;
  }
  def getBlankIntPairUtil : QCUtility[(Int,Int)] = new blankIntPairQCUtility;
  
  
  //def getBlankUtil[T] : QCUtility[T] = {
    
  //}
  /*
  class blankQCUtility extends QCUtility[Nothing] {
   def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : Nothing = {
     return Nothing;
   }
   def writeOutput(outfile : String, summaryWriter : WriterUtil){
     //do nothing!
   }
   def getUtilityName : String = "blankQCUtility";
  }*/
  
  
  
  class blankIntPairQCUtility extends QCUtility[(Int,Int)] {
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : (Int,Int) = {
      (-1, -1);
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
      //do nothing!
    }
    def getUtilityName : String = "blankIntPairQCUtility";
  }
   
  class blankIntQCUtility extends QCUtility[Int] {
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : Int = {
      -1;
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
      //do nothing!
    }
    def getUtilityName : String = "blankIntQCUtility";
  }
  
  class blankStringQCUtility extends QCUtility[String] {
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : String = {
      //do nothing!
      ""
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
      //do nothing!
    }
    def getUtilityName : String = "blankStringQCUtility";
  }
  
  class blankUnitQCUtility extends QCUtility[Unit] {
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
      //do nothing!
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
      //do nothing!
    }
    def getUtilityName : String = "blankUnitQCUtility";
  }
  
  class blankBooleanQCUtility extends QCUtility[Boolean] {
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int) : Boolean = {
      return true;
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
      //do nothing!
    }
    def getUtilityName : String = "blankUnitQCUtility";
  }
}