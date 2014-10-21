package internalUtils

import internalUtils.Reporter._
//import scala.reflect.runtime.universe._
 

object optionHolder {
  //val QC_Options : scala.collection.mutable.Map[String,
  /*
  private var paramMap = Map.empty[String,(Type, Any)];
  
  def registerGlobalParam[T](name : String, item : T)(implicit t : TypeTag[T]){
    paramMap = paramMap.updated(name, typeOf[T] -> item);
  }
  def getGlobalParam[T](key : String)(implicit m : TypeTag[T]): Option[T] = {
    paramMap.get(key).flatMap{
      case (t, s) => if (t <:< typeOf[T]) Some(s.asInstanceOf[T]) else {error("FATAL INTERNAL ERROR: PARAMETER " + key + " OF WRONG TYPE!"); None};
    }
  }*/
  

  
  var OPTION_noGzipOutput = false;
  
  //val DEFAULT_debugMode = true;
  
  var OPTION_debugMode = true;
  
  
  

 /*private var paramMap = Map.empty[String,(Manifest[_], Any)];
  
  def registerParam[T](name : String, item : T)(implicit m : Manifest[T]){
    paramMap = paramMap.updated(name, m -> item);
  }
  def get[T](key : String)(implicit m : Manifest[T]): Option[T] = {
    paramMap.get(key).flatMap{
      case (om, s) => if (om <:< m) Some(s.asInstanceOf[T]) else None;
    }
  }*/
   
  /*
  private val stringParams : scala.collection.mutable.Map[String,String] = scala.collection.mutable.Map[String,String]();
  private val booleanParams : scala.collection.mutable.Map[String,Boolean] = scala.collection.mutable.Map[String,Boolean]();
  private val intParams : scala.collection.mutable.Map[String,Int] = scala.collection.mutable.Map[String,Int]();
  private val doubleParams : scala.collection.mutable.Map[String,Double] = scala.collection.mutable.Map[String,Double]();
  private val stringSeqParams : scala.collection.mutable.Map[String,Seq[String]] = scala.collection.mutable.Map[String,Seq[String]]();
  private val intSeqParams : scala.collection.mutable.Map[String,Seq[Int]] = scala.collection.mutable.Map[String,Seq[Int]]();
  private val doubleSeqParams : scala.collection.mutable.Map[String,Seq[Double]] = scala.collection.mutable.Map[String,Seq[Double]]();
  private val booleanSeqParams : scala.collection.mutable.Map[String,Seq[Boolean]] = scala.collection.mutable.Map[String,Seq[Boolean]]();


  //def getStringParam(name : String)
  def setParam(name : String, value : Any){
    val n = name;
    value match {
      case v : Seq[Double] => doubleSeqParams(n) = v;
      case v : Seq[Int] => intSeqParams(n) = v;
      case v : Seq[String] => stringSeqParams(n) = v;
      case v : Seq[Boolean] => booleanSeqParams(n) = v;
      case v : String => stringParams(n) = v;
      case v : Boolean => booleanParams(n) = v;
      case v : Double => doubleParams(n) = v;
      case v : Int => intParams(n) = v;
      case v : Any => {
        error("UNKNOWN PARAMETER TYPE: [" + n.toString + "] = [" + n.toString() +"]");
      }
    }
  }
  
  def getParam[T](n : String)(implicit m : Manifest[T]) : T = {
    
  }
  
  def getStringParam(n : String) : Option[String] = stringParams.get(n);
  def stringParam(n : String) : String = stringParams(n);
  
  def getBooleanParam(n : String) : Option[Boolean] = booleanParams.get(n);
  def booleanParam(n : String) : Boolean = booleanParams(n);
  
  def getIntParam(n : String) : Option[Int] = intParams.get(n);
  def intParam(n : String) : Int = intParams(n);
  
  def getDoubleParam(n : String) : Option[Double] = doubleParams.get(n);
  def doubleParam(n : String) : Double = doubleParams(n);
  
  //INCOMPLETE: finish for all possible types!
  
  
  
  
  */
  
  
    /*
  private val stringParams : scala.collection.mutable.Map[String,QC_StringParam] = scala.collection.mutable.Map[String,QC_StringParam]();
  private val booleanParams : scala.collection.mutable.Map[String,QC_BooleanParam] = scala.collection.mutable.Map[String,QC_BooleanParam]();
  private val intParams : scala.collection.mutable.Map[String,QC_IntParam] = scala.collection.mutable.Map[String,QC_IntParam]();
  private val doubleParams : scala.collection.mutable.Map[String,QC_DoubleParam] = scala.collection.mutable.Map[String,QC_DoubleParam]();
  private val stringSeqParams : scala.collection.mutable.Map[String,QC_StringSeqParam] = scala.collection.mutable.Map[String,QC_StringSeqParam]();
  private val intSeqParams : scala.collection.mutable.Map[String,QC_IntSeqParam] = scala.collection.mutable.Map[String,QC_IntSeqParam]();
  private val doubleSeqParams : scala.collection.mutable.Map[String,QC_DoubleSeqParam] = scala.collection.mutable.Map[String,QC_DoubleSeqParam]();
  private val booleanSeqParams : scala.collection.mutable.Map[String,QC_BooleanSeqParam] = scala.collection.mutable.Map[String,QC_BooleanSeqParam]();
  
  private abstract class QC_Param[B]{
    def name : String;
    def value : B;
  }
  private case class QC_StringParam(in_name : String, in_value : String) extends QC_Param[String]{
    def name = in_name;
    def value = in_value;
  }
  private case class QC_BooleanParam(in_name : String, in_value : Boolean) extends QC_Param[Boolean]{
    def name = in_name;
    def value = in_value;
  }
  private case class QC_IntParam(in_name : String, in_value : Int) extends QC_Param[Int]{
    def name = in_name;
    def value = in_value;
  }
  private case class QC_DoubleParam(in_name : String, in_value : Double) extends QC_Param[Double]{
    def name = in_name;
    def value = in_value;
  }
  private case class QC_StringSeqParam(in_name : String, in_value : Seq[String]) extends QC_Param[Seq[String]]{
    def name = in_name;
    def value = in_value;
  }
  private case class QC_IntSeqParam(in_name : String, in_value : Seq[Int]) extends QC_Param[Seq[Int]]{
    def name = in_name;
    def value = in_value;
  }
  private case class QC_DoubleSeqParam(in_name : String, in_value : Seq[Double]) extends QC_Param[Seq[Double]]{
    def name = in_name;
    def value = in_value;
  }
  private case class QC_BooleanSeqParam(in_name : String, in_value : Seq[Boolean]) extends QC_Param[Seq[Boolean]]{
    def name = in_name;
    def value = in_value;
  }
  */
  
} 