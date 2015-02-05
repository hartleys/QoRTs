package internalUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;

//import scala.reflect.runtime.universe._


object commandLineUI {

  final val HELP_COMMAND_LIST : List[String] = List("help", "-help", "--help");
  final val MANUAL_COMMAND_LIST : List[String] = List("?","'?'", "\"?\"","man","-man","--man");
  
  var CLUI_CONSOLE_LINE_WIDTH = 68;
  
  final val ALWAYS_DEBUG_MODE = true;
  final val DEBUG_MODE_FLAG = "--debug";
  
  val DEFAULT_QUICK_SYNOPSIS = "<Synopsis not written>";
  val DEFAULT_SYNOPSIS = "<Synopsis not written>";
  val DEFAULT_DESCRIPTION = "<Description not written>";
  val DEFAULT_AUTHOR = List("Stephen W. Hartley, Ph.D. <stephen.hartley@nih.gov>");
  val DEFAULT_LEGAL =  
List(
     " This software is \"United States Government Work\" under the terms of the United States Copyright "+
     " Act.  It was written as part of the authors' official duties for the United States Government and "+
     " thus cannot be copyrighted.  This software is freely available to the public for use without a "+
     " copyright notice.  Restrictions cannot be placed on its present or future use.",

     " Although all reasonable efforts have been taken to ensure the accuracy and reliability of the "+
     " software and data, the National Human Genome Research Institute (NHGRI) and the U.S. Government "+
     " does not and cannot warrant the performance or results that may be obtained by using this software "+
     " or data.  NHGRI and the U.S. Government disclaims all warranties as to performance, merchantability "+
     " or fitness for any particular purpose.",

     " In any work or product derived from this material, proper attribution of the authors as the source "+
     " of the software or data should be made, using \"NHGRI Genome Technology Branch\" as the citation.",

     " NOTE: This package includes (internally) the sam-1.113.jar library from picard tools."+
     " That package uses the MIT license, which can be accessed using the command:",
     
     " java -jar thisjarfile.jar help samjdkinfo");
   
  //type Command = ((Array[String]) => Unit);

  //WARNING: STRIPQUOTES IS CURRENTLY NONFUNCTIONAL!
  
  abstract class CommandLineRunUtil {
    val parser : CommandLineArgParser;
    
    def run(args : Array[String]);
  }
  

  
  /*final val helpCommandList : Map[String, () => CommandLineRunUtil] = 
    Map(
           ("'?'" -> (() => new helpDocs())),
           ("\"?\"" -> (() => new helpDocs())),
           ("?" -> (() => new helpDocs())),
           ("help" -> (() => new helpDocs())),
           ("man" -> (() => new helpDocs())),
           ("--help" -> (() => new helpDocs())),
           ("--man" -> (() => new helpDocs())),
           ("-help" -> (() => new helpDocs())),
           ("-man" -> (() => new helpDocs()))
       );*/
  
  
  //
  //command: the name of the command used to execute the utility
  //quickSynopsis: a short (less than 48-char) description of what the utility does.
  //synopsis: DEPRECIATED. The synopsis is now generated automatically.
  //description: a multi-line or even multi-paragraph description of the utility.
  
  val CLUI_UNIVERSAL_ARGS : List[Argument[Any]] = 
                    new UnaryArgument( name = "verbose",
                                         arg = List("--verbose"), // name of value
                                         argDesc = "Flag to indicate that debugging information and extra progress information should be sent to stderr." // description
                                       ) ::
                    new UnaryArgument( name = "quiet",
                                         arg = List("--quiet","-s"), // name of value
                                         argDesc = "Flag to indicate that only errors and warnings should be sent to stderr." // description
                                       ) :: 
                                       List();
  
  class CommandLineArgParser(command : String, quickSynopsis : String, synopsis : String, description : String,  argList : List[Argument[Any]], authors : List[String] = DEFAULT_AUTHOR, legal : List[String] = DEFAULT_LEGAL) {

    //private var argMap = Map.empty[String,(Type, Any)];
    //private def registerArg[T](name : String, item : T)(implicit t : TypeTag[T]){
    //  argMap = argMap.updated(name, typeOf[T] -> item);
    //}
    
    //private var argMap = Map.empty[String,(Type, Any)];
    //private def registerArg[T](name : String, item : T)(implicit t : TypeTag[T]){
    //  argMap = argMap.updated(name, typeOf[T] -> item);
    //}
    
    //private 
    
    
    private def filterOutFinalArgs[T](arg : Argument[T]) : Boolean = {
      arg match {
        case a : FinalArgument[T] => true;
        case a : Any => false;
      }
    }
    val finalArgList = argList.filter(filterOutFinalArgs);
    val paramArgList = argList.filterNot(filterOutFinalArgs);
    
    def get[T](key : String) : T = {
      argList.find(arg => arg.getName() == key) match {
        case Some(arg) => {
          val s = arg.getValue();
          return s.asInstanceOf[T];
        }
        case None => {
          error("FATAL INTERNAL ERROR: IMPOSSIBLE STATE! parameter "+key +" does not exist!");
          return None.get;
        }
      }
    }
    //def get[T](key : String)(implicit m : TypeTag[T]): T = {
    //  val (t,s) = argMap(key)
    //  s.asInstanceOf[T];
    //  //getArgOption[T](key).get
    //}
    
    def getQuickSynopsis : String = quickSynopsis;
    def getDescription : String = description;
    def getSynopsis : String = description;
    
    //def getArgOption[T](key : String)(implicit m : TypeTag[T]): Option[T] = {
    //  argMap.get(key).flatMap{
    //    case (t, s) => if (t <:< typeOf[T]) Some(s.asInstanceOf[T]) else {error("FATAL INTERNAL ERROR: PARAMETER " + key + " OF WRONG TYPE! Found: " + t + " expected: " + typeOf[T]); None};
    //  }
    //}
    
    //case (t, s) => if (t <:< typeOf[T]) Some(s.asInstanceOf[T]) else {error("FATAL INTERNAL ERROR: PARAMETER " + key + " OF WRONG TYPE!"); None};
    def parseArguments(args : List[String], debugMode : Boolean = internalUtils.optionHolder.OPTION_debugMode) : Boolean = {
      try {
      if(args.length < 1){
        reportShortHelp();
        false;
      } else if(MANUAL_COMMAND_LIST.contains(args(0))){
        reportManual();
        false;
      } else if(HELP_COMMAND_LIST.contains(args(0))){
        reportShortHelp();
        false;
      } else {
        parseArgs_master(args.toList, debugMode);
        true;
      }
      } catch {
        case e : Exception => {
          reportln("Syntax Error? Syntax must be:","warn");
          reportShortHelp();
          reportln("For more information, use option --man","warn");
          reportln("Error is:","warn");
          throw e;
        }
      }
    }
    
    def reportManual(verb : String = "output") {
      report(getManual(),verb);
    }
    def getManual() : String = {
      val sb = new StringBuilder("");
      
      sb.append("NAME\n");
      sb.append("	" + command + "\n");
      sb.append("	Version: " + runner.runner.QORTS_VERSION + "\n");
      sb.append("\n");
      sb.append("USAGE\n"); 
      sb.append(wrapLinesWithIndent(getShortHelp(),CLUI_CONSOLE_LINE_WIDTH,"        ", false).substring(4)+"\n");
      sb.append("\n");
      sb.append("DESCRIPTION:\n");
      sb.append(wrapLinesWithIndent(description,CLUI_CONSOLE_LINE_WIDTH,"    ",false)+"\n");
      //lineseq2string(wrapLinesWithIndent(DESCRIPTION, internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH, "    ", false))
      sb.append("\n");
      sb.append("REQUIRED ARGUMENTS:\n");
      for(arg <- argList.filter(_.argIsMandatory)){
        sb.append(arg.getFullDescription+"\n");
      }
      sb.append("\n");
      
      sb.append("OPTIONS:\n");
      for(arg <- argList.filter(! _.argIsMandatory)){
        sb.append(arg.getFullDescription+"\n");
        sb.append("\n");
      }
      sb.append("AUTHORS:\n");
      sb.append(lineseq2string(wrapLinesWithIndent(authors,CLUI_CONSOLE_LINE_WIDTH,"    ",false)) + "\n");
      
      sb.append("LEGAL:\n");
      sb.append(lineseq2string(wrapLinesWithIndent(legal,CLUI_CONSOLE_LINE_WIDTH,"    ",false)) + "\n");
      
      return sb.toString;
    }
    
    private def sectionFormat(s : String) : String = {
      return(s);
    }
    private def sectionFormat(ss : Seq[String]) : String = {
      return(ss.foldLeft("")((soFar,s) => soFar + "\n" + s));
    }
    
    def reportShortHelp(verb : String = "output"){
      report(getShortHelp(),verb);
    } 
    def getShortHelp() : String = {
      "java [Java Options] -jar " + runner.runner.Runner_ThisProgramsExecutableJarFileName + " " + command + " [options] " +
         argList.foldLeft[String]("")((soFar : String, curr : Argument[Any]) => soFar +" "+ curr.getShortSyntax()) +"\n"
         //getForMoreHelp();
    }
    
    def getForMoreHelp() : String = {
      "For more info, use:\njava -jar " + runner.runner.Runner_ThisProgramsExecutableJarFileName + " --man " + command + ""
    }
    
    def getMandatoryArgumentsHelp() : String = {
      val sb = new StringBuilder("");
      for(arg <- argList){
        arg match {
          case (a : BinaryArgument[Any]) => {
            sb.append("");
          }
          case (a : FinalArgument[Any]) =>{
            
          }
        }
      }
      return sb.toString();
    }
    
    private def parseArgs_master(inputArguments : List[String], debugMode : Boolean = internalUtils.optionHolder.OPTION_debugMode){
      if(inputArguments.length < finalArgList.length){
         throwSyntaxErrorMessage("Not enough arguments: Require at least " + finalArgList.length + " arguments!\nRequired syntax is:\n" + getShortHelp());
      }
      
      val (inputParamArgs, inputFinalArgs) = inputArguments.splitAt( inputArguments.length - finalArgList.length );
      
      reportln("INPUT_COMMAND("+command+")","note");
      
      for((p,arg) <- inputFinalArgs.zip(finalArgList)){
        arg.parse(List(p));
        if(debugMode) reportln("  INPUT_ARG("+arg.getName()+")="+arg(),"note");
      }
      parseParamArgs(inputParamArgs, debugMode);
      
      //Now check to make sure all mandatory parameters are set:
      argList.find(! _.isReady()) match {
        case Some(unreadyArg) =>{
          throwSyntaxErrorMessage("Mandatory argument is not set: " + unreadyArg.getShortSyntax());
        } 
        case None =>{
          //do nothing
        }
      }
      //////for(a <- argList){
      //////  registerArg(a.getName,a.getValue());
      //////}
    }
    
    private def parseParamArgs(inputArguments : List[String], debugMode : Boolean){
      if(inputArguments.length != 0){
        //if(inputArguments.exists(isArgument(_))){
          //if(isArgument(inputArguments.head)){
            argList.find(_.isNamed(inputArguments.head)) match {
              case Some(arg) => {
                val remainder = arg.parse(inputArguments);
                if(debugMode) reportln("  INPUT_ARG("+arg.getName()+")="+arg(),"note");
                parseParamArgs(remainder, debugMode);
              }
              case None => {
                throwSyntaxErrorMessage("Unrecognized Argument: " + inputArguments.head);
              }
            }
         // } else {
          // throwSyntaxErrorMessage("Unexpected string (not recognized as an option name or argument): " + inputArguments.head);
         // }
        //} else {
        //  throwSyntaxErrorMessage("Unexpected/unrecognized commands/options/arguments: " + inputArguments);
        //}
      }
    }
  }
  
  private def throwSyntaxErrorMessage(s : String){
     error("SYNTAX ERROR! "+ s);
  }
  private def looksLikeArgument(arg : String) : Boolean = {
      if(arg.length == 0) false;
      else if(arg.charAt(0) == '-'){
        if(arg.length == 1) false;
        else true;
      } else false;
  }
  
  abstract class Argument[+T] {
    
    def getName() : String;
    def getValue() : T;
    //def getType()(m : TypeTag[T]) : reflect.runtime.universe.Type;
    def isNamed(argName : String) : Boolean;
    //def setValue(t : Any) : Unit;
    def parse(args : List[String]) : List[String];
    def isReady() : Boolean;
    
    def describe() : String;
    def getShortSyntax() : String;
    def getFullSyntax() : String;
    def argMasterType : String;
    def argSubType : String;
    def argType : String;
    
    def apply() : T = getValue();
    def getFullDescription() : String = {
      "    "+getFullSyntax()+"\n"+wrapLineWithIndent(describe(),CLUI_CONSOLE_LINE_WIDTH,8)+"\n        ("+argType+")";
    }
    
    def argIsMandatory : Boolean;
  }

  case class UnaryArgument(name : String, arg : List[String], argDesc : String, defaultValue : Boolean = false, isImportant : Boolean = false) extends Argument[Boolean] {
    def argIsMandatory = false;
    var value : Boolean = defaultValue;
    def getName = name;
    //def getType()(implicit m : TypeTag[Boolean]) = typeOf[Boolean];
    def describe() : String = {
      argDesc;
    }
    def getShortSyntax() : String = {
      if(isImportant) "["+arg(0)+"]";
      else "";
    }
    def getFullSyntax() : String = {
      arg(0);
    }
    def argMasterType() : String = "flag";
    def argSubType() : String = "";
    def argType() : String = "flag";
    
    def getValue() : Boolean = value;
    def isNamed(an : String) = {
      if(arg.exists(_ == an)) true;
      else {
        if(an.length > 1){
          val substr = an.substring(0,2);
          arg.exists(_ == substr);
        } else false;
      }
    }
    def setValue(t : Boolean) {
      value = t;
    }
    def isReady : Boolean = true;
    def parse(args : List[String]) : List[String] = {
      if(args.head.charAt(0) == '-'){
        if(args.head.charAt(1) == '-'){
          setValue(true);
          args.tail;
        } else {
          if(args.head.length == 2){
            setValue(true);
            args.tail;
          } else if(args.head.length < 2){
            throwSyntaxErrorMessage("?ErrorMessageNotAdded?");
            List();
          } else {            
            val newArg = "-" + args.head.substring(2);
            newArg :: args.tail;
          }
        }
      } else {
        throwSyntaxErrorMessage("?ErrorMessageNotAdded?");
        List();
      }
    }
  }
  case class BinaryOptionArgument[T](name : String, arg: List[String], valueName : String,  argDesc : String, defaultValue : Option[T] = None, isMandatory : Boolean = false, isImportant : Boolean = false, stripQuotes : Boolean = false)(implicit stringParser : StringParser[T]) extends Argument[Option[T]] {
    def argIsMandatory = isMandatory;
    def getName = name;
    var value : Option[T] = defaultValue;
    
    def getShortSyntax() : String = {
      if(isMandatory) "<"+arg(0)+" "+valueName+">";
      else if(isImportant) "["+arg(0)+" "+valueName+"]";
      else "";
    }
    def getFullSyntax() : String = {
      if(isMandatory) ""+arg(0)+" "+valueName+"";
      else ""+arg(0)+" "+valueName+"";
    }
    def argMasterType() : String = "monadic";
    def argSubType() : String = stringParser.argType;
    def argType() : String = argSubType();
    
    def describe() : String = {
      argDesc;
    }
    def getValue() : Option[T] = value;
    def isNamed(an : String) = arg.exists(_ == an);
    def setValue(t : T){
      value = Some(t);
    }
    def isReady : Boolean = ! ((value.isEmpty) && (isMandatory));
    def parse(args : List[String]) : List[String] = {
      if(args.length < 2) throwSyntaxErrorMessage("Variable " + arg(0) + " not set to anything!");
      
      val valueString = args.tail.head;
      setValue(stringParser.parse(valueString));
      args.tail.tail;
    }
  }
  case class BinaryArgument[T](name : String, arg: List[String], valueName : String,  argDesc : String, defaultValue : Option[T] = None, isMandatory : Boolean = false, isImportant : Boolean = false, stripQuotes : Boolean = false)(implicit stringParser : StringParser[T]) extends Argument[T] {
    def argIsMandatory = isMandatory;
    def getName = name;
    //def getType()(implicit m : TypeTag[T]) : reflect.runtime.universe.Type = typeOf[T];
    var isSet = false;
    var hasDefault = false;
    var value : T = defaultValue match {
      case Some(t) => {isSet = true; hasDefault = true; t;} 
      case None => stringParser.unsetValue;
    }
    
    def getShortSyntax() : String = {
      if(isMandatory) "<"+arg(0)+" "+valueName+">";
      else if(isImportant) "["+arg(0)+" "+valueName+"]";
      else "";
    }
    def getFullSyntax() : String = {
      ""+arg(0)+" "+valueName+"";
    }
    def argMasterType() : String = "monadic";
    def argSubType() : String = stringParser.argType;
    def argType() : String = argSubType();
    
    def describe() : String = {
      argDesc;
    }
    def getValue() : T = value;
    def isNamed(an : String) = arg.exists(_ == an);
    def setValue(t : T){
      isSet = true;
      value = t;
    }
    def isReady : Boolean = (isSet || hasDefault) && (isSet || (! isMandatory));
    def parse(args : List[String]) : List[String] = {
      if(args.length < 2) throwSyntaxErrorMessage("Variable " + arg(0) + " not set to anything!");
      
      val valueString = args.tail.head;
      setValue(stringParser.parse(valueString));
      args.tail.tail;
    }
  }
  case class ListArgument[T](name : String, arg: List[String], valueName : String, argDesc : String, defaultValue : Option[List[T]] = None, isMandatory : Boolean = false, isImportant : Boolean = false, stripQuotes : Boolean = false)(implicit stringParser : StringParser[T]) extends Argument[List[T]] {
    def argIsMandatory = isMandatory;
    def getName = name;
    var isSet = false;
    var hasDefault = false;
    var value : List[T] = defaultValue match {
      case Some(t) => {isSet = true; hasDefault = true; t;} 
      case None => List();
    }
    
    def getShortSyntax() : String = {
      if(isMandatory) "<"+arg(0)+" "+valueName+" ..."+">";
      else if(isImportant) "["+arg(0)+" "+valueName+" ..."+"]";
      else "";
    }
    def getFullSyntax() : String = {
      if(isMandatory) "<"+arg(0)+" "+valueName+" ..."+">";
      else "["+arg(0)+" "+valueName+" ..."+"]";
    }
    def argMasterType() : String = "WhiteSpaceDelimitedList";
    def argSubType() : String = stringParser.argType;
    def argType() : String = "WhiteSpaceDelimitedList of "+argSubType();
    
    def describe() : String = {
      argDesc;
    }
    def getValue() : List[T] = value;
    def setValue(t : List[T]) {value = t;}
    def addToValue(t : T) {
      value = value ++ List(t);
    }
    def isNamed(an : String) = arg.exists(_ == an);
    def isReady : Boolean = (isSet || hasDefault) && (isSet || (! isMandatory));
    
    def parse(args : List[String]) : List[String] = {
      if(args.length < 2) throwSyntaxErrorMessage("Variable " + arg(0) + " not set to anything!");
      
      val stringList = args.takeWhile(! looksLikeArgument(_));
      val theRestList = args.takeRight(args.length - stringList.length);
      
      val valueList = stringList.map(stringParser.parse(_));
      
      setValue(valueList);
      
      theRestList;
    }
  }
  
  case class FinalArgument[T](name : String, valueName : String,  argDesc : String, isImportant : Boolean = false, stripQuotes : Boolean = false)(implicit stringParser : StringParser[T]) extends Argument[T] {
    def argIsMandatory = isMandatory;
    def getName = name;
    val defaultValue : Option[T] = None;
    val isMandatory : Boolean = true;
    var isSet = false;
    var hasDefault = false;
    var value : T = stringParser.unsetValue;
    
    def getShortSyntax() : String = {
      "_"+valueName+"_";
    }
    def getFullSyntax() : String = {
      ""+valueName+"";
    }
    def argMasterType() : String = "trailingMonadic";
    def argSubType() : String = stringParser.argType;
    def argType() : String = argSubType();
    
    def describe() : String = argDesc;
    def getValue() : T = if(isSet) value; else{ error("Syntax error!"); value; }
    def setValue(t : T) { isSet = true; value = t;}
    def isNamed(an : String) : Boolean = ! isSet;
    def isReady : Boolean = isSet;
    def parse(args : List[String]) : List[String] = {
      val valueString = args.head;
      setValue(stringParser.parse(valueString));
      args.tail;
    }
  }
  
  
  abstract class StringParser[T]{
    def parse(s : String) : T;
    def argType : String;
    def unsetValue : T;
  }
  implicit object stringStringParser extends StringParser[String]{
    def parse(s : String) : String = s;
    def argType : String = "String";
    def unsetValue : String = "";
  }
  implicit object intStringParser extends StringParser[Int]{
    def parse(s : String) : Int = string2int(s);
    def argType : String = "Int";
    def unsetValue : Int = -1;
  }
  implicit object doubleStringParser extends StringParser[Double]{
    def parse(s : String) : Double = string2double(s);
    def argType : String = "Double";
    def unsetValue : Double = -1;
  }
  implicit object floatStringParser extends StringParser[Float]{
    def parse(s : String) : Float = string2float(s);
    def argType : String = "Float";
    def unsetValue : Float = -1;
  }
  implicit object commaListStringParser extends StringParser[List[String]]{
    def parse(s : String) : List[String] = s.split(",").toList;
    def argType : String = "CommaDelimitedListOfStrings";
    def unsetValue : List[String] = List();
  }
  implicit object commaListDoubleParser extends StringParser[List[Double]]{
    def parse(s : String) : List[Double] = s.split(",").map(string2double(_)).toList;
    def argType : String = "CommaDelimitedListOfDoubles";
    def unsetValue : List[Double] = List();
  }
  implicit object commaListIntParser extends StringParser[List[Int]]{
    def parse(s : String) : List[Int] = s.split(",").map(string2int(_)).toList;
    def argType : String = "CommaDelimitedListOfDoubles";
    def unsetValue : List[Int] = List();
  }
  implicit object commaListFloatParser extends StringParser[List[Float]]{
    def parse(s : String) : List[Float] = s.split(",").map(string2float(_)).toList;
    def argType : String = "CommaDelimitedListOfDoubles";
    def unsetValue : List[Float] = List();
  }  
  
}