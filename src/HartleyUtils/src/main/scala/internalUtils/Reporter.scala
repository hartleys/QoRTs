package internalUtils

/*
 * This is a standard 
 * 
 * I use this in most of my scala programs.
 */

object Reporter {

  
  private abstract class ReportLogger(){
    
    
    def reportln(str : String, verb : String){
      if(isWorthy(verb)){
        freshenLine;
        report_basefunction(str + "\n");
      }
    } 
     
    def startReport(str : String, verb : String){
      if(isWorthy(verb)){
        freshenLine;
        report_basefunction(str);
        isMidline = str.last != '\n';
      }
    }
    
    /*def freshen(verb : String){
      if(isWorthy(verb)){
        freshenLine;
        isMidline = false;
      }
    }*/
    def report(str : String, verb : String){
      if(isWorthy(verb)){
        report_basefunction(str);
        isMidline = str.last != '\n';
      }
    }
    
    /*
     * Do not use the below functions:
     */
    
    final val verbosityNames = List("output","error","warn","report","note","progress","debug","deepdebug");
    var isMidline = false;
    def isWorthy(verb : String) : Boolean = {
      if(verbosityNames.indexOf(verb) == -1) true;
      else verbositySetting(verbosityNames.indexOf(verb))
    }
    def freshenLine() {
      if(isMidline) report_basefunction("\n");
      isMidline = false;
    }
    
    /*
     * These must be defined in the inheriting classes:
     */
    //def init(filename : String);
    val verbositySetting: Array[Boolean]
    //def open()
    def report_basefunction(str : String)
    def close()
  }
 // case class FileReportLogger(var verbositySetting: Array[Boolean])(file : fileUtils.WriterUtil) extends ReportLogger {
    
 // }
  
  private case class FileReportLogger(vset : Array[Boolean], outfile : String) extends ReportLogger {
    val verbositySetting = vset;
    
    val writer : fileUtils.WriterUtil = fileUtils.openWriter(outfile);
    
    def report_basefunction(str : String){
      writer.write(str);
      writer.flush;
    }
    def close() {
      fileUtils.close(writer);
    }
  }
  
  private case class ConsoleReportLogger(vset : Array[Boolean]) extends ReportLogger {
    val verbositySetting = vset;
    
    def report_basefunction(str : String){
      print(str);
    }
    def close() {
      //do nothing.
    }
  }
  
  private case class ErrConsoleReportLogger(vset : Array[Boolean]) extends ReportLogger {
    val verbositySetting = vset;
    
    def report_basefunction(str : String){
      Console.err.print(str);
    }
    def close() {
      //do nothing.
    }
  }
  
  private case class StringReportLogger(vset : Array[Boolean]) extends ReportLogger {
    val verbositySetting = vset;
    val sb : StringBuilder = new StringBuilder();
    
    def report_basefunction(str : String){
      sb.append(str);
    }
    def close() {
      //do nothing.
    }
    def getLogString() : String = return sb.toString();
  }
  
  private case class AddedFileReportLogger(logfile : String, srl : StringReportLogger) extends ReportLogger {
    val verbositySetting = srl.verbositySetting;
    val writer = fileUtils.openWriter(logfile);
    writer.write(srl.getLogString());
    
    def report_basefunction(str : String){
      writer.write(str);
      writer.flush();
    }
    def close() {
      writer.close();
      /*try {
        
        writer.write(srl.getLogString());
        
      } catch {
        case e : Exception => {
          //do nothing!
        }
      }*/
    }
  }
  /*
   * 
   *     } catch {
      case e : Exception => {
        internalUtils.Reporter.reportln("============================FATAL_ERROR============================\n"+
                                        "QoRTs encountered a FATAL ERROR. For general help, use command:\n"+
                                        "          java -jar path/to/jar/QoRTs.jar --help\n"+
                                        "============================FATAL_ERROR============================\n"+
                                        "Error info:","note");
        throw e;
      }
    }
   */
 
  private var loggers : List[ReportLogger] = List[ReportLogger]();

  
/*
 * USAGE METHODS:
 */
  
  //final val verbosityNames = List("output","error","warn","report","note","progress","debug","deepdebug");
  val DEFAULT_CONSOLE_VERBOSITY = Array(false,true,true,true,true,true,false,false);
  val QUIET_CONSOLE_VERBOSITY = Array(false,true,true,false,false,false,false,false);
  val VERBOSE_CONSOLE_VERBOSITY = Array(false,true,true,true,true,true,true,true);
  val OUTPUT_VERBOSITY = Array(true, false, false, false, false,false,false,false);
  
  val logVerbositySetting = Array(false,true,true,true,true,true,true,true);
  val debugLogVerbositySetting = Array(false,true,true,true,true,true,true,true);
  val warningLogVerbositySetting = Array(false,true,true,false,false,false,false,false);
  
  //default internal logger:
  private val internalLog : StringReportLogger = StringReportLogger(logVerbositySetting);
  private val warningLog : StringReportLogger = StringReportLogger(warningLogVerbositySetting);
  private val outputLog : ConsoleReportLogger = ConsoleReportLogger(OUTPUT_VERBOSITY);
  
  def getWarnings : String = warningLog.getLogString();
  /*
   * Initializers:
   */
  def init_full(logDir : String) {
    val logfile = logDir + "log.log";
    val debugLogfile = logDir + "debugLog.log";
    
    val fileLogger = FileReportLogger(logVerbositySetting, logfile);
    val debugFileLogger = FileReportLogger(debugLogVerbositySetting, debugLogfile);
    val consoleLogger = ConsoleReportLogger(DEFAULT_CONSOLE_VERBOSITY);
    
    loggers = fileLogger :: debugFileLogger :: consoleLogger :: loggers;
  }
  
  def init_simple(logDir : String) {
    val logfile = logDir + "log.log";
    //val debugLogfile = logDir + "debugLog.log";
    
    val fileLogger = FileReportLogger(logVerbositySetting, logfile);
    val consoleLogger = ConsoleReportLogger(DEFAULT_CONSOLE_VERBOSITY);
    
    loggers = fileLogger :: consoleLogger :: loggers;
  }
  
  def init_logfilefree {
    val consoleLogger = ConsoleReportLogger(DEFAULT_CONSOLE_VERBOSITY);
    loggers = consoleLogger :: loggers;
  }
  
  def init_stderrOnly(verbositySetting : Array[Boolean] = DEFAULT_CONSOLE_VERBOSITY) {
    val errLogger = ErrConsoleReportLogger(verbositySetting);
    loggers = errLogger :: loggers;
  }
  
  def init_completeLogFile(logfile : String) {
    val fileLogger = AddedFileReportLogger(logfile, internalLog)
    loggers = fileLogger :: loggers;
  }
  def init_warningLogFile(logfile : String) {
    val fileLogger = AddedFileReportLogger(logfile, warningLog)
    loggers = fileLogger :: loggers;
  }
  
  def init_base(){
    loggers = internalLog :: loggers;
    loggers = warningLog :: loggers;
    loggers = outputLog :: loggers;
    
  }
  
  /*
   * Reporting options:
   */
  
  var anyWarning : Boolean = false;
  
  def hasWarningOccurred() : Boolean = anyWarning;
  
  def reportln(str : String, verb : String) {
    if(verb == "warn"){
      anyWarning = true;
    }
    
    loggers.map((logger) => logger.reportln(str,verb))
  }
  
  def report(str : String, verb : String){
    if(verb == "warn"){
      anyWarning = true;
    }
    loggers.map((logger) => logger.report(str,verb))
  }
  
  def startReport(str : String, verb : String){
    if(verb == "warn"){
      anyWarning = true;
    }
    loggers.map((logger) => logger.startReport(str,verb))
  }
  
  def error(str : String){
    reportln("<====== FATAL ERROR! ======>","error");
    reportln("----------------------------","error");
    reportln("     Error message: \"" + str + "\"","error");
    reportln("     Stack Trace:","error");
    val stackTrace = Thread.currentThread.getStackTrace;
    stackTrace.map((ste) => reportln("        " + ste.toString, "error"))
    
    reportln("<==========================>","error");
    closeLogs;
    throw new Exception(str);
  }
  def error(e : Exception){
    reportln("<====== FATAL ERROR! ======>","error");
    reportln("----------------------------","error");
    reportln("     Exception message: \"" + e.toString + "\"","error");
    reportln("     Stack Trace:","error");
    val stackTrace = e.getStackTrace;
    stackTrace.map((ste) => reportln("        " + ste.toString, "error"))
    
    reportln("<==========================>","error");
    closeLogs;
    throw e;
  }
  
  def closeLogs() {
    loggers.map((logger) => logger.close);
    loggers = List();
  }
}
