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
        isMidline = true;
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
        isMidline = true;
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
 
  private var loggers : List[ReportLogger] = List[ReportLogger]();
  
/*
 * USAGE METHODS:
 */
  
  val consoleVerbositySetting = Array(true,true,true,true,true,true,true,true);
  val logVerbositySetting = Array(true,true,true,true,true,true,true,true);
  val debugLogVerbositySetting = Array(true,true,true,true,true,true,true,true);
  
  /*
   * Initializers:
   */
  def init_full(logDir : String) {
    val logfile = logDir + "log.log";
    val debugLogfile = logDir + "debugLog.log";
    
    val fileLogger = FileReportLogger(logVerbositySetting, logfile);
    val debugFileLogger = FileReportLogger(debugLogVerbositySetting, debugLogfile);
    val consoleLogger = ConsoleReportLogger(consoleVerbositySetting);
    
    loggers = fileLogger :: debugFileLogger :: consoleLogger :: loggers;
  }
  
  def init_simple(logDir : String) {
    val logfile = logDir + "log.log";
    //val debugLogfile = logDir + "debugLog.log";
    
    val fileLogger = FileReportLogger(logVerbositySetting, logfile);
    val consoleLogger = ConsoleReportLogger(consoleVerbositySetting);
    
    loggers = fileLogger :: consoleLogger :: loggers;
  }
  
  def init_logfilefree {
    val consoleLogger = ConsoleReportLogger(consoleVerbositySetting);
    loggers = consoleLogger :: loggers;
  }
  
  def init_stderrOnly {
    val errLogger = ErrConsoleReportLogger(consoleVerbositySetting);
    loggers = errLogger :: loggers;
  }
  
  /*
   * Reporting options:
   */
  
  def reportln(str : String, verb : String) {
    loggers.map((logger) => logger.reportln(str,verb))
  }
  
  def report(str : String, verb : String){
    loggers.map((logger) => logger.report(str,verb))
  }
  
  def startReport(str : String, verb : String){
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
    loggers.map((logger) => logger.close)
  }
}
