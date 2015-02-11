package internalTests

import net.sf.samtools._

object testing {;import org.scalaide.worksheet.runtime.library.WorksheetSupport._; def main(args: Array[String])=$execute{;$skip(109); 
  println("Welcome to the Scala worksheet");$skip(146); 
  
  
  def escapeToMarkdown(s : String) : String = {
    escapifyString(s, List("`","\\*","_","{","}","[","]","(",")","#","+","-",".","!"));
  };System.out.println("""escapeToMarkdown: (s: String)String""");$skip(228); 
  def escapifyString(s : String, escapifyStrings : Seq[String], escapeString : String = "\\\\") : String = {
    escapifyStrings.foldLeft[String](s)((soFar, curr) => {
      soFar.replaceAll(curr,escapeString+curr);
    });
  };System.out.println("""escapifyString: (s: String, escapifyStrings: Seq[String], escapeString: String)String""");$skip(47); 
  
  val testStr = "Hello_My_Name_Is_\\Steve!";System.out.println("""testStr  : String = """ + $show(testStr ));$skip(111); val res$0 = 
  
  escapifyString(testStr,List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"));System.out.println("""res0: String = """ + $show(res$0))}
    
}
