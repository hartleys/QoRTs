package internalTests

import net.sf.samtools._

object testing {
  println("Welcome to the Scala worksheet")       //> Welcome to the Scala worksheet
  
  
  def escapeToMarkdown(s : String) : String = {
    escapifyString(s, List("`","\\*","_","{","}","[","]","(",")","#","+","-",".","!"));
  }                                               //> escapeToMarkdown: (s: String)String
  def escapifyString(s : String, escapifyStrings : Seq[String], escapeString : String = "\\\\") : String = {
    escapifyStrings.foldLeft[String](s)((soFar, curr) => {
      soFar.replaceAll(curr,escapeString+curr);
    });
  }                                               //> escapifyString: (s: String, escapifyStrings: Seq[String], escapeString: Stri
                                                  //| ng)String
  
  val testStr = "Hello_My_Name_Is_\\Steve!"       //> testStr  : String = Hello_My_Name_Is_\Steve!
  
  escapifyString(testStr,List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"))
                                                  //> res0: String = Hello\_My\_Name\_Is\_\Steve\!
    
}