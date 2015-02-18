package internalTests

import net.sf.samtools._

object testing {
  println("Welcome to the Scala worksheet")       //> Welcome to the Scala worksheet
  
  /*
  def escapeToMarkdown(s : String) : String = {
    escapifyString(s, List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"));
  }
  def escapifyString(s : String, escapifyStrings : Seq[String], escapeString : String = "\\\\") : String = {
    escapifyStrings.foldLeft[String](s)((soFar, curr) => {
      soFar.replaceAll(curr,escapeString+curr);
    });
  }
  
  */

  
  def escapeToMarkdown(s : String) : String = {
    escapifyString(s, List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"));
  }                                               //> escapeToMarkdown: (s: String)String
  def escapifyString(s : String, escapifyStrings : Seq[String], escapeString : String = "\\\\") : String = {
    escapifyStrings.foldLeft[String](s)((soFar, curr) => {
      soFar.replaceAll(curr,escapeString+curr);
    });
  }                                               //> escapifyString: (s: String, escapifyStrings: Seq[String], escapeString: Stri
                                                  //| ng)String
  
  val testStr = "--stranded_fr_secondstrand"      //> testStr  : String = --stranded_fr_secondstrand
  
  val out = escapifyString(testStr,List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"))
                                                  //> out  : String = \-\-stranded\_fr\_secondstrand
  
  val out2 = escapeToMarkdown(testStr)            //> out2  : String = \-\-stranded\_fr\_secondstrand
  
  //val writer = internalUtils.fileUtils.openWriter("test.txt")
  //writer.write(out2)
  //writer.close()
  
  new java.io.File( "." ).getCanonicalPath        //> res0: String = C:\eclipseScalav211\eclipse
  
  "test\\_string"                                 //> res1: String("test\\_string") = test\_string
  
  "test_string".replaceAll("_","\\\\_")           //> res2: String = test\_string
  
  
  CigarOperator.SOFT_CLIP                         //> res3: net.sf.samtools.CigarOperator = S
  
  CigarOperator.SOFT_CLIP.consumesReadBases()     //> res4: Boolean = true
  
  CigarOperator.SOFT_CLIP.consumesReferenceBases()//> res5: Boolean = false
  
  
  
  
  
  
}