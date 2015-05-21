package internalTests

import net.sf.samtools._

object testing {;import org.scalaide.worksheet.runtime.library.WorksheetSupport._; def main(args: Array[String])=$execute{;$skip(109); 
  println("Welcome to the Scala worksheet");$skip(121); 
  

  val x = scala.collection.mutable.AnyRefMap[internalUtils.commonSeqUtils.GenomicInterval,Int]().withDefault(k => 0);System.out.println("""x  : scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInterval,Int] = """ + $show(x ));$skip(79); 
  
  val giv = internalUtils.commonSeqUtils.GenomicInterval("chrX",'+',10,100);System.out.println("""giv  : internalUtils.commonSeqUtils.GenomicInterval = """ + $show(giv ));$skip(17); 
  
  x(giv) += 1;$skip(7); val res$0 = 
  
  x;System.out.println("""res0: scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInterval,Int] = """ + $show(res$0));$skip(14); val res$1 = 
  
  x.keySet;System.out.println("""res1: scala.collection.Set[internalUtils.commonSeqUtils.GenomicInterval] = """ + $show(res$1));$skip(429); val res$2 = 
  
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
  
(10 until 10).toVector;System.out.println("""res2: Vector[Int] = """ + $show(res$2));$skip(27); val res$3 = 

Range(10, 1, -1).toVector;System.out.println("""res3: Vector[Int] = """ + $show(res$3));$skip(33); 
val op = CigarOperator.SOFT_CLIP;System.out.println("""op  : net.sf.samtools.CigarOperator = """ + $show(op ));$skip(24); val res$4 = 

op.consumesReadBases();System.out.println("""res4: Boolean = """ + $show(res$4));$skip(28); val res$5 = 
op.consumesReferenceBases();System.out.println("""res5: Boolean = """ + $show(res$5));$skip(43); val res$6 = 

System.getProperty("sun.arch.data.model");System.out.println("""res6: String = """ + $show(res$6));$skip(161); 
  
  def escapeToMarkdown(s : String) : String = {
    escapifyString(s, List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"));
  };System.out.println("""escapeToMarkdown: (s: String)String""");$skip(228); 
  def escapifyString(s : String, escapifyStrings : Seq[String], escapeString : String = "\\\\") : String = {
    escapifyStrings.foldLeft[String](s)((soFar, curr) => {
      soFar.replaceAll(curr,escapeString+curr);
    });
  };System.out.println("""escapifyString: (s: String, escapifyStrings: Seq[String], escapeString: String)String""");$skip(49); 
  
  val testStr = "--stranded_fr_secondstrand";System.out.println("""testStr  : String = """ + $show(testStr ));$skip(121); 
  
  val out = escapifyString(testStr,List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"));System.out.println("""out  : String = """ + $show(out ));$skip(42); 
  
  val out2 = escapeToMarkdown(testStr);System.out.println("""out2  : String = """ + $show(out2 ));$skip(155); val res$7 = 
  
  //val writer = internalUtils.fileUtils.openWriter("test.txt")
  //writer.write(out2)
  //writer.close()
  
  new java.io.File( "." ).getCanonicalPath;System.out.println("""res7: String = """ + $show(res$7));$skip(21); val res$8 = 
  
  "test\\_string";System.out.println("""res8: String("test\\_string") = """ + $show(res$8));$skip(43); val res$9 = 
  
  "test_string".replaceAll("_","\\\\_");System.out.println("""res9: String = """ + $show(res$9));$skip(32); val res$10 = 
  
  
  CigarOperator.SOFT_CLIP;System.out.println("""res10: net.sf.samtools.CigarOperator = """ + $show(res$10));$skip(49); val res$11 = 
  
  CigarOperator.SOFT_CLIP.consumesReadBases();System.out.println("""res11: Boolean = """ + $show(res$11));$skip(54); val res$12 = 
  
  CigarOperator.SOFT_CLIP.consumesReferenceBases();System.out.println("""res12: Boolean = """ + $show(res$12));$skip(22); 
  
  
  
  val X = 10;System.out.println("""X  : Int = """ + $show(X ));$skip(12); 
  val Y = 0;System.out.println("""Y  : Int = """ + $show(Y ));$skip(29); val res$13 = 
  
  X.toDouble / Y.toDouble;System.out.println("""res13: Double = """ + $show(res$13))}
  
  
  
  
  
}
