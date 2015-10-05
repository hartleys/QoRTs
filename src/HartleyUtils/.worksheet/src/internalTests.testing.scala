package internalTests

import net.sf.samtools._

object testing {;import org.scalaide.worksheet.runtime.library.WorksheetSupport._; def main(args: Array[String])=$execute{;$skip(109); 
  println("Welcome to the Scala worksheet");$skip(25); val res$0 = 
  "Hello".substring(0,2);System.out.println("""res0: String = """ + $show(res$0));$skip(334); 

  //val reader =  new SAMFileReader(new java.io.File("C:\\Users\\hartleys\\work\\nihwork\\home_copy\\projects\\ZZZ-ExampleDataset\\TestSets\\readLength50\\inputData\\fastq\\SAMP1_RG1.50bp.1.fq.gz"));
  
     

  def zeroPad(i : Int, cols : Int) : String = {
    val s = i.toString;
    return 0.toString * (cols - s.length) + s;
  };System.out.println("""zeroPad: (i: Int, cols: Int)String""");$skip(20); val res$1 = 

   0.toString * -2;System.out.println("""res1: String = """ + $show(res$1));$skip(23); val res$2 = 
   
   zeroPad(3000,3);System.out.println("""res2: String = """ + $show(res$2));$skip(125); val res$3 = 

  internalUtils.stdUtils.zipIteratorWithCount(List("A","B","C","D","E","F","G").iterator).filter(x => x._2 % 2 == 0).toList;System.out.println("""res3: List[(String, Int)] = """ + $show(res$3));$skip(118); 

  val x = scala.collection.mutable.AnyRefMap[internalUtils.commonSeqUtils.GenomicInterval,Int]().withDefault(k => 0);System.out.println("""x  : scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInterval,Int] = """ + $show(x ));$skip(79); 
  
  val giv = internalUtils.commonSeqUtils.GenomicInterval("chrX",'+',10,100);System.out.println("""giv  : internalUtils.commonSeqUtils.GenomicInterval = """ + $show(giv ));$skip(17); 
  
  x(giv) += 1;$skip(7); val res$4 = 
  
  x;System.out.println("""res4: scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInterval,Int] = """ + $show(res$4));$skip(14); val res$5 = 
  
  x.keySet;System.out.println("""res5: scala.collection.Set[internalUtils.commonSeqUtils.GenomicInterval] = """ + $show(res$5));$skip(22); 
  
  val line = "x	y";System.out.println("""line  : String = """ + $show(line ));$skip(436); val res$6 = 
  

  
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
  
  
(10 until 10).toVector;System.out.println("""res6: Vector[Int] = """ + $show(res$6));$skip(27); val res$7 = 

Range(10, 1, -1).toVector;System.out.println("""res7: Vector[Int] = """ + $show(res$7));$skip(33); 
val op = CigarOperator.SOFT_CLIP;System.out.println("""op  : net.sf.samtools.CigarOperator = """ + $show(op ));$skip(24); val res$8 = 

op.consumesReadBases();System.out.println("""res8: Boolean = """ + $show(res$8));$skip(28); val res$9 = 
op.consumesReferenceBases();System.out.println("""res9: Boolean = """ + $show(res$9));$skip(43); val res$10 = 

System.getProperty("sun.arch.data.model");System.out.println("""res10: String = """ + $show(res$10));$skip(161); 
  
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
  
  val out2 = escapeToMarkdown(testStr);System.out.println("""out2  : String = """ + $show(out2 ));$skip(155); val res$11 = 
  
  //val writer = internalUtils.fileUtils.openWriter("test.txt")
  //writer.write(out2)
  //writer.close()
  
  new java.io.File( "." ).getCanonicalPath;System.out.println("""res11: String = """ + $show(res$11));$skip(21); val res$12 = 
  
  "test\\_string";System.out.println("""res12: String("test\\_string") = """ + $show(res$12));$skip(43); val res$13 = 
  
  "test_string".replaceAll("_","\\\\_");System.out.println("""res13: String = """ + $show(res$13));$skip(32); val res$14 = 
  
  
  CigarOperator.SOFT_CLIP;System.out.println("""res14: net.sf.samtools.CigarOperator = """ + $show(res$14));$skip(49); val res$15 = 
  
  CigarOperator.SOFT_CLIP.consumesReadBases();System.out.println("""res15: Boolean = """ + $show(res$15));$skip(54); val res$16 = 
  
  CigarOperator.SOFT_CLIP.consumesReferenceBases();System.out.println("""res16: Boolean = """ + $show(res$16));$skip(22); 
  
  
  
  val X = 10;System.out.println("""X  : Int = """ + $show(X ));$skip(12); 
  val Y = 0;System.out.println("""Y  : Int = """ + $show(Y ));$skip(29); val res$17 = 
  
  X.toDouble / Y.toDouble;System.out.println("""res17: Double = """ + $show(res$17))}
  
  
  
  
  
}
