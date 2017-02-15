package internalTests

import net.sf.samtools._

object testing {
  println("Welcome to the Scala worksheet")
  "Hello".substring(0,2)
  "Hello".substring(1)
    
  scala.math.log(10.0)
  
  //val reader =  new SAMFileReader(new java.io.File("C:\\Users\\hartleys\\work\\nihwork\\home_copy\\projects\\ZZZ-ExampleDataset\\TestSets\\readLength50\\inputData\\fastq\\SAMP1_RG1.50bp.1.fq.gz"));
  
  "Hello+World".split("\\+")
  
  Range(0,10).toVector
  
  Vector("hi").tail

  scala.math.log(0).isInfinite


  print(Vector("Hello","World","How","Are","You"));
  print(Some("Hi"))
  def zeroPad(i : Int, cols : Int) : String = {
    val s = i.toString;
    return 0.toString * (cols - s.length) + s;
  }
   
   0.toDouble / 0.toDouble
   
   
   
   0.toString * -2
   
   zeroPad(3000,3)

  internalUtils.stdUtils.zipIteratorWithCount(List("A","B","C","D","E","F","G").iterator).filter(x => x._2 % 2 == 0).toList

  val x = scala.collection.mutable.AnyRefMap[internalUtils.commonSeqUtils.GenomicInterval,Int]().withDefault(k => 0)
  
  val giv = internalUtils.commonSeqUtils.GenomicInterval("chrX",'+',10,100)
  
  x(giv) += 1
  
  x
  
  x.keySet
  
  val line = "x	y"
  

  
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
  
  
(10 until 10).toVector

Range(10, 1, -1).toVector
val op = CigarOperator.SOFT_CLIP

op.consumesReadBases()
op.consumesReferenceBases()

System.getProperty("sun.arch.data.model")
  
  def escapeToMarkdown(s : String) : String = {
    escapifyString(s, List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"));
  }
  def escapifyString(s : String, escapifyStrings : Seq[String], escapeString : String = "\\\\") : String = {
    escapifyStrings.foldLeft[String](s)((soFar, curr) => {
      soFar.replaceAll(curr,escapeString+curr);
    });
  }
  
  val testStr = "--stranded_fr_secondstrand"
  
  val out = escapifyString(testStr,List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"))
  
  val out2 = escapeToMarkdown(testStr)
  
  //val writer = internalUtils.fileUtils.openWriter("test.txt")
  //writer.write(out2)
  //writer.close()
  
  new java.io.File( "." ).getCanonicalPath
  
  "test\\_string"
  
  "test_string".replaceAll("_","\\\\_")
  
  
  CigarOperator.SOFT_CLIP
  
  CigarOperator.SOFT_CLIP.consumesReadBases()
  
  CigarOperator.SOFT_CLIP.consumesReferenceBases()
  
  
  
  val X = 10
  val Y = 0
  
  X.toDouble / Y.toDouble
  
  
   Seq(1,2,3,4,5).filter(_ < 3)
  
  
  val splitstr = "Gene = \"HRG10;2\"; TX = \"Blah\"; exon = 1; variant = \"3\"; Gene = \"Say \\\"Hello\\\"\";"
  
  internalUtils.stdUtils.splitRespectQuote(splitstr,";")

   
  internalUtils.stdUtils.parseTokens(splitstr,';').toList.foreach((s : String) => println("'"+s+"'"))
  
  
  splitstr.split(";(?=([^\"]*\"[^\"]*\")*[^\"]*$)").toList.foreach((s : String) => println("'"+s+"'"))

  
  
  
  
}