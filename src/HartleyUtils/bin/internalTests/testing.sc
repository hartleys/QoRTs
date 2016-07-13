package internalTests

import net.sf.samtools._

object testing {
  println("Welcome to the Scala worksheet")       //> Welcome to the Scala worksheet
  "Hello".substring(0,2)                          //> res0: String = He
  "Hello".substring(1)                            //> res1: String = ello
  
  scala.math.log(10.0)                            //> res2: Double = 2.302585092994046
  
  //val reader =  new SAMFileReader(new java.io.File("C:\\Users\\hartleys\\work\\nihwork\\home_copy\\projects\\ZZZ-ExampleDataset\\TestSets\\readLength50\\inputData\\fastq\\SAMP1_RG1.50bp.1.fq.gz"));
  
  "Hello+World".split("\\+")                      //> res3: Array[String] = Array(Hello, World)
  
  Range(0,10).toVector                            //> res4: Vector[Int] = Vector(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
  
  Vector("hi").tail                               //> res5: scala.collection.immutable.Vector[String] = Vector()

  scala.math.log(0).isInfinite                    //> res6: Boolean = true


  print(Vector("Hello","World","How","Are","You"));
                                                  //> Vector(Hello, World, How, Are, You)
  print(Some("Hi"))                               //> Some(Hi)
  def zeroPad(i : Int, cols : Int) : String = {
    val s = i.toString;
    return 0.toString * (cols - s.length) + s;
  }                                               //> zeroPad: (i: Int, cols: Int)String
   
   0.toDouble / 0.toDouble                        //> res7: Double = NaN
   
   
   
   0.toString * -2                                //> res8: String = ""
   
   zeroPad(3000,3)                                //> res9: String = 3000

  internalUtils.stdUtils.zipIteratorWithCount(List("A","B","C","D","E","F","G").iterator).filter(x => x._2 % 2 == 0).toList
                                                  //> res10: List[(String, Int)] = List((A,0), (C,2), (E,4), (G,6))

  val x = scala.collection.mutable.AnyRefMap[internalUtils.commonSeqUtils.GenomicInterval,Int]().withDefault(k => 0)
                                                  //> x  : scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInter
                                                  //| val,Int] = Map()
  
  val giv = internalUtils.commonSeqUtils.GenomicInterval("chrX",'+',10,100)
                                                  //> giv  : internalUtils.commonSeqUtils.GenomicInterval = GenomicInterval(chrX,
                                                  //| +,10,100)
  
  x(giv) += 1
  
  x                                               //> res11: scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInt
                                                  //| erval,Int] = Map(GenomicInterval(chrX,+,10,100) -> 1)
  
  x.keySet                                        //> res12: scala.collection.Set[internalUtils.commonSeqUtils.GenomicInterval] =
                                                  //|  Set(GenomicInterval(chrX,+,10,100))
  
  val line = "x	y"                                //> line  : String = x	y
  

  
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
  
  
(10 until 10).toVector                            //> res13: Vector[Int] = Vector()

Range(10, 1, -1).toVector                         //> res14: Vector[Int] = Vector(10, 9, 8, 7, 6, 5, 4, 3, 2)
val op = CigarOperator.SOFT_CLIP                  //> op  : net.sf.samtools.CigarOperator = S

op.consumesReadBases()                            //> res15: Boolean = true
op.consumesReferenceBases()                       //> res16: Boolean = false

System.getProperty("sun.arch.data.model")         //> res17: String = 64
  
  def escapeToMarkdown(s : String) : String = {
    escapifyString(s, List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"));
  }                                               //> escapeToMarkdown: (s: String)String
  def escapifyString(s : String, escapifyStrings : Seq[String], escapeString : String = "\\\\") : String = {
    escapifyStrings.foldLeft[String](s)((soFar, curr) => {
      soFar.replaceAll(curr,escapeString+curr);
    });
  }                                               //> escapifyString: (s: String, escapifyStrings: Seq[String], escapeString: Str
                                                  //| ing)String
  
  val testStr = "--stranded_fr_secondstrand"      //> testStr  : String = --stranded_fr_secondstrand
  
  val out = escapifyString(testStr,List("`","\\*","_","\\{","\\}","\\[","\\]","\\(","\\)","\\#","\\+","-","\\.","!"))
                                                  //> out  : String = \-\-stranded\_fr\_secondstrand
  
  val out2 = escapeToMarkdown(testStr)            //> out2  : String = \-\-stranded\_fr\_secondstrand
  
  //val writer = internalUtils.fileUtils.openWriter("test.txt")
  //writer.write(out2)
  //writer.close()
  
  new java.io.File( "." ).getCanonicalPath        //> res18: String = C:\eclipseScalav400
  
  "test\\_string"                                 //> res19: String("test\\_string") = test\_string
  
  "test_string".replaceAll("_","\\\\_")           //> res20: String = test\_string
  
  
  CigarOperator.SOFT_CLIP                         //> res21: net.sf.samtools.CigarOperator = S
  
  CigarOperator.SOFT_CLIP.consumesReadBases()     //> res22: Boolean = true
  
  CigarOperator.SOFT_CLIP.consumesReferenceBases()//> res23: Boolean = false
  
  
  
  val X = 10                                      //> X  : Int = 10
  val Y = 0                                       //> Y  : Int = 0
  
  X.toDouble / Y.toDouble                         //> res24: Double = Infinity
  
  
   Seq(1,2,3,4,5).filter(_ < 3)                   //> res25: Seq[Int] = List(1, 2)
  
  
  val splitstr = "Gene = \"HRG10;2\"; TX = \"Blah\"; exon = 1; variant = \"3\"; Gene = \"Say \\\"Hello\\\"\";"
                                                  //> splitstr  : String = Gene = "HRG10;2"; TX = "Blah"; exon = 1; variant = "3"
                                                  //| ; Gene = "Say \"Hello\"";
  
  internalUtils.stdUtils.splitRespectQuote(splitstr,";")
                                                  //> res26: Array[String] = Array(Gene = "HRG10;2", " TX = "Blah"", " exon = 1",
                                                  //|  " variant = "3"", " Gene = "Say \"Hello\""")

   
  internalUtils.stdUtils.parseTokens(splitstr,';').toList.foreach((s : String) => println("'"+s+"'"))
                                                  //> 'Gene = "HRG10;2"'
                                                  //| ' TX = "Blah"'
                                                  //| ' exon = 1'
                                                  //| ' variant = "3"'
                                                  //| ' Gene = "Say \"Hello\""'
  
  
  splitstr.split(";(?=([^\"]*\"[^\"]*\")*[^\"]*$)").toList.foreach((s : String) => println("'"+s+"'"))
                                                  //> 'Gene = "HRG10;2"'
                                                  //| ' TX = "Blah"'
                                                  //| ' exon = 1'
                                                  //| ' variant = "3"'
                                                  //| ' Gene = "Say \"Hello\""'

  
  
  
  
}