package internalTests

import net.sf.samtools._

object testing {
  println("Welcome to the Scala worksheet")       //> Welcome to the Scala worksheet
  "Hello".substring(0,2)                          //> res0: String = He

  //val reader =  new SAMFileReader(new java.io.File("C:\\Users\\hartleys\\work\\nihwork\\home_copy\\projects\\ZZZ-ExampleDataset\\TestSets\\readLength50\\inputData\\fastq\\SAMP1_RG1.50bp.1.fq.gz"));
  
  "Hello+World".split("\\+")                      //> res1: Array[String] = Array(Hello, World)
  

  Vector("hi").tail                               //> res2: scala.collection.immutable.Vector[String] = Vector()

  print(Vector("Hello","World","How","Are","You"));
                                                  //> Vector(Hello, World, How, Are, You)
  print(Some("Hi"))                               //> Some(Hi)
  def zeroPad(i : Int, cols : Int) : String = {
    val s = i.toString;
    return 0.toString * (cols - s.length) + s;
  }                                               //> zeroPad: (i: Int, cols: Int)String
   
   0.toDouble / 0.toDouble                        //> res3: Double = NaN
   
   
   
   0.toString * -2                                //> res4: String = ""
   
   zeroPad(3000,3)                                //> res5: String = 3000

  internalUtils.stdUtils.zipIteratorWithCount(List("A","B","C","D","E","F","G").iterator).filter(x => x._2 % 2 == 0).toList
                                                  //> res6: List[(String, Int)] = List((A,0), (C,2), (E,4), (G,6))

  val x = scala.collection.mutable.AnyRefMap[internalUtils.commonSeqUtils.GenomicInterval,Int]().withDefault(k => 0)
                                                  //> x  : scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInterv
                                                  //| al,Int] = Map()
  
  val giv = internalUtils.commonSeqUtils.GenomicInterval("chrX",'+',10,100)
                                                  //> giv  : internalUtils.commonSeqUtils.GenomicInterval = GenomicInterval(chrX,+
                                                  //| ,10,100)
  
  x(giv) += 1
  
  x                                               //> res7: scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInte
                                                  //| rval,Int] = Map(GenomicInterval(chrX,+,10,100) -> 1)
  
  x.keySet                                        //> res8: scala.collection.Set[internalUtils.commonSeqUtils.GenomicInterval] = 
                                                  //| Set(GenomicInterval(chrX,+,10,100))
  
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
  
  
(10 until 10).toVector                            //> res9: Vector[Int] = Vector()

Range(10, 1, -1).toVector                         //> res10: Vector[Int] = Vector(10, 9, 8, 7, 6, 5, 4, 3, 2)
val op = CigarOperator.SOFT_CLIP                  //> op  : net.sf.samtools.CigarOperator = S

op.consumesReadBases()                            //> res11: Boolean = true
op.consumesReferenceBases()                       //> res12: Boolean = false

System.getProperty("sun.arch.data.model")         //> res13: String = 64
  
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
  
  new java.io.File( "." ).getCanonicalPath        //> res14: String = C:\eclipseScalav400
  
  "test\\_string"                                 //> res15: String("test\\_string") = test\_string
  
  "test_string".replaceAll("_","\\\\_")           //> res16: String = test\_string
  
  
  CigarOperator.SOFT_CLIP                         //> res17: net.sf.samtools.CigarOperator = S
  
  CigarOperator.SOFT_CLIP.consumesReadBases()     //> res18: Boolean = true
  
  CigarOperator.SOFT_CLIP.consumesReferenceBases()//> res19: Boolean = false
  
  
  
  val X = 10                                      //> X  : Int = 10
  val Y = 0                                       //> Y  : Int = 0
  
  X.toDouble / Y.toDouble                         //> res20: Double = Infinity
  
  
  
  
  
}