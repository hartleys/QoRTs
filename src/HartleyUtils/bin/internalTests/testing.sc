package internalTests

import net.sf.samtools._

object testing {
  println("Welcome to the Scala worksheet")       //> Welcome to the Scala worksheet
  

  val x = scala.collection.mutable.AnyRefMap[internalUtils.commonSeqUtils.GenomicInterval,Int]().withDefault(k => 0)
                                                  //> x  : scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInterv
                                                  //| al,Int] = Map()
  
  val giv = internalUtils.commonSeqUtils.GenomicInterval("chrX",'+',10,100)
                                                  //> giv  : internalUtils.commonSeqUtils.GenomicInterval = GenomicInterval(chrX,+
                                                  //| ,10,100)
  
  x(giv) += 1
  
  x                                               //> res0: scala.collection.mutable.Map[internalUtils.commonSeqUtils.GenomicInter
                                                  //| val,Int] = Map(GenomicInterval(chrX,+,10,100) -> 1)
  
  x.keySet                                        //> res1: scala.collection.Set[internalUtils.commonSeqUtils.GenomicInterval] = S
                                                  //| et(GenomicInterval(chrX,+,10,100))
  
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
  
(10 until 10).toVector                            //> res2: Vector[Int] = Vector()

Range(10, 1, -1).toVector                         //> res3: Vector[Int] = Vector(10, 9, 8, 7, 6, 5, 4, 3, 2)
val op = CigarOperator.SOFT_CLIP                  //> op  : net.sf.samtools.CigarOperator = S

op.consumesReadBases()                            //> res4: Boolean = true
op.consumesReferenceBases()                       //> res5: Boolean = false

System.getProperty("sun.arch.data.model")         //> res6: String = 64
  
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
  
  new java.io.File( "." ).getCanonicalPath        //> res7: String = C:\eclipseScalav400
  
  "test\\_string"                                 //> res8: String("test\\_string") = test\_string
  
  "test_string".replaceAll("_","\\\\_")           //> res9: String = test\_string
  
  
  CigarOperator.SOFT_CLIP                         //> res10: net.sf.samtools.CigarOperator = S
  
  CigarOperator.SOFT_CLIP.consumesReadBases()     //> res11: Boolean = true
  
  CigarOperator.SOFT_CLIP.consumesReferenceBases()//> res12: Boolean = false
  
  
  
  val X = 10                                      //> X  : Int = 10
  val Y = 0                                       //> Y  : Int = 0
  
  X.toDouble / Y.toDouble                         //> res13: Double = Infinity
  
  
  
  
  
}