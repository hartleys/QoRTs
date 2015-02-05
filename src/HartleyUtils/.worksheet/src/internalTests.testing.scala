package internalTests

import net.sf.samtools._

object testing {;import org.scalaide.worksheet.runtime.library.WorksheetSupport._; def main(args: Array[String])=$execute{;$skip(109); 
  println("Welcome to the Scala worksheet");$skip(46); 
  
  
  val s = "Hello\nMy\nName is Steve\n";System.out.println("""s  : String = """ + $show(s ));$skip(74); ;

  val L1 = List("A","B","C","D","E","F","G","H","I","J","K","L","M","N");System.out.println("""L1  : List[String] = """ + $show(L1 ));$skip(27); 

  val iterL = L1.iterator;System.out.println("""iterL  : Iterator[String] = """ + $show(iterL ));$skip(34); 
 
  val giterL = iterL.grouped(3);System.out.println("""giterL  : internalTests.testing.iterL.GroupedIterator[String] = """ + $show(giterL ));$skip(17); val res$0 = 
  
  giterL.next;System.out.println("""res0: List[String] = """ + $show(res$0));$skip(37); 
  
  val t = CigarOperator.INSERTION;System.out.println("""t  : net.sf.samtools.CigarOperator = """ + $show(t ));$skip(55); 
  val t2 = CigarOperator.SOFT_CLIP.consumesReadBases();System.out.println("""t2  : Boolean = """ + $show(t2 ));$skip(60); 
  val t3 = CigarOperator.SOFT_CLIP.consumesReferenceBases();System.out.println("""t3  : Boolean = """ + $show(t3 ));$skip(24); 

  val t4 = Range(0,10);System.out.println("""t4  : scala.collection.immutable.Range = """ + $show(t4 ));$skip(33); 

  val i1 = Iterator[Int](1,2,3);System.out.println("""i1  : Iterator[Int] = """ + $show(i1 ));$skip(21); 

  val i11 = i1.next


  class TestClass() {
    val myArray = Array[Int](0,1,2,3,4,5);
    def apply(i : Int) : Int = {
      return myArray(i);
    }
    def update(i : Int, j : Int){
      myArray(i) = j;
    }
  };System.out.println("""i11  : Int = """ + $show(i11 ));$skip(228); 
  
  val tc = new TestClass();System.out.println("""tc  : internalTests.testing.TestClass = """ + $show(tc ));$skip(23); ;
  
  val tct1 = tc(1);System.out.println("""tct1  : Int = """ + $show(tct1 ));$skip(20); ;
  val tct3 = tc(2);System.out.println("""tct3  : Int = """ + $show(tct3 ));$skip(14); ;
  tc(2) += 1;$skip(21); ;
  
val tct2 = tc(2);System.out.println("""tct2  : Int = """ + $show(tct2 ));$skip(34); ;

 val str = "%1.3f".format(115.0);System.out.println("""str  : String = """ + $show(str ));$skip(97); 
 
 val tttsts = Seq('A'.toInt ,'B'.toInt, 'C'.toInt,'G'.toInt,'N'.toInt,'S'.toInt,'H'.toInt).max;System.out.println("""tttsts  : Int = """ + $show(tttsts ));$skip(18); 
 
  val x1 = 200;System.out.println("""x1  : Int = """ + $show(x1 ));$skip(29); ;
  val x1s = 200.toHexString;System.out.println("""x1s  : String = """ + $show(x1s ));$skip(38); ;
  val s1x = Integer.parseInt("c8",16);System.out.println("""s1x  : Int = """ + $show(s1x ));$skip(55); 
  val srsg = "sample.ID	size.factorx".substring(0,21);System.out.println("""srsg  : String = """ + $show(srsg ));$skip(139); ;
                                                  
                                                  
  val r1b = Vector((386666,386766));System.out.println("""r1b  : scala.collection.immutable.Vector[(Int, Int)] = """ + $show(r1b ));$skip(37); ;
  val r2b = Vector((386793,386893));System.out.println("""r2b  : scala.collection.immutable.Vector[(Int, Int)] = """ + $show(r2b ));$skip(38); ;
  
  val merged = (r1b ++ r2b).sorted;System.out.println("""merged  : scala.collection.immutable.Vector[(Int, Int)] = """ + $show(merged ));$skip(17); val res$1 = 
  
  merged.tail;System.out.println("""res1: scala.collection.immutable.Vector[(Int, Int)] = """ + $show(res$1));$skip(14); val res$2 = 
  merged.head;System.out.println("""res2: (Int, Int) = """ + $show(res$2));$skip(24); val res$3 = 
  
  "A".compareTo("B");System.out.println("""res3: Int = """ + $show(res$3));$skip(21); val res$4 = 
  "B".compareTo("A");System.out.println("""res4: Int = """ + $show(res$4));$skip(21); val res$5 = 
  "B".compareTo("B");System.out.println("""res5: Int = """ + $show(res$5));$skip(424); val res$6 = 
  
  merged.tail.foldLeft(Vector(merged.head))( (soFar,curr) =>{
      if(curr._1 <= soFar.last._2){
        print("soFar = " + soFar+"\n");
        print("curr = " + curr+"\n");
        print("soFar.last._2 = " + soFar.last._2+"\n");
        print("curr._1 = " + curr._1+"\n");
        soFar.updated(soFar.length - 1, (soFar.last._1, math.max(curr._2, soFar.last._2)));
      } else {
        soFar :+ curr;
      }
    });System.out.println("""res6: scala.collection.immutable.Vector[(Int, Int)] = """ + $show(res$6))}
    

    
}
