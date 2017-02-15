package internalUtils

import java.io.BufferedReader;
import java.io.BufferedInputStream;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.Writer;
import java.io.File;
import java.util.zip.GZIPOutputStream;
import java.util.zip.GZIPInputStream;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;

import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import scala.collection.JavaConversions._
import internalUtils.optionHolder._;
import internalUtils.Reporter._;


object fileUtils {
  
  /*
   * Bottom-level output commands:
   */
  
  type WriterUtil = Writer;
  
   class DocWriterUtil(writer : WriterUtil) {
     //close()
     //flush()
     //write(char[] cbuf, int off, int len)

     def flush(){
       writer.flush();
     }
     def write(cbuf : Array[Char],off : Int, len : Int){
       writer.write(cbuf,off,len);
     }
     
     var fileMap : Map[String,(String,List[(String,String,String)])] = Map[String,(String,List[(String,String,String)])]();
     
     def registerFile(fileSuffix : String, fileSummary : String, tableLines : (String,String,String)*){
       fileMap = fileMap + ((fileSuffix -> (fileSummary,tableLines.toList)));
     }
     
     def close(){
       fileMap.keys.toList.sorted(ord = internalUtils.stdUtils.AlphabetOrdering).foreach{case (fileSuffix : String) => {
         val (fileSummary,tableLines) = fileMap(fileSuffix); 
         writer.write("\nFILE\t"+fileSuffix+"\t"+fileSummary+"\n");
         tableLines.foreach{case (a,b,c) => {
           writer.write(a+"\t"+b+"\t"+c+"\n");
         }}
       }}
       
       writer.close();
     }
   }
  
  def openWriter(filename : String) : WriterUtil = {
    new BufferedWriter(new FileWriter(filename));
  }
  
  def openGzipWriter(filename : String) : WriterUtil = {
    new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(filename)))));
  }
  def openZipWriter(filename : String) : WriterUtil = {
    new BufferedWriter(new OutputStreamWriter(new ZipOutputStream(new FileOutputStream(new File(filename)))));
  }
  
  class SummaryWriter(summaryFile : String, includeDesc : Boolean) {
    val writer : WriterUtil = openWriter(summaryFile);
    if(includeDesc){
      writer.write("FIELD\tCOUNT\tDESC\n");
    } else {
      writer.write("FIELD\tCOUNT\n");
    }
    def write(field : String, count : String,desc : String = ""){
      if(includeDesc){
        writer.write(field+"\t"+count+"\t"+desc+"\n");
      } else {
        writer.write(field+"\t"+count+"\n");
      }
    }
    def close(){
      writer.close();
    }
  }
  
  def write2Tab(writer : internalUtils.fileUtils.WriterUtil, tableLines : (String,String)*){
     if(writer == null){
       //do nothing
     } else {
       for((a,b) <- tableLines){
         writer.write(a+"\t"+b+"\n");
       }
     }
  }
  
  def write3Tab(writer : internalUtils.fileUtils.WriterUtil, tableLines : (String,String,String)*){
     if(writer == null){
       //do nothing
     } else {
       for((a,b,c) <- tableLines){
         writer.write(a+"\t"+b+"\t"+c+"\n");
       }
     }
  }
  
  def openWriterSmart(filename : String, allowStdout : Boolean = false) : WriterUtil = {
    if(allowStdout && filename == "-"){
      return new BufferedWriter(new OutputStreamWriter(System.out));
    }
    
    if(filename.length >= 3){
      if(filename.takeRight(3) == ".gz"){
        return openGzipWriter(filename);
      }
      if(filename.length >= 4){
        if(filename.takeRight(4) == ".zip"){
          return openZipWriter(filename);
        }
      }
    }
    return openWriter(filename);
  } 
  def openWriterSmart_viaGlobalParam(filename : String) : WriterUtil = {
    if(! internalUtils.optionHolder.OPTION_noGzipOutput){
      return openGzipWriter(filename + ".gz");
    } else {
      return openWriter(filename);
    }
    
    //getGlobalParam[Boolean]("noGzipOutput") match {
    //  case Some(b) =>{
    //    if(! b){
     //     return openGzipWriter(filename + ".gz");
    //    }
    //  }
    //  case None =>{
    //    Reporter.report("warning: globalParam noGzipOutput not found!","warn");
     // }
    //}
    //return openWriter(filename);
  }
  
  
  def writeln(str : String, writer: WriterUtil){
    writer.write(str + "\n");
  }
  def write(str : String, writer : WriterUtil){
    writer.write(str);
  }
  def close(writer : WriterUtil){
    writer.close;
  } 
  
  /*
   * Advanced output commands:
   */
  
  def skipCommentedInitialLines(iter : Iterator[String], commentString : String = "#") : Iterator[String] = {
    
    while(iter.hasNext){
      val next = iter.next;
      if(! next.startsWith(commentString)){
        return( Iterator[String](next) ++ iter );
      }
      //else keep going.
    }
    error("ERROR opening file! No un-commented lines!")
    
    return(null);
    
  }
  
  
  def prepList(l : scala.collection.GenTraversableOnce[String], delim : String = "	") : String = {
    l.mkString(delim);
    //val list = l.toList;
    //if(list.length == 0){
    //  ""
    //} 
    //else {
    //  l.toList.tail.foldLeft(l.toList.head)((sofar,curr) => sofar + delim + curr)
    //}
  }
  
  //def getWriter(outfile: String) : WriterUtil = {
   //  new BufferedWriter(new FileWriter(outfile));
  //}
  
  def writeListLine(l : scala.collection.GenTraversableOnce[String], delim : String = "	", writer : WriterUtil) = {
    writeln(prepList(l,delim), writer);
  }
  
  def fileExists(filename : String) : Boolean = {
    java.nio.file.Files.exists(java.nio.file.Paths.get(filename));
  }
  
  /***
   *** Basic input commands:
   ***/
  
  def getLinesSmartSearchUnzip(infile : String) : Iterator[String] = {
    val plainFile = new File(infile);
    val gzFile = new File(infile + ".gz");
    val zipFile = new File(infile + ".zip");
    
    if(gzFile.exists){
      return scala.io.Source.fromInputStream(new GZIPInputStream(new BufferedInputStream(new FileInputStream(gzFile)))).getLines;
    } else if(zipFile.exists){
      return scala.io.Source.fromInputStream(new ZipInputStream(new BufferedInputStream(new FileInputStream(zipFile)))).getLines;
    } else if(plainFile.exists){
      return getLines(infile);
    } else {
      error("File: "+infile +" does not exist. Checked for .gz and .zip variants.");
      return null;
    }
  }
  
  def getLinesSmartUnzipFromSeq(infiles : Seq[String], allowStdin : Boolean = false) : Iterator[String] = {
    return infiles.tail.foldLeft(getLinesSmartUnzip(infiles.head,allowStdin))((soFar,inf) => {
      soFar ++ getLinesSmartUnzip(inf,false);
    });
  }
  
  def getLinesSmartUnzip(infile : String, allowStdin : Boolean = false) : Iterator[String] = {
    if(allowStdin && infile == "-"){
      return io.Source.stdin.getLines;
    }
    
    if(infile.length > 3 && infile.takeRight(3) == ".gz"){
        return scala.io.Source.fromInputStream(new GZIPInputStream(new BufferedInputStream(new FileInputStream(infile)))).getLines;
    }
    if(infile.length > 4 && infile.takeRight(4) == ".zip"){
        return scala.io.Source.fromInputStream(new ZipInputStream(new BufferedInputStream(new FileInputStream(infile)))).getLines;
    }
    return getLines(infile);
  }
  
  def getLines(infile: String) : Iterator[String] = {
    scala.io.Source.fromFile(new File(infile)).getLines();
  }
  
  def getCellIterator(infile : String, delim : String, smartUnzip : Boolean = false) : Iterator[Array[String]] = {
    val lines = if(smartUnzip) getLinesSmartUnzip(infile) else getLines(infile);
    return new Iterator[Array[String]]{
      def hasNext : Boolean = lines.hasNext;
      def next : Array[String] = {
        lines.next.split(delim);
      }
    }
  }
  
  /*
   * Advanced input commands:
   */
  
  /*DEPRECIATED
   * def getCellsFromFileCols(infile : String, delim : String) : List[Array[String]] = {
    val lines = getLines(infile);
    (for(line <- lines) yield {
      line.split(delim); 
    }).toList
  }*/
  
  /*
   * DEPRECIATED
   def getMatrixFromFileCols(infile : String, cols : List[Int], delim : String) : List[List[String]] = {
    val lines = getLines(infile);
    (for(line <- lines) yield {
      val cells = line.split(delim)
      (for(i <- cols) yield {
        cells(i);
      }).toList      
    }).toList
  }*/
  
  /*
   * This just gives you a list of the values found in the given column of a file.
   
  def getListFromFileCols(infile : String, col : Int, delim : String) : List[String] = {
    val lines = getLines(infile);
    (for(line <- lines) yield {line.split(delim)(col)}).toList
  }*/
  
  /*
   * If you want a map for a file, in which the key is taken from a specific column in the file
   * and the value is taken from another column in the file.
   * Note: This does not deal with collisions. Make sure the keys are unique!
   */
  def getMapFromFileCols(infile : String, keyCol : Int, valCol : Int, delim : String) : Map[String,String] = {
    val lines = getLines(infile);
    val pairs : Seq[(String,String)] = 
      (for(line <- lines) yield {
        val cells = line.split(delim)
        (cells(keyCol), cells(valCol))
      }).toSeq
    Map() ++ pairs;
  }
  
  /*
   * If you want a map for a file, in which the key is taken from a specific column in the file
   * and the value for each key is the entire line from the file.
   * Note: This does not deal with collisions. Make sure the keys are unique!
   */
  def getLineMapFromFileCols(infile : String, keyCol : Int, delim : String) : Map[String,String] = {
    val lines = getLines(infile);
    val pairs : Seq[(String,String)] = 
      (for(line <- lines) yield {
        val cells = line.split(delim)
        (cells(keyCol), line)
      }).toSeq
    Map() ++ pairs;
  }
  

  
  
  
}