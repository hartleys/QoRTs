package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;

import java.io.FileInputStream;

import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
import java.io.OutputStream
import java.io.FileOutputStream
import java.io.InputStream
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.File
import scala.collection.JavaConversions._
import java.util.zip.GZIPOutputStream;
import java.io._;

object convertSoftToHardClipping {

  val cigarCol : Int = 5;
  val seqCol : Int = 9;
  val qualCol : Int = 10;
  
  class convertSoftToHardClipping_runner extends CommandLineRunUtil {
    val parser : CommandLineArgParser = 
      new CommandLineArgParser(
          command = "",
          quickSynopsis = "",
          synopsis = "",
          description = "",
          argList = 
            new FinalArgument[String](
                                         name = "infile",
                                         valueName = "infile",
                                         argDesc = "" // description
                                        ) ::
            new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "" // description
            ) :: List()
      );
    
    def run(args : Array[String]){
      convertSoftToHardClipping.run(parser.get[String]("infile"),
                          parser.get[String]("outfile"));
    }
  }
  
  def run(infile :String, outfile :String) {
    val in = if(infile == "-"){
      io.Source.stdin.getLines;
    } else {
      getLines(infile);
    }
    
    val out = if(outfile == "-"){
      new BufferedWriter(new OutputStreamWriter(System.out));
    } else {
      new BufferedWriter(new FileWriter(outfile));
    }
    
    for(line <- in){
      try{
      if(line.startsWith("@")){
        out.write(line+"\n");
      } else {
        val cells = line.split("\\s+");
        if(cells(cigarCol).length() == 1){
          out.write(line+"\n");
        } else {
          val cigops = cells(cigarCol).split("\\d+").tail;
          val cigcts = cells(cigarCol).split("[MIDNSHPX=]");
          
          if(cigops.length != cigcts.length){
            scala.Console.err.println("Error! cigops.length ("+cigops.length+") != cigcts.length ("+cigcts.length+")\nCIGOPS: "+cigops.foldLeft("")((b,a) => b+", \""+a+"\"")+"\nCIGCTS: "+cigcts.foldLeft("")((b,a) => b+", \""+a+"\"")+"\nFOR LINE: \n\""+line+"\"");
            error("");
          }
          
          if(cigops.head == "S" || cigops.last == "S"){
            if(cigops.head == "S"){
              val ct = string2int(cigcts.head)
              cells(seqCol) = cells(seqCol).substring(ct);
              cells(qualCol) = cells(qualCol).substring(ct);
              cigops(0) = "H";
              //cigcts(0) = "";
            }
            if(cigops.last == "S"){
              val ct = string2int(cigcts.last);
              cells(seqCol) = cells(seqCol).substring(0,cells(seqCol).length - ct);
              cells(qualCol) = cells(qualCol).substring(0,cells(qualCol).length - ct);
              cigops(cigops.length - 1) = "H";
              //cigcts(cigcts.length - 1) = "";
            }
          
            cells(cigarCol) = cigcts.zip(cigops).foldLeft("")((b : String,a : (String,String))  =>   b+a._1+a._2   );
            
            val outLine = cells.tail.foldLeft(cells.head)((a : String,b : String) => a + "	"+ b);
            out.write(outLine + "\n");
            
          } else {
            out.write(line+"\n");
          }
        }
      }
      } catch {
        case e : Exception => {
          scala.Console.err.println("Error encountered at line: \n\""+line+"\"");
          out.close();
          scala.Console.err.println("output stream closed.\nThe Original Exception:");
          throw e;
        }
      }
    }
    out.close();
    //for (line <- io.Source.stdin.getLines){
      
    //}
  }
  
}