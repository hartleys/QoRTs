package fileConversionUtils


import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.File;
import internalUtils.stdUtils._;
import internalUtils.Reporter._;
import internalUtils.commandLineUI._;

object makeBedFromGtf {

class converter extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "makeBedFromGtf", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This simple utility converts a gtf transcript annotation file into a .bed transcript annotation file."+
                        "Note that this conversion may theoretically be \"lossy\", as it is possible for gtf files to contain information that "+
                        "cannot be perfectly translated into the bed format."+
                        ""+
                        "WARNING: THIS SUB-UTILITY IS BETA! NOT FOR GENERAL USE!"+
                        ""+
                        ""+
                        ""+
                        "",   
          argList = 
                    new BinaryOptionArgument[String](
                                         name = "rgb", 
                                         arg = List("--rgb"), 
                                         valueName = "r,g,b",  
                                         argDesc = "The rgb color for all the bed file lines."
                                        ) ::
                    new FinalArgument[String](
                                         name = "gtffile",
                                         valueName = "annotation.gtf.gz",
                                         argDesc = "The gtf file, or '-' to read from stdin. If the filename ends with \".gz\" or \".zip\" then the file will be decompressed using the appropriate method." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output bed file. If the filename ends with \".gz\" or \".zip\" then the file will be compressed using the appropriate method." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         makeBedFromGtf.run(
             parser.get[String]("gtffile"),
             parser.get[String]("outfile"),
             parser.get[Option[String]]("rgb")
           );
         }
     }
   }
  
 // def test() {
  def run(infile : String, outfile : String, rgb : Option[String]){

	val extraAnnoBed =  false;
     
    val lines : Iterator[String] = internalUtils.fileUtils.getLinesSmartUnzip(infile, true);
    
   // var line : String = reader.readLine();
    var lnct : Int = 0;
    var featurect : Int = 0;
    var transct : Int = 0;
    var cdsct : Int = 0;
    var mdMap : Map[String,ReadyMetadata] = Map[String,ReadyMetadata]();
        
    for( line <- lines){
      val cells : Array[String] = line.split("	");
      
      if(cells(featureIndex) == exonString || cells(featureIndex) == cdsString || cells(featureIndex) == startString || cells(featureIndex) == endString){
        try {
          //println("starting line: " + lnct);
          val newMd = getMetaData(cells);
          val name = newMd.name;
                    
          if(mdMap.contains(newMd.name)){
            mdMap = (mdMap - name) + ((name,mergeMetadata(mdMap(name),newMd, lnct))) ;
          } else {
            mdMap = mdMap + ((name,newMd));
          }
          featurect = featurect + 1;
        } catch {
          case e : Throwable => {
            println("Error on line: " + lnct);
            println("Line is: \n" + line);
            throw e;
          }
        }
      }
      
     if(lnct % 100000 == 0) reportln(lnct + " gtf lines read.","debug")
     //if(lnct % 100000 == 0){if(lnct%500000==0) print("("+exonct+"/"+genect+":"+md._1+","+md._2+")\n"); else print("("+exonct+"/"+genect+":"+md._1+","+md._2+")");}
     //if(lnct % 100000 == 0){if(lnct%500000==0) print("("+featurect+"/"+transct+")\n"); else print("("+featurect+"/"+transct+")");}
     
      lnct = lnct + 1;
      //line = reader.readLine();
    }
    
    println("Found " + mdMap.size + " transcripts, over " + lnct + " lines.");
    
    println("done reading. Writing...");
    
    val writer : BufferedWriter = new BufferedWriter(new FileWriter(outfile));

    for(md <- mdMap) {
      writeLine(md._2 , writer, extraAnnoBed );
    }
    //reader.close;
    writer.close;
    println("done");
  }

  
  final def chrIndex = 0;
  final def sourceIndex = 1;
  final def featureIndex = 2;
  final def startIndex = 3;
  final def endIndex = 4;
  final def scoreIndex = 5;
  final def strandIndex = 6;
  
  final def attrIndex = 8;
  
  final def geneIdKey = "gene_id";
  final def transIdKey = "transcript_id";
  
  final val exonString = "exon";
  final val cdsString = "CDS";
  final val startString = "start_codon";
  final val endString = "stop_codon";
  final val useLines : Set[String] = Set(exonString,cdsString,startString,endString)

  def getAttribute(attrString : String, key : String) : String = {
    def splitted : List[Array[String]] = attrString.trim().split(";").map(x => x.trim().split(" ")).toList;
    (splitted.find(x => x(0) == key)) match {
      case Some(x : Array[String]) => x(1);
      case None => {
        error("FATAL ERROR: could not find key \""+key+"\" in attribute string:\n\""+attrString+"\"\nWhich was split into: \n"+
            splitted.foldLeft[String]("|")((a : String, x : Array[String]) => 
              a + "|" + x.foldLeft("")((b : String, y : String) =>  b+"/"+y))  );
        ""
      }
    }
  }
  def getTxName(cells : Array[String]) : String = {
    getAttribute(cells(attrIndex),transIdKey);
  }

  abstract class Metadata;
 // case class EmptyMetadata() extends Metadata;
  case class ReadyMetadata(chr : String, start : IntOption, end : IntOption, name : String, geneID : String, score : IntOption, strand : Char, thickStart : IntOption, thickEnd : IntOption, itemRgb : String, blockCount : Int, blocks : List[(Int,Int)]) extends Metadata;

  def writeBedLine(md : Metadata, delim : String, extendedBed : Boolean): String = md match{
    //case EmptyMetadata() => throw new Exception("ERROR: Impossible state! Something is wrong!")
    case ReadyMetadata(chr, start , end , name , geneID, score, strand , thickStart , thickEnd , itemRgb, blockCount, blocks) => 
      
      chr +"	"+ start +"	"+ end +"	"+ name +"	"+ score.toStringWithDefault(IntOpt(0)) +"	"+ strand +
      "	"+ thickStart.toStringWithDefault(start) +"	"+ thickEnd.toStringWithDefault(end) +"	"+ itemRgb +"	"+ blockCount +"	" + 
      blocksToString(blocks,start.getInt) + (if(extendedBed) "	" + geneID + blocksToString_extendedVer(blocks) else "");
  }
  
  def blocksToString(blocks : List[(Int,Int)], start : Int) : String = {
    val blocksSorted = blocks.map((b) => (b._1,b._2 - start)).sortBy((b) => b._2);
    
    blocksSorted.tail.foldLeft(blocksSorted.head._1.toString)((s,i) => s + "," + i._1) + "	" +
    blocksSorted.tail.foldLeft(blocksSorted.head._2.toString)((s,i) => s + "," + i._2);
  }
  
  def blocksToString_extendedVer(blocks : List[(Int,Int)]) : String = {
    val blocksSorted = blocks.map((b) => (b._2,b._1+b._2)).sortBy((b) => b._1);
    
    blocksSorted.foldLeft("")((a,b) => a + "	"+b._1+"	"+b._2);
  }
  
  class IntOption {
    override def toString() : String = this match {
      case IntOpt(value) => value.toString
      case IntNone() => "."
    }
    def toStringWithDefault(defInt : IntOption) : String = this match {
      case IntOpt(value) => value.toString;
      case IntNone() => defInt.toString;
    }
    
    def merge(that : IntOption)(f : (Int, Int) => Int) : IntOption = this match {
      case IntOpt(value) => that match {
        case IntOpt(thatVal) => IntOpt(f(value,thatVal));
        case IntNone() => this;
      }
      case IntNone() => that;
    }
    def plusOne() : IntOption = this match {
      case IntOpt(ivalue) => IntOpt(ivalue + 1)
      case IntNone() => IntNone();
    }
    def minusOne() : IntOption = this match {
      case IntOpt(ivalue) => IntOpt(ivalue - 1)
      case IntNone() => IntNone();
    }
    def getInt() : Int = this match {
      case IntOpt(value) => value;
      case IntNone() => throw new Exception("GTF Formatting Error: Some that should be set for every transcript is not set for at least one transcript. Probably the transcript start or end.");
    }
  }
  case class IntOpt(value : Int) extends IntOption;
  case class IntNone() extends IntOption;
  
  def getIntOption(s : String) : IntOption = {
    if(s == ".") IntNone();
    else IntOpt(augmentString(s).toInt);
  }

  
  def writeLine(md : Metadata, writer : BufferedWriter, extraAnnoBed : Boolean) = {
    writer.write(writeBedLine(md,"	",extraAnnoBed) + "\n");
  }
  
  def getMetaData(cells : Array[String]) : ReadyMetadata = {
    ReadyMetadata(cells(chrIndex), 
        getIntOption(cells(startIndex)).minusOne ,
        getIntOption(cells(endIndex)) , 
        getAttribute(cells(attrIndex),transIdKey) , 
        getAttribute(cells(attrIndex),geneIdKey) , 
        getIntOption(cells(scoreIndex)) , 
        cells(strandIndex).charAt(0) , 
        if(cells(featureIndex) == cdsString) getIntOption(cells(startIndex)).minusOne else IntNone() , 
        if(cells(featureIndex) == cdsString) getIntOption(cells(endIndex)) else IntNone() , 
        "0", 
        if(cells(featureIndex) != exonString) 0 else 1, 
        if(cells(featureIndex) != exonString) List() else List((augmentString(cells(endIndex)).toInt + 1 - augmentString(cells(startIndex)).toInt, augmentString(cells(startIndex)).toInt - 1) ));
  } 
  
  
  def mergeMetadata(md1 : Metadata, md2 : ReadyMetadata, lnct : Int) : ReadyMetadata = md1 match {
   // case EmptyMetadata() => md2;
    case ReadyMetadata(chr, start , end , name , geneID, score, strand , thickStart , thickEnd , itemRgb, blockCount, blocks) => {
      ReadyMetadata(
          if(chr == md2.chr) chr else {throw new Exception("Gtf Format Error on line " + lnct + ". Transcript " + name + " found to span multiple chromosomes!")}, 
          start.merge(md2.start)(math.min(_,_)) , 
          end.merge(md2.end)(math.max(_,_)) ,
          if(name == md2.name) name else {throw new Exception("Error: Impossible state! There is a bug in my code! On gtf line " + lnct + " with transcript " + name)} , 
          if(geneID == md2.geneID) geneID else {throw new Exception("Error: Impossible state! There is a bug in my code! On gtf line " + lnct + " with transcript " + name)} ,
          score.merge(md2.score)((x,y) => x), 
          if(strand == md2.strand) strand else throw new Exception("Gtf Format Error on line " + lnct + ". Transcript " + name + " has inconsistant strand!"), 
          thickStart.merge(md2.thickStart)(math.min(_,_)) , 
          thickEnd.merge(md2.thickEnd)(math.max(_,_)) ,
          itemRgb, 
          blockCount + md2.blockCount, 
          blocks ++ md2.blocks);
    }
  }
   

}
