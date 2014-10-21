package fileConversionUtils

import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
import java.io.OutputStream
import java.io.FileOutputStream
import java.io.InputStream
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.File
import scala.collection.JavaConversions._

import net.sf.samtools._

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;

import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.JavaConversions._;
import internalUtils.optionHolder._;

import internalUtils.genomicUtils._;
import qcUtils.QCUtility;

object bamToWiggle { 
  
  class wiggleMaker extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "bamTowiggle", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "Calculates depth-of-coverage for each equally-sized window across the entire genome. Notes: This counts paired reads seperately. Unlike mpileup, it does not count areas as covered when they span a read's splice junction.",   
          argList = 
                    new BinaryArgument[Double](   name = "sizefactor",
                                                        arg = List("--sizefactor"),  
                                                        valueName = "float", 
                                                        argDesc = "The size factor, for normalization. Defaults to 1. If set, then the count of every cell will be divided by the given value.", 
                                                        defaultValue = Some(1.0)
                                                        ) ::
                    new UnaryArgument(    name = "noGzipOutput",
                                         arg = List("--noGzipOutput"), // name of value
                                         argDesc = "Flag to indicate that output should NOT be gzipped." // description
                                       ) ::
                    new UnaryArgument(    name = "negativeReverseStrand",
                                         arg = List("--negativeReverseStrand"), // name of value
                                         argDesc = "Flag to indicate that the reverse strand should be counted in negative values. Can be useful for plotting stranded data on a single multiwig track, via trackhubs." // description
                                       ) ::
                    new UnaryArgument(    name = "testRun",
                                         arg = List("--testRun","-t"), // name of value
                                         argDesc = "Flag to indicate that only the first 100k reads should be read in, used for test runs." // description
                                       ) ::
                    new UnaryArgument(   name = "stranded", 
                                         arg = List("--stranded"), // name of value
                                         argDesc = "The stranded flag. If this is set, then two wiggle files are created, one for each strand." // description
                                       ) ::
                    new UnaryArgument(   name = "fr_secondstrand", 
                                         arg = List("--stranded_fr_secondstrand"), // name of value
                                         argDesc = "If this option is set, the data is assumed to be a fr_secondstrand type stranded library. By default the assumed library type for stranded data is fr_firststrand." // description
                                       ) ::
                    new UnaryArgument(   name = "isSingleEnd", 
                                         arg = List("--isSingleEnd"), // name of value
                                         argDesc = "Flag for single-end data. CURRENTLY NONFUNCTIONAL: THIS UTILITY IS INTENDED ONLY FOR USE WITH PAIRED-END DATA!" 
                                       ) ::
                    new UnaryArgument(   name = "noTruncate", 
                                         arg = List("--noTruncate"), // name of value
                                         argDesc = "The UCSC tool wigToBigWig only allows wiggle files in which every window is of equal size. "+
                                                   "This means that if the chromosome size is not divisible by the window size, a few bases are not "+
                                                   "counted on the end of the chromosome. Using this flag will cause this utility to not truncate off"+
                                                   " the last odd-sized window. However, be aware that this will mean that you cannot use wigToBigWig to"+
                                                   " convert the wiggle file to a (smaller and more efficient) bigWig file." 
                                       ) ::
                    new BinaryArgument[Int](   name = "windowSize",
                                                        arg = List("--windowSize"),  
                                                        valueName = "num", 
                                                        argDesc = "The length, in base-pairs, for each counting bin, or \"window\". Note: if this is set low the utility will take longer to run and will consume more memory. The default window size is 100bp.", 
                                                        defaultValue = Some(100)
                                                        ) ::
                    new FinalArgument[List[String]](
                                         name = "infiles",
                                         valueName = "infile.bam[,infile2.bam,...]",
                                         argDesc = "The input sam or bam file or files, or '-' to read from stdin. Note: if you include more than one bam file, the list must be comma delimited with no whitespace!" // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "trackName",
                                         valueName = "trackName",
                                         argDesc = "The name for the wiggle track." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "chromLengthFile",
                                         valueName = "chromLengthFile",
                                         argDesc = "The chrom length file. This should be a simple tab-delimited text file that includes "+
                                                   "each chromosome in the first column and the chromosome's length (in base-pairs) in the second "+
                                                   "column. If the wiggle file is intended for use with a standard genome on the UCSC genome browser,"+
                                                   " then the UCSC utility \"fetchChromSizes\" should be used to generate this file, as the UCSC genome browser "+
                                                   "and related utilities may have difficulties if the chrom lengths do not match the expected values."
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfilePrefix",
                                         valueName = "outfilePrefix",
                                         argDesc = "The output file prefix. If unstranded, this will produce one file named \"outfilePrefix.unstranded.wig.gz\". If stranded, this will produce two files: \"outfilePrefix.fwd.wig.gz\" and \"outfilePrefix.rev.wig.gz\". If the \"--noGzipOutput\" flag is raised then the output files will not have the \".gz\" extension at the end. IMPORTANT NOTE: if the window size is set to any size other than the default of 100, the window size will be added to the filename. The filename will be  \"outfilePrefix.win50.unstranded.wig.gz\" for unstranded wiggles with a 50bp window, and so on." // description
                                        ) :: List() );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         bamToWiggle.run(
             parser.get[List[String]]("infiles"),
             parser.get[String]("outfilePrefix"),
             parser.get[String]("trackName"),
             parser.get[String]("chromLengthFile"),
             parser.get[Boolean]("noTruncate"),
             parser.get[Int]("windowSize"),
             parser.get[Boolean]("isSingleEnd"),
             parser.get[Boolean]("stranded"),
             parser.get[Boolean]("fr_secondstrand"),
             parser.get[Double]("sizefactor"),
             parser.get[Boolean]("testRun"),
             parser.get[Boolean]("noGzipOutput"),
             parser.get[Boolean]("negativeReverseStrand")
           );
         }
     }
   }
  //createWiggle(inputSamOrBamFile : File, 
  //trackName : String, 
  //countMethod : Int, 
  //stranded : Boolean, 
  //noFilter : Boolean, 
  //truncate : Boolean, 
  //singleEnd : Boolean, 
  //span : Int,
  //chromLengthFile : String, 
  //outfile : String, 
  //adjustmentFactor : Option[Double], 
  //negativeStrand : Boolean, 
  //getStrandFromFirstRead : Boolean)
  
  def run(infiles : List[String], outfilePrefix : String, trackName : String, chromLengthFile : String, noTruncate : Boolean, windowSize : Int, isSingleEnd : Boolean, stranded : Boolean, fr_secondStrand : Boolean, sizeFactor : Double, testRun : Boolean, noGzipOutput : Boolean, negativeReverseStrand : Boolean){
    
    //registerGlobalParam[Boolean]("noGzipOutput", noGzipOutput);
    internalUtils.optionHolder.OPTION_noGzipOutput = noGzipOutput;
    
    val initialTimeStamp = TimeStampUtil();
    standardStatusReport(initialTimeStamp);
    
    val qcBTW : QcBamToWig = new QcBamToWig(trackName , chromLengthFile , noTruncate , windowSize , isSingleEnd, stranded , fr_secondStrand, sizeFactor, negativeReverseStrand );
    
    for(infile <- infiles){
      runOnFile(infile , qcBTW , testRun );
    }
    
    val postRunStamp = TimeStampUtil();
    standardStatusReport(postRunStamp);
    
    qcBTW.writeOutput(outfilePrefix, null);
    
    val finalStamp = TimeStampUtil();
    standardStatusReport(finalStamp);
  }
    
  def runOnFile(infile : String, qcBTW : QcBamToWig, testRun : Boolean){
    val (samFileAttributes, recordIter) = initSamRecordIterator(infile);
    val pairedIter : Iterator[(SAMRecord,SAMRecord)] = if(testRun) samRecordPairIterator_withMulti(recordIter, true, 200000) else samRecordPairIterator_withMulti(recordIter);
    
    val coda : Array[Int] = internalUtils.commonSeqUtils.getNewCauseOfDropArray;
    val coda_options : Array[Boolean] = internalUtils.commonSeqUtils.CODA_DEFAULT_OPTIONS.toArray;
    
    var readNum = 0;
    val samIterationTimeStamp = TimeStampUtil();
    for(pair <- pairedIter){
    //for((pair,readNum) <- numberedIter){
      val (r1,r2) = pair;
      readNum += 1;
      if(internalUtils.commonSeqUtils.useReadPair(r1,r2,coda, coda_options, Set[String](), None)){
        qcBTW.runOnReadPair(r1,r2,readNum);
      }
    }
  }
  
  
  class QcBamToWig(trackName : String, chromLengthFile : String, noTruncate : Boolean, windowSize : Int, isSingleEnd: Boolean, stranded : Boolean, fr_secondStrand : Boolean, sizeFactor : Double, negativeReverseStrand : Boolean) extends QCUtility[Unit] {
    val chromMap : Map[(String,Char),Chrom] = genChrom(chromLengthFile, windowSize, stranded, ! noTruncate);
    var unknownChromSet : Set[String] = Set[String]();
    
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
      runOnRead(r1);
      runOnRead(r2);
    }
    def runOnRead(r : SAMRecord){
      val chromName = r.getReferenceName();
      val strand = internalUtils.commonSeqUtils.getStrand(r , stranded , fr_secondStrand );
      chromMap.get((chromName,strand)) match {
        case Some(chrom) => {
          chrom.countSamRecord(r);
        }
        case None => {
          // error("Chromosome " + chromName + " not found! Strand: \'" + strand+ "\'. Something is wrong with the bam file read!");
          //  report("!["+chromName+","+strand+"]","deepDebug");
          if(! unknownChromSet.contains(chromName)){
            reportln("WARNING: Chromosome not found in the chromLengthFile: ["+chromName+","+strand+"]","deepDebug");
            unknownChromSet += chromName;
          }
        }
      }
    }
    def writeOutput(outfile : String, summaryWriter : WriterUtil){
      val windowString = if(windowSize == 100) "" else (".Win"+windowSize.toString());
      
      if(stranded){
        val writerF = openWriterSmart_viaGlobalParam(outfile + windowString + ".fwd.wig");
        val writerR = openWriterSmart_viaGlobalParam(outfile + windowString + ".rev.wig");
        
        val sortedKeyList : Vector[(String,Char)] = chromMap.keySet.toVector.sorted
        for(chromPairs <- sortedKeyList){
          	val chrom = chromMap(chromPairs);
	        if(chrom.chromStrand == '+'){
	          writeChrom(chrom,writerF,sizeFactor,false);
	        } else {
	          writeChrom(chrom,writerR,sizeFactor,negativeReverseStrand);
	        }
        }
        close(writerF);
        close(writerR);
      } else {
        val writer = openWriterSmart_viaGlobalParam(outfile + windowString + ".unstranded.wig");
        //writer.write("track type=wiggle_0 name="+trackName+" visibility=full\n");
        
        val sortedKeyList : Vector[(String,Char)] = chromMap.keySet.toVector.sorted
        for(chromPairs <- sortedKeyList){
          val chrom = chromMap(chromPairs);
	      writeChrom(chrom,writer,sizeFactor,false);
        }
        close(writer);
      }
    }
    
    def getUtilityName : String = "bamToWig";

  }
  
  case class Chrom(chromName :  String, chromStrand : Char, windowCounts : Array[Int], span : Int, truncate : Boolean){
    def countSamRecord(samRecord : SAMRecord) {
        val blocks = samRecord.getAlignmentBlocks();
      
        var i = 0;
        while(i < blocks.size()){
          val block : AlignmentBlock = blocks.get(i);
          val start = block.getReferenceStart();
          val length = block.getLength();
          val end = start + length;
          
          if(! truncate && (end - 1 / span) >=  windowCounts.length){
            error("ERROR: read extends outside chromosome length!");
          }
          
          var j = start;
          while(j < end && (j / span) < windowCounts.length){
            windowCounts(j / span) = windowCounts(j / span) + 1;
            j = j + 1;
          }
          
          i = i + 1;
        }
    }
  }
  

  def writeChrom(chrom : Chrom, writer : WriterUtil, sizeFactor : Double, negativeStrand : Boolean){
    writer.write("fixedStep  chrom="+chrom.chromName+"  start=1  step="+chrom.span+"\n");
    val adjustmentFactor = if(negativeStrand) ((-1).toDouble * sizeFactor * chrom.span.toDouble) else (sizeFactor * chrom.span.toDouble);        
    chrom.windowCounts.foreach(
      (ct : Int) => {
        writer.write((ct.toDouble / adjustmentFactor)  +"\n");
      }
    );
  }
  
  def genChrom(chromLengthFile : String, span : Int, stranded : Boolean, truncate : Boolean) : Map[(String,Char),Chrom] = {
    var chromMap = Map[(String,Char),Chrom]();
    
    val lines = getLines(chromLengthFile);
    
    for(line <- lines){
      val cells = line.split("\\s+");
      val chromName = cells(0);
      val chromLength = string2int(cells(1));
      val spanCt =  if(chromLength % span == 0 || truncate){
                       chromLength / span;
                    } else {
                       (chromLength / span) + 1;
                    }
      
      if(stranded){
        val chromP = new Chrom(chromName,'+',Array.ofDim[Int](spanCt),span, truncate);
        val chromM = new Chrom(chromName,'-',Array.ofDim[Int](spanCt),span, truncate);
        chromMap = chromMap + (((chromName,'+'), chromP));
        chromMap = chromMap + (((chromName,'-'), chromM));        
      } else {
        val chrom = new Chrom(chromName,'.',Array.ofDim[Int](spanCt),span, truncate);
        chromMap = chromMap + (((chromName,'.'), chrom));
      }
    }
    chromMap;
  }
}