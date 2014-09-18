package qcUtils

import java.util.zip.GZIPInputStream

import java.util.zip.GZIPOutputStream
import java.io.OutputStream
import java.io.FileOutputStream
import java.io.InputStream
import java.io.ByteArrayInputStream
import java.io.FileInputStream
import java.io.File
import scala.collection.JavaConversions._
import scala.collection.mutable.HashMap;

import net.sf.samtools._

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;

import scala.collection.parallel.immutable.ParVector;

object runAllQC {
  
  final val QC_DEFAULT_ON_FUNCTION_LIST  : Set[String] = Set("InsertSize",
                                                             "NVC",
                                                             "CigarOpDistribution",
                                                             "QualityScoreDistribution",
                                                             "GCDistribution",
                                                             "GeneCalcs",
                                                             "JunctionCalcs", 
                                                             "StrandCheck", 
                                                             "writeKnownSplices",
                                                             "writeNovelSplices",
                                                             "writeSpliceExon", 
                                                             "writeDESeq",
                                                             "writeDEXSeq",
                                                             "writeGenewiseGeneBody", 
                                                             "writeGeneCounts", 
                                                             "writeClippedNVC", 
                                                             "chromCounts");
  final val QC_DEFAULT_OFF_FUNCTION_LIST : Set[String] = Set("FPKM","makeAllBrowserTracks", "cigarMatch");
  final val QC_FUNCTION_LIST : Set[String] = QC_DEFAULT_ON_FUNCTION_LIST ++ QC_DEFAULT_OFF_FUNCTION_LIST;

  
  //"InsertSize","NVC","CigarOpDistribution","QualityScoreDistribution","GCDistribution","GeneCalcs",
  //"JunctionCounts", "StrandCheck", "writeKnownSplices","writeNovelSplices","writeSpliceExon", 
  //"writeDESeq","writeDEXSeq","writeGenewiseGeneBody", "writeGeneCounts"
  final val QC_FUNCTION_DEPENDANCIES : Map[String,String] = Map[String,String](
      ("writeDESeq" -> "GeneCalcs"),
      ("writeGeneCounts" -> "GeneCalcs"),
      ("writeGenewiseGeneBody" -> "GeneCalcs"),
      ("writeDEXSeq" -> "JunctionCalcs"),
      ("writeSpliceExon" -> "JunctionCalcs"),
      ("writeKnownSplices" -> "JunctionCalcs"),
      ("writeNovelSplices" -> "JunctionCalcs"),
      ("writeClippedNVC" -> "NVC")
      );
  
  
  //def run(args : Array[String]){
  
  class allQC_runner extends CommandLineRunUtil {
    
    val parser : CommandLineArgParser = 
      new CommandLineArgParser(
          command = "QC", 
          quickSynopsis = "Runs a battery of QC tools", 
          synopsis = "", 
          description = "This utility runs a large battery of QC tools on a single given sam or bam file."+
                        "This is the primary function of the QoRT utility. "+
                        ""+
                        ""+
                        ""+
                        "",   
          argList = 
                    new BinaryOptionArgument[String](
                                         name = "flatgfffile", 
                                         arg = List("--flatgff"), 
                                         valueName = "flattenedGffFile.gff.gz",  
                                         argDesc = "A \"flattened\" gtf file that matches the standard gtf file. Optional, "+
                                                   "but using this option may make the utility run faster and consume less memory."+
                                                   " It may also be useful for downstream analyses, as it assigns unique identifiers to all exons and splice "+
                                                   "junctions. The flattened gtf file can be generated using "+
                                                   "the \"makeFlatGtf\" command.\n"+
                                                   "If the filename ends with \".gz\" or \".zip\", the file will be parsed using the appropriate decompression method."
                                        ) ::
                    new UnaryArgument(   name = "isSingleEnd", 
                                         arg = List("--isSingleEnd","-e"), // name of value
                                         argDesc = "Flag to indicate that reads are single end. WARNING: UNIMPLEMENTED and UNSUPPORTED AT THIS TIME! This utility is not designed for, and will not function properly on single-strand data!" // description
                                       ) ::
                    new UnaryArgument( name = "stranded",
                                         arg = List("--stranded","-s"), // name of value
                                         argDesc = "Flag to indicate that data is stranded." // description
                                       ) ::
                    new UnaryArgument( name = "fr_secondStrand",
                                         arg = List("--stranded_fr_secondstrand","-a"), // name of value
                                         argDesc = "Flag to indicate that reads are from a fr_secondstrand type of stranded library (equivalent to the \"stranded = yes\" option in HTSeq). "+
                                                   "If your data is stranded, you must know the library type in order to analyze it properly. This utility uses the same "+
                                                   "definitions as cufflinks to define strandedness type. By default, the fr_firststrand "+
                                                   "library type is assumed for all stranded data (equivalent to the \"stranded = reverse\" option in HTSeq)." // description
                                       ) ::
                    new UnaryArgument(    name = "testRun",
                                         arg = List("--testRun","-t"), // name of value
                                         argDesc = "Flag to indicate that only the first 100k reads should be read in. Used for testing." // description
                                       ) ::
                    new UnaryArgument( name = "noMultiMapped",
                                         arg = List("--fileContainsNoMultiMappedReads"), // name of value
                                         argDesc = "Flag to indicate that the input sam/bam file contains only primary alignments (ie, no multi-mapped reads). This flag is ALWAYS OPTIONAL, but when applicable this utility will run (slightly) faster using this argument." // description
                                       ) ::
                    new UnaryArgument( name = "keepMultiMapped",
                                         arg = List("--keepMultiMapped"), // name of value
                                         argDesc = "Flag to indicate that the tool should NOT filter out multi-mapped reads." // description
                                       ) ::
                    new UnaryArgument( name = "noGzipOutput",
                                         arg = List("--noGzipOutput"), // name of value
                                         argDesc = "Flag to indicate that output files should NOT be compressed into the gzip format. By default almost all output files are compressed to save space." // description
                                       ) ::
                    new BinaryOptionArgument[String](
                                         name = "readGroup", 
                                         arg = List("--readGroup"), 
                                         valueName = "readGroupName",  
                                         argDesc =  "If this option is set, all analyses will be restricted to reads that are tagged with the given "+
                                                    "readGroupName using an RG tag. This can be used if multiple read groups have already been combined "+
                                                    "into a single bam file."
                                        ) ::
                    new BinaryArgument[Int](name = "numThreads",
                                                        arg = List("--numThreads"),  
                                                        valueName = "num", 
                                                        argDesc = "The number of threads to allow. By default this utility will only allow one thread to be used.", 
                                                        defaultValue = Some(1)
                                                        ) :: 
                    new BinaryArgument[List[String]](   name = "dropChromList",
                                                        arg = List("--dropChrom"),  
                                                        valueName = "dropChromosomes", 
                                                        argDesc = "A comma-delimited list of chromosomes to ignore and exclude from all analyses. Important: no whitespace!", 
                                                        defaultValue = Some(List[String]())
                                                        ) ::
                    new BinaryArgument[List[String]](   name = "skipFunctions",
                                                        arg = List("--skipFunctions"),  
                                                        valueName = "func1,func2,...", 
                                                        argDesc = "A comma-delimited list of functions to skip. Important: No whitespace! The default-on functions are: "+QC_DEFAULT_ON_FUNCTION_LIST.mkString(", "), 
                                                        defaultValue = Some(List[String]())
                                                        ) ::
                    new BinaryArgument[List[String]](name = "addFunctions",
                                                        arg = List("--addFunctions"),  
                                                        valueName = "func1,func2,...", 
                                                        argDesc = "A list of functions to add. This can be used to add functions that are off by default. Followed by a comma delimited list, with no internal whitespace. The default-off functions are: "+QC_DEFAULT_OFF_FUNCTION_LIST.mkString(", "), 
                                                        defaultValue = Some(List[String]())
                                                        ) :: 
                    new BinaryArgument[List[String]](name = "runFunctions",
                                                        arg = List("--runFunctions"),  
                                                        valueName = "func1,func2,...", 
                                                        argDesc = "The complete list of functions to run. Setting this option turns off ALL functions EXCEPT for the ones explicitly requested here. Followed by a comma delimited list, with no internal whitespace. Allowed options are: "+QC_FUNCTION_LIST.mkString(", "), 
                                                        defaultValue = Some(List[String]())
                                                        ) :: 
                    new FinalArgument[String]( 
                                         name = "infile",
                                         valueName = "infile",
                                         argDesc = "The input .bam or .sam file of aligned sequencing reads. Or \'-\' to read from stdin."
                                        ) :: 
                    new FinalArgument[String](
                                         name = "gtffile",
                                         valueName = "gtffile.gtf",
                                         argDesc = "The gtf annotation file. This tool was designed to use the standard gtf annotations provided by Ensembl, but other annotations can be used as well.\n" +
                                                   "If the filename ends with \".gz\" or \".zip\", the file will be parsed using the appropriate decompression method."
                                        ) :: 
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfiledir",
                                         argDesc = "The output file directory." // description
                                        ) :: List()
      );
    
    def run(args : Array[String]){
      val out = parser.parseArguments(args.toList.tail);
      
      if(out){
      runAllQC.run(parser.get[String]("infile") ,
          parser.get[String]("outfile") + "/" + "QC",
          parser.get[String]("gtffile"),
          parser.get[Option[String]]("flatgfffile"),
          parser.get[List[String]]("dropChromList"),
          parser.get[Boolean]("isSingleEnd"),
          parser.get[Boolean]("stranded"),
          parser.get[Boolean]("fr_secondStrand"),
          parser.get[Boolean]("testRun"),
          parser.get[List[String]]("runFunctions"),
          parser.get[List[String]]("addFunctions"),
          parser.get[List[String]]("skipFunctions"),
          parser.get[Boolean]("noGzipOutput"),
          parser.get[Boolean]("noMultiMapped"),
          parser.get[Boolean]("keepMultiMapped"),
          parser.get[Int]("numThreads"),
          parser.get[Option[String]]("readGroup")
      );
      }
    }
  }
  
  def run(infile : String, 
          outfile : String, 
          gtffile : String, 
          flatgtffile : Option[String], 
          dropChromList : List[String], 
          isSingleEnd : Boolean, 
          stranded : Boolean, 
          fr_secondStrand : Boolean, 
          testRun : Boolean, 
          runFunctions : List[String], 
          addFunctions : List[String], 
          dropFunctions : List[String], 
          noGzipOutput : Boolean, 
          noMultiMapped : Boolean,
          keepMultiMapped : Boolean,
          numThreads : Int,
          readGroup : Option[String]){
    reportln("Starting ALLQC:","note");
    val initialTimeStamp = TimeStampUtil();
    standardStatusReport(initialTimeStamp);
    
    val stdGtfCodes = new internalUtils.GtfTool.GtfCodes();
    val flatGtfCodes = new internalUtils.GtfTool.GtfCodes();
    
    if(isSingleEnd) error("FATAL ERROR: Single-end option SET. Single-end data is not currently supported!");
    
    val dropChrom = dropChromList.toSet;
    
    val runFunc_initial = if(runFunctions.isEmpty){
      (QC_DEFAULT_ON_FUNCTION_LIST ++ addFunctions.toSet) -- dropFunctions.toSet;
    } else {
      runFunctions.toSet -- dropFunctions.toSet;
    }
    
    val runFunc = runFunc_initial.foldLeft(runFunc_initial)((soFar,currFunc) => {
      QC_FUNCTION_DEPENDANCIES.get(currFunc) match {
        case Some(reqFunc) => {
          if(soFar.contains(reqFunc)){
            soFar;
          } else {
            reportln("Function \"" + currFunc +"\" requires top-level function \"" + reqFunc+ "\". Adding the required function to the active function list.","warning");
            soFar + reqFunc;
          }
        }
        case None => {
          soFar;
        }
      }
    });
    
    val notFoundFunc = runFunc.find(  ! QC_FUNCTION_LIST.contains(_) ) ;
    if( ! notFoundFunc.isEmpty ){
      error("ERROR: Function not found: \"" + notFoundFunc.get +"\"");
    }
    
    reportln("Running functions: " + runFunc.mkString(","),"progress");
    
    //reportln("infile: " + infile , "note");
    //reportln("outfile: " + outfile , "note");
    //reportln("gtffile: " + gtffile , "note");
    //reportln("flatgtffile: " + flatgtffile , "note");
    //reportln("stranded: " + stranded , "note");
    //reportln("fr_secondStrand: " + fr_secondStrand , "note");
    //reportln("testRun: " + testRun , "note");
    //reportln("dropChrom: "+dropChromList.mkString(";"),"note");
    //reportln("run_functions: " + runFunc.mkString(";"),"note");
    
    setQcOptions(noGzipOutput);
    val anno_holder = new qcGtfAnnotationBuilder(gtffile , flatgtffile , stranded , stdGtfCodes, flatGtfCodes);

    runOnSeqFile(initialTimeStamp = initialTimeStamp, infile = infile, outfile = outfile, anno_holder = anno_holder, testRun = testRun, runFunc = runFunc, stranded = stranded, fr_secondStrand = fr_secondStrand, dropChrom = dropChrom, keepMultiMapped = keepMultiMapped, noMultiMapped = noMultiMapped, numThreads = numThreads, readGroup )
  }
  
  def runOnSeqFile(initialTimeStamp : TimeStampUtil, 
                   infile : String, 
                   outfile : String, 
                   anno_holder : qcGtfAnnotationBuilder, 
                   testRun  :Boolean, 
                   runFunc : Set[String], 
                   stranded : Boolean, 
                   fr_secondStrand : Boolean, 
                   dropChrom : Set[String], 
                   keepMultiMapped : Boolean, 
                   noMultiMapped : Boolean, 
                   numThreads : Int,
                   readGroup : Option[String]){
    
    
    val (samFileAttributes, recordIter) = initSamRecordIterator(infile);
    val pairedIter : Iterator[(SAMRecord,SAMRecord)] = 
      if(noMultiMapped){
        if(testRun) samRecordPairIterator(recordIter, true, 200000) else samRecordPairIterator(recordIter)
      } else {
        if(testRun) samRecordPairIterator_withMulti(recordIter, true, 200000) else samRecordPairIterator_withMulti(recordIter)
      }
    
    val readLength = samFileAttributes.readLength;
    val isSortedByName = samFileAttributes.isSortedByName;
    val isSortedByPosition = samFileAttributes.isSortedByPosition;
    val isDefinitelyPairedEnd = samFileAttributes.isDefinitelyPairedEnd;
    val minReadLength = samFileAttributes.minReadLength;
    
    if(readLength != minReadLength){reportln("Warning: Read length is not consistent! In the first 1000 reads, read length varies from "+minReadLength+" to " +readLength+"!\nThis may cause odd things to happen. In general, it is STRONGLY recommended that you always avoid hard clipping of reads.","warning")}
    if(! isDefinitelyPairedEnd){ reportln("Warning: Have not found any matched read pairs in the first 1000 reads. Is data paired end? WARNING: This utility is only designed for use on paired end data!","warning"); }
    if(isSortedByPosition){ reportln("SAM/BAM file looks like it might be sorted by position. If so: this mode is not currently supported!","warning"); }
    if(! isSortedByName) error("FATAL ERROR: SAM/BAM file is not sorted by name! Sort the file by name!");
    reportln("SAMRecord Reader Generated. Based on the first 1000 reads, the reads appear to be of length: "+readLength+".","note");
    standardStatusReport(initialTimeStamp);

    val coda : Array[Int] = internalUtils.commonSeqUtils.getNewCauseOfDropArray;
    val coda_options : Array[Boolean] = internalUtils.commonSeqUtils.CODA_DEFAULT_OPTIONS.toArray;
    if(keepMultiMapped) coda_options(internalUtils.commonSeqUtils.CODA_NOT_UNIQUE_ALIGNMENT) = false;
    if(! readGroup.isEmpty) coda_options(internalUtils.commonSeqUtils.CODA_NOT_MARKED_RG) = true;
    
    //"writeKnownSplices","writeNovelSplices","writeSpliceExon", "writeDESeq","writeDEXSeq","writeGenewiseGeneBody"
    //  final val QC_FUNCTION_LIST : Seq[String] = Seq("InsertSize","NVC","CigarOpDistribution","QualityScoreDistribution","GCDistribution","GeneCounts","JunctionCounts");
    val qcGGC:  QCUtility[String] =   if(runFunc.contains("GeneCalcs"))                 new qcGetGeneCounts(stranded,fr_secondStrand,anno_holder,coda,coda_options,40, runFunc.contains("FPKM"), runFunc.contains("writeGenewiseGeneBody"), runFunc.contains("writeDESeq"), runFunc.contains("writeGeneCounts")) else QCUtility.getBlankStringUtil;
    val qcIS :  QCUtility[Unit]   =   if(runFunc.contains("InsertSize"))                new qcInnerDistance(anno_holder, stranded, fr_secondStrand, readLength)        else QCUtility.getBlankUnitUtil;
    val qcCS :  QCUtility[Unit]   =   if(runFunc.contains("NVC"))                       new qcNVC(readLength, runFunc.contains("writeClippedNVC"))                     else QCUtility.getBlankUnitUtil;
    val qcJD :  QCUtility[Unit]   =   if(runFunc.contains("CigarOpDistribution"))       new qcCigarDistribution(readLength)                                            else QCUtility.getBlankUnitUtil;
    val qcQSC : QCUtility[Unit]   =   if(runFunc.contains("QualityScoreDistribution"))  new qcQualityScoreCounter(readLength, qcQualityScoreCounter.MAX_QUALITY_SCORE) else QCUtility.getBlankUnitUtil;
    val qcGC :  QCUtility[Unit]   =   if(runFunc.contains("GCDistribution"))            new qcGCContentCount(readLength)                                               else QCUtility.getBlankUnitUtil;
    val qcJC :  QCUtility[Unit]   =   if(runFunc.contains("JunctionCalcs"))             new qcJunctionCounts(anno_holder, stranded, fr_secondStrand, runFunc.contains("writeDEXSeq"), runFunc.contains("writeSpliceExon"), runFunc.contains("writeKnownSplices"), runFunc.contains("writeNovelSplices"))                   else QCUtility.getBlankUnitUtil;
    val qcST :  QCUtility[Unit]   =   if(runFunc.contains("StrandCheck"))               new qcStrandTest(anno_holder, stranded, fr_secondStrand)                       else QCUtility.getBlankUnitUtil;
    val qcCC :  QCUtility[Unit]   =   if(runFunc.contains("chromCounts"))               new qcChromCount( fr_secondStrand)                                             else QCUtility.getBlankUnitUtil;
    val qcCM :  QCUtility[Unit]   =   if(runFunc.contains("cigarMatch"))                new qcCigarMatch(readLength)                                                   else QCUtility.getBlankUnitUtil;
    
    val qcALL = parConvert(Vector(qcGGC, qcIS, qcCS, qcJD, qcQSC, qcGC, qcJC, qcST, qcCC, qcCM), numThreads);
    
    reportln("QC Utilities Generated!","note");
    standardStatusReport(initialTimeStamp);
    //GenomicArrayOfSets.printGenomicArrayToFile("TEST.OUT.gtf",geneArray);
    var readNum = 0;
    val samIterationTimeStamp = TimeStampUtil();
    for(pair <- pairedIter){
    //for((pair,readNum) <- numberedIter){
      val (r1,r2) = pair;
      readNum += 1;
      
      if(internalUtils.commonSeqUtils.useReadPair(r1,r2,coda, coda_options, dropChrom, readGroup)){
        
        qcALL.foreach( _.runOnReadPair(r1,r2,readNum) );
        //qcGGC.runOnReadPair(r1,r2,readNum);
        //qcCS.runOnReadPair(r1,r2,readNum);
        //qcJD.runOnReadPair(r1,r2,readNum);
        //qcQSC.runOnReadPair(r1,r2,readNum);
        //qcIS.runOnReadPair(r1,r2,readNum);
        //qcGC.runOnReadPair(r1,r2,readNum);
        //qcJC.runOnReadPair(r1,r2,readNum);
        //qcST.runOnReadPair(r1,r2,readNum);
        //qcCC.runOnReadPair(r1,r2,readNum);
      }
    }
    reportln("Finished reading SAM. Read: " + readNum + " read-pairs","note");
    standardStatusReport(initialTimeStamp);
    
    val outputIterationTimeStamp = TimeStampUtil();
    reportln("Read Stats:\n" + internalUtils.commonSeqUtils.causeOfDropArrayToString(coda, coda_options),"note");
    
    val summaryWriter = openWriter(outfile + "summary.txt");
    val strandedCode = if(! stranded){ 0 } else {if(fr_secondStrand) 2; else 1;}

    summaryWriter.write("FIELD	COUNT\n");
    summaryWriter.write("Stranded_Rule_Code	"+strandedCode+"\n");
    
    val iterationMinutes = (outputIterationTimeStamp.compareTo(samIterationTimeStamp) / 1000).toDouble / 60.toDouble;
    val minutesPerMillion = iterationMinutes / (readNum.toDouble / 1000000.toDouble);
    val minutesPerMillionPF = iterationMinutes / ((coda(internalUtils.commonSeqUtils.CODA_READ_PAIR_OK)).toDouble / 1000000.toDouble);
    
    summaryWriter.write("BENCHMARK_MinutesOnSamIteration	" + "%1.2f".format(iterationMinutes) + "\n");
    summaryWriter.write("BENCHMARK_MinutesPerMillionReads	" + "%1.2f".format(minutesPerMillion) + "\n");
    summaryWriter.write("BENCHMARK_MinutesPerMillionGoodReads	" + "%1.2f".format(minutesPerMillionPF) + "\n");
    
    reportln("Writing Output...","note");
    qcALL.seq.foreach( _.writeOutput(outfile, summaryWriter) );
    //qcGGC.writeOutput(outfile, summaryWriter);
    //qcCS.writeOutput(outfile, summaryWriter);
    //qcJD.writeOutput(outfile, summaryWriter);
    //qcQSC.writeOutput(outfile, summaryWriter);
    //qcIS.writeOutput(outfile, summaryWriter);
    //qcGC.writeOutput(outfile, summaryWriter);
    //qcJC.writeOutput(outfile, summaryWriter);
    //qcST.writeOutput(outfile, summaryWriter);
    //qcCC.writeOutput(outfile, summaryWriter);
    reportln("Done.","note");
    
    summaryWriter.write("COMPLETED_WITHOUT_ERROR	1\n");
    
    close(summaryWriter);
    
    if(runFunc.contains("makeAllBrowserTracks")){
      reportln("Making (optional) browser tracks...","note");
      //TO DO!
      reportln("Done making browser tracks.","note");
    }
    
    
    val finalTimeStamp = TimeStampUtil();
    reportln("Time spent on setup:           " + TimeStampUtil.timeDifferenceFormatter(samIterationTimeStamp.compareTo(initialTimeStamp)),"note");
    reportln("Time spent on SAM iteration:   " + TimeStampUtil.timeDifferenceFormatter(outputIterationTimeStamp.compareTo(samIterationTimeStamp)),"note");
    reportln("                               (" + minutesPerMillion + " minutes per million read-pairs)","note");
    reportln("                               (" + minutesPerMillionPF + " minutes per million read-pairs used)","note");
    reportln("Time spent on file output:     " + TimeStampUtil.timeDifferenceFormatter(finalTimeStamp.compareTo(outputIterationTimeStamp)),"note");
    reportln("Total runtime:                 " + TimeStampUtil.timeDifferenceFormatter(finalTimeStamp.compareTo(initialTimeStamp)),"note");
  }
  /*
   * SETTING GLOBAL OPTIONS:
   */
  def setQcOptions(noGzipOutput : Boolean){
    //registerGlobalParam[Boolean]("noGzipOutput", noGzipOutput)
    internalUtils.optionHolder.OPTION_noGzipOutput = noGzipOutput;
  }
  
}















