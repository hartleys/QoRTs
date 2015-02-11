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
  
  final val QC_DEFAULT_ON_FUNCTION_LIST  : scala.collection.immutable.Set[String] = scala.collection.immutable.Set("InsertSize",
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
  final val QC_DEFAULT_OFF_FUNCTION_LIST : scala.collection.immutable.Set[String] = scala.collection.immutable.Set("FPKM","makeWiggles","makeJunctionBed","makeAllBrowserTracks", "cigarMatch");
  final val QC_FUNCTION_LIST : scala.collection.immutable.Set[String] = QC_DEFAULT_ON_FUNCTION_LIST ++ QC_DEFAULT_OFF_FUNCTION_LIST;
  final val COMPLETED_OK_FILENAME = ".QORTS_COMPLETED_OK";
  final val COMPLETED_WARN_FILENAME = ".QORTS_COMPLETED_WARN";
  final val MASTERLEVEL_FUNCTION_LIST = List[String]("GeneCalcs", "InsertSize","NVC","CigarOpDistribution","QualityScoreDistribution","GCDistribution","JunctionCalcs","StrandCheck","chromCounts","cigarMatch","makeWiggles");
  
  final val QC_INCOMPATIBLE_WITH_SINGLE_END_FUNCTION_LIST : scala.collection.immutable.Set[String] = scala.collection.immutable.Set(
      "InsertSize","cigarMatch"
  );
  
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
      ("writeClippedNVC" -> "NVC"),
      ("makeAllBrowserTracks" -> "makeWiggles"),
      ("makeAllBrowserTracks" -> "makeJunctionBed")
  );
  
  
  //def run(args : Array[String]){
  
  class allQC_runner extends CommandLineRunUtil {
    override def priority = 1;
    val parser : CommandLineArgParser = 
      new CommandLineArgParser(
          command = "QC", 
          quickSynopsis = "Runs a battery of QC tools", 
          synopsis = "", 
          description = "This utility runs a large battery of QC / data processing tools on a single given sam or bam file."+
                        "This is the primary function of the QoRT utility."+
                        ""+
                        ""+
                        ""+
                        "",
          argList = 


                    new UnaryArgument(   name = "singleEnded", 
                                         arg = List("--singleEnded","-e"), // name of value
                                         argDesc = "Flag to indicate that reads are single end." // description
                                       ) ::
                    new UnaryArgument( name = "coordSorted",
                                         arg = List("--coordSorted"), // name of value
                                         argDesc = "Flag to indicate that input bam file is coordinate-sorted, rather than name-sorted. "+
                                                   "Note that QoRTs will take longer to run and use more memory in this mode. To improve performance, sort the data by name prior to using of QoRTs. "+
                                                   "In addition, if an (extremely) large fraction of the read-pairs are "+
                                                   "mapped to extremely distant loci (or different chromosomes), then memory issues may arise. However, this should not be a problem with most datasets. "+
                                                   "Technically this function will also allow QoRTs to work on unsorted bam files, but this is STRONGLY not recommended, as memory usage will by greatly increased." // description
                                       ) ::
                    new UnaryArgument( name = "stranded",
                                         arg = List("--stranded","-s"), // name of value
                                         argDesc = "Flag to indicate that data is stranded." // description
                                       ) ::
                    new UnaryArgument( name = "fr_secondStrand",
                                         arg = List("--stranded_fr_secondstrand","-a"), // name of value
                                         argDesc = "Flag to indicate that reads are from a fr_secondstrand type of stranded library (equivalent to the \"stranded = yes\" option in HTSeq or the \"fr_secondStrand\" library-type option in TopHat/CuffLinks). "+
                                                   "If your data is stranded, you must know the library type in order to analyze it properly. This utility uses the same "+
                                                   "definitions as cufflinks to define strandedness type. By default, the fr_firststrand "+
                                                   "library type is assumed for all stranded data (equivalent to the \"stranded = reverse\" option in HTSeq)." // description
                                       ) ::
                    new BinaryOptionArgument[Int](
                                         name = "maxReadLength", 
                                         arg = List("--maxReadLength"), 
                                         valueName = "len",
                                         argDesc =  "Sets the maximum read length. For unclipped datasets this option is not necessary since the read length can be determined from the data. "+
                                                    "By default, QoRTs will attempt to determine the max read length by examining the first 1000 reads. "+
                                                    "If your data is hard-clipped prior to alignment, then it is strongly recommended that this option be included, or else an error may occur. "+
                                                    "Note that hard-clipping data prior to alignment is generally not recommended, because this makes it difficult (or impossible) "+
                                                    "to determine the sequencer read-cycle of each nucleotide base. This may obfuscate cycle-specific artifacts, trends, or errors, the detection of which is one of the primary purposes of QoRTs! "+
                                                    "In addition, hard clipping (whether before or after alignment) removes quality score data, and thus quality score metrics may be misleadingly optimistic. "+
                                                    "A MUCH preferable method of removing undesired sequence is to replace such sequence with N's, which preserves the quality score and the sequencer cycle information while still removing undesired sequence. "+
                                                    ""+
                                                    ""
                                        ) ::

                    new UnaryArgument(    name = "testRun",
                                         arg = List("--testRun","-t"), // name of value
                                         argDesc = "Flag to indicate that only the first 100k reads should be read in. Used for testing." // description
                                       ) ::

                    new UnaryArgument( name = "keepMultiMapped",
                                         arg = List("--keepMultiMapped"), // name of value
                                         argDesc = "Flag to indicate that the tool should NOT filter out multi-mapped reads. Note that even with this flag raised this utility will still only "+
                                                    "use the 'primary' alignment location for each read. By default any reads that are marked as multi-mapped will be ignored entirely."+
                                                    " Most aligners use the MAPQ value to mark multi-mapped reads. Any read with MAPQ < 255 is assumed to be non-uniquely mapped. "+
                                                    " Thus: this option is equivalent to setting --minMAPQ to 0."// description
                                       ) ::

                    //new UnaryArgument( name = "neverPartiallyUnpaired",
                    //                     arg = List("--neverPartiallyUnpaired"), // name of value
                    //                     argDesc = "Flag to indicate that when a read's mate is unmapped, it always appears in the input file. This option is always optional, but will improve performance." // description
                    //                  ) ::
                    new UnaryArgument( name = "noGzipOutput",
                                         arg = List("--noGzipOutput"), // name of value
                                         argDesc = "Flag to indicate that output files should NOT be compressed into the gzip format. By default almost all output files are compressed to save space." // description
                                       ) ::

                    new BinaryOptionArgument[String](
                                         name = "readGroup", 
                                         arg = List("--readGroup"), 
                                         valueName = "readGroupName",  
                                         argDesc =  "If this option is set, all analyses will be restricted to ONLY reads that are tagged with the given "+
                                                    "readGroupName (using an RG tag). This can be used if multiple read-groups have already been combined "+
                                                    "into a single bam file, but you want to summarize each read-group separately."
                                        ) ::

                    new BinaryArgument[Int](name = "minMAPQ",
                                           arg = List("--minMAPQ"),  
                                           valueName = "num", 
                                           argDesc = "Filter out reads with less than the given MAPQ. Set to 0 to turn off mapq filtering.", 
                                           defaultValue = Some(255)
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
                                                        argDesc = "The complete list of functions to run. Setting this option turns off ALL functions EXCEPT for the ones explicitly requested here. Some functions require other functions. If these functions are requested, all functions it is dependent on will also run. Followed by a comma delimited list, with no internal whitespace. Allowed options are: "+QC_FUNCTION_LIST.mkString(", "), 
                                                        defaultValue = Some(List[String]())
                                                        ) :: 
                    new BinaryOptionArgument[Int](
                                         name = "seqReadCt", 
                                         arg = List("--seqReadCt"), 
                                         valueName = "val",  
                                         argDesc = "(Optional) The number of reads for the replicate, prior to alignment."+
                                                   "This will be passed on into the QC.summary.txt file."+
                                                   ""+
                                                   ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "rawfastq", 
                                         arg = List("--rawfastq"), 
                                         valueName = "myfastq.fq.gz",  
                                         argDesc = "(Optional) The raw fastq, prior to alignment. This is used ONLY to calculate the number of pre-alignment reads(or read-pairs) simply by counting the number of lines and dividing by 4. "+
                                                   "The number of pre-alignment read-pairs can be included explicitly via the --seqReadCt option, or added in the "+
                                                   "plotting / cross-comparison step by including the input.read.pair.count column in the replicate decoder."+
                                                   "In general, the --seqReadCt option is recommended when possible.\n"+
                                                   "If the filename ends with \".gz\" or \".zip\", the file will be parsed using the appropriate decompression method."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "chromSizes", 
                                         arg = List("--chromSizes"), 
                                         valueName = "chrom.sizes.txt",  
                                         argDesc = "A chrom.sizes file. The first (tab-delimited) column must contain all chromosomes found in the dataset. "+
                                                   "The second column must contain chromosome sizes (in base-pairs). If a standard genome is being used, it is strongly recommended that this be generated by "+
                                                   "the UCSC utility 'fetchChromSizes'.\n"+
                                                   "This file is ONLY needed to produce wiggle files. If this is provided, then by default QoRTs will produce 100-bp-window wiggle files (and junction '.bed' files) for the supplied data."+
                                                   "In order to produce wiggle files, this parameter is REQUIRED."
                                        ) ::
                    new BinaryArgument[String](name = "trackTitlePrefix",
                                           arg = List("--trackTitlePrefix"),  
                                           valueName = "titlePrefix", 
                                           argDesc = "The prefix used for the track name in the track definition line of any browser tracks ('.wig' or '.bed' files) generated by this utility."+
                                                     "Note that no browser tracks will be created by default, unless the '--chromSizes' option is set. Bed files can also be generated using the option '--addFunction makeJunctionBed'", 
                                           defaultValue = Some("UntitledTrack")
                                           ) :: 
                    new BinaryOptionArgument[String](
                                         name = "flatgfffile", 
                                         arg = List("--flatgff"), 
                                         valueName = "flattenedGffFile.gff.gz",  
                                         argDesc = "A \"flattened\" gtf file that matches the standard gtf file. Optional."+
                                                   "It may also be useful for downstream analyses, as it assigns unique identifiers to all exons and splice "+
                                                   "junctions. The flattened gtf file can be generated using "+
                                                   "the \"makeFlatGff\" command. Note that the command must be run with the same strandedness code.\n"+
                                                   "If the filename ends with \".gz\" or \".zip\", the file will be parsed using the appropriate decompression method."
                                        ) ::

                      new BinaryOptionArgument[String](
                                         name = "restrictToGeneList", 
                                         arg = List("--restrictToGeneList"), 
                                         valueName = "geneList.txt",  
                                         argDesc =  "If this option is set, almost all analyses will be restricted to reads that are found on genes named in the "+
                                                    "supplied gene list file. The file should contain a gene ID on each line and nothing else. "+
                                                    "The only functions that will be run on the full set of all reads will be the functions that calculate "+
                                                    "the gene mapping itself. NOTE: if you want to include ambiguous reads, include a line with the text: '_ambiguous'. "+
                                                    "If you want to include reads that do not map to any known feature, include a line with the text: '_no_feature'. "+
                                                    "WARNING: this is not intended for default use. It is intended to be used when re-running QoRTs, with the intention of "+
                                                    "examining artifacts that can be caused in various plots by a small number of genes with extremely high coverage. For example, "+
                                                    "GC content plots sometimes contain visible spikes caused by small mitochondrial genes with extremely high expression."+
                                                    "ADDITIONAL WARNING: This feature is in BETA, and is not yet fully tested."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "dropGeneList", 
                                         arg = List("--dropGeneList"), 
                                         valueName = "geneList.txt",  
                                         argDesc =  "If this option is set, almost all analyses will be restricted to reads that are NOT found on genes named in the "+
                                                    "supplied gene list file. The file should contain a gene ID on each line and nothing else. "+
                                                    "The only functions that will be run on the full set of all reads will be the functions that calculate "+
                                                    "the gene mapping itself. NOTE: if you want to EXCLUDE ambiguous reads, include a line with the text: '_ambiguous'. "+
                                                    "If you want to EXCLUDE reads that do not map to any known feature, include a line with the text: '_no_feature'. "+
                                                    "WARNING: this is not intended for default use. It is intended to be used when re-running QoRTs, with the intention of "+
                                                    "examining artifacts that can be caused by certain individual 'problem genes'. For example, "+
                                                    "GC content plots sometimes contain visible spikes caused by small transcripts / RNA's with extremely high expression levels."+
                                                    "ADDITIONAL WARNING: This feature is in BETA, and is not yet fully tested."
                                        ) ::                                      
                                        //DEPRECIATED OPTIONS:
                    new UnaryArgument( name = "noMultiMapped",
                                         arg = List("--fileContainsNoMultiMappedReads"), // name of value
                                         argDesc = "Flag to indicate that the input sam/bam file contains only primary alignments (ie, no multi-mapped reads). This flag is ALWAYS OPTIONAL, but when applicable this utility will run (slightly) faster when using this argument. (DEPRECIATED! The performance improvement was marginal)" // description
                                       ) ::
                    new UnaryArgument( name = "parallelFileRead",
                                         arg = List("--parallelFileRead"), // name of value
                                         argDesc = "DEPRECIATED: DO NOT USE. Flag to indicate that bam file reading should be run in paralell for increased speed. Note that in this mode you CANNOT read from stdin. Also note that for this to do anything useful, the numThreads option must be set to some number greater than 1. Also note that additional threads above 9 will have no appreciable affect on speed." // description
                                       ) ::    
                    new BinaryArgument[Int](name = "numThreads",
                                                        arg = List("--numThreads"),  
                                                        valueName = "num", 
                                                        argDesc = "DEPRECIATED, nonfunctional.", 
                                                        defaultValue = Some(1)
                                                        ) :: 
//MANDATORY OPTIONS:
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
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS
      );
    
    def run(args : Array[String]){
      val out = parser.parseArguments(args.toList.tail);
      
      if(out){
      runAllQC.run(parser.get[String]("infile") ,
          parser.get[String]("outfile") + "/" + "QC",
          parser.get[String]("gtffile"),
          parser.get[Option[String]]("flatgfffile"),
          parser.get[List[String]]("dropChromList"),
          parser.get[Boolean]("singleEnded"),
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
          parser.get[Option[String]]("readGroup"),
          parser.get[Boolean]("parallelFileRead"),
          parser.get[Int]("minMAPQ"),
          parser.get[Option[String]]("restrictToGeneList"),
          parser.get[Option[String]]("dropGeneList"),
          parser.get[Boolean]("coordSorted"),
          parser.get[Option[Int]]("maxReadLength"),
          parser.get[Option[Int]]("seqReadCt"),
          parser.get[Option[String]]("rawfastq"),
          parser.get[Option[String]]("chromSizes"),
          parser.get[String]("trackTitlePrefix")
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
          readGroup : Option[String],
          parallelFileRead : Boolean,
          minMAPQ: Int,
          restrictToGeneList : Option[String],
          dropGeneList : Option[String],
          unsorted : Boolean,
          maxReadLength : Option[Int],
          seqReadCt : Option[Int],
          rawfastq : Option[String],
          chromSizes : Option[String],
          trackTitlePrefix : String){

    internalUtils.Reporter.init_completeLogFile(outfile + ".log");
    
    val gigsMaxMem = getMaxMemoryXmxInGigs;
    if(gigsMaxMem < 4){
      reportln("NOTE: maximum allocation memory = " + gigsMaxMem + " gigabytes.\n"+
               "    This might be ok, or might cause OutOfMemoryExceptions later.\n"+
               "    For most large datasets/genomes at least 4 gb is recommended.\n"+
               "    (Actual required memory may be less than this.)\n"+
               "    To increase the memory maximum, include the parameter -Xmx\n"+
               "    in between the java command and the -jar parameter.\n"+
               "    For example: to increase the memory maximum to 4 gigabytes:\n"+
               "        java -Xmx4G -jar /path/to/jar/QoRTs.jar QC ...","note");
    }
    
    reportln("Starting QC","note");
    val initialTimeStamp = TimeStampUtil();
    standardStatusReport(initialTimeStamp);

    val COMPLETED_OK_FILEPATH = outfile + COMPLETED_OK_FILENAME;
    val COMPLETED_OK_OLDFILE = new java.io.File(COMPLETED_OK_FILEPATH);
    if(COMPLETED_OK_OLDFILE.exists()){
      report("Deleting old \"QORTS_COMPLETED_OK\" file.","note");
      COMPLETED_OK_OLDFILE.delete();
    }
    
    val stdGtfCodes = new internalUtils.GtfTool.GtfCodes();
    val flatGtfCodes = new internalUtils.GtfTool.GtfCodes();
    
    if(isSingleEnd) reportln("QoRTs is Running in single-end mode.","note");
    else reportln("QoRTs is Running in paired-end mode.","note");
    
    val dropChrom = dropChromList.toSet;
    
    val defaultFunctonList = if(chromSizes.isEmpty){
      QC_DEFAULT_ON_FUNCTION_LIST;
    } else {
      QC_DEFAULT_ON_FUNCTION_LIST ++ Set("makeWiggles");
    }
    
    val runFunc_initial = if(runFunctions.isEmpty){
      (defaultFunctonList ++ addFunctions.toSet) -- dropFunctions.toSet;
    } else {
      runFunctions.toSet -- dropFunctions.toSet;
    }
    
    val runFuncTEMP : scala.collection.immutable.Set[String] = (if(restrictToGeneList.isEmpty & dropGeneList.isEmpty) scala.collection.immutable.Set[String]() else scala.collection.immutable.Set[String]("GeneCalcs"));
    val runFuncTEMP2 = (runFunc_initial ++ runFuncTEMP).foldLeft(runFunc_initial ++ runFuncTEMP)((soFar,currFunc) => {
      QC_FUNCTION_DEPENDANCIES.get(currFunc) match {
        case Some(reqFunc) => {
          if(soFar.contains(reqFunc)){
            soFar;
          } else {
            reportln("Function \"" + currFunc +"\" requires top-level function \"" + reqFunc+ "\". Adding the required function to the active function list.","note");
            soFar + reqFunc;
          }
        }
        case None => {
          soFar;
        }
      }
    });
    val runFunc = if(! isSingleEnd) runFuncTEMP2 else {
      runFuncTEMP2 -- QC_INCOMPATIBLE_WITH_SINGLE_END_FUNCTION_LIST;
    }
    
    val notFoundFunc = runFunc.find(  ! QC_FUNCTION_LIST.contains(_) ) ;
    if( ! notFoundFunc.isEmpty ){
      error("ERROR: Function not found: \"" + notFoundFunc.get +"\"");
    }
    
    if(runFunc.contains("makeWiggles") && chromSizes.isEmpty){
      error("Error: function makeWiggles REQUIRES the chromSizes parameter! Set parameter '--chromSizes' or turn off function 'makeWiggles'.");
    }
    
    
    //wrapSimpleLineWithIndent_staggered(line : String, width : Int, indent : String, firstLineIndent : String)
    reportln("Running functions: " + wrapSimpleLineWithIndent_staggered(runFunc.mkString(", "), 68, "        ", ""),"progress");
    
    //reportln("infile: " + infile , "note");
    //reportln("outfile: " + outfile , "note");
    //reportln("gtffile: " + gtffile , "note");
    //reportln("flatgtffile: " + flatgtffile , "note");
    //reportln("stranded: " + stranded , "note");
    //reportln("fr_secondStrand: " + fr_secondStrand , "note");
    //reportln("testRun: " + testRun , "note");
    //reportln("dropChrom: "+dropChromList.mkString(";"),"note");
    //reportln("run_functions: " + runFunc.mkString(";"),"note");
    
    val geneListKeep : Option[Set[String]] = restrictToGeneList match {
      case Some(glkf) => {
        Some(internalUtils.fileUtils.getLinesSmartUnzip(glkf).toSet[String]);
      }
      case None => None;
    }
    val geneListDrop : Option[Set[String]] = dropGeneList match {
      case Some(glkf) => {
        Some(internalUtils.fileUtils.getLinesSmartUnzip(glkf).toSet[String]);
      }
      case None => None;
    }
    val geneKeepFunc : (String => Boolean) = if(geneListKeep.isEmpty && geneListDrop.isEmpty){
      (g : String) => true;
    } else {
      (g : String) => {
        (geneListKeep.isEmpty || (  geneListKeep.get.contains(g))) && 
        (geneListDrop.isEmpty || (! geneListDrop.get.contains(g)))
      }
    }
    
    setQcOptions(noGzipOutput);
    val anno_holder = new qcGtfAnnotationBuilder(gtffile , flatgtffile , stranded , stdGtfCodes, flatGtfCodes);
    
    if(parallelFileRead){
      reportln("ERROR ERROR ERROR: parallell file read is NOT IMPLEMENTED AT THIS TIME!","warn");
      //runOnSeqFile_PAR(initialTimeStamp = initialTimeStamp, infile = infile, outfile = outfile, anno_holder = anno_holder, testRun = testRun, runFunc = runFunc, stranded = stranded, fr_secondStrand = fr_secondStrand, dropChrom = dropChrom, keepMultiMapped = keepMultiMapped, noMultiMapped = noMultiMapped, numThreads = numThreads, readGroup )
    } else {
      runOnSeqFile(initialTimeStamp = initialTimeStamp, infile = infile, outfile = outfile, anno_holder = anno_holder, testRun = testRun, 
          runFunc = runFunc, stranded = stranded, fr_secondStrand = fr_secondStrand, dropChrom = dropChrom, keepMultiMapped = keepMultiMapped, 
          noMultiMapped = noMultiMapped, numThreads = numThreads, readGroup , minMAPQ = minMAPQ, geneKeepFunc = geneKeepFunc, 
          isSingleEnd = isSingleEnd, unsorted = unsorted, maxReadLength = maxReadLength,
          seqReadCt = seqReadCt, rawfastq = rawfastq,
          chromSizes = chromSizes, trackTitlePrefix = trackTitlePrefix)
    }
  }
  
  def runOnSeqFile(initialTimeStamp : TimeStampUtil, 
                   infile : String, 
                   outfile : String, 
                   anno_holder : qcGtfAnnotationBuilder, 
                   testRun  : Boolean, 
                   runFunc : Set[String], 
                   stranded : Boolean, 
                   fr_secondStrand : Boolean, 
                   dropChrom : Set[String], 
                   keepMultiMapped : Boolean, 
                   noMultiMapped : Boolean, 
                   numThreads : Int,
                   readGroup : Option[String],
                   minMAPQ : Int,
                   geneKeepFunc : (String => Boolean),
                   isSingleEnd : Boolean,
                   unsorted : Boolean,
                   maxReadLength : Option[Int],
                   seqReadCt : Option[Int],
                   rawfastq : Option[String],
                   chromSizes : Option[String],
                   trackTitlePrefix : String){
    
    val inputReadCt : Option[Int] = seqReadCt match {
      case Some(ct) => Some(ct);
      case None => {
        rawfastq match {
          case Some(fqfile) => Some(internalUtils.fileUtils.getLinesSmartUnzip(fqfile).length / 4);
          case None => None;
        }
      }
    }
      
      
     /* if(! rawfastq.isEmpty && (seqReadCt.isEmpty)){
      reportln("> Reading fastq file to determine unaligned read count","note");
      val fastqLineCt = internalUtils.fileUtils.getLinesSmartUnzip(rawfastq.get).length;
      
      reportln("> done.","note");
    }*/
    
    
    val peekCt = 2000;
    val COMPLETED_OK_FILEPATH = outfile + COMPLETED_OK_FILENAME;
    val COMPLETED_WARN_FILEPATH = outfile + COMPLETED_WARN_FILENAME;
    val (samFileAttributes, recordIter) = initSamRecordIterator(infile, peekCt);
    
    val pairedIter : Iterator[(SAMRecord,SAMRecord)] = 
      if(isSingleEnd){
        if(testRun) samRecordPairIterator_withMulti_singleEnd(recordIter, true, 200000) else samRecordPairIterator_withMulti_singleEnd(recordIter);
      } else {
        if(unsorted){
          if(testRun) samRecordPairIterator_unsorted(recordIter, true, 200000) else samRecordPairIterator_unsorted(recordIter)
        // Faster noMultiMapped running is DEPRECIATED!
        //} else if(noMultiMapped){
        //  if(testRun) samRecordPairIterator(recordIter, true, 200000) else samRecordPairIterator(recordIter)
        } else {
          if(testRun) samRecordPairIterator_withMulti(recordIter, true, 200000) else samRecordPairIterator_withMulti(recordIter)
        }
      } 
    
    val maxObservedReadLength = samFileAttributes.readLength;
    val readLength = if(maxReadLength.isEmpty) maxObservedReadLength else maxReadLength.get;
    val isSortedByName = samFileAttributes.isSortedByName;
    val isSortedByPosition = samFileAttributes.isSortedByPosition;
    val isDefinitelyPairedEnd = samFileAttributes.isDefinitelyPairedEnd;
    val minReadLength = samFileAttributes.minReadLength;
    
    if(readLength != minReadLength){reportln("Note: Read length is not consistent. "+
                                             "In the first "+peekCt+" reads, read length varies from "+minReadLength+" to " +maxObservedReadLength+"!\n"+
                                             "Note that using data that is hard-clipped prior to alignment is NOT recommended, because this makes it difficult (or impossible) "+
                                             "to determine the sequencer read-cycle of each nucleotide base. This may obfuscate cycle-specific artifacts, trends, or errors, the detection of which is one of the primary purposes of QoRTs!"+
                                             "In addition, hard clipping (whether before or after alignment) removes quality score data, and thus quality score metrics may be misleadingly optimistic. "+
                                             "A MUCH preferable method of removing undesired sequence is to replace such sequence with N's, which preserves the quality score and the sequencer cycle information.","note")}
    if((readLength != minReadLength) & (maxReadLength.isEmpty)){
      reportln("WARNING WARNING WARNING: Read length is not consistent, AND \"--maxReadLength\" option is not set!\n"+
               "QoRTs has ATTEMPTED to determine the maximum read length ("+readLength+").\n"+
               "It is STRONGLY recommended that you use the --maxReadLength option \n"+
               "to set the maximum possible read length, or else errors may occur if/when \n"+
               "reads longer than "+readLength+ " appear.","warn")
    }
     
    if(samFileAttributes.allReadsMarkedPaired & isSingleEnd) reportln("WARNING WARNING WARNING! Running in single-end mode, but reads appear to be paired-end! Errors may follow.\nStrongly recommend removing the '--isSingleEnd' option!","warn");
    if(samFileAttributes.allReadsMarkedSingle & (! isSingleEnd)) reportln("WARNING WARNING WARNING! Running in paired-end mode, but reads appear to be single-end! Errors may follow.\nStrongly recommend using the '--isSingleEnd' option","warn");
    if(samFileAttributes.mixedSingleAndPaired) reportln("WARNING WARNING WARNING! Data appears to be a mixture of single-end and paired-end reads! QoRTs was not designed to function under these conditions. Errors may follow!","warn");
    
    if(! isSingleEnd){
      if( (! isDefinitelyPairedEnd)){ reportln("Warning: Have not found any matched read pairs in the first "+peekCt+" reads. Is data paired-end? Use option --singleEnd for single-end data.","warn"); }
      if( isSortedByPosition & (! unsorted )){ reportln("Based on the first "+peekCt+" reads, SAM/BAM file appears to be sorted by read position. If this is so, you should probably use the \"--coordSorted\" option.","warn"); }
      if( ((! isSortedByPosition) & ( unsorted ))){ reportln("WARNING: You are using the \"--coordSorted\" option, but data does NOT appear to be sorted by read position (based on the first "+peekCt+" reads)! This is technically ok, but may cause QoRTs to use too much memory!","warn"); }
      if( ((! isSortedByName) & (! unsorted ))) error("FATAL ERROR: SAM/BAM file is not sorted by name (based on the first "+peekCt+" reads)! Either sort the file by name, or sort by read position and use the \"--coordSorted\" option.");
    }
    
    reportln("SAMRecord Reader Generated. Read length: "+readLength+".","note");
    standardStatusReport(initialTimeStamp);

    val coda : Array[Int] = internalUtils.commonSeqUtils.getNewCauseOfDropArray;
    val coda_options : Array[Boolean] = internalUtils.commonSeqUtils.CODA_DEFAULT_OPTIONS.toArray;
    if(isSingleEnd) CODA_SINGLE_END_OFF_OPTIONS.foreach( coda_options(_) = false );
    if(keepMultiMapped) coda_options(internalUtils.commonSeqUtils.CODA_NOT_UNIQUE_ALIGNMENT) = false;
    if(! readGroup.isEmpty) coda_options(internalUtils.commonSeqUtils.CODA_NOT_MARKED_RG) = true;
    
    //"writeKnownSplices","writeNovelSplices","writeSpliceExon", "writeDESeq","writeDEXSeq","writeGenewiseGeneBody"
    //  final val QC_FUNCTION_LIST : Seq[String] = Seq("InsertSize","NVC","CigarOpDistribution","QualityScoreDistribution","GCDistribution","GeneCounts","JunctionCounts");
    val qcGGC:  QCUtility[String] =   if(runFunc.contains("GeneCalcs"))                 new qcGetGeneCounts(stranded,fr_secondStrand,anno_holder,coda,coda_options,40, runFunc.contains("FPKM"), runFunc.contains("writeGenewiseGeneBody"), runFunc.contains("writeDESeq"), runFunc.contains("writeGeneCounts"), geneKeepFunc) else QCUtility.getBlankStringUtil;
    val qcIS :  QCUtility[Unit]   =   if(runFunc.contains("InsertSize"))                new qcInnerDistance(anno_holder, stranded, fr_secondStrand, readLength)        else QCUtility.getBlankUnitUtil;
    val qcCS :  QCUtility[Unit]   =   if(runFunc.contains("NVC"))                       new qcNVC(isSingleEnd, readLength, runFunc.contains("writeClippedNVC"))                     else QCUtility.getBlankUnitUtil;
    val qcJD :  QCUtility[Unit]   =   if(runFunc.contains("CigarOpDistribution"))       new qcCigarDistribution(isSingleEnd, readLength)                                            else QCUtility.getBlankUnitUtil;
    val qcQSC : QCUtility[Unit]   =   if(runFunc.contains("QualityScoreDistribution"))  new qcQualityScoreCounter(isSingleEnd, readLength, qcQualityScoreCounter.MAX_QUALITY_SCORE) else QCUtility.getBlankUnitUtil;
    val qcGC :  QCUtility[Unit]   =   if(runFunc.contains("GCDistribution"))            new qcGCContentCount(isSingleEnd, readLength)                                               else QCUtility.getBlankUnitUtil;
    val qcJC :  QCUtility[Unit]   =   if(runFunc.contains("JunctionCalcs"))             new qcJunctionCounts(anno_holder, stranded, fr_secondStrand, runFunc.contains("writeDEXSeq"), runFunc.contains("writeSpliceExon"), runFunc.contains("writeKnownSplices"), runFunc.contains("writeNovelSplices"))                   else QCUtility.getBlankUnitUtil;
    val qcST :  QCUtility[Unit]   =   if(runFunc.contains("StrandCheck"))               new qcStrandTest(isSingleEnd, anno_holder, stranded, fr_secondStrand)                       else QCUtility.getBlankUnitUtil;
    val qcCC :  QCUtility[Unit]   =   if(runFunc.contains("chromCounts"))               new qcChromCount(isSingleEnd, fr_secondStrand)                                              else QCUtility.getBlankUnitUtil;
    val qcCM :  QCUtility[Unit]   =   if(runFunc.contains("cigarMatch"))                new qcCigarMatch(readLength)                                                   else QCUtility.getBlankUnitUtil;
    val qcWIG : QCUtility[Unit]   =   if(runFunc.contains("makeWiggles"))               new fileConversionUtils.bamToWiggle.QcBamToWig(trackTitlePrefix,
                                                                                                                   chromSizes.get,false,100,
                                                                                                                   isSingleEnd, stranded, fr_secondStrand, 
                                                                                                                   1.0, true, true, true, None, "") else QCUtility.getBlankUnitUtil;
    
    val qcALL = parConvert(Vector(qcGGC, qcIS, qcCS, qcJD, qcQSC, qcGC, qcJC, qcST, qcCC, qcCM), numThreads);
    
    reportln("QC Utilities Generated!","note");
    standardStatusReport(initialTimeStamp);
    //GenomicArrayOfSets.printGenomicArrayToFile("TEST.OUT.gtf",geneArray);
    var readNum = 0;
    var keptMultiMappedCt = 0;
    val samIterationTimeStamp = TimeStampUtil();
    for(pair <- pairedIter){
    //for((pair,readNum) <- numberedIter){
      val (r1,r2) = pair;
      readNum += 1;
      
      if(internalUtils.commonSeqUtils.useReadPair(r1,r2,coda, coda_options, dropChrom, readGroup, minMAPQ)){
          val gene = qcGGC.runOnReadPair(r1,r2,readNum);
          if( geneKeepFunc(gene) ){
          //if( geneListKeep.get.contains(gene) ){
            qcIS.runOnReadPair(r1,r2,readNum);
            qcCS.runOnReadPair(r1,r2,readNum);
            qcJD.runOnReadPair(r1,r2,readNum);
            qcQSC.runOnReadPair(r1,r2,readNum);
            qcGC.runOnReadPair(r1,r2,readNum);
            qcJC.runOnReadPair(r1,r2,readNum);
            qcST.runOnReadPair(r1,r2,readNum);
            qcCC.runOnReadPair(r1,r2,readNum);
            qcCM.runOnReadPair(r1,r2,readNum);
            qcWIG.runOnReadPair(r1,r2,readNum);
            
            if(internalUtils.commonSeqUtils.isReadMultiMapped(r1) || internalUtils.commonSeqUtils.isReadMultiMapped(r2)){
              keptMultiMappedCt += 1;
            }
          }
      }
    }
    
    reportln("Finished reading SAM. Read: " + readNum + " read-pairs","note");
    standardStatusReport(initialTimeStamp);
    
    val outputIterationTimeStamp = TimeStampUtil();
    report("> Read Stats:\n" + stripFinalNewline(indentifyLines(internalUtils.commonSeqUtils.causeOfDropArrayToString(coda, coda_options),">   ")),"note");
    
    val summaryWriter = openWriter(outfile + ".summary.txt");
    val strandedCode = if(! stranded){ 0 } else {if(fr_secondStrand) 2; else 1;}

    summaryWriter.write("FIELD	COUNT\n");
    summaryWriter.write("Stranded_Rule_Code	"+strandedCode+"\n");
    summaryWriter.write(internalUtils.commonSeqUtils.causeOfDropArrayToStringTabbed(coda, coda_options));
    
    summaryWriter.write("KEPT_NOT_UNIQUE_ALIGNMENT	"+keptMultiMappedCt+"\n");
    
    if(isSingleEnd){
      summaryWriter.write("IS_SINGLE_END	1\n");
    } else {
      summaryWriter.write("IS_SINGLE_END	0\n");
    }
    
    if(inputReadCt.isEmpty){
      summaryWriter.write("PREALIGNMENT_READ_CT	-1\n");
    } else {
      summaryWriter.write("PREALIGNMENT_READ_CT	"+inputReadCt.get+"\n");
    }
    
    
    
    val iterationMinutes = (outputIterationTimeStamp.compareTo(samIterationTimeStamp) / 1000).toDouble / 60.toDouble;
    val minutesPerMillion = iterationMinutes / (readNum.toDouble / 1000000.toDouble);
    val minutesPerMillionPF = iterationMinutes / ((coda(internalUtils.commonSeqUtils.CODA_READ_PAIR_OK)).toDouble / 1000000.toDouble);
    
    summaryWriter.write("BENCHMARK_MinutesOnSamIteration	" + "%1.2f".format(iterationMinutes) + "\n");
    summaryWriter.write("BENCHMARK_MinutesPerMillionReads	" + "%1.2f".format(minutesPerMillion) + "\n");
    summaryWriter.write("BENCHMARK_MinutesPerMillionGoodReads	" + "%1.2f".format(minutesPerMillionPF) + "\n");
    
    reportln("Writing Output...","note");
    qcALL.seq.foreach( _.writeOutput(outfile, summaryWriter) );
    
    qcWIG.writeOutput(outfile + "wiggle.", summaryWriter);
    //qcGGC.writeOutput(outfile, summaryWriter);
    //qcCS.writeOutput(outfile, summaryWriter);
    //qcJD.writeOutput(outfile, summaryWriter);
    //qcQSC.writeOutput(outfile, summaryWriter);
    //qcIS.writeOutput(outfile, summaryWriter);
    //qcGC.writeOutput(outfile, summaryWriter);
    //qcJC.writeOutput(outfile, summaryWriter);
    //qcST.writeOutput(outfile, summaryWriter);
    //qcCC.writeOutput(outfile, summaryWriter);
    
    if(runFunc.contains("makeJunctionBed")){
      reportln("Making '.bed' junction count tracks... ","note");
      
      val bedfile = if(internalUtils.optionHolder.OPTION_noGzipOutput){
        Some(List(outfile + ".spliceJunctionAndExonCounts.forJunctionSeq.txt"))
      } else {
        Some(List(outfile + ".spliceJunctionAndExonCounts.forJunctionSeq.txt.gz"))
      }
      
      fileConversionUtils.makeSpliceJunctionBed.run(
           sizeFactorFile = None, 
           sizeFactors = None, 
           filenames = bedfile,
           sampleList = None,
           title = Some(trackTitlePrefix),
           ignoreSizeFactors = true,
           outfile = outfile + ".junctionBed.known.bed.gz",
           infilePrefix = "",
           infileSuffix = "",
           gffIterator = Some(anno_holder.makeFlatReader()),
           gff = "NONEXISTANTGFF.gff",
           stranded = stranded,
           digits = 2,
           includeFullSpliceNames = false,
           calcMean = false,
           nonflatgtf = false,
           rgb = None,
           trackTitle = trackTitlePrefix,
           additionalTrackOptions = "",
           skipAnnotatedJunctions = false, skipNovelJunctions = false
      );
      reportln("Done making browser tracks.","note");
    }
    
    if(internalUtils.Reporter.hasWarningOccurred()){
      summaryWriter.write("COMPLETED_WITHOUT_WARNING	0\n");
      reportln("QoRTs completed WITH WARNINGS! See log for details.","warn");
      val completedWarnWriter = openWriter(COMPLETED_WARN_FILEPATH);
      completedWarnWriter.write("# Note: if this file EXISTS, then QoRTs QC completed WITH WARNINGS. Warning messages follow:\n");
      completedWarnWriter.write(internalUtils.Reporter.getWarnings+"\n");
      completedWarnWriter.close();
    } else {
      summaryWriter.write("COMPLETED_WITHOUT_WARNING	1\n");
      reportln("QoRTs QC complete with no problems.","note");
    }
    
    summaryWriter.write("COMPLETED_WITHOUT_ERROR	1\n");
    
    close(summaryWriter);
    
    reportln("Done.","note");
    
    val finalTimeStamp = TimeStampUtil();
    reportln("Time spent on setup:           " + TimeStampUtil.timeDifferenceFormatter(samIterationTimeStamp.compareTo(initialTimeStamp)),"note");
    reportln("Time spent on SAM iteration:   " + TimeStampUtil.timeDifferenceFormatter(outputIterationTimeStamp.compareTo(samIterationTimeStamp)),"note");
    reportln("                               (" + minutesPerMillion + " minutes per million read-pairs)","note");
    reportln("                               (" + minutesPerMillionPF + " minutes per million read-pairs used)","note");
    reportln("Time spent on file output:     " + TimeStampUtil.timeDifferenceFormatter(finalTimeStamp.compareTo(outputIterationTimeStamp)),"note");
    reportln("Total runtime:                 " + TimeStampUtil.timeDifferenceFormatter(finalTimeStamp.compareTo(initialTimeStamp)),"note");
    
    val completedOkWriter = openWriter(COMPLETED_OK_FILEPATH);
    completedOkWriter.write("# Note: if this file EXISTS, then QoRTs completed without ERRORS.\n#If there were any warnings, then a file \""+COMPLETED_WARN_FILEPATH+"\" will also exist.\n#See QC.log for details.");
    completedOkWriter.close();
  }
  /*
  def runOnSeqFile_PAR(initialTimeStamp : TimeStampUtil, 
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
    
    val COMPLETED_OK_FILEPATH = outfile + COMPLETED_OK_FILENAME;
    
    val runMasterLevelFunctions = runFunc.filter(MASTERLEVEL_FUNCTION_LIST.contains(_)).toVector;

    val samFileAttributes = peekSamRecordIterator(infile);
    
    if(infile == "-"){
      error("FATAL ERROR: Cannot perform multithreaded file reading when reading from standard input! Set infile to a file name, rather than '-'!");
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
    

    val qcAllVector : Vector[QCUtility[Any]] = runMasterLevelFunctions.map((funcName : String) => {
      if(funcName == "GeneCalcs") new qcGetGeneCounts(stranded,fr_secondStrand,anno_holder,coda,coda_options,40, runFunc.contains("FPKM"), runFunc.contains("writeGenewiseGeneBody"), runFunc.contains("writeDESeq"), runFunc.contains("writeGeneCounts"), );
      else if(funcName == "InsertSize") new qcInnerDistance(anno_holder, stranded, fr_secondStrand, readLength);
      else if(funcName == "NVC") new qcNVC(readLength, runFunc.contains("writeClippedNVC"));
      else if(funcName == "CigarOpDistribution") new qcCigarDistribution(readLength) ;
      else if(funcName == "QualityScoreDistribution") new qcQualityScoreCounter(readLength, qcQualityScoreCounter.MAX_QUALITY_SCORE);
      else if(funcName == "GCDistribution") new qcGCContentCount(readLength)  ;
      else if(funcName == "JunctionCalcs") new qcJunctionCounts(anno_holder, stranded, fr_secondStrand, runFunc.contains("writeDEXSeq"), runFunc.contains("writeSpliceExon"), runFunc.contains("writeKnownSplices"), runFunc.contains("writeNovelSplices")) ;
      else if(funcName == "StrandCheck") new qcStrandTest(anno_holder, stranded, fr_secondStrand) ;
      else if(funcName == "chromCounts") new qcChromCount( fr_secondStrand);
      else if(funcName == "cigarMatch") new qcCigarMatch(readLength);
      else QCUtility.getBlankUnitUtil;
    })
    
    val qcALL = parConvert(qcAllVector, numThreads);
    
    reportln("QC Utilities Generated! ("+numThreads+" Threads)","note");
    standardStatusReport(initialTimeStamp);
    val samIterationTimeStamp = TimeStampUtil();
    
    qcALL.foreach((qcu : QCUtility[Any]) => {
       val coda_buffer : Array[Int] = if(qcu.getUtilityName == "GeneCalcs") {
         coda;
       } else {
         internalUtils.commonSeqUtils.getNewCauseOfDropArray;
       }
       val recordIter : Iterator[SAMRecord] = (new SAMFileReader(new File(infile))).iterator;
       val pairedIter : Iterator[(SAMRecord,SAMRecord)] = if(noMultiMapped){
           if(testRun) samRecordPairIterator(recordIter, false, 200000) else samRecordPairIterator(recordIter)
         } else {
           if(testRun) samRecordPairIterator_withMulti(recordIter, false, 200000) else samRecordPairIterator_withMulti(recordIter)
         }
       var readNum = 0;
       for(pair <- pairedIter){
         val (r1,r2) = pair;
         readNum += 1;
         if(internalUtils.commonSeqUtils.useReadPair(r1,r2,coda_buffer, coda_options, dropChrom, readGroup)){
           qcu.runOnReadPair(r1, r2, readNum);
         }
       }
    })
    
    standardStatusReport(initialTimeStamp);
    
    val outputIterationTimeStamp = TimeStampUtil();
    reportln("Read Stats:\n" + internalUtils.commonSeqUtils.causeOfDropArrayToString(coda, coda_options),"note");
    
    val summaryWriter = openWriter(outfile + "summary.txt");
    val strandedCode = if(! stranded){ 0 } else {if(fr_secondStrand) 2; else 1;}

    summaryWriter.write("FIELD	COUNT\n");
    summaryWriter.write("Stranded_Rule_Code	"+strandedCode+"\n");
    
    val readNum = coda(internalUtils.commonSeqUtils.CODA_TOTAL_READ_PAIRS);
    val iterationMinutes = (outputIterationTimeStamp.compareTo(samIterationTimeStamp) / 1000).toDouble / 60.toDouble;
    val minutesPerMillion = iterationMinutes / (readNum.toDouble / 1000000.toDouble);
    val minutesPerMillionPF = iterationMinutes / ((coda(internalUtils.commonSeqUtils.CODA_READ_PAIR_OK)).toDouble / 1000000.toDouble);
    
    summaryWriter.write("BENCHMARK_MinutesOnSamIteration	" + "%1.2f".format(iterationMinutes) + "\n");
    summaryWriter.write("BENCHMARK_MinutesPerMillionReads	" + "%1.2f".format(minutesPerMillion) + "\n");
    summaryWriter.write("BENCHMARK_MinutesPerMillionGoodReads	" + "%1.2f".format(minutesPerMillionPF) + "\n");
    
    reportln("Writing Output...","note");
    qcALL.seq.foreach( _.writeOutput(outfile, summaryWriter) );
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
    
    val completedOkWriter = openWriter(COMPLETED_OK_FILEPATH);
    completedOkWriter.write("# Note: if this file EXISTS, then QoRTs completed without errors.");
    completedOkWriter.close();
  }*/
  
  
  /*
   * SETTING GLOBAL OPTIONS:
   */
  def setQcOptions(noGzipOutput : Boolean){
    //registerGlobalParam[Boolean]("noGzipOutput", noGzipOutput)
    internalUtils.optionHolder.OPTION_noGzipOutput = noGzipOutput;
  }
  
}















