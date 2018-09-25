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
   
  final val QC_DEFAULT_ON_FUNCTION_MAP : scala.collection.immutable.Map[String,String] = Map(
                      ("InsertSize"                -> "Insert size distribution (paired-end data only)."),
                      ("NVC"                       -> "Nucleotide-vs-Cycle counts."),
                      ("CigarOpDistribution"       -> "Cigar operation rates by cycle and cigar operation length rates (deletions, insertions, splicing, clipping, etc)."),
                      ("QualityScoreDistribution"  -> "Calculate quality scores by cycle."),
                      ("GCDistribution"            -> "Calculate GC content distribution."),
                      ("GeneCalcs"                 -> "Find gene assignment and gene body calculations."),
                      ("JunctionCalcs"             -> "Find splice junctions (both novel and annotated)."),
                      ("StrandCheck"               -> "Check the strandedness of the data. Note that if the stranded option is set incorrectly, this tool will automatically print a warning to that effect."),
                      ("writeKnownSplices"         -> "Write known splice junction counts."),
                      ("writeNovelSplices"         -> "Write novel splice junction counts."),
                      ("writeJunctionSeqCounts"    -> "Write counts file designed for use with JunctionSeq (contains known splice junctions, gene counts, and exon counts)."),
                      ("writeDESeq"                -> "Write gene-level read/read-pair counts file, suitable for use with DESeq, EdgeR, etc."),
                      ("writeDEXSeq"               -> "Write exon-level read/read-pair counts file, designed for use with DEXSeq."),
                      ("writeGeneBody"             -> "Write gene-body distribution file."),
                      ("writeGenewiseGeneBody"     -> "Write file containing gene-body distributions for each (non-overlapping) gene."),
                      ("writeGeneCounts"           -> "Write extended gene-level read/read-pair counts file (includes counts for CDS/UTR, ambiguous regions, etc)."),
                      ("writeClippedNVC"           -> "Write NVC file containing clipped sequences."),
                      ("writeBiotypeCounts"        -> "Write a table listing read counts for each gene BioType (uses the optional \"gene_biotype\" GTF attribute)."),
                      ("chromCounts"               -> "Calculate chromosome counts"),
                      ("readLengthDistro"          -> "Tabulates the distribution of read lengths. "),
                      ("writeSpliceExon"           -> "Synonym for function \"writeJunctionSeqCounts\" (for backwards-compatibility)"),
                      ("cigarLocusCounts"          -> "BETA: This function is still undergoing basic testing. It is not intended for production use at this time."),
                      ("overlapMatch"              -> "BETA: This function calculates the matching of overlapping sections of paired reads.")
  );
  
  final val QC_DEFAULT_OFF_FUNCTION_MAP : scala.collection.immutable.Map[String,String] = Map(
                      ("calcOnTarget"              -> "BETA: requires --targetRegionBed parameter. This function calculates the rates at which reads intersect with the On-Target area. Intended for whole exome sequencing data. Make sure to use the --targetRegionBed parameter or else this function will deactivate! (Default: ON iff targetRegionBed param is found)"),
                      ("referenceMatch"            -> "BETA: requires --genomeFA parameter. This function calculates the matching against the reference. Requires the specification of genome fasta file(s). REQUIRES COORDINATE-SORTED BAM FILES! REQUIRES THAT FA AND BAM HAVE THE SAME CHROMOSOME ORDERING! (Default: ON iff genomeFA param is found)"),
                      ("fastqUtils"                -> "BETA: requires --rawfastq parameter.  Adds additional tests that use the supplied raw fastq file. Requires that one (or two) fastq files be supplied. (Default: ON iff rawfastq param is found)"),

                      ("writeDocs"                 -> "Writes a QC.documentation.txt file that documents all output files."),
                      

                      ("FPKM"                      -> "Write FPKM values. Note: FPKMs are generally NOT the recommended normalization method. We recommend using a more advanced normalization as provided by DESeq, edgeR, CuffLinks, or similar (default: OFF)"),
                      ("annotatedSpliceExonCounts" -> "Write counts for exons, known-splice-junctions, and genes, with annotation columns indicating chromosome, etc (default: OFF)"),
                      ("makeWiggles"               -> "Write \"wiggle\" coverage files with 100-bp window size. Note: this REQUIRES that the --chromSizes parameter be included! (default: OFF)"),
                      ("makeJunctionBed"           -> "Write splice-junction count \"bed\" files. (default: OFF)"),
                      ("makeAllBrowserTracks"      -> "Write both the \"wiggle\" and the splice-junction bed files (default: OFF)"),
                      ("calcDetailedGeneCounts"    -> "Calculate more detailed read counts for each gene, counting the number of reads that cover introns, cross-strand, etc (default: OFF)"),
                      ("cigarMatch"                -> "Work-In-Progress: this function is a placeholder for future functionality, and is not intended for use at this time. (default: OFF)"),
                      ("writeGeneBodyIv"           -> "Writes an optional additional file detailing the intervals used in the gene-body coverage calculations (\"QC.geneBodyCoverage.DEBUG.intervals.txt.gz\") (default: OFF)"),
                      ("testDataDump"              -> "EXPERIMENTAL: This function dumps a bunch of information for internal testing purposes. NOT FOR GENERAL USE! (default: OFF)"),
                      ("mismatchEngine"            -> "Internal module that runs overlap/reference mismatch calculations. Automatically included on any runs that include these functions.")
  );
  
  final val QC_HIDDEN_BETA_FUNCTION_MAP : scala.collection.immutable.Map[String,String] = Map(
                      ("outputFastqFile" -> "UTILITY METHOD! Simply outputs fastq file(s)"),
                      ("outputBAMFile"   -> "UTILITY METHOD! Simply outputs BAM file with all PF reads.")
                      //("featureComboCount" -> "Write a table including all combinations of exon/junction features, and the counts. (BETA FUNCTION! NOT FOR PRODUCTION USE!)")
      );
  
  final val QC_DEFAULT_ON_FUNCTION_LIST  : scala.collection.immutable.Set[String] = QC_DEFAULT_ON_FUNCTION_MAP.keySet;
  final val QC_DEFAULT_OFF_FUNCTION_LIST  : scala.collection.immutable.Set[String] = QC_DEFAULT_OFF_FUNCTION_MAP.keySet;
  
  /*final val QC_DEFAULT_ON_FUNCTION_LIST  : scala.collection.immutable.Set[String] = scala.collection.immutable.Set("InsertSize",
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
                                                             "writeGeneBody",
                                                             "writeGenewiseGeneBody", 
                                                             "writeGeneCounts", 
                                                             "writeClippedNVC", 
                                                             "chromCounts");
  final val QC_DEFAULT_OFF_FUNCTION_LIST : scala.collection.immutable.Set[String] = scala.collection.immutable.Set(
                                                             "FPKM",
                                                             "annotatedSpliceExonCounts",
                                                             "makeWiggles",
                                                             "makeJunctionBed",
                                                             "makeAllBrowserTracks", 
                                                             "cigarMatch",
                                                             "cigarLocusCounts");*/
  final val QC_FUNCTION_LIST : scala.collection.immutable.Set[String] = QC_DEFAULT_ON_FUNCTION_LIST ++ QC_DEFAULT_OFF_FUNCTION_LIST;
  final val COMPLETED_OK_FILENAME = ".QORTS_COMPLETED_OK";
  final val IN_PROGRESS_FILENAME = ".QORTS_RUNNING";
  final val COMPLETED_WARN_FILENAME = ".QORTS_COMPLETED_WARN";
  //final val MASTERLEVEL_FUNCTION_LIST = List[String]("GeneCalcs", "InsertSize","NVC","CigarOpDistribution","QualityScoreDistribution","GCDistribution","JunctionCalcs","StrandCheck","chromCounts","cigarMatch","makeWiggles");
  
  
  final val QC_INCOMPATIBLE_WITH_SINGLE_END_FUNCTION_LIST : scala.collection.immutable.Set[String] = scala.collection.immutable.Set(
      "InsertSize","cigarMatch","overlapMatch"
  );
  
  //"InsertSize","NVC","CigarOpDistribution","QualityScoreDistribution","GCDistribution","GeneCalcs",
  //"JunctionCounts", "StrandCheck", "writeKnownSplices","writeNovelSplices","writeSpliceExon", 
  //"writeDESeq","writeDEXSeq","writeGenewiseGeneBody", "writeGeneCounts"
  final val QC_FUNCTION_DEPENDANCIES : Map[String,Set[String]] = Map[String,Set[String]](
      ("writeDESeq" -> Set("GeneCalcs")),
      ("writeGeneCounts" -> Set("GeneCalcs")),
      ("writeGenewiseGeneBody" -> Set("writeGeneBody")),
      ("writeGeneBody" -> Set("GeneCalcs")),
      ("writeDEXSeq" -> Set("JunctionCalcs")),
      ("writeJunctionSeqCounts" -> Set("writeSpliceExon")), //backwards compatible synonym "writeSpliceExon"
      ("writeSpliceExon" -> Set("JunctionCalcs")),
      ("writeKnownSplices" -> Set("JunctionCalcs")),
      ("writeNovelSplices" -> Set("JunctionCalcs")),
      ("annotatedSpliceExonCounts" -> Set("JunctionCalcs")),
      ("writeClippedNVC" -> Set("NVC")),
      ("makeAllBrowserTracks" -> Set("makeJunctionBed","makeWiggles")),
      ("writeBiotypeCounts" -> Set("GeneCalcs")),
      ("referenceMatch" -> Set("mismatchEngine")),
      ("overlapMatch" -> Set("mismatchEngine")),
      ("writeGeneBodyIv" -> Set("writeGeneBody"))
  );
  
  final val SAM_PEEK_LINECT = 10000;
  final val SAM_TESTRUN_LINECT = 100000;
  
  //def run(args : Array[String]){
  
  class allQC_runner extends CommandLineRunUtil {
    override def priority = 1;
    val parser : CommandLineArgParser = 
      new CommandLineArgParser(
          command = "QC", 
          quickSynopsis = "Runs a battery of QC tools", 
          synopsis = "", 
          description = "This utility runs a large battery of QC / data processing tools on a single given sam or bam file. "+
                        "This is the primary function of the QoRTs utility. All analyses are run via a single pass through the sam/bam file."+
                        ""+
                        ""+
                        ""+
                        "",
          argList = 


                    new UnaryArgument(   name = "singleEnded", 
                                         arg = List("--singleEnded","-e"), // name of value
                                         argDesc = "Flag to indicate that reads are single end." // description
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
                    new BinaryArgument[Int](name = "minMAPQ",
                                           arg = List("--minMAPQ"),  
                                           valueName = "num", 
                                           argDesc = "Filter out reads with less than the given MAPQ. Most RNA-Seq aligners use the MAPQ field "+
                                                     "to differentiate uniquely-mapped and multi-mapped reads. However, different aligners "+
                                                     "use a different MAPQ conventions. By default, all reads with a MAPQ "+
                                                     "of less than 255 will be excluded, as this is the MAPQ associated with "+
                                                     "uniquely-aligned reads generated by the RNA-STAR aligner. For use with "+
                                                     "TopHat2 you should set this to 50. The MAPQ behavior for GSNAP is not well "+
                                                     "documented, but it appears that a filtering threshold of 30 should be adequate. "+
                                                     "Set this to 0 to turn off mapq filtering.", 
                                           defaultValue = Some(255)
                                           ) :: 
                    new UnaryArgument( name = "generatePlots",
                                         arg = List("--generatePlots"), // name of value
                                         argDesc = "Generate all single-replicate QC plots. Equivalent to the combination of: --generateMultiPlot --generateSeparatePlots and --generatePdfReport. "+
                                                   "This option will cause QoRTs to make an Rscript system call, loading the R package QoRTs. "+
                                                   ""+
                                                   "(Note: this requires that R be installed and in the PATH, and that QoRTs be installed on that R installation)"
                                       ) ::

                    new UnaryArgument(    name = "testRun",
                                         arg = List("--testRun","-t"), // name of value
                                         argDesc = "Flag to indicate that only the first 100k reads should be read in. Used for testing." // description
                                       ) ::
                                       

                    new UnaryArgument( name = "keepMultiMapped",
                                         arg = List("--keepMultiMapped"), // name of value
                                         argDesc = "Flag to indicate that the tool should NOT filter out multi-mapped reads. Note that even with this flag raised this utility will still only "+
                                                    "use the 'primary' alignment location for each read. By default any reads that are marked as multi-mapped will be ignored entirely."+
                                                    " Most aligners use the MAPQ value to mark multi-mapped reads. Any read with MAPQ < 255 is assumed to be non-uniquely mapped (this is the standard used by RNA-STAR and TopHat/TopHat2). "+
                                                    " This option is equivalent to \"--minMAPQ 0\"."// description
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

                    new BinaryArgument[List[String]](   name = "dropChromList",
                                                        arg = List("--dropChrom"),  
                                                        valueName = "dropChromosomes", 
                                                        argDesc = "A comma-delimited list of chromosomes to ignore and exclude from all analyses. Important: no whitespace!", 
                                                        defaultValue = Some(List[String]())
                                                        ) ::
                    new BinaryArgument[List[String]](   name = "skipFunctions",
                                                        arg = List("--skipFunctions"),  
                                                        valueName = "func1,func2,...", 
                                                        argDesc = "A list of functions to skip (comma-delimited, no whitespace). See the sub-functions list, below. The default-on functions are: "+QC_DEFAULT_ON_FUNCTION_LIST.mkString(", "), 
                                                        defaultValue = Some(List[String]())
                                                        ) ::
                    new BinaryArgument[List[String]](name = "addFunctions",
                                                        arg = List("--addFunctions"),  
                                                        valueName = "func1,func2,...", 
                                                        argDesc = "A list of functions to add (comma-delimited, no whitespace). This can be used to add functions that are off by default. Followed by a comma delimited list, with no internal whitespace. See the sub-functions list, below. The default-off functions are: "+QC_DEFAULT_OFF_FUNCTION_LIST.mkString(", "), 
                                                        defaultValue = Some(List[String]())
                                                        ) :: 
                    new BinaryArgument[List[String]](name = "runFunctions",
                                                        arg = List("--runFunctions"),  
                                                        valueName = "func1,func2,...", 
                                                        argDesc = "The complete list of functions to run  (comma-delimited, no whitespace). Setting this option runs ONLY for the functions explicitly requested here (along with any functions upon which the assigned functions are dependent). See the sub-functions list, below. Allowed options are: "+QC_FUNCTION_LIST.mkString(", "), 
                                                        defaultValue = Some(List[String]())
                                                        ) :: 
                    new BinaryOptionArgument[Int](
                                         name = "seqReadCt", 
                                         arg = List("--seqReadCt"), 
                                         valueName = "val",  
                                         argDesc = "(Optional) The total number of reads (or read-pairs, for paired-end data) generated by the sequencer for this sample, prior to alignment. "+
                                                   "This will be passed on into the QC.summary.txt file and used to calculate mapping rate."+
                                                   ""+
                                                   ""
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "rawfastq", 
                                         arg = List("--rawfastq"), 
                                         valueName = "myfastq.1.fq.gz,myfastq.2.fq.gz",  
                                         argDesc = "(Optional) The raw fastq, prior to alignment. In normal operation, this is used ONLY to calculate the number of pre-alignment reads (or read-pairs) simply by counting the number of lines and dividing by 4. "+
                                                   "Alternatively, the number of pre-alignment read-pairs can be included explicitly via the --seqReadCt option, or added in the "+
                                                   "plotting / cross-comparison step by including the input.read.pair.count column in the replicate decoder."+
                                                   "In general, the --seqReadCt option is recommended when available.\n"+
                                                   "Certain optional QC functions are also available that utilize the raw Fastq file in other ways. "+
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
                    new BinaryArgument[String](name = "title",
                                           arg = List("--title"),  
                                           valueName = "myTitle", 
                                           argDesc = "The title of the replicate. Used for the track name in the track definition line of "+
                                                      "any browser tracks ('.wig' or '.bed' files) generated by this utility. "+
                                                      "Also may be used in the figure text, if figures are being generated."+
                                                     "Note that no browser tracks will be created by default, unless the '--chromSizes' option is set. Bed files can also be generated using the option '--addFunction makeJunctionBed'", 
                                           defaultValue = Some("Untitled")
                                           ) :: 
                    new BinaryOptionArgument[String](
                                         name = "flatgfffile", 
                                         arg = List("--flatgff"), 
                                         valueName = "flattenedGffFile.gff.gz",  
                                         argDesc = "A \"flattened\" gff file that matches the standard gtf file. Optional. "+
                                                   "The \"flattened\" gff file assigns unique identifiers for all exons, splice junctions, and aggregate-genes. "+
                                                   "This is used for the junction counts and exon counts (for DEXSeq). "+
                                                   "The flattened gtf file can be generated using "+
                                                   "the \"makeFlatGff\" command. Flattened GFF files containing novel splice junctions can be generated using the \"mergeNovelSplices\" function. "+
                                                   "Note that (for most purposes) the command should be run with the same strandedness code as found in the dataset. "+
                                                   "Running a flattened gff that was generated using a different strandedness mode may be useful for certain purposes, but is generally not supported "+
                                                   "and is for advanced users only."+
                                                   "See the documentation for makeFlatGff for more information. "+
                                                   "\n"+
                                                   "If the filename ends with \".gz\" or \".zip\", the file will be parsed using the appropriate decompression method."
                                        ) ::


                    new UnaryArgument( name = "generateMultiPlot",
                                         arg = List("--generateMultiPlot"), // name of value
                                         argDesc = "Generate a multi-frame figure, containing a visual summary of all QC stats. "+
                                                   "(Note: this requires that R be installed and in the PATH, and that QoRTs be installed on that R installation)"+
                                                   ""+
                                                   ""
                                       ) ::
                    new UnaryArgument( name = "generateSeparatePlots",
                                         arg = List("--generateSeparatePlots"), // name of value
                                         argDesc = "Generate seperate plots for each QC stat, rather than only one big multiplot. "+
                                                   "(Note: this requires that R be installed and in the PATH, and that QoRTs be installed on that R installation) "+
                                                   ""+
                                                   ""+
                                                   "" 
                                       ) ::
                    new UnaryArgument( name = "generatePdfReport",
                                       arg = List("--generatePdfReport"), // name of value
                                       argDesc = "Generate a pdf report. "+
                                                 "(Note: this requires that R be installed and in the PATH, and that QoRTs be installed on that R installation)"+
                                                 ""+
                                                 ""+
                                                 ""
                                       ) ::
                    new UnaryArgument( name = "prefilterImproperPairs",
                                       arg = List("--prefilterImproperPairs"), // name of value
                                       argDesc = "If this parameter is used, QoRTs will ignore all reads that lack the 0x2 flag, which "+
                                                 "indicates that the read is properly paired. It will not attempt to match such reads into pairs. "+
                                                 "This option may be necessary with certain aligners (including BWA). Note that if this filter is used, "+
                                                 "QoRTs will not tally improperly paired reads."+
                                                 ""
                                       ) ::
                    new BinaryArgument[Int](   name = "adjustPhredScore",
                                                        arg = List("--adjustPhredScore"),  
                                                        valueName = "val", 
                                                        argDesc = "QoRTs expects input files to conform to the SAM format specification, "+
                                                                  "which requires all Phred scores to be in Phred+33 encoding. However "+
                                                                  "some older tools produce SAM files with nonstandard encodings. To read such data, "+
                                                                  "you can set this parameter to subtract from the apparent (phred+33) phred score. Thus, to read "+
                                                                  "Phred+64 data (produced by Illumina v1.3-1.7), set this parameter to 31. Note: QoRTs does not support "+
                                                                  "negative Phred scores. NOTE: THIS OPTION IS EXPERIMENTAL!", 
                                                        defaultValue = Some(0)
                                                        ) ::
                    new BinaryArgument[Int](   name = "maxPhredScore",
                                                        arg = List("--maxPhredScore"),  
                                                        valueName = "val", 
                                                        argDesc = "According to the standard FASTQ and SAM format specification, Phred "+
                                                                  "quality scores are supposed to range from 0 to 41. However, certain sequencing machines "+
                                                                  "such as the HiSeq4000 supposedly produce occasional quality scores as high as 45. If your dataset "+
                                                                  "contains quality scores in excess of 41, then you must use this option to set the maximum legal quality "+
                                                                  "score. Otherwise, QoRTs will throw an error.", 
                                                        defaultValue = Some(41)
                                                        ) ::
                   //BETA OPTIONS:
                    /* Impossible: requires foreknowledge of top genes...
                     * new BinaryOptionArgument[Double](
                                         name = "ignoreTopQuantileGenes", 
                                         arg = List("--ignoreTopQuantileGenes"), 
                                         valueName = "val",  
                                         argDesc =  "If this option is set, almost all analyses will ignore reads that appear to map to genes in the upper "+
                                                    "quantile. The supplied value specifies the proportion of genes that should be ignored. This option overrides the \"--restrictToGeneList\" and \"--dropGeneList\" options.\n"+
                                                    "WARNING: This feature is still EXPERIMENTAL, and is not yet fully tested or ready for production use."+
                                                    ""+
                                                    ""+
                                                    ""+
                                                    ""+
                                                    ""+
                                                    ""
                                        ) ::*/
                    new BinaryArgument[String](   name = "summaryFileName",
                                               arg = List("--summaryFileSuffix"),  
                                               valueName = ".summary.txt", 
                                               argDesc = "The suffix of the 'summary' file. This is useful to set if you want to run multiple QC runs in parallel "+
                                                         "to reduce runtime, without overwriting one another's summary files."+
                                                         "In particular, the NVC metrics often take a long time to run, so splitting those "+
                                                         "off using the --runFunctions parameter might speed things up considerably. "+
                                                         "Note that 'QC' will be appended in the actual filename. "+
                                                         "THIS OPTION IS BETA!", 
                                               defaultValue = Some(".summary.txt")
                                           ) ::
                                        
                                        
                   new BinaryOptionArgument[String](
                                         name = "extractReadsByMetric", 
                                         arg = List("--extractReadsByMetric"), 
                                         valueName = "metric=value",  
                                         argDesc =  "THIS OPTIONAL PARAMETER IS STILL UNDER BETA TESTING. This parameter allows "+
                                                    "the user to extract anomalous reads that showed up in previous QoRTs runs. "+
                                                    "Currently reads can be extracted based on the following metrics: StrandTestStatus, InsertSize and GCcount."+
                                                    ""+
                                                    ""+
                                                    ""
                                        ) ::
                    new UnaryArgument( name = "keepOnTarget",
                                         arg = List("--keepOnlyOnTarget"), // name of value
                                         argDesc = "Experimental flag. Ignores reads that DO NOT fall within the target region (specified by the required bedfile using the --targetRegionBed parameter)."+
                                                   "" // description
                                       ) ::
                    new UnaryArgument( name = "dropOnTarget",
                                         arg = List("--dropOnTarget"), // name of value
                                         argDesc = "Experimental flag. Ignores reads that DO fall within the target region (specified by the required bedfile using the --targetRegionBed parameter)."+
                                                   "" // description
                                       ) ::
                    new BinaryOptionArgument[Double](
                                         name = "randomSubsample", 
                                         arg = List("--randomSubsample"), 
                                         valueName = "1.00",  
                                         argDesc =  "If this option is set, QoRTs will ignore a random fraction of the input read pairs. This can drastically reduce runtime, though it may reduce the accuracy of the output QC metrics."
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
                    new UnaryArgument( name = "WholeExomeSeq",
                                         arg = List("--DNA"), // name of value
                                         argDesc = "BETA: This flag makes various changes to allow QoRTs to "+
                                                   "run on whole-exome or whole-genome DNA-Seq data."+
                                                   "" // description
                                       ) :: 
                    new UnaryArgument( name = "isRNASeq",
                                         arg = List("--RNA"), // name of value
                                         argDesc = "Indicates that the data is RNA-Seq (this is the default: flag does nothing)."+
                                                   "" // description
                                       ) :: 

                    new BinaryOptionArgument[List[String]](
                                         name = "genomeFA", 
                                         arg = List("--genomeFA"), 
                                         valueName = "chr.fa.gz[,chr2.fa,...]",  
                                         argDesc = "Reference genome sequence. This can either be a single FASTA file with all the chromosomes included, or a comma-delimited list of fasta files with 1 chromosome each. "+
                                                   "Note: IF multiple fasta files are specificed, each must contain ONLY ONE chromosome. "+
                                                   "If a single multi-chromosome fasta file is specificed, performance will be improved if the chromosomes are in the same order as they are found in the BAM file, "+
                                                   "however, this is not required. "+
                                                   "The genomic sequence is used by certain experimental sub-utilities (currently only the referenceMatch utility). "+
                                                   "Comma delimited, no spaces. "+
                                                   "Fasta files can be in plaintext, gzipped or zipped."+
                                                   ""
                                        ) ::                              
                                        
                    new BinaryArgument[Int](   name = "genomeBufferSize",
                                                        arg = List("--genomeBufferSize"),  
                                                        valueName = "val", 
                                                        argDesc = "The size of the genome fasta buffer. Tuning this parameter may improve performance."+
                                                                  ""+
                                                                  ""+
                                                                  ""+
                                                                  "", 
                                                        defaultValue = Some(10000)
                                                        ) ::

                    new BinaryArgument[String](
                                         name = "outfilePrefix", 
                                         arg = List("--outfilePrefix"), 
                                         valueName = "sampID",  
                                         argDesc = "Prefix to be prepended to all output files. If this is set, all output files will use the format: \"outfiledir/<prefix>QC.qcfilename.txt.gz\""+
                                                   "",
                                          defaultValue = Some("")
                                        ) ::      
                                        
//DEPRECIATED or BETA OPTIONS:
                                       
                    new UnaryArgument( name = "nameSorted",
                                         arg = List("--nameSorted"), // name of value
                                         argDesc = "DEPRECATED: Relevant for paired-end reads only. \n"+
                                                   "This flag is used to run QoRTs in \"name-sorted\" mode. This flag is optional, as under the "+"default mode QoRTs will accept BAM files sorted by either name OR position. "+
                                                   "However, "+
                                                   "The only actual requirement in this mode is that "+
                                                   "read pairs be adjacent. \n"+
                                                   "Errors may occur if the SAM flags are inconsistent: for example, if orphaned reads appear with the \"mate mapped\" SAM flag set."+
                                                   ""+
                                                   ""+
                                                   "" // description
                                       ) ::    
                                       
                    new UnaryArgument( name = "coordSorted",
                                         arg = List("--coordSorted"), // name of value
                                         argDesc = ""+
                                                   "DEPRECATED: this mode is now subsumed by the default mode and as such this parameter is now nonfunctional.\n"+
                                                   "Note that, in the default mode, for paired-end data QoRTs will accept "+
                                                   "EITHER coordinate-sorted OR name-sorted bam files. In \"--nameSorted\" mode, QoRTs ONLY accepts "+
                                                   "name-sorted bam files.\n"+
                                                   "If a large fraction of the read-pairs are "+
                                                   "mapped to extremely distant loci (or to different chromosomes), then memory "+
                                                   "issues may arise. However, this should not be a problem with most datasets. "+
                                                   "Technically by default QoRTs can run on arbitrarily-ordered bam files, "+
                                                   "but this is STRONGLY not recommended, as memory usage will by greatly increased." // description
                                       ) ::
                    new UnaryArgument( name = "noMultiMapped",
                                         arg = List("--fileContainsNoMultiMappedReads"), // name of value
                                         argDesc = "DEPRECATED. Flag to indicate that the input sam/bam file contains only primary alignments (ie, no multi-mapped reads). This flag is ALWAYS OPTIONAL, but when applicable this utility will run (slightly) faster when using this argument. (DEPRECIATED! The performance improvement was marginal)" // description
                                       ) ::
                    new UnaryArgument( name = "parallelFileRead",
                                         arg = List("--parallelFileRead"), // name of value
                                         argDesc = "DEPRECATED: DO NOT USE. Flag to indicate that bam file reading should be run in paralell for increased speed. Note that in this mode you CANNOT read from stdin. Also note that for this to do anything useful, the numThreads option must be set to some number greater than 1. Also note that additional threads above 9 will have no appreciable affect on speed." // description
                                       ) ::    
                    new BinaryArgument[Int](name = "numThreads",
                                                        arg = List("--numThreads"),  
                                                        valueName = "num", 
                                                        argDesc = "DEPRECIATED, nonfunctional.", 
                                                        defaultValue = Some(1)
                                                        ) :: 
                    new UnaryArgument( name = "checkForAlignmentBlocks",
                                         arg = List("--checkForAlignmentBlocks"), // name of value
                                         argDesc = "Certain aligners will mark reads 'aligned' even though they have no aligned "+
                                                   "bases. This option will automatically check for some reads and ignore them, "+
                                                   "rather than throwing an error."// description
                                       ) ::
                    new BinaryOptionArgument[String](
                                         name = "targetRegionBed", 
                                         arg = List("--targetRegionBed"), 
                                         valueName = "targetRegion.bed",  
                                         argDesc =  "For whole exome sequencing, this specifies the exome target regions."
                                        ) :: 
                    new BinaryOptionArgument[Int](
                                         name = "stopAfterNReads", 
                                         arg = List("--stopAfterNReads"), 
                                         valueName = "n",  
                                         argDesc =  "Stop after reading in n reads or read-pairs."
                                        ) :: 
                    new BinaryOptionArgument[Long](
                                         name = "randomSeed", 
                                         arg = List("--randomSeed"), 
                                         valueName = "n",  
                                         argDesc =  "Use specified random seed."
                                        ) :: 
                    new UnaryArgument(
                                         name = "parseIlluminaStyleReadIDs", 
                                         arg = List("--parseIlluminaStyleReadIDs"), 
                                         argDesc =  "Specifies that the read-names are in the illumina style. CURRENTLY NONFUNCTIONAL!"
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
                                         name = "outdir",
                                         valueName = "qcDataDir",
                                         argDesc = "The output directory." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS,
         manualExtras = QC_DEFAULT_ON_FUNCTION_MAP.foldLeft("DEFAULT SUB-FUNCTIONS\n")((soFar,curr) => {
           val depString = QC_FUNCTION_DEPENDANCIES.get(curr._1) match {
             case Some(depset) => " [Depends: "+depset.toList.mkString(", ")+"]";
             case None => "";
           }
           soFar + "    "+curr._1+"\n"+wrapLineWithIndent(curr._2 + depString,internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH,8)+"\n";
         }) + QC_DEFAULT_OFF_FUNCTION_MAP.foldLeft("NON-DEFAULT SUB-FUNCTIONS\n")((soFar,curr) => {
           val depString = QC_FUNCTION_DEPENDANCIES.get(curr._1) match {
             case Some(depset) => " [Depends: "+depset.toList.mkString(", ")+"]";
             case None => "";
           }
           soFar + "    "+curr._1+"\n"+wrapLineWithIndent(curr._2+depString,internalUtils.commandLineUI.CLUI_CONSOLE_LINE_WIDTH,8)+"\n";
         }),
         markdownManualExtras = QC_DEFAULT_ON_FUNCTION_MAP.foldLeft("## DEFAULT SUB-FUNCTIONS:\n")((soFar,curr) => {
           val depString = QC_FUNCTION_DEPENDANCIES.get(curr._1) match {
             case Some(depset) => " [Depends: "+depset.toList.mkString(", ")+"]";
             case None => "";
           }
           //"### "+(getFullSyntax()).replaceAll("_","\\\\_")+":\n\n> "+(describe()).replaceAll("_","\\\\_")+ (" ("+argType+")\n\n").replaceAll("_","\\\\_");
           soFar + "* " + curr._1 + ": " + curr._2 + depString + "\n\n";
         }) + QC_DEFAULT_OFF_FUNCTION_MAP.foldLeft("## NON-DEFAULT SUB-FUNCTIONS:\n")((soFar,curr) => {
           val depString = QC_FUNCTION_DEPENDANCIES.get(curr._1) match {
             case Some(depset) => " [Depends: "+depset.toList.mkString(", ")+"]";
             case None => "";
           }
           soFar + "* " + curr._1 + ": " + curr._2 + depString + "\n\n";
         })
      );
    
    def run(args : Array[String]){
      val out = parser.parseArguments(args.toList.tail);

      
      if(out){
      runAllQC.run(parser.get[String]("infile"),
          parser.get[String]("outdir"),
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
          ! parser.get[Boolean]("nameSorted"),
          parser.get[Option[Int]]("maxReadLength"),
          parser.get[Option[Int]]("seqReadCt"),
          parser.get[Option[List[String]]]("rawfastq"),
          parser.get[Option[String]]("chromSizes"),
          parser.get[String]("title"),
          parser.get[Boolean]("generatePlots") | parser.get[Boolean]("generateMultiPlot"),
          parser.get[Boolean]("generatePlots") | parser.get[Boolean]("generateSeparatePlots"),
          parser.get[Boolean]("generatePlots") | parser.get[Boolean]("generatePdfReport"),
          parser.get[Option[String]]("extractReadsByMetric"),
          parser.get[Int]("adjustPhredScore"),
          parser.get[Int]("maxPhredScore"),
          parser.get[Boolean]("WholeExomeSeq"),
          parser.get[Boolean]("checkForAlignmentBlocks"),
          parser.get[String]("summaryFileName"),
          parser.get[Option[String]]("targetRegionBed"),
          parser.get[Option[Int]]("stopAfterNReads"),
          parser.get[Option[List[String]]]("genomeFA"),
          parser.get[Boolean]("dropOnTarget"),
          parser.get[Boolean]("keepOnTarget"),
          parser.get[Option[Double]]("randomSubsample"),
          parser.get[Option[Long]]("randomSeed"),
          parser.get[String]("outfilePrefix"),
          parser.get[Int]("genomeBufferSize"),
          prefilterImproperPairs = parser.get[Boolean]("prefilterImproperPairs")
      );
      }
    }
  }
  
  def run(infile : String, 
          outdir : String, 
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
          rawfastq : Option[List[String]],
          chromSizes : Option[String],
          trackTitlePrefix : String,
          generateMultiPlot : Boolean,
          generateSeparatePlots : Boolean,
          generatePdfReport : Boolean,
          extractReadsByMetric : Option[String],
          adjustPhredScore : Int,
          maxPhredScore : Int,
          isWES : Boolean,
          checkForAlignmentBlocks : Boolean,
          summaryFileName : String,
          targetRegionBed : Option[String],
          stopAfterNReads : Option[Int],
          genomeFA : Option[List[String]],
          dropOnTarget : Boolean,
          keepOnTarget : Boolean,
          randomSubsample : Option[Double],
          randomSeed : Option[Long],
          outfilePrefix : String,
          genomeBufferSize : Int,
          prefilterImproperPairs : Boolean){
    
    randomSeed match {
      case Some(seed) => scala.util.Random.setSeed(seed);
      case None => {
        //do nothing
      }
    }
    
    val outDirFile = new File(outdir);
    if(! outDirFile.exists()){
      reportln("Creating Directory: "+ outdir,"note");
      outDirFile.mkdir();
    }
    val outfile = outdir + "/" + outfilePrefix + "QC";
    
    val logfile = outfile +"."+ internalUtils.stdUtils.getRandomString(12) + ".log";
    internalUtils.Reporter.init_completeLogFile(logfile);
    reportln("Created Log File: " + logfile,"note");
    
    internalUtils.Reporter.init_inProgressFile(outfile + IN_PROGRESS_FILENAME);
    
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
    
    if(isSingleEnd){
      reportln("QoRTs is Running in single-end mode.","note");
      reportln("Note: read-sorting is irrelevant in single-ended mode.","note");
    } else {
      reportln("QoRTs is Running in paired-end mode.","note");
      if(unsorted) reportln("QoRTs is Running in any-sorted mode.","note");
      else         reportln("QoRTs is Running in name-sorted mode.","note");
    }
    
    val dropChrom = dropChromList.toSet;
    
    val defaultFunctonList : Set[String] = (if(chromSizes.isEmpty){
      QC_DEFAULT_ON_FUNCTION_LIST;
    } else {
      reportln("Chromosome size file added. Adding target wiggle plot generation.","debug");
      QC_DEFAULT_ON_FUNCTION_LIST ++ Set("makeWiggles");
    }) ++ (if(targetRegionBed.isEmpty){
      Set[String]();
    } else {
      reportln("Target BED file specified. Adding target bed file testing.","debug");
      Set[String]("calcOnTarget");
    }) ++ (if(rawfastq.isEmpty){
      Set[String]();
    } else {
      reportln("Raw fastq files specified. Adding fastq testing.","debug");
      Set[String]("fastqUtils");
    }) ++ (if(genomeFA.isEmpty){
      reportln("No parameter --genomeFA found.","debug");
      Set[String]();
    } else {
      reportln("Parameter --genomeFA found. Adding reference mismatch testing.","debug");
      Set[String]("referenceMatch");
    })
    
    //reportln("Default functions: " + wrapSimpleLineWithIndent_staggered(defaultFunctonList.toList.sorted.mkString(", "), 68, "        ", ""),"debug");
    
    val runFunc_initial = (if(runFunctions.isEmpty){
      (defaultFunctonList ++ addFunctions.toSet) -- dropFunctions.toSet;
    } else {
      runFunctions.toSet -- dropFunctions.toSet;
    }) -- ( if(! isSingleEnd) Set() else QC_INCOMPATIBLE_WITH_SINGLE_END_FUNCTION_LIST )
    
    val runFuncTEMP : scala.collection.immutable.Set[String] = (if(restrictToGeneList.isEmpty & dropGeneList.isEmpty) scala.collection.immutable.Set[String]() else {
      reportln("NOTE: Options \"--restrictToGeneList\" and \"--dropGeneList\" require function \"GeneCalcs\". Adding \"GeneCalcs\" to the active function list...","note")
      scala.collection.immutable.Set[String]("GeneCalcs");
    });
        
    def addDependencies(currFunc : String, soFar : scala.collection.immutable.Set[String]) : scala.collection.immutable.Set[String] = {
      QC_FUNCTION_DEPENDANCIES.get(currFunc) match {
        case Some(reqFuncSet) => {
          val addFunc = reqFuncSet -- soFar;
          addFunc.foreach((rf) => { reportln("NOTE: Function \"" + currFunc +"\" requires function \"" + rf + "\". Adding \""+rf+"\" to the active function list...","note") });
          val newFuncSet = soFar ++ addFunc;
          
          addFunc.foldLeft(newFuncSet)((sf, cf) => {
            addDependencies(cf,sf);
          });
        }
        case None => {
          soFar;
        }
      }
    }
    
    val runFuncTEMP2 = (runFunc_initial ++ runFuncTEMP).foldLeft(runFunc_initial ++ runFuncTEMP)((soFar,currFunc) => {
      addDependencies(currFunc,soFar);
    });
    val runFunc = if(! isSingleEnd) runFuncTEMP2 else {
      runFuncTEMP2 -- QC_INCOMPATIBLE_WITH_SINGLE_END_FUNCTION_LIST;
    }
    
    val notFoundFunc = runFunc.find(  ! (QC_FUNCTION_LIST ++ QC_HIDDEN_BETA_FUNCTION_MAP.keySet).contains(_) ) ;
    if( ! notFoundFunc.isEmpty ){
      error("ERROR: Function not found: \"" + notFoundFunc.get +"\"");
    }
    
    if(runFunc.contains("makeWiggles") && chromSizes.isEmpty){
      error("Error: function makeWiggles REQUIRES the chromSizes parameter! Set parameter '--chromSizes' or turn off function 'makeWiggles'.");
    }
    
    if(runFunc.contains("referenceMatch") && genomeFA.isEmpty){
       error("Error: function referenceMatch REQUIRES the genomeFA parameter! Set parameter '--genomeFA' or turn off function 'referenceMatch'.");
    }
    
    //wrapSimpleLineWithIndent_staggered(line : String, width : Int, indent : String, firstLineIndent : String)
    reportln("Running functions: " + wrapSimpleLineWithIndent_staggered(runFunc.toList.sorted(internalUtils.stdUtils.AlphabetOrdering).mkString(", "), 68, "        ", ""),"progress");
    
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
    val anno_holder = new qcGtfAnnotationBuilder(gtffile , flatgtffile , stranded , stdGtfCodes, flatGtfCodes, targetRegionBed);
    
    //Add extract genes by metric function:
    //StrandTestStatus, InsertSize and GCcount
    val extractReadsFunction : Option[((String, Int, (Int,Int)) => Boolean)] = extractReadsByMetric match {
      case Some(erbm) => {
        val f = erbm.split(",").map((metricString : String) => {
          val metricArray = metricString.split("=");
          if(metricArray.length != 2) {
            error("Badly formatted extractReadsByMetric parameter");
          }
          if(metricArray(0) == "StrandTestStatus"){
            (st : String, insertSize : Int, gc : (Int,Int)) => {
              st == metricArray(1);
            }; 
          } else if(metricArray(0) == "InsertSize"){
            val targetSize = string2int(metricArray(1));
            (st : String, insertSize : Int, gc : (Int,Int)) => {
              insertSize == targetSize;
            }; 
          } else if(metricArray(0) == "GCcount"){
            val target = string2int(metricArray(1));
            (st : String, insertSize : Int, gc : (Int,Int)) => {
              target == gc._1 + gc._2;
            }; 
          } else {
            error("Badly formatted extractReadsByMetric parameter. Unreciognized metric: " + metricArray(0));
            (st : String, insertSize : Int, gc : (Int,Int)) => {
              false;
            }; 
          }
        });
        
        Some(f.tail.foldLeft(f.head)((soFar,curr) => {
          (st : String, insertSize : Int, gc : (Int,Int)) => {
            soFar(st,insertSize,gc) && curr(st,insertSize,gc);
          };
        }));
      }
      case None => {
        None;
      }
    }
    
    
    
    if(parallelFileRead){
      reportln("ERROR ERROR ERROR: parallel file read is NOT IMPLEMENTED AT THIS TIME!","warn");
      //runOnSeqFile_PAR(initialTimeStamp = initialTimeStamp, infile = infile, outfile = outfile, anno_holder = anno_holder, testRun = testRun, runFunc = runFunc, stranded = stranded, fr_secondStrand = fr_secondStrand, dropChrom = dropChrom, keepMultiMapped = keepMultiMapped, noMultiMapped = noMultiMapped, numThreads = numThreads, readGroup )
    } else {
      runOnSeqFile(initialTimeStamp = initialTimeStamp, infile = infile, outfile = outfile, anno_holder = anno_holder, 
          testRun = testRun || (! stopAfterNReads.isEmpty), 
          runFunc = runFunc, stranded = stranded, fr_secondStrand = fr_secondStrand, dropChrom = dropChrom, keepMultiMapped = keepMultiMapped, 
          noMultiMapped = noMultiMapped, numThreads = numThreads, readGroup , minMAPQ = minMAPQ, geneKeepFunc = geneKeepFunc, 
          isSingleEnd = isSingleEnd, unsorted = unsorted, maxReadLength = maxReadLength,
          seqReadCt = seqReadCt, rawfastq = rawfastq,
          chromSizes = chromSizes, trackTitlePrefix = trackTitlePrefix,
          outdir = outdir,
          generateMultiPlot = generateMultiPlot,
          generateSeparatePlots=generateSeparatePlots,
          generatePdfReport=generatePdfReport,
          extractReadsFunction = extractReadsFunction,
          adjustPhredScore = adjustPhredScore,
          maxPhredScore = maxPhredScore,
          calcDetailedGeneCounts = runFunc.contains("calcDetailedGeneCounts"),
          isWES = isWES,checkForAlignmentBlocks=checkForAlignmentBlocks,
          summaryFileName =summaryFileName,
          targetRegionBed=targetRegionBed,
          stopAfterNReads=stopAfterNReads,
          genomeFA = genomeFA,dropOnTarget=dropOnTarget,keepOnTarget=keepOnTarget,randomSubsample=randomSubsample,
          genomeBufferSize = genomeBufferSize,
          outfilePrefix = outfilePrefix,
          prefilterImproperPairs = prefilterImproperPairs)
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
                   rawfastq : Option[List[String]],
                   chromSizes : Option[String],
                   trackTitlePrefix : String,
                   outdir : String,
                   generateMultiPlot : Boolean,
                   generateSeparatePlots : Boolean,
                   generatePdfReport : Boolean,
                   extractReadsFunction : Option[((String, Int, (Int,Int)) => Boolean)],
                   adjustPhredScore : Int,
                   maxPhredScore : Int,
                   calcDetailedGeneCounts : Boolean,
                   isWES : Boolean,
                   checkForAlignmentBlocks : Boolean,
                   summaryFileName : String,
                   targetRegionBed : Option[String],
                   stopAfterNReads : Option[Int],
                   genomeFA : Option[List[String]],
                   dropOnTarget : Boolean,
                   keepOnTarget : Boolean,
                   randomSubsample : Option[Double],
                   genomeBufferSize : Int,
                   outfilePrefix : String,
                   prefilterImproperPairs : Boolean
                   ){

    
    reportln("Checking first " + SAM_PEEK_LINECT + " reads. Checking SAM file for formatting errors...","note");
    val reader : SAMFileReader = if(infile == "-"){
      new SAMFileReader(System.in);
    } else {
      new SAMFileReader(new File(infile));
    }
    
    val peekCt = SAM_PEEK_LINECT;
    val testRunLineCt = stopAfterNReads.getOrElse(SAM_TESTRUN_LINECT);

    val COMPLETED_OK_FILEPATH = outfile + COMPLETED_OK_FILENAME;
    val COMPLETED_WARN_FILEPATH = outfile + COMPLETED_WARN_FILENAME;
    val (samFileAttributes, recordIter) = initSamRecordIterator(reader, peekCt);
    
    val maxObservedReadLength = samFileAttributes.readLength;
    val readLength = if(maxReadLength.isEmpty) maxObservedReadLength else maxReadLength.get;
    val isSortedByNameLexicographically = samFileAttributes.isSortedByNameLexicographically;
    val isSortedByPosition = samFileAttributes.isSortedByPosition;
    val isDefinitelyPairedEnd = samFileAttributes.isDefinitelyPairedEnd;
    val minReadLength = samFileAttributes.minReadLength;
    
    //Debugging info:
    reportln("   Stats on the first "+peekCt+ " reads:","debug");
    reportln("        Num Reads Primary Map:    " + samFileAttributes.numPeekReadsMapped,"debug");
    reportln("        Num Reads Paired-ended:   " + samFileAttributes.numPeekReadsPaired,"debug");
    reportln("        Num Reads mapped pair:    " + samFileAttributes.numPeekReadsPairMapped,"debug");
    reportln("        Num Pair names found:     " + samFileAttributes.numPeekPairs,"debug");
    reportln("        Num Pairs matched:        " + samFileAttributes.numPeekPairsMatched,"debug");
    reportln("        Read Seq length:          " + samFileAttributes.simpleMinReadLength + " to " + samFileAttributes.simpleMaxReadLength,"debug");
    reportln("        Unclipped Read length:    " + samFileAttributes.minReadLength + " to " + samFileAttributes.readLength,"debug");
    reportln("        Final maxReadLength:      " + readLength,"debug");
    reportln("        maxPhredScore:            " + samFileAttributes.maxObservedQual,"debug");
    reportln("        minPhredScore:            " + samFileAttributes.minObservedQual,"debug");
    
    if(samFileAttributes.maxObservedQual > maxPhredScore){
            reportln("WARNING WARNING WARNING: \n"+
               "   SAM format check:\n"+
               "      Phred Qual > "+maxPhredScore+"!\n"+
               "      You will need to set either --adjustPhredScores or --maxPhredScores\n"+
               "      in order to compute Phred quality metrics! QoRTs WILL throw an error\n"+
               "      if quality metrics are attempted!","warn");
    }
    
    if(readLength != minReadLength){ reportln("NOTE: Read length is not consistent.\n"+
                                             "   In the first "+peekCt+" reads, read length varies from "+minReadLength+" to " +maxObservedReadLength+" (param maxReadLength="+readLength+")\n"+
                                             "Note that using data that is hard-clipped prior to alignment is NOT recommended, because this makes it difficult (or impossible) "+
                                             "to determine the sequencer read-cycle of each nucleotide base. This may obfuscate cycle-specific artifacts, trends, or errors, the detection of which is one of the primary purposes of QoRTs!"+
                                             "In addition, hard clipping (whether before or after alignment) removes quality score data, and thus quality score metrics may be misleadingly optimistic. "+
                                             "A MUCH preferable method of removing undesired sequence is to replace such sequence with N's, which preserves the quality score and the sequencer cycle information.","note")}
    if((readLength != minReadLength) & (maxReadLength.isEmpty)){
      reportln("   WARNING WARNING WARNING: Read length is not consistent, AND \"--maxReadLength\" option is not set!\n"+
               "      QoRTs has ATTEMPTED to determine the maximum read length ("+readLength+").\n"+
               "      It is STRONGLY recommended that you use the --maxReadLength option \n"+
               "      to set the maximum possible read length, or else errors may occur if/when \n"+
               "      reads longer than "+readLength+ " appear.","warn")
    }
    
    val variableReadLen = readLength != minReadLength || (! maxReadLength.isEmpty);
     
    if(samFileAttributes.allReadsMarkedPaired & isSingleEnd) reportln("   WARNING WARNING WARNING! Running in single-end mode, but reads appear to be paired-end! Errors may follow.\n"+
                                                                      "           Strongly recommend removing the '--singleEnded' option!","warn");
    if(samFileAttributes.allReadsMarkedSingle & (! isSingleEnd)) reportln("   WARNING WARNING WARNING! Running in paired-end mode, but reads appear to be single-end! Errors may follow.\n"+
                                                                          "           Strongly recommend using the '--singleEnded' option","warn");
    if(samFileAttributes.mixedSingleAndPaired) reportln("   WARNING WARNING WARNING! Data appears to be a mixture of single-end and paired-end reads!\n"+
                                                        "           QoRTs was not designed to function under these conditions. Errors may follow!","warn");
    
    
    /*
     * Start Fastq Read:
     */
    //val FQSU = if(! rawfastq.isEmpty) new fqcStdUtils(readLength,isSingleEnd) else new NullFqcUtility();
    
    
    
    if(! isSingleEnd){ 
      reportln("   Note: Data appears to be paired-ended.","debug");
      
      var sortWarning = false;
      
      if( ! (samFileAttributes.perfectPairing || isSortedByPosition) ){
        reportln("   WARNING: Reads do not appear to be sorted by coordinate or by name. Sorting input data is STRONGLY recommended, but not technically required.","warn");
        sortWarning = true;
      }
      if(samFileAttributes.numPeekPairsMatched == 0){
        reportln("   Warning: Have not found any matched read-pairs in the first "+peekCt+" reads. Is data paired-end? Is data sorted?","warn"); 
        if(samFileAttributes.malformedPairNameCt > 0){
          reportln("   WARNING: No read-pairs found, but there are reads that match exactly\n"+
                   "            except for the last character, which is \"1\" in one read \n"+
                   "            and \"2\" in the other. This may indicate a malformed SAM \n"+
                   "            file in which the read-pairs are named with their readID \n"+
                   "            rather than read-pair ID. In standard SAM files, paired \n"+
                   "            reads MUST have the EXACT SAME first column.",
                   "warn");
        }
        sortWarning = true;
      }
      if( (! isDefinitelyPairedEnd)){ 
        reportln("   Warning: Have not found any matched read pairs in the first "+peekCt+" reads. Is data paired-end? Use option --singleEnd for single-end data.","warn");
        sortWarning = true;
      }
      if( isSortedByPosition & (! unsorted )){ 
        reportln("   Warning: SAM/BAM file appears to be sorted by read position, but you are running in --nameSorted mode.\n"+
                 "            If this is so, you should probably omit the '--nameSorted' option, as errors may follow.","warn"); 
        sortWarning = true;
      }
      
      val isOkNote = if(! sortWarning){"(This is OK)."} else {""}
      if(samFileAttributes.perfectPairing){
        reportln("   Sorting Note: Reads appear to be grouped by read-pair, probably sorted by name"+isOkNote,"note");
      } else {
        reportln("   Sorting Note: Reads are not sorted by name "+isOkNote,"note");
      }
      if(isSortedByPosition){
        reportln("   Sorting Note: Reads are sorted by position "+isOkNote,"note");
      } else {
        if(runFunc.contains("referenceMatch")){
          reportln("   Sorting Note: Reads are NOT sorted by position! Position-sorting is REQUIRED for the referenceMismatch module to function properly. Errors may follow.","note");
        } else {
          reportln("   Sorting Note: Reads are not sorted by position "+isOkNote,"note");
        }
      }
      
      //Samtools sorts in an odd way! Delete name sort check:
      //if( ((! isSortedByNameLexicographically) & (! unsorted ))) reportln("Note: SAM/BAM file does not appear to be sorted lexicographically by name (based on the first "+peekCt+" reads). It is (hopefully) sorted by read name using the samtools method.","debug");
    }
    if(internalUtils.Reporter.hasWarningOccurred()){
      reportln("Done checking first " + SAM_PEEK_LINECT + " reads. WARNINGS FOUND!","note");
    } else{
      reportln("Done checking first " + SAM_PEEK_LINECT + " reads. No major problems detected!","note");
    }
    
    val pairedIter : Iterator[(SAMRecord,SAMRecord)] = 
      if(isSingleEnd){
        if(testRun) samRecordPairIterator_withMulti_singleEnd(recordIter, true, testRunLineCt) else samRecordPairIterator_withMulti_singleEnd(recordIter);
      } else {
        if(unsorted){
          if(runFunc.contains("referenceMatch")){
             if(testRun) samRecordPairIterator_resorted(recordIter, true, testRunLineCt, prefilterImproperPairs = prefilterImproperPairs) else samRecordPairIterator_resorted(recordIter, prefilterImproperPairs = prefilterImproperPairs)
          } else {
             if(testRun) samRecordPairIterator_unsorted(recordIter, true, testRunLineCt, prefilterImproperPairs = prefilterImproperPairs) else samRecordPairIterator_unsorted(recordIter, prefilterImproperPairs = prefilterImproperPairs)
             //if(testRun) samRecordPairIterator_resorted(recordIter, true, testRunLineCt, prefilterImproperPairs = prefilterImproperPairs) else samRecordPairIterator_resorted(recordIter, prefilterImproperPairs = prefilterImproperPairs)
          }
        // Faster noMultiMapped running is DEPRECIATED!
        //} else if(noMultiMapped){
        //  if(testRun) samRecordPairIterator(recordIter, true, 200000) else samRecordPairIterator(recordIter)
        } else {
          if(testRun) samRecordPairIterator_withMulti(recordIter, true, testRunLineCt, prefilterImproperPairs = prefilterImproperPairs) else samRecordPairIterator_withMulti(recordIter, prefilterImproperPairs = prefilterImproperPairs)
        }
      }
    
    reportln("SAMRecord Reader Generated. Read length: "+readLength+".","note");
    standardStatusReport(initialTimeStamp);
    
    val coda : Array[Int] = internalUtils.commonSeqUtils.getNewCauseOfDropArray;
    val coda_options : Array[Boolean] = internalUtils.commonSeqUtils.CODA_DEFAULT_OPTIONS.toArray;
    if(isSingleEnd) CODA_SINGLE_END_OFF_OPTIONS.foreach( coda_options(_) = false );
    if(keepMultiMapped) coda_options(internalUtils.commonSeqUtils.CODA_NOT_UNIQUE_ALIGNMENT) = false;
    if(! readGroup.isEmpty) coda_options(internalUtils.commonSeqUtils.CODA_NOT_MARKED_RG) = true;
    if(checkForAlignmentBlocks) coda_options(internalUtils.commonSeqUtils.CODA_NO_ALN_BLOCKS) = true;
    
    //"writeKnownSplices","writeNovelSplices","writeSpliceExon", "writeDESeq","writeDEXSeq","writeGenewiseGeneBody"
    //  final val QC_FUNCTION_LIST : Seq[String] = Seq("InsertSize","NVC","CigarOpDistribution","QualityScoreDistribution","GCDistribution","GeneCounts","JunctionCounts");
    
    //FASTQ READTHROUGH:
    
    val (fqLineCount,fqUtilSeq) = rawfastq match {
       case Some(fqfile) => {
         reportln("Starting fastq readthrough.","debug");
         
         val fq1 = new FastqReader(fqfile(0), maxPhredScore = maxPhredScore, adjustPhredScore = adjustPhredScore, strict  = true, allowStdin = false);
         val fq2 = if((! isSingleEnd) && fqfile.length > 1) new FastqReader(fqfile(1), maxPhredScore = maxPhredScore, adjustPhredScore = adjustPhredScore, strict  = true, allowStdin = false) else new EmptyFastqReader;
         val SE = isSingleEnd || fqfile.length < 2
         
         val fqLines = if(testRun) presetProgressReporters.wrapIterator_readPairs(fq1.zip(fq2),verbose=true,cutoff=testRunLineCt) 
                       else        presetProgressReporters.wrapIterator_readPairs(fq1.zip(fq2),verbose=true);
         
         val fqcGC =    if(runFunc.contains("GCDistribution"))           new fqcGC(readLength, SE) else fqcUtility.NullFqcUtility;
         val fqcQS =    if(runFunc.contains("QualityScoreDistribution")) new fqcQualScores(readLength, SE, maxPhredScore) else fqcUtility.NullFqcUtility;
         val fqcNVC =   if(runFunc.contains("QualityScoreDistribution")) new fqcNVC(readLength, SE) else fqcUtility.NullFqcUtility;
         val fqUtilSeq : Vector[fqcUtility[Any]] = if(runFunc.contains("fastqUtils")) Vector(fqcGC,fqcQS,fqcNVC) else Vector();
         
         var lnct = 0;
         for((f1,f2) <- fqLines){
           fqUtilSeq.foreach(u => { u.runOnReadPair(f1,f2,lnct) });
           lnct += 1;
           //if(lnct % 100000 == 0){ if(lnct % 500000 == 0) { if(lnct % 1000000 == 0) report(".["+lnct+" FQ blocks processed] [Time: "+getDateAndTimeString+"] \n","progress"); else report(". ","progress"); } else { report(".","progress");}}
         }
         reportln("Finished fastq readthrough.","debug");
         (Some(lnct),fqUtilSeq);
       }
       case None => {
         (None,Seq[qcUtils.fqcUtility[Nothing]]())
       }
    }
    
    val inputReadCt : Option[Int] = seqReadCt match {
      case Some(ct) => Some(ct);
      case None => {
        fqLineCount match {
          case Some(ct) => Some(ct);
          case None => None;
        }
      }
    }

    val dumpfile = if(runFunc.contains("testDataDump")) Some(outfile + ".testDataDump.txt") else None;
    
    val qcGGC:  QCUtility[String] =   if(runFunc.contains("GeneCalcs"))                 new qcGetGeneCounts(stranded,fr_secondStrand,anno_holder,coda,coda_options,40, runFunc.contains("FPKM"), runFunc.contains("writeGenewiseGeneBody"), runFunc.contains("writeDESeq"), runFunc.contains("writeGeneCounts"), runFunc.contains("writeGeneBody"), runFunc.contains("writeBiotypeCounts"), geneKeepFunc, calcDetailedGeneCounts= calcDetailedGeneCounts, writeGeneBodyIv=runFunc.contains("writeGeneBodyIv")) else QCUtility.getBlankStringUtil;
    val qcIS :  QCUtility[Int]    =   if(runFunc.contains("InsertSize"))                new qcInnerDistance(anno_holder, stranded, fr_secondStrand, readLength,variableReadLen=variableReadLen, jumpSplices = ! isWES)        else QCUtility.getBlankIntUtil;
    val qcCS :  QCUtility[Unit]   =   if(runFunc.contains("NVC"))                       new qcNVC(isSingleEnd, readLength, runFunc.contains("writeClippedNVC"))                     else QCUtility.getBlankUnitUtil;
    val qcJD :  QCUtility[Unit]   =   if(runFunc.contains("CigarOpDistribution"))       new qcCigarDistribution(isSingleEnd, readLen = readLength,variableReadLen=variableReadLen)                                            else QCUtility.getBlankUnitUtil;
    val qcQSC : QCUtility[Unit]   =   if(runFunc.contains("QualityScoreDistribution"))  new qcQualityScoreCounter(isSingleEnd, readLength, maxPhredScore, adjustPhredScore) else QCUtility.getBlankUnitUtil;
    val qcGC :  QCUtility[(Int,Int)]= if(runFunc.contains("GCDistribution"))      new qcGCContentCount(isSingleEnd, readLength)                                               else QCUtility.getBlankIntPairUtil;
    val qcJC :  QCUtility[Unit]   =   if(runFunc.contains("JunctionCalcs"))             new qcJunctionCounts(anno_holder, stranded, fr_secondStrand, runFunc.contains("writeDEXSeq"), runFunc.contains("writeSpliceExon"), runFunc.contains("writeKnownSplices"), runFunc.contains("writeNovelSplices"), runFunc.contains("annotatedSpliceExonCounts"))                   else QCUtility.getBlankUnitUtil;
    val qcST :  QCUtility[String] =   if(runFunc.contains("StrandCheck"))               new qcStrandTest(isSingleEnd, anno_holder, stranded, fr_secondStrand)                       else QCUtility.getBlankStringUtil;
    val qcCC :  QCUtility[String] =   if(runFunc.contains("chromCounts"))               new qcChromCount(isSingleEnd, fr_secondStrand)                                              else QCUtility.getBlankStringUtil;
    val qcCM :  QCUtility[Unit]   =   if(runFunc.contains("cigarMatch"))                new qcCigarMatch(readLength)                                                   else QCUtility.getBlankUnitUtil;
    val qcCLC : QCUtility[(Int,Int)]= if(runFunc.contains("cigarLocusCounts"))          new qcCigarLocusCounts(stranded,fr_secondStrand,true,true,5)                                    else QCUtility.getBlankIntPairUtil;
    val qcWIG : QCUtility[Unit]   =  if(runFunc.contains("makeWiggles"))                 new fileConversionUtils.bamToWiggle.QcBamToWig(trackTitlePrefix,
                                                                                                                   chromSizes.get,false,100,
                                                                                                                   isSingleEnd, stranded, fr_secondStrand, 
                                                                                                                   1.0, true, true, true, None, "") else QCUtility.getBlankUnitUtil;
    val qcMAT : QCUtility[Unit] =     if(runFunc.contains("mismatchEngine"))  new qcOverlapMatch(readLength,dumpfile,
                                                                                               maxQualScore=maxPhredScore, adjustPhredScore=adjustPhredScore,isSingleEnd=isSingleEnd,
                                                                                               calcReferenceMatch=runFunc.contains("referenceMatch"),
                                                                                               calcOverlapMismatch = runFunc.contains("overlapMatch"),
                                                                                               genomeFa=genomeFA,genomeBufferSize=genomeBufferSize)  else QCUtility.getBlankUnitUtil;
    
    val outputFastqFile = if(runFunc.contains("outputFastqFile")) Some(outfile) else None;
    val outputBAMFile   = if(runFunc.contains("outputBAMFile")) Some(outfile) else None;

    val qcMU : QCUtility[Unit] = if(true) new qcMinorUtils(isSingleEnd =isSingleEnd, stranded =stranded, fr_secondStrand_bool =fr_secondStrand,readLen =readLength,
                                              writeChromList =true,
                                              writeReadLengthDistro = runFunc.contains("readLengthDistro"),
                                              outputFastqFile = outputFastqFile,
                                              outputBAMFile = outputBAMFile) else QCUtility.getBlankUnitUtil;
    
    if(runFunc.contains("calcOnTarget") && targetRegionBed.isEmpty) error("FATAL ERROR: CANNOT CALCULATE ON-TARGET REGION WITHOUT TARGET BED FILE! Use parameter --targetRegionBed");
    if(targetRegionBed.isEmpty && (dropOnTarget || keepOnTarget)) error("FATAL ERROR: CANNOT USE dropOnTarget or keepOnTarget without specified target region bed file. Use parameter --targetRegionBed");
     
    //val iv1 = getGenomicIntervalsFromRead(r1, false,false).toVector;
    //val iv2 = getGenomicIntervalsFromRead(r2, false,false).toVector;
    //val ivB = mergeGenomicIntervalSets(iv1 ++ iv2);
    //val targets = getTargetSet(ivB);
    
    val qcOT : QCUtility[Boolean] = if(runFunc.contains("calcOnTarget"))                new qcOnTargetRegion(anno_holder , isSingleEnd ) else new QCUtility.blankBooleanQCUtility;
    
    //Special super secret utility: 
    //val qcFCC :  QCUtility[Unit]   =   if(runFunc.contains("featureComboCount"))        new qcFeatureComboCt(anno_holder, stranded, fr_secondStrand,  outfile = outfile)                   else QCUtility.getBlankUnitUtil;
        
    val extractWriter : Option[WriterUtil] = extractReadsFunction match {
      case Some(f) => Some(openWriter(outfile + ".extractedReads.sam"));
      case None => None;
    }
    val writeExtractedReads = ((r1 : SAMRecord, r2 : SAMRecord) => {
      extractWriter match {
        case Some(w) => {
          w.write(r1.getSAMString+"\n");
          w.write(r2.getSAMString+"\n");
        }
        case None => {
          reportln("Warning: impossible state.\n","warn");
          //Do nothing? This should never happen.
        }
      }
    });
    
    
    //val qcALL = parConvert(Vector(qcGGC, qcIS, qcCS, qcJD, qcQSC, qcGC, qcJC, qcST, qcCC, qcCM), numThreads);
    val qcALL = Vector(qcGGC, qcIS, qcCS, qcJD, qcQSC, qcGC, qcJC, qcST, qcCC, qcCM, qcCLC,qcOT,qcMAT,qcMU);
    
    reportln("QC Utilities Generated!","note");
    standardStatusReport(initialTimeStamp);
    //GenomicArrayOfSets.printGenomicArrayToFile("TEST.OUT.gtf",geneArray);
    var readNum = 0;
    var useReadNum = 0;
    var keptMultiMappedCt = 0;
    var minObsReadLength = readLength;
    var maxObsReadLength = readLength;
    val samIterationTimeStamp = TimeStampUtil();
    //while(pairedIter.hasNext){
    for(pair <- pairedIter){
    //for((pair,readNum) <- numberedIter){
    //  val pair = pairedIter.next;      
      val (r1,r2) = pair;
      readNum += 1;
      try{
        if(internalUtils.commonSeqUtils.useReadPair(r1,r2,coda, coda_options, dropChrom, readGroup, minMAPQ, dropOnTarget=dropOnTarget,keepOnTarget=keepOnTarget,randomSubsample=randomSubsample,anno_holder=anno_holder)){
            val gene = qcGGC.runOnReadPair(r1,r2,readNum);
            
            if( geneKeepFunc(gene) ){
            //if( geneListKeep.get.contains(gene) ){
              minObsReadLength = math.min(minObsReadLength, math.min(r1.getReadLength(), r2.getReadLength()));
              maxObsReadLength = math.max(maxObsReadLength, math.max(r1.getReadLength(), r2.getReadLength()));
              
              val onTarget = qcOT.runOnReadPair(r1,r2,readNum);
              val ins = qcIS.runOnReadPair(r1,r2,readNum);
              qcMU.runOnReadPair(r1,r2,readNum);
              qcCS.runOnReadPair(r1,r2,readNum);
              qcJD.runOnReadPair(r1,r2,readNum);
              qcQSC.runOnReadPair(r1,r2,readNum);
              val gc = qcGC.runOnReadPair(r1,r2,readNum);
              qcJC.runOnReadPair(r1,r2,readNum);
              val st = qcST.runOnReadPair(r1,r2,readNum);
              qcCC.runOnReadPair(r1,r2,readNum);
              qcCM.runOnReadPair(r1,r2,readNum);
              qcCLC.runOnReadPair(r1,r2,readNum);
              qcWIG.runOnReadPair(r1,r2,readNum);
              qcMAT.runOnReadPair(r1,r2,readNum);
              //qcFCC.runOnReadPair(r1,r2,readNum);
              useReadNum += 1;
              
              val extract = extractReadsFunction match {
                case Some(erf) => {
                  erf(st,ins,gc);
                }
                case None => {
                  false;
                }
              }
              if(extract) {
                writeExtractedReads(r1,r2);
              }
              
              if(internalUtils.commonSeqUtils.isReadMultiMapped(r1, minMAPQ) || internalUtils.commonSeqUtils.isReadMultiMapped(r2, minMAPQ)){
                keptMultiMappedCt += 1;
              }
            }
        }
      } catch {
        case e : Exception => {
          internalUtils.Reporter.reportln("Fatal error thrown for read: "+r1.getReadName(),"note");
          throw e;
        }
      }
    }
    
    extractWriter match {
      case Some(w) => {
        w.close();
      }
      case None => {
        //Do nothing.
      }
    }
    
    reportln("Finished reading SAM. Read: " + readNum + " reads/read-pairs.","note");
    reportln("Finished reading SAM. Used: " + useReadNum + " reads/read-pairs.","note");
    standardStatusReport(initialTimeStamp);
    
    val outputIterationTimeStamp = TimeStampUtil();
    report("> Read Stats:\n" + stripFinalNewline(indentifyLines(internalUtils.commonSeqUtils.causeOfDropArrayToString(coda, coda_options),">   ")),"note");
    
    val summaryWriter = openWriter(outfile + ".summary.txt");
    val docWriter = if(runFunc.contains("writeDocs")) new DocWriterUtil(openWriter(outfile + ".documentation.txt")) else null;
    //val summaryWriter = SummaryWriter(outfile + summaryFileName);
    
    val strandedCode = if(! stranded){ 0 } else {if(fr_secondStrand) 2; else 1;}

    summaryWriter.write("FIELD\tCOUNT\tDESC\n");
    summaryWriter.write("Stranded_Rule_Code	"+strandedCode+"\t"+"Code for the strandedness rule used. 0 if data is unstranded, 1 if data is fr_firstStrand, 2 if data is fr_secondStrand."+"\n");
    summaryWriter.write(internalUtils.commonSeqUtils.causeOfDropArrayToStringTabbed(coda, coda_options));
    
    summaryWriter.write("KEPT_NOT_UNIQUE_ALIGNMENT	"+keptMultiMappedCt+"\t"+"Number of reads or read-pairs kept despite not being uniquely aligned."+"\n");
    summaryWriter.write("minObservedReadLength\t"+minObsReadLength +"\t"+"The base-pair length of the smallest read"+ "\n");
    summaryWriter.write("maxObservedReadLength\t"+maxObsReadLength +"\t"+"The base-pair length of the largest read"+ "\n");
    summaryWriter.write("maxLegalPhredScore\t"+maxPhredScore +"\t"+"The largest observed PHRED score."+ "\n");
    
    if(isSingleEnd){
      summaryWriter.write("IS_SINGLE_END	1"+"\t"+"0 if data is paired-ended, 1 if data is single-ended"+"\n");
    } else {
      summaryWriter.write("IS_SINGLE_END	0"+"\t"+"0 if data is paired-ended, 1 if data is single-ended"+"\n");
    }
    
    if(inputReadCt.isEmpty){
      reportln("Pre-alignment read count unknown (Set --seqReadCt or --rawfastq)","note");
      summaryWriter.write("PREALIGNMENT_READ_CT\t-1"+"\t"+"The number of reads found pre-alignment. Can be set using --seqReadCt or --rawfastq. -1 if unknown."+"\n");
    } else {
      reportln("Pre-alignment read count: "+inputReadCt.get,"note");
      summaryWriter.write("PREALIGNMENT_READ_CT\t"+inputReadCt.get+"\t"+"The number of reads found pre-alignment. Can be set using --seqReadCt or --rawfastq. -1 if unknown."+"\n");
    }
    
    val iterationMinutes = (outputIterationTimeStamp.compareTo(samIterationTimeStamp) / 1000).toDouble / 60.toDouble;
    val minutesPerMillion = iterationMinutes / (readNum.toDouble / 1000000.toDouble);
    val minutesPerMillionPF = iterationMinutes / ((coda(internalUtils.commonSeqUtils.CODA_READ_PAIR_OK)).toDouble / 1000000.toDouble);
    
    summaryWriter.write("BENCHMARK_MinutesOnSamIteration	" + "%1.2f".format(iterationMinutes) +"\t"+"The number of minutes spent on the SAM iteration step."+ "\n");
    summaryWriter.write("BENCHMARK_MinutesPerMillionReads	" + "%1.2f".format(minutesPerMillion) +"\t"+"The number of minutes per million reads spent on the SAM iteration step"+ "\n");
    summaryWriter.write("BENCHMARK_MinutesPerMillionGoodReads	" + "%1.2f".format(minutesPerMillionPF) +"\t"+"The number of minutes per million reads that passed the initial filtering step."+ "\n");
    
    if(useReadNum == 0) reportln("WARNING WARNING WARNING: Zero \"usable\" reads found! This could be due to a number of factors: \n"+
          "If the reads were not aligned via one of the standard RNA-Seq aligners such as RNA-STAR or TopHat/TopHat2, then "+
          "the alignments might not use the common convention of using MAPQ to indicate multi-mapping status. \n"+
          "RNA-STAR and TopHat both mark multi-mapped reads by "+
          "assigning them a MAPQ score of less than 255. By default QoRTs ignores these multi-mapped reads. "+
          "You can deactivate this filtering step using the \"--keepMultiMapped\" option.\n"+
          "Note: Alignment via BowTie, BowTie2, or other non-spliced aligners is NOT RECOMMENDED for RNA-Seq data. \n"+
          "If the data was aligned using such methods, it is strongly recommended that it be realigned using a splice-aware aligner.\n"+
          "\n"+
          "Continuing with output execution. Errors will likely follow..."+
          "\n","warn");

    reportln("Writing Output...","note");
    val qcAllOut = qcALL ++ fqUtilSeq;
    
    
    qcAllOut.seq.foreach( _.writeOutput(outfile, summaryWriter, docWriter = docWriter) );
    
    qcWIG.writeOutput(outfile + ".wiggle", summaryWriter, docWriter = docWriter);
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
    
    summaryWriter.write("READ_LENGTH	"+readLength+"\t"+"The read length."+"\n");
    
    if(internalUtils.Reporter.hasWarningOccurred()){
      summaryWriter.write("COMPLETED_WITHOUT_WARNING\t0\t"+"0 if complete without throwing any warnings. 1 if warnings were thrown."+"\n");
      reportln("QoRTs completed WITH WARNINGS! See log for details.","warn");
      val completedWarnWriter = openWriter(COMPLETED_WARN_FILEPATH);
      completedWarnWriter.write("# Note: if this file EXISTS, then QoRTs QC completed WITH WARNINGS. Warning messages follow:\n");
      completedWarnWriter.write(internalUtils.Reporter.getWarnings+"\n");
      completedWarnWriter.close();
    } else {
      summaryWriter.write("COMPLETED_WITHOUT_WARNING	1\t"+"0 if complete without throwing any warnings. 1 if warnings were thrown."+"\n");
      reportln("QoRTs QC complete with no problems.","note");
    }
    
    summaryWriter.write("QoRTs_initTimeStamp\t"+initialTimeStamp.ts+"\t"+"Time stamp for when QoRTs QC began."+"\n");
    summaryWriter.write("QoRTs_samDoneTimeStamp\t"+outputIterationTimeStamp.ts+"\t"+"Time stamp for when QoRTs QC finished SAM iteration."+"\n");
    summaryWriter.write("QoRTs_majorVer\t"+runner.runner.QORTS_MAJOR_VERSION+"\t"+"QoRTs major version number"+"\n");
    summaryWriter.write("QoRTs_minorVer\t"+runner.runner.QORTS_MINOR_VERSION+"\t"+"QoRTs minor version number"+"\n");
    summaryWriter.write("QoRTs_patchVer\t"+runner.runner.QORTS_PATCH_VERSION+"\t"+"QoRTs patch version number"+"\n");
    summaryWriter.write("QoRTs_compileTimeStamp\t"+runner.runner.QORTS_COMPILE_TIME+"\t"+"The timestamp for when the version of QoRTs was built."+"\n");
    
    summaryWriter.write("COMPLETED_WITHOUT_ERROR\t1\t1 if QoRTs completed without errors. If QoRTs encountered an error, this file should not exist.\n");
    
    
    if(docWriter != null) docWriter.close();
    close(summaryWriter);
    
    if(generateMultiPlot | generateSeparatePlots | generatePdfReport){
      report("Generating plots...","progress")
      fileConversionUtils.generatePlotsWithR.generateSimplePlots(
             qcdir = outdir,
             qcprefix = outfilePrefix,
             uniqueID = trackTitlePrefix,
             makePng = generateMultiPlot,
             makePdf = generatePdfReport,
             makeSeparatePngs = generateSeparatePlots
           );
    }
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















