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
     override def priority = 24;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "bamToWiggle", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "Generates '.wig' wiggle files, suitable for use with the UCSC genome browser or "+
                        "similar tools. Wiggle files contain depth-of-coverage counts for all equally-sized "+
                        "windows across the entire genome. These depth-of-coverage counts can be either be by read (the default) or "+
                        "by read-pair."+
                        ""+
                        ""+
                        "",
          argList = 
           
                    new BinaryArgument[Double](   name = "sizefactor",
                                                        arg = List("--sizefactor"),  
                                                        valueName = "float", 
                                                        argDesc = "The size factor, for normalization. Defaults to 1. If set, then the count of every cell will be divided by the given value.", 
                                                        defaultValue = Some(1.0)
                                                        ) ::
                    /*new UnaryArgument( name = "verbose",
                                         arg = List("--verbose"), // name of value
                                         argDesc = "Flag to indicate that debugging information and extra progress information should be sent to stderr." // description
                                       ) ::
                    new UnaryArgument( name = "quiet",
                                         arg = List("--quiet","-s"), // name of value
                                         argDesc = "Flag to indicate that only errors and warnings should be sent to stderr." // description
                                       ) ::*/
                    new UnaryArgument( name = "coordSorted",
                                         arg = List("--coordSorted"), // name of value
                                         argDesc = "Flag to indicate that input bam file is coordinate-sorted, rather than name-sorted. "+
                                                   "Note that QoRTs will take longer to run and use more memory in this mode. To improve performance, sort the data by name prior to using of QoRTs. "+
                                                   "In addition, if an (extremely) large fraction of the read-pairs are "+
                                                   "mapped to extremely distant loci (or different chromosomes), then memory issues may arise. However, this should not be a problem with most datasets. "+
                                                   "Technically this function will also allow QoRTs to work on unsorted bam files, but this is STRONGLY not recommended, as memory usage will by greatly increased." // description
                                       ) ::
                    new BinaryArgument[Int](name = "minMAPQ",
                                           arg = List("--minMAPQ"),  
                                           valueName = "num", 
                                           argDesc = "Filter reads for the given minimum MAPQ. Set to 0 to turn off mapq filtering.", 
                                           defaultValue = Some(255)
                                           ) :: 
                    new UnaryArgument(    name = "noGzipOutput",
                                         arg = List("--noGzipOutput"), // name of value
                                         argDesc = "Flag to indicate that output should NOT be gzipped." // description
                                       ) ::
                    new UnaryArgument(    name = "negativeReverseStrand",
                                         arg = List("--negativeReverseStrand"), // name of value
                                         argDesc = "Flag to indicate that the reverse strand should be counted in negative values. Can be useful for plotting stranded data on a single multiwig track, via trackhubs." // description
                                       ) ::
                    new UnaryArgument(    name = "countByReadPair",
                                         arg = List("--countByReadPair"), // name of value
                                         argDesc = "DEPRECIATED, nonfunctional." // description
                                       ) ::
                    new UnaryArgument(    name = "simpleCountByRead",
                                         arg = List("--simpleCountByRead"), // name of value
                                         argDesc = "Flag to indicate that this utility should count 'read coverage' rather than 'read-pair coverage'. "+
                                                   "If this flag is raised, then regions where a read-pair overlaps will be counted "+
                                                   "TWICE, one for each read. By default each read is counted separately, so a read-pair where the reads overlap "+
                                                   "will be double-counted. By default, overlapping regions of read-pairs are counted only once." // description
                                       ) ::
                    new UnaryArgument(    name = "testRun",
                                         arg = List("--testRun","-t"), // name of value
                                         argDesc = "Flag to indicate that only the first 100k reads should be read in, used for test runs." // description
                                       ) ::
                    new UnaryArgument(   name = "stranded", 
                                         arg = List("--stranded"), // name of value
                                         argDesc = "The stranded flag, used for strand-specific data. If this is set, then two wiggle files are created, one for each strand." // description
                                       ) ::
                    new UnaryArgument(   name = "fr_secondstrand", 
                                         arg = List("--stranded_fr_secondstrand"), // name of value
                                         argDesc = "If this option is set, the data is assumed to be a fr_secondstrand type stranded library. By default the assumed library type for stranded data is fr_firststrand." // description
                                       ) ::
                    new UnaryArgument(   name = "singleEnded", 
                                         arg = List("--singleEnded"), // name of value
                                         argDesc = "Flag for single-end data. Note that many other options do not apply in this case (for example: option --countPairsTogether does nothing in single-end mode)" 
                                       ) ::
                    new UnaryArgument( name = "keepMultiMapped",
                                         arg = List("--keepMultiMapped"), // name of value
                                         argDesc = "Flag to indicate that the tool should NOT filter out multi-mapped reads. Note that even with this flag raised this utility will still only use the 'primary' alignment location for each read. By default any reads that are marked as multi-mapped will be ignored entirely." // description
                                       ) ::
                    new UnaryArgument( name = "includeTrackDefLine",
                                         arg = List("--includeTrackDefLine"), // name of value
                                         argDesc = "DEPRECIATED, nonfunctional" // description
                                       ) ::
                    new UnaryArgument( name = "omitTrackDefLine",
                                         arg = List("--omitTrackDefLine"), // name of value
                                         argDesc = "Flag to indicate that a track definition line should NOT be included." // description
                                       ) ::
                    new BinaryArgument[String](
                                         name = "additionalTrackOptions", 
                                         arg = List("--additionalTrackOptions"), 
                                         valueName = "\"track options\"",  
                                         argDesc =  "More options for the track. For more information refer to the wiggle track definition on the UCSC genome browser website.",
                                         defaultValue = Some("")
                                        ) ::                                        
                    new BinaryOptionArgument[String](
                                         name = "readGroup", 
                                         arg = List("--readGroup"), 
                                         valueName = "readGroupName",  
                                         argDesc =  "If this option is set, all analyses will be restricted to ONLY reads that are tagged with the given "+
                                                    "readGroupName (using an RG tag). This can be used if multiple read-groups have already been combined "+
                                                    "into a single bam file, but you want to summarize each read-group separately."
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "rgbColor", 
                                         arg = List("--rgbColor"), 
                                         valueName = "rval,gval,bval",  
                                         argDesc =  "A comma delimited list of integers from 0-255, indicating the rgb color and alt-color for the wiggle file."+
                                                    "This is equivalent to using --additionalTrackOptions \"color=r,g,b altColor=r,g,b itemRgb=On\""
                                        ) ::
                    new UnaryArgument(   name = "noTruncate", 
                                         arg = List("--noTruncate"), // name of value
                                         argDesc = "The UCSC tool wigToBigWig only allows wiggle files in which every window is of equal size. "+
                                                   "This means that if the chromosome size is not divisible by the window size, a few bases are not "+
                                                   "counted on the end of the chromosome. Using this flag will cause this utility to not truncate off"+
                                                   " the last odd-sized window. However, be aware that this will mean that you cannot use the UCSC 'wigToBigWig' utility to"+
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
                                         name = "trackTitlePrefix",
                                         valueName = "trackTitlePrefix",
                                         argDesc = "The prefix for the title for the wiggle track. Track titles will be " // description
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
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
       
       if(out){
         bamToWiggle.run(
             parser.get[List[String]]("infiles"),
             parser.get[String]("outfilePrefix"),
             parser.get[String]("trackTitlePrefix"),
             parser.get[String]("chromLengthFile"),
             parser.get[Boolean]("noTruncate"),
             parser.get[Int]("windowSize"),
             parser.get[Boolean]("singleEnded"),
             parser.get[Boolean]("stranded"),
             parser.get[Boolean]("fr_secondstrand"),
             parser.get[Double]("sizefactor"),
             parser.get[Boolean]("testRun"),
             parser.get[Boolean]("noGzipOutput"),
             parser.get[Boolean]("negativeReverseStrand"),
             ! parser.get[Boolean]("simpleCountByRead"),
             parser.get[Boolean]("keepMultiMapped"),
             parser.get[Option[String]]("readGroup"),
             ! parser.get[Boolean]("omitTrackDefLine"),
             parser.get[Option[String]]("rgbColor"),
             parser.get[String]("additionalTrackOptions"),
             parser.get[Boolean]("coordSorted"),
             parser.get[Int]("minMAPQ")
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
  
  def run(infiles : List[String], outfilePrefix : String, trackName : String, chromLengthFile : String, noTruncate : Boolean, 
      windowSize : Int, isSingleEnd : Boolean, stranded : Boolean, fr_secondStrand : Boolean, sizeFactor : Double, testRun : Boolean, noGzipOutput : Boolean, 
      negativeReverseStrand : Boolean, countPairsTogether : Boolean,
      keepMultiMapped : Boolean, readGroup : Option[String],
      includeTrackDef : Boolean,
      rgbColor : Option[String],
      additionalTrackOptions : String,
      coordSorted : Boolean, minMAPQ : Int){
    
    //registerGlobalParam[Boolean]("noGzipOutput", noGzipOutput);
    internalUtils.optionHolder.OPTION_noGzipOutput = noGzipOutput;
    
    val initialTimeStamp = TimeStampUtil();
    
    if(isSingleEnd) reportln("QoRTs is Running in single-end mode.","note");
    else reportln("QoRTs is Running in paired-end mode.","note");
    
    if(coordSorted) reportln("QoRTs is Running in coordinate-sorted mode.","note");
    else reportln("QoRTs is Running in name-sorted mode.","note");
    
    val qcBTW : QcBamToWig = new QcBamToWig(trackName , chromLengthFile , noTruncate , windowSize , isSingleEnd, stranded , fr_secondStrand, sizeFactor, negativeReverseStrand, countPairsTogether , includeTrackDef, rgbColor, additionalTrackOptions);
    
    standardStatusReport(initialTimeStamp);
    val postSetupStamp = TimeStampUtil();
    reportln("> Time spent on setup:           " + TimeStampUtil.timeDifferenceFormatter(postSetupStamp.compareTo(initialTimeStamp)) +" ]","note");
    
    for(infile <- infiles){
      //runOnFile(infile , qcBTW , testRun, isSingleEnd, keepMultiMapped, readGroup);
      runOnFile2(infile, qcBTW, testRun, keepMultiMapped, readGroup, minMAPQ, isSingleEnd, coordSorted);
    }
    
    val postRunStamp = TimeStampUtil();
    //standardStatusReport(initialTimeStamp);
    
    qcBTW.writeOutput(outfilePrefix, null);
    
    val finalStamp = TimeStampUtil();
    standardStatusReport(postRunStamp);
    reportln("> Time spent on file output:     " + TimeStampUtil.timeDifferenceFormatter(finalStamp.compareTo(postRunStamp)) +" ]","note");
    reportln("> Total runtime:                 " + TimeStampUtil.timeDifferenceFormatter(finalStamp.compareTo(initialTimeStamp)) +" ]","note");
  }
    
  def runOnFile2(  
                   infile : String, 
                   qcBTW : QcBamToWig,
                   testRun  :Boolean, 
                   keepMultiMapped : Boolean, 
                   readGroup : Option[String],
                   minMAPQ : Int,
                   isSingleEnd : Boolean,
                   unsorted : Boolean){
    
    
    val peekCt = 2000;
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
    val readLength = maxObservedReadLength 
    val isSortedByNameLexicographically = samFileAttributes.isSortedByNameLexicographically;
    val isSortedByPosition = samFileAttributes.isSortedByPosition;
    val isDefinitelyPairedEnd = samFileAttributes.isDefinitelyPairedEnd;
    val minReadLength = samFileAttributes.minReadLength;
    
    //if(readLength != minReadLength){reportln("Note: Read length is not consistent. "+
    //                                         "In the first "+peekCt+" reads, read length varies from "+minReadLength+" to " +maxObservedReadLength+"!\n"+
    //                                         "Note that using data that is hard-clipped prior to alignment is NOT recommended, because this makes it difficult (or impossible) "+
    //                                         "to determine the sequencer read-cycle of each nucleotide base. This may obfuscate cycle-specific artifacts, trends, or errors, the detection of which is one of the primary purposes of QoRTs!"+
    //                                         "In addition, hard clipping (whether before or after alignment) removes quality score data, and thus quality score metrics may be misleadingly optimistic. "+
    //                                         "A MUCH preferable method of removing undesired sequence is to replace such sequence with N's, which preserves the quality score and the sequencer cycle information.","note")}
    //if((readLength != minReadLength) & (maxReadLength.isEmpty)){
    //  reportln("WARNING WARNING WARNING: Read length is not consistent, AND \"--maxReadLength\" option is not set!\n"+
    //          "QoRTs has ATTEMPTED to determine the maximum read length ("+readLength+").\n"+
    //           "It is STRONGLY recommended that you use the --maxReadLength option \n"+
    //           "to set the maximum possible read length, or else errors may occur if/when \n"+
    //           "reads longer than "+readLength+ " appear.","note")
    //}
    
    if(samFileAttributes.allReadsMarkedPaired & isSingleEnd) reportln("WARNING WARNING WARNING! Running in single-end mode, but reads appear to be paired-end! Errors may follow.\nStrongly recommend removing the '--isSingleEnd' option!","warn");
    if(samFileAttributes.allReadsMarkedSingle & (! isSingleEnd)) reportln("WARNING WARNING WARNING! Running in paired-end mode, but reads appear to be single-end! Errors may follow.\nStrongly recommend using the '--isSingleEnd' option","warn");
    if(samFileAttributes.mixedSingleAndPaired) reportln("WARNING WARNING WARNING! Data appears to be a mixture of single-end and paired-end reads! QoRTs was not designed to function under these conditions. Errors may follow!","warn");
    
    if(! isSingleEnd){
      if( (! isDefinitelyPairedEnd)){ reportln("Warning: Have not found any matched read pairs in the first "+peekCt+" reads. Is data paired-end? Use option --singleEnd for single-end data.","warn"); }
      if( isSortedByPosition & (! unsorted )){ reportln("Based on the first "+peekCt+" reads, SAM/BAM file appears to be sorted by read position. If this is so, you should probably use the \"--coordSorted\" option.","warn"); }
      if( ((! isSortedByPosition) & ( unsorted ))){ reportln("WARNING: You are using the \"--coordSorted\" option, but data does NOT appear to be sorted by read position (based on the first "+peekCt+" reads)! This is technically ok, but may cause QoRTs to use too much memory!","warn"); }
      //if( ((! isSortedByNameLexicographically) & (! unsorted ))) error("FATAL ERROR: SAM/BAM file is not sorted by name (based on the first "+peekCt+" reads)! Either sort the file by name, or sort by read position and use the \"--coordSorted\" option.");
    }
   
    
    //reportln("SAMRecord Reader Generated. Read length: "+readLength+".","note");

    val coda : Array[Int] = internalUtils.commonSeqUtils.getNewCauseOfDropArray;
    val coda_options : Array[Boolean] = internalUtils.commonSeqUtils.CODA_DEFAULT_OPTIONS.toArray;
    if(isSingleEnd) CODA_SINGLE_END_OFF_OPTIONS.foreach( coda_options(_) = false );
    if(keepMultiMapped) coda_options(internalUtils.commonSeqUtils.CODA_NOT_UNIQUE_ALIGNMENT) = false;
    if(! readGroup.isEmpty) coda_options(internalUtils.commonSeqUtils.CODA_NOT_MARKED_RG) = true;
    
    var readNum = 0;
    val samIterationTimeStamp = TimeStampUtil();
    
    
    for(pair <- pairedIter){
    //for((pair,readNum) <- numberedIter){
      val (r1,r2) = pair;
      readNum += 1;
      if(internalUtils.commonSeqUtils.useReadPair(r1,r2,coda, coda_options, Set[String](), readGroup, minMAPQ = minMAPQ)){
        qcBTW.runOnReadPair(r1,r2,readNum);
      }
    }
    
    val postRunStamp = TimeStampUtil();
    
    val iterationMinutes = (postRunStamp.compareTo(samIterationTimeStamp) / 1000).toDouble / 60.toDouble;
    val minutesPerMillion = iterationMinutes / (readNum.toDouble / 1000000.toDouble);
    val minutesPerMillionPF = iterationMinutes / ((coda(internalUtils.commonSeqUtils.CODA_READ_PAIR_OK)).toDouble / 1000000.toDouble);

    reportln("> Time spent on SAM iteration:   " + TimeStampUtil.timeDifferenceFormatter(postRunStamp.compareTo(samIterationTimeStamp)) + " ","note");
    reportln(">                                (" + minutesPerMillion + " minutes per million read-pairs)"+ " ","note");
    reportln(">                                (" + minutesPerMillionPF + " minutes per million read-pairs used)"+ " ","note");
  } 
  
  
  def runOnFile(infile : String, qcBTW : QcBamToWig, testRun : Boolean, isSingleEnd : Boolean, keepMultiMapped : Boolean, readGroup : Option[String]){
    val (samFileAttributes, recordIter) = initSamRecordIterator(infile);
    val pairedIter : Iterator[(SAMRecord,SAMRecord)] = 
      if(isSingleEnd){
        if(testRun) samRecordPairIterator_withMulti_singleEnd(recordIter, true, 200000) else samRecordPairIterator_withMulti_singleEnd(recordIter);
      } else {
        if(testRun) samRecordPairIterator_withMulti(recordIter, true, 200000) else samRecordPairIterator_withMulti(recordIter)
      }
    
    val coda : Array[Int] = internalUtils.commonSeqUtils.getNewCauseOfDropArray;
    val coda_options : Array[Boolean] = internalUtils.commonSeqUtils.CODA_DEFAULT_OPTIONS.toArray;
    if(isSingleEnd) CODA_SINGLE_END_OFF_OPTIONS.foreach( coda_options(_) = false );
    if(keepMultiMapped) coda_options(internalUtils.commonSeqUtils.CODA_NOT_UNIQUE_ALIGNMENT) = false;
    if(! readGroup.isEmpty) coda_options(internalUtils.commonSeqUtils.CODA_NOT_MARKED_RG) = true;
    
    var readNum = 0;
    val samIterationTimeStamp = TimeStampUtil();
    for(pair <- pairedIter){
    //for((pair,readNum) <- numberedIter){
      val (r1,r2) = pair;
      readNum += 1;
      if(internalUtils.commonSeqUtils.useReadPair(r1,r2,coda, coda_options, Set[String](), readGroup)){
        qcBTW.runOnReadPair(r1,r2,readNum);
      }
    }
  }
  
  def getReadBlocks(r : SAMRecord) : Vector[(Int,Int)] = {
    r.getAlignmentBlocks().toVector.map((block) => {
      (block.getReferenceStart() - 1, block.getReferenceStart() - 1 + block.getLength());
    });
  }
  
  def getOverlappedPairBlocks(r1 : SAMRecord, r2 : SAMRecord) : Vector[(Int,Int)] = {
    val r1b = getReadBlocks(r1);
    val r2b = getReadBlocks(r2);
    
    //def blocksOverlap(b1 : (Int,Int), b2 : (Int,Int)) : Boolean = {
    //  b1._1 <= b2._2 && b2._1 <= b1._2;
    //}
    //def mergeBlocks(b1 : (Int,Int), b2 : (Int,Int)) : (Int,Int) = {
    //  (math.min(b1._1,b2._1), math.max(b1._2, b2._2))
    //}
    //val r2bOverlap : Seq[(Int,Int)] = r2b.filter((b2) => r1b.exists(blocksOverlap(_,b2)));
    //val r2bNonOverlap : Seq[(Int,Int)] = r2b.filterNot((b2) => r1b.exists(blocksOverlap(_,b2)));
    val merged = (r1b ++ r2b).sorted
    merged.tail.foldLeft(Vector(merged.head))( (soFar,curr) =>{
      if(curr._1 <= soFar.last._2){
        soFar.updated(soFar.length - 1, (soFar.last._1, math.max(curr._2, soFar.last._2)));
      } else {
        soFar :+ curr;
      }
    })
  }
  
  
  class QcBamToWig(trackName : String, chromLengthFile : String, 
                   noTruncate : Boolean, windowSize : Int, 
                   isSingleEnd: Boolean, stranded : Boolean,
                   fr_secondStrand : Boolean, sizeFactor : Double, 
                   negativeReverseStrand : Boolean, countPairsTogether : Boolean, 
                   includeTrackDef : Boolean, rgbColor : Option[String],
                   additionalTrackOptions : String) extends QCUtility[Unit] {
    
    val chromMap : Map[(String,Char),Chrom] = genChrom(chromLengthFile, windowSize, stranded, ! noTruncate);
    var unknownChromSet : Set[String] = Set[String]();
    
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
      if(isSingleEnd){
        runOnRead(r1);
      } else if(! countPairsTogether){
        runOnRead(r1);
        runOnRead(r2);
      } else {
        val chromName = r1.getReferenceName();
        val strand = internalUtils.commonSeqUtils.getStrand(r1 , stranded , fr_secondStrand );
        val blocks = getOverlappedPairBlocks(r1,r2);
        chromMap.get((chromName,strand)) match {
          case Some(chrom) => {
            chrom.countBlocks(blocks);
          }
          case None => {
            if(! unknownChromSet.contains(chromName)){
              reportln("WARNING: Chromosome not found in the chromLengthFile: ["+chromName+","+strand+"]","deepDebug");
              unknownChromSet += chromName;
            }
          }
        }
      }
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
        val rgbColorString = if(rgbColor.isEmpty) " " else " color="+rgbColor.get+" altColor="+rgbColor.get +" itemRgb=On";
        val additionalOptionsString = " "+additionalTrackOptions+" ";
        if(includeTrackDef){
          writerF.write("track name="+trackName+"_FWD type=wiggle_0 visibility=full"+rgbColorString+additionalOptionsString+"\n");
          writerR.write("track name="+trackName+"_REV type=wiggle_0 visibility=full"+rgbColorString+additionalOptionsString+"\n");
        }
        
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
        val rgbColorString = if(rgbColor.isEmpty) " " else " color="+rgbColor.get+" altColor="+rgbColor.get +" ";
        val additionalOptionsString = " "+additionalTrackOptions+" ";
        if(includeTrackDef){
          writer.write("track name="+trackName+" type=wiggle_0 "+rgbColorString+" "+additionalOptionsString+"\n");
        }
        
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
  
  case class Chrom(chromName :  String, chromStrand : Char, windowCounts : Array[Long], span : Int, truncate : Boolean){ 
    def countSamRecord(samRecord : SAMRecord) {
        //val blocks = samRecord.getAlignmentBlocks();
        val blocks = getReadBlocks(samRecord);
        
        countBlocks(blocks);
        /*var i = 0; 
        for((start,end) <- blocks){
        //while(i < blocks.size()){
          //val block : AlignmentBlock = blocks.get(i);
          //val start = block.getReferenceStart();
          //val length = block.getLength();
          //val end = start + length;
          
          if(! truncate && (end - 1 / span) >=  windowCounts.length){
            error("ERROR: read extends outside chromosome length!");
          }
          
          var j = start;
          while(j < end && (j / span) < windowCounts.length){
            windowCounts(j / span) = windowCounts(j / span) + 1;
            j = j + 1;
          }
          
          i = i + 1;
        }*/
    }
    
    def countBlocks(blocks : Seq[(Int,Int)]){
      for((start,end) <- blocks){
        if(! truncate && (end - 1 / span) >=  windowCounts.length){
          error("ERROR: read extends outside chromosome length!");
        }
        var j = start;
        while(j < end && (j / span) < windowCounts.length){
          windowCounts(j / span) = windowCounts(j / span) + 1;
          j = j + 1;
        }
      }
    }
  }
  

  def writeChrom(chrom : Chrom, writer : WriterUtil, sizeFactor : Double, negativeStrand : Boolean){
    writer.write("fixedStep  chrom="+chrom.chromName+"  start=1  step="+chrom.span+" span="+chrom.span+"\n");
    val adjustmentFactor = if(negativeStrand) ((-1).toDouble * sizeFactor * chrom.span.toDouble) else (sizeFactor * chrom.span.toDouble);        
    chrom.windowCounts.foreach(
      (ct : Long) => {
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
        val chromP = new Chrom(chromName,'+',Array.ofDim[Long](spanCt),span, truncate);
        val chromM = new Chrom(chromName,'-',Array.ofDim[Long](spanCt),span, truncate);
        chromMap = chromMap + (((chromName,'+'), chromP));
        chromMap = chromMap + (((chromName,'-'), chromM));        
      } else {
        val chrom = new Chrom(chromName,'.',Array.ofDim[Long](spanCt),span, truncate);
        chromMap = chromMap + (((chromName,'.'), chrom));
      }
    }
    chromMap;
  }
}