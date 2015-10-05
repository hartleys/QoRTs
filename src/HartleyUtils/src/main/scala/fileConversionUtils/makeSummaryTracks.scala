package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;

import java.io.File;

object makeSummaryTracks {

  /*
   * INCOMPLETE! UNTESTED!
   */
  
  val DEFAULT_MERGE_FILE_LIST : List[String] = List("UnstrandedWiggleTrack","StrandedWiggleTrack","KnownOnlyJunction","WithNovelJunction");
  
  class summaryTrackMerge extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "summaryTrackMerge", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "",   
          argList = 
                    new BinaryOptionArgument[String](
                                         name = "sampleID", 
                                         arg = List("--sampleID"), 
                                         valueName = "sampid",  
                                         argDesc = "Optional: the id of the specific sample that is to be merged. By default, this utility will merge all samples found in the decoder. With this option selected, it will ONLY merge the one sample named here."
                                        ) ::
                    new BinaryArgument[Int](
                                         name = "wiggleWindow", 
                                         arg = List("--wiggleWindow"), 
                                         valueName = "val",  
                                         argDesc = "The window size of the alternate-size wiggle track, if applicable." ,
                                         defaultValue = Some(100)
                                        ) ::
                    new BinaryArgument[List[String]](
                                         name = "mergeFiles", 
                                         arg = List("--mergeFiles"), 
                                         valueName = "file1[,file2,...]",  
                                         argDesc = "A comma-delimited list of strings, indicating which file types to attempt to merge. By default, this utility autodetects the presence of all mergable qc files and merges all standard files. Valid codes are:" + DEFAULT_MERGE_FILE_LIST.mkString(", ") ,
                                         defaultValue = Some(DEFAULT_MERGE_FILE_LIST)
                                        ) ::
                    new FinalArgument[String](
                                         name = "infileDir",
                                         valueName = "infileDir",
                                         argDesc = "The top-level directory in which all the QC output can be found." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "decoderFile",
                                         valueName = "decoderFile",
                                         argDesc = "The decoder file, which must conform to the requirements of the QoRT decoder specification. "+
                                                   "In particular it MUST have two specific columns: \n"+
                                                   "\"sample.ID\" This utility will merge the count data output from all bamfiles that have the same sample.ID"+
                                                   "\nand\n"+
                                                   "\"group.ID\" This must be the name of the experimental condition group." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output file directory." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS
     );
     
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         mergeQcOutput.run(
             parser.get[String]("infileDirs"),
             parser.get[String]("outfile"),
             parser.get[List[String]]("mergeFiles"),
             parser.get[Int]("wiggleWindow"),
             "",
             "Untitled"
           );
         }
     }
  }
  
  
  def runAll(infileDir : String, outfileDir : String, decoderFile : String, sizeFactorFile : Option[String], 
             mergeGroups : Option[List[String]], ignoreSizeFactors : Boolean, 
             groupByColName : String = "group.ID", mergeList : List[String]){
    val decoderRaw = getLinesSmartUnzip(decoderFile, true).toVector;
    val decoderHeadLine = decoderRaw.head;
    val decoderRawTail = decoderRaw.tail;
    val decoderVals = decoderRawTail.map(_.split("\\s+").toVector);
    
    val decoderColNames = decoderHeadLine.split("\\s+");

    val sampleCol = if(! decoderColNames.contains("sample.ID")){
      if(decoderColNames.contains("unique.ID")){
        decoderColNames.indexOf("unique.ID")
      } else if(decoderColNames.contains("lanebam.ID")){
        decoderColNames.indexOf("lanebam.ID")
      } else {
        error("FATAL ERROR: decoder must contain sample.ID or unique.ID column.");
        -1;
      }
    } else {
      decoderColNames.indexOf("sample.ID")
    }

    val groupCol = if(decoderColNames.contains(groupByColName)){
      decoderColNames.indexOf(groupByColName);
    } else {
      error("FATAL ERROR: Decoder must contain a column named: \""+groupByColName+"\"");
      -1;
    }
    
    val sizeFactorCol = if(decoderColNames.contains("size.factor")){
      decoderColNames.indexOf("size.factor");
    } else {
      -1;
    }
    
    if(! decoderVals.forall( _.length == decoderColNames.length )){
      val minLength = decoderVals.map(_.length).min;
      val minLengthLineNum = decoderVals.map(_.length).indexOf(minLength);
      if(minLength > sampleCol && minLength > groupCol && minLength > sizeFactorCol) {
        reportln("WARNING: not all decoder lines of same length! Title line is of length " + decoderColNames.length + " whereas line " + minLengthLineNum + " is of length " + minLength,"warn");
      } else {
        error("ERROR: not all decoder lines of same length! Title line is of length " + decoderColNames.length + " whereas line " + minLengthLineNum + " is of length " + minLength + ", which means that either sample.ID or group.ID is absent from this line.");
      }
    }
    
    val sampleIDs = decoderVals.map( _.apply(sampleCol) );
    val groupIDs = decoderVals.map( _.apply(groupCol) );
    
    val sizeFactors : Vector[Double] = if(ignoreSizeFactors){
      repToSeq(1.0,sampleIDs.length).toVector;
    } else {
      if(sizeFactorFile.isEmpty){
        if(sizeFactorCol == -1){
          repToSeq(1.0,sampleIDs.length).toVector;
        } else {
          decoderVals.map( (v) => string2double(v(sizeFactorCol)) );
        }
      } else {
        val sizeFactorMap = readSizeFactorMap(sizeFactorFile.get);
        val sfvector = sampleIDs.map((id) => {
          sizeFactorMap.get(id) match {
            case Some(sf) => sf;
            case None => {
              error("ERROR: sample.ID "+id+" not found in size factor file. Size factor file only contains sample.IDs: [" + sizeFactorMap.keys.mkString(", ") + "]");
              -1.0;
            }
          }
        })
        sfvector;
      }
    }
    
    val decoder : Vector[((String, Double),String)] = sampleIDs.zip(sizeFactors).zip(groupIDs);
    
    val groups : Vector[String] = if(mergeGroups.isEmpty){
      groupIDs.toSet.toVector;
    } else {
      groupIDs.toSet.toVector.filter(mergeGroups.get.contains(_));
    }
    
    for(grp <- groups){
      reportln("> Starting summary merger of "+groupByColName+"=\""+grp+"\"","progress");
      
      val subsetDecoder = decoder.filter(t => t._2 == grp);
      
      //print testing details:
      reportln("> Merging Samples:","progress");
      for(i <- 0 until subsetDecoder.length){
        reportln(">   sample.ID=\""+subsetDecoder(i)._1._1+"\", group.ID=\""+subsetDecoder(i)._2+"\", size.factor="+subsetDecoder(i)._1._2,"note");
      }
      
      val subsetPairs : Vector[(String,Double)] = subsetDecoder.map(_._1);
      
      
    }
    
  }
  
  
  
  def makeSummaryTracksFromSet(infilePathPrefix : String, infilePrefix : String = "/",
                              decoder : Seq[(String,Double)], mergeList : List[String] = DEFAULT_MERGE_FILE_LIST,
                              outfilePrefix : String,
                              
                              
                              junctionBedWithNovelSuffix : String = "QC.junctionBed.withNovel.bed.gz",
                              junctionBedKnownSuffix : String = "QC.junctionBed.known.bed.gz",
                              unstrandedWiggleSuffix : String = "QC.wiggle.unstranded.wig.gz",
                              fwdWiggleSuffix : String = "QC.wiggle.fwd.wig.gz",
                              revWiggleSuffix : String = "QC.wiggle.rev.wig.gz"
                              ){
    
    //List("UnstrandedWiggleTrack","StrandedWiggleTrack","KnownJunction","NovelJunction");
    if(mergeList.contains("UnstrandedWiggleTrack")){
      val firstFile = infilePathPrefix + decoder.head._1 + infilePrefix + unstrandedWiggleSuffix;
      if((new File(firstFile)).exists()){
        //runHelper2(pairlist : Seq[(String,Double)], outfile : String, quiet : Boolean, trackDefLine : Option[String])
        val fileDecoder = decoder.map(p => (infilePathPrefix + p._1 + infilePrefix + unstrandedWiggleSuffix, p._2) );
        fileConversionUtils.SumWigglesFast.runHelper2(fileDecoder, outfilePrefix, None);
      } else {
        reportln("WARNING: file "+firstFile+" not found. Skipping UnstrandedWiggleTrack summary.","warn");
      }
    }
    
    if(mergeList.contains("StrandedWiggleTrack")){
      val firstFile = infilePathPrefix + decoder.head._1 + infilePrefix + fwdWiggleSuffix;
      if((new File(firstFile)).exists()){
        val fileDecoderFwd = decoder.map(p => (infilePathPrefix + p._1 + infilePrefix + fwdWiggleSuffix, p._2) );
        fileConversionUtils.SumWigglesFast.runHelper2(fileDecoderFwd, outfilePrefix, None);
        
        val fileDecoderRev = decoder.map(p => (infilePathPrefix + p._1 + infilePrefix + revWiggleSuffix, p._2) );
        fileConversionUtils.SumWigglesFast.runHelper2(fileDecoderRev, outfilePrefix, None);
      } else {
        reportln("WARNING: file "+firstFile+" not found. Skipping StrandedWiggleTrack summary.","warn");
      }
    }
    
    if(mergeList.contains("KnownOnlyJunction")){
      val firstFile = infilePathPrefix + decoder.head._1 + infilePrefix + junctionBedKnownSuffix;
      if((new File(firstFile)).exists()){
        
      } else {
        reportln("WARNING: file "+firstFile+" not found. Skipping KnownOnlyJunction summary.","warn");
      }
    }
    
    if(mergeList.contains("WithNovelJunction")){
      val firstFile = infilePathPrefix + decoder.head._1 + infilePrefix + junctionBedWithNovelSuffix;
      if((new File(firstFile)).exists()){
        
      } else {
        reportln("WARNING: file "+firstFile+" not found. Skipping WithNovelJunction summary.","warn");
      }
    }
    
    
  }
  
  private def readSizeFactorMap(sfFile : String) : Map[String,Double] = {
    val rawsf = getLinesSmartUnzip(sfFile).toVector;
    
    (if(rawsf.head.substring(0,9) == "sample.ID") {
        if(rawsf.head.substring(0,21) != "sample.ID	size.factor"){
          error("Error: first two columns of table must be sample.ID and size.factor. The table must be tab-delimited.");
        }
        rawsf.tail.map(line => { 
          val cells = line.split("	");
          (cells(0),string2double(cells(1)));
        }).toVector;
      } else {
        rawsf.map(line => { 
          val cells = line.split("	");
          (cells(0),string2double(cells(1)));
        }).toVector;
    }).toMap
    
    /*val out = sf.map((s) => {
      val cells = s.split("\\s+");
      if(cells.length < 2) error("ERROR: size factor line has fewer than 2 columns.");
      (cells(0), string2double(cells(1)))
    }).toMap;
    
    out;*/
  }
  
  
}











