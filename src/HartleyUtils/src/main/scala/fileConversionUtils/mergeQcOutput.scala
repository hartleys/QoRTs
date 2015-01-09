package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;

import java.io.File;

object mergeQcOutput {
  val mergeFileList : List[String] = List("DESeq","DEXSeq","JunctionSeq","NovelSplice","KnownSplice","WiggleTrack","WiggleTrackAltWin");
  
  class merger extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "mergeCounts", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This utility merges count data from multiple QoRTs QC runs."+
                        ""+
                        ""+
                        ""+
                        ""+
                        ""+
                        "",   
          argList = 
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
                                         argDesc = "A comma-delimited list of strings, indicating which file types to attempt to merge. By default, this utility autodetects the presence of all mergable qc files and merges all standard files. Valid codes are:" + mergeFileList.mkString(",") ,
                                         defaultValue = Some(mergeFileList)
                                        ) ::
                    new FinalArgument[String](
                                         name = "infileDirs",
                                         valueName = "infileDirs",
                                         argDesc = "The input files' directories, as a comma-delimited list with no whitespace." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfilePrefix",
                                         valueName = "outfilePrefix",
                                         argDesc = "The output file prefix (or directory)" // description
                                        ) :: List() );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         mergeQcOutput.run(
             parser.get[String]("infileDirs"),
             parser.get[String]("outfile"),
             parser.get[List[String]]("mergeFiles"),
             parser.get[Int]("wiggleWindow")
           );
         }
     }
   }
  
  class multiMerger extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "mergeAllCounts", 
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
                                         argDesc = "A comma-delimited list of strings, indicating which file types to attempt to merge. By default, this utility autodetects the presence of all mergable qc files and merges all standard files. Valid codes are:" + mergeFileList.mkString(",") ,
                                         defaultValue = Some(mergeFileList)
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
                                                   "\"qc.data.dir\" This must be the file path to the output data directory, from the infileDir file location." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output file directory." // description
                                        ) :: List() );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         mergeQcOutput.multirun(
             parser.get[String]("infileDir"),
             parser.get[String]("decoderFile"),
             parser.get[String]("outfile"),
             parser.get[List[String]]("mergeFiles"),
             parser.get[Int]("wiggleWindow"),
             parser.get[Option[String]]("sampleID")
           );
         }
     }
  }
  
  def multirun(infileDir : String, decoderFile : String, outfile : String, mergeFiles : List[String], wiggleWindow : Int, sampID : Option[String]){
    val decoder : Map[String,Set[String]] = decode(decoderFile,infileDir);
    
    val initialTimeStamp = TimeStampUtil();
    
    if(sampID.isEmpty){
      for((sampleID, infiles) <- decoder){
        merge(infiles.toSeq, outfile + "/"+sampleID+"/",  mergeFiles, wiggleWindow);
      }
    } else {
      val infiles = decoder(sampID.get);
      merge(infiles.toSeq, outfile + "/"+sampID.get+"/",  mergeFiles, wiggleWindow);
    }
    
    standardStatusReport(initialTimeStamp);
  }
  def decode(decoder : String, infileDir : String) : Map[String,Set[String]] = {
    val lines = getLinesSmartUnzip(decoder);
    val header = lines.next.split("\\s+");
    val qcDirCol = header.indexOf("qc.data.dir");
    val idCol = header.indexOf("sample.ID");
    
    if(idCol == -1) error("Fatal error! Decoder has no column named sample.ID!");
    if(qcDirCol == -1) error("Fatal error! Decoder has no column named qc.data.dir!");
    
    lines.foldLeft(Map[String,Set[String]]().withDefault(x => Set[String]()))((soFar, line) =>{
      val cells = line.split("\\s+");
      val qcDir = cells(qcDirCol) + "/";
      val id = cells(idCol);
      soFar.updated(id,soFar(id) + (infileDir + qcDir));
    });
  }
   
   def run(infile : String, outfile : String, mergeFiles : List[String], wiggleWindow : Int){
     val infiles = infile.split(",");
     merge(infiles.toSeq, outfile, mergeFiles, wiggleWindow);
   }
   
   def merge(infiles : Seq[String], outfile : String, mergeFiles : List[String], wiggleWindow : Int){
     
     var currTimeStamp = TimeStampUtil();
     
     val DESeq_suffix = "QC.geneCounts.formatted.for.DESeq.txt.gz";
     val DEXSeq_suffix = "QC.exonCounts.formatted.for.DEXSeq.txt.gz";
     val Splice_suffix = "QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz";
     val NovelSplice_suffix = "QC.spliceJunctionCounts.novelSplices.txt.gz";
     val KnownSplice_suffix = "QC.spliceJunctionCounts.knownSplices.txt.gz";
     
     val unstranded_wiggle_suffix = "QC.wiggle.unstranded.wig.gz";
     val stranded_wiggle_fwd_suffix = "QC.wiggle.fwd.wig.gz";
     val stranded_wiggle_rev_suffix = "QC.wiggle.rev.wig.gz";
     
     val unstranded_alt_wiggle_suffix = "QC.wiggle.Win"+wiggleWindow.toString()+".unstranded.wig.gz";
     val stranded_alt_wiggle_fwd_suffix = "QC.wiggle.Win"+wiggleWindow.toString()+".fwd.wig.gz";
     val stranded_alt_wiggle_rev_suffix = "QC.wiggle.Win"+wiggleWindow.toString()+".rev.wig.gz";
     
     
     //if(outfile.last == '/' || outfile.last == '\\'){
       val outDir = new File(outfile);
       if(! outDir.exists()){
         reportln("Creating Directory: "+ outDir,"note");
         outDir.mkdir();
         reportln("Successfully Created Directory: " + outDir, "note");
       }
     //}
     
     if(! mergeFiles.contains("DESeq")){
       reportln("Skipping DESeq merge" ,"note");
     } else if((new File(infiles.head +DESeq_suffix)).exists()){
       report("Merging DESeq data...","note");
       mergeSimpleData(infiles, DESeq_suffix, outfile + DESeq_suffix);
       report("done\n","note");
     } else {
       reportln("DESeq data not found at: " + infiles.head +DESeq_suffix + "\n	Skipping DESeq merge" ,"note");
     }
     
     standardStatusReport(currTimeStamp, ">    ");
     currTimeStamp = TimeStampUtil();
     
     if(! mergeFiles.contains("DEXSeq")){
       reportln("Skipping DEXSeq merge" ,"note");
     } else if((new File(infiles.head +DEXSeq_suffix)).exists()){
       report("Merging DEXSeq data...","note");
       mergeSimpleData(infiles, DEXSeq_suffix, outfile + DEXSeq_suffix);
       report("done\n","note");
     } else {
       reportln("DEXSeq data not found at: " + infiles.head +DEXSeq_suffix + "\n	Skipping DEXSeq merge" ,"note");
     }
     
     standardStatusReport(currTimeStamp, ">    ");
     currTimeStamp = TimeStampUtil();
     
     if(! mergeFiles.contains("JunctionSeq")){
       reportln("Skipping JunctionSeq merge" ,"note");
     } else if((new File(infiles.head +Splice_suffix)).exists()){
       report("Merging JunctionSeq data...","note");
       mergeSimpleData(infiles, Splice_suffix, outfile + Splice_suffix);
       report("done\n","note");
     } else {
       reportln("Splice data not found at: " + infiles.head +Splice_suffix + "\n	Skipping SpliceSeq merge" ,"note");
     }
     
     standardStatusReport(currTimeStamp, ">    ");
     currTimeStamp = TimeStampUtil();
     
     if(! mergeFiles.contains("NovelSplice")){
       reportln("Skipping NovelSplice merge" ,"note");
     } else if((new File(infiles.head +NovelSplice_suffix)).exists()){
       report("Merging novel splice data...","note");
       mergeNovelSpliceData(infiles, NovelSplice_suffix, outfile + NovelSplice_suffix);
       report("done\n","note");
     } else {
       reportln("Splice data not found at: " + infiles.head +NovelSplice_suffix + "\n	Skipping Novel Splice merge" ,"note");
     }
     
     standardStatusReport(currTimeStamp, ">    ");
     currTimeStamp = TimeStampUtil();
     
     if(! mergeFiles.contains("KnownSplice")){
       reportln("Skipping KnownSplice merge" ,"note");
     } else if((new File(infiles.head +KnownSplice_suffix)).exists()){
       report("Merging known splice data...","note");
       mergeComplexData(infiles, KnownSplice_suffix, outfile + KnownSplice_suffix);
       report("done\n","note");
     } else {
       reportln("Splice data not found at: " + infiles.head +KnownSplice_suffix + "\n	Skipping Known Splice merge" ,"note");
     }
     
     standardStatusReport(currTimeStamp, ">    ");
     currTimeStamp = TimeStampUtil();
     
     if(! mergeFiles.contains("WiggleTrack")){
       reportln("Skipping WiggleTrack merge" ,"note");
     } else {
       //TODO!!!
       //reportln("NOTE: makeAllBrowserTracks is UNIMPLEMENTED at this time!","debug");
       
       if( ((new File(infiles.head + unstranded_wiggle_suffix)).exists()) | ((new File(infiles.head + stranded_wiggle_fwd_suffix)).exists())){
         //Do nothing?
       } else {
         reportln("WARNING: Cannot find wiggle files! Attempted to find:","note");
         reportln("          \""+infiles.head + unstranded_wiggle_suffix+"\"","note");
         reportln("          OR","note");
         reportln("          \""+infiles.head + stranded_wiggle_fwd_suffix+"\"","note");
       }
       
       if((new File(infiles.head + unstranded_wiggle_suffix)).exists()){
         report("Merging unstranded wiggle data...","note");
         val pairlist = infiles.map(infile => (infile + unstranded_wiggle_suffix, 1.0));
         SumWigglesFast.runHelper2(pairlist, outfile + unstranded_wiggle_suffix, false);
         report("done\n","note");
       }
       if((new File(infiles.head + stranded_wiggle_fwd_suffix)).exists()){
         report("Merging stranded wiggle data...","note");
         val fwdpairlist = infiles.map(infile => (infile + stranded_wiggle_fwd_suffix, 1.0));
         SumWigglesFast.runHelper2(fwdpairlist, outfile + stranded_wiggle_fwd_suffix, false);
         
         val revpairlist = infiles.map(infile => (infile + stranded_wiggle_rev_suffix, 1.0));
         SumWigglesFast.runHelper2(revpairlist, outfile + stranded_wiggle_rev_suffix, false);
         report("done\n","note");
       }
        
       //SumWigglesFast.run(filelist , outfile , false , false , false , true , None)
     }
     
     standardStatusReport(currTimeStamp, ">    ");
     currTimeStamp = TimeStampUtil();
     
     if(! mergeFiles.contains("WiggleTrackAltWin")){
       reportln("Skipping WiggleTrack (with non-default window size) merge" ,"note");
     } else {
       val unstrandedWigExists = (new File(infiles.head + unstranded_alt_wiggle_suffix)).exists();
       val strandedWigExists = (new File(infiles.head + stranded_alt_wiggle_fwd_suffix)).exists();
       
       if(! (unstrandedWigExists || strandedWigExists) ){
         report("No wiggle track with alternate window size found. This is probably normal.","note");
       }
       
       if(unstrandedWigExists){
         report("Merging unstranded alt-window wiggle data...","note");
         val pairlist = infiles.map(infile => (infile + unstranded_alt_wiggle_suffix, 1.0));
         SumWigglesFast.runHelper2(pairlist, outfile + unstranded_alt_wiggle_suffix, false);
         report("done\n","note");
       }
       if(strandedWigExists){
         report("Merging stranded alt-window wiggle data...","note");
         val fwdpairlist = infiles.map(infile => (infile + stranded_alt_wiggle_fwd_suffix, 1.0));
         SumWigglesFast.runHelper2(fwdpairlist, outfile + stranded_alt_wiggle_fwd_suffix, false);
         
         val revpairlist = infiles.map(infile => (infile + stranded_alt_wiggle_rev_suffix, 1.0));
         SumWigglesFast.runHelper2(revpairlist, outfile + stranded_alt_wiggle_rev_suffix, false);
         report("done\n","note");
       }
       //SumWigglesFast.run(filelist , outfile , false , false , false , true , None)
     }
     
     standardStatusReport(currTimeStamp, ">    ");
     currTimeStamp = TimeStampUtil();
     
   }
   
   def mergeNovelSpliceData(infilePrefixes : Seq[String], infileSuffix : String, outfile : String){
     //val iterSeq : Seq[Iterator[String]] = infilePrefixes.map(prefix => getLinesSmartUnzip(prefix + infileSuffix));
     
     val outmap : Map[GenomicInterval,Int] = infilePrefixes.foldLeft(Map[GenomicInterval,Int]().withDefault(x => 0))((soFar,currFilePrefix) =>{
       val lines = getLinesSmartUnzip(currFilePrefix + infileSuffix);
       lines.next;
       lines.foldLeft(soFar)((soFar2, line) => {
         val cells = line.split("	");
         val iv = GenomicInterval(cells(0), cells(1).charAt(0), string2int(cells(2)), string2int(cells(3)));
         val ct = string2int(cells(4)) + soFar2(iv);
         soFar2.updated(iv,ct);
       })
     });
     
     val writer = openWriterSmart(outfile);
     writer.write("chrom	strand	start	end	CT\n");
     for(iv <- outmap.keys.toVector.sorted){
       val ct = outmap(iv);
       writer.write(iv.chromName+"	"+iv.strand+"	"+iv.start+"	"+iv.end+"	"+ct+"\n");
     }
     close(writer);
   }
   
   def mergeComplexData(infilePrefixes : Seq[String], infileSuffix : String, outfile : String){
     
     val iterSeq : Seq[Iterator[String]] = infilePrefixes.map(prefix => getLinesSmartUnzip(prefix + infileSuffix));
     iterSeq.foreach(_.next);
     
     val sumIterator : Iterator[String] = new Iterator[String](){
       def hasNext : Boolean = iterSeq.head.hasNext;
       def next : String = {
         if(iterSeq.forall( _.hasNext )){
           val lines = iterSeq.map(_.next);
           val pairs = lines.map(getPairComplex(_));
           val name = pairs.head._1;
           if(! pairs.forall(_._1 == name)) error("FATAL ERROR: files do not have identical feature names and ordering!");
           val sum = pairs.foldLeft(0)((soFar, curr) => soFar + curr._2);
           name +"	"+sum;
         } else {
           error("FATAL ERROR: not all files of the same length!");
           null;
         }
       }
     }
     
     val writer = openWriterSmart(outfile);
     for(line <- sumIterator){
       writer.write(line +"\n");
     }
     close(writer);
   }
   
   def getPairComplex(str : String) : (String,Int) = {
     val index = str.lastIndexOf('	');
     val split = str.splitAt(index);
     return (split._1, string2int(split._2.tail));
   }
   
   
   def mergeSimpleData(infilePrefixes : Seq[String], infileSuffix : String, outfile : String){
     
     val iterSeq : Seq[Iterator[String]] = infilePrefixes.map(prefix => getLinesSmartUnzip(prefix + infileSuffix));
     val sumIterator : Iterator[String] = new Iterator[String](){
       def hasNext : Boolean = iterSeq.head.hasNext;
       def next : String = {
         if(iterSeq.forall( _.hasNext )){
           val lines = iterSeq.map(_.next);
           val pairs = lines.map(getPair(_));
           val name = pairs.head._1;
           if(! pairs.forall(_._1 == name)) error("FATAL ERROR: files do not have identical feature names and ordering!");
           val sum = pairs.foldLeft(0)((soFar, curr) => soFar + curr._2);
           name +"	"+sum;
         } else {
           error("FATAL ERROR: not all files of the same length!");
           null;
         }
       }
     }
     
     val writer = openWriterSmart(outfile);
     for(line <- sumIterator){
       writer.write(line +"\n");
     }
     close(writer);
   }
   
   private def getPair(in : String) : (String, Int) = {
     val cells = in.split("\\s+");
     if(cells.length != 2){
       error("FATAL ERROR: file misformatted!");
     }
     if(cells(1) == "NA"){
       return (cells(0),-1);
     } else {
       return ((cells(0), string2int(cells(1)) ));
     }
   }
}

















