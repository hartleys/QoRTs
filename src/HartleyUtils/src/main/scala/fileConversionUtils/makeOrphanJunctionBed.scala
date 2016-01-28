package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;
import internalUtils.GtfTool._;

object makeOrphanJunctionBed {

  class converter extends CommandLineRunUtil { 
     override def priority = 28;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "makeOrphanJunctionTrack", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This utility takes the 'orphan' splice junction count files created by the QoRTs QC utility "+
                        "(optionally across multiple samples) "+
                        "and creates a single merged splice junction 'bed' file that lists each splice junction along with the "+
                        "read-pair coverage counts. It can optionally calculate the mean counts, and/or normalize the counts "+
                        "using the supplied normalization size factors."+
                        "The output splice junction bed file "+
                        "can be used to visualize splice junction counts using the UCSC genome browser, IGV, "+
                        "or other similar utilities."+
                        "Note: Either the '--filenames' or the '--sampleList' option MUST be set! The sampleList option is "+
                        "generally used with the --infilePrefix and --infileSuffix options to determine the input filenames."+
                        ""+
                        ""+
                        ""+
                        ""+
                        "",
          argList =             
                    new BinaryOptionArgument[List[String]](
                                         name = "filenames", 
                                         arg = List("--filenames"), 
                                         valueName = "file1,file2,file3,...",  
                                         argDesc = "A comma-delimited list of wiggle files to merge. "+
                                                   "Allows input files to be specified manually. Alternatively, filenames can be inferred from --infilePrefix, --infileSuffix, and the --sampleList, if those options are specified. "+
                                                   "Either this parameter OR --sampleList MUST BE SPECIFIED."+
                                                   ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "sampleList", 
                                         arg = List("--sampleList"), 
                                         valueName = "[sampleList.txt | - | samp1,samp2,samp3,...]",  
                                         argDesc = "Either a comma-delimited list of sample id's or a file containing a list of sample id's. "+
                                                   "The file must either contain no title line, or contain a title line that includes a \"sample.ID\" column. "+
                                                   "Either this option OR --filenames MUST BE SPECIFIED. Note if the sample list is a file then it must end with the extension '.txt'"+
                                                   ""
                                        ) ::
                    new BinaryArgument[String](   name = "infilePrefix",
                                                        arg = List("--infilePrefix"),  
                                                        valueName = "infilePrefix", 
                                                        argDesc = "A file prefix for all input junction count files. By default the full file path should be specified by the infile parameter.", 
                                                        defaultValue = Some("")
                                                        ) ::
                    new BinaryArgument[String](   name = "infileSuffix",
                                                        arg = List("--infileSuffix"),  
                                                        valueName = "infileSuffix", 
                                                        argDesc = "A file suffix for all input junction count files. By default the full file path should be specified by the infile parameter.", 
                                                        defaultValue = Some("")
                                                        ) :: 
                    new BinaryOptionArgument[String](
                                         name = "sizeFactorFile", 
                                         arg = List("--sizeFactorFile"), 
                                         valueName = "val",  
                                         argDesc = "A file containing (at least) two columns: a list of sample ID's and their double-precision floating-point size factors. "+
                                                   "The first line must include at least two columns: \"sample.ID\" and \"size.factor\""+
                                                   "If this option is set, all counts will be divided by the given normalization factors. Sample ID's will be matched from the --sampleList parameter."+
                                                   ""
                                        ) ::
                    new BinaryOptionArgument[List[Double]](
                                         name = "sizeFactors", 
                                         arg = List("--sizeFactors"), 
                                         valueName = "val",  
                                         argDesc = "A list of double-precision floating-point values. "+
                                                   "If this or any size factor option is set,"+
                                                   " all counts will be divided by the given normalization factors."+
                                                   " The length must be the same as the number of files to merge. Must have the same length and ordering as the --sampleList or --filenames parameter."
                                        ) ::
                    new UnaryArgument(name = "calcMean",
                              arg = List("--calcMean"), // name of value
                              argDesc = "Flag to indicate that the splice junction counts should be averaged, rather than added up." // description
                              ) ::
                    new UnaryArgument(name = "stranded",
                              arg = List("--stranded"), // name of value
                              argDesc = "Flag to indicate that data is stranded." // description
                              ) :: 
                    new BinaryArgument[List[String]](   name = "junctionTypes",
                                                        arg = List("--junctionTypes"),  
                                                        valueName = "junctionTypes", 
                                                        argDesc = "Whether to include ambiguous junctions, orphaned junctions, or both. Comma-delimited list (no spaces!).", 
                                                        defaultValue = Some(List("ALL","ambig","orphan"))
                                                        ) :: 
                              
                   new BinaryOptionArgument[String](
                                         name = "rgb", 
                                         arg = List("--rgb"), 
                                         valueName = "r,g,b",  
                                         argDesc = "The rgb color for all the bed file lines."
                                        ) ::
                    new BinaryArgument[String](
                                         name = "trackTitle", 
                                         arg = List("--trackTitle"), 
                                         valueName = "title",  
                                         argDesc = "The title of the track. By default this will be the same as the title parameter.",
                                         defaultValue = Some("UntitledBed")
                                        ) ::
                    new BinaryArgument[String](         name = "additionalTrackOptions",
                                                        arg = List("--additionalTrackOptions"),  
                                                        valueName = "options", 
                                                        argDesc = "Additional track definition options, added to the track definition line. See the UCSC documentation for more information.", 
                                                        defaultValue = Some("")
                                                        ) ::
                    //new UnaryArgument(name = "includeFullSpliceNames",
                    //          arg = List("--includeFullSpliceNames"), // name of value
                    //          argDesc = "Flag to indicate that full splice names, including gene ID, should be used." // description
                    //          ) :: 

                    new UnaryArgument(name = "ignoreSizeFactors",
                              arg = List("--ignoreSizeFactors"), // name of value
                              argDesc = "Flag to indicate that this utility should ignore size factors even if they are found in the input listFile." // description
                              ) ::

                    //new UnaryArgument(name = "nonflatgtf",
                    //          arg = List("--nonflatgtf"), // name of value
                    //         argDesc = "Flag to indicate that instead of a \"flattened\" gff file, a standard-format gtf file has been specified. If this flag is raised, it will automatically create the standard flattened gtf information, in memory. It will not be written to disk" // description
                    //          ) :: 
                    //new UnaryArgument(name = "",
                    //          arg = List("--skipAnnotatedJunctions"), // name of value
                    //          argDesc = "If this option is used, annotated splice junctions will not be included in the output file. "+
                    //                    "Note: this only works if there are novel junctions in the input file."+
                    //                    ""+
                    //                    ""
                    //          ) :: 
                    //new UnaryArgument(name = "skipNovelJunctions",
                    //          arg = List("--skipNovelJunctions"), // name of value
                    //          argDesc = "If this option is used, novel splice junctions will not be included in the output file."+
                    //                    ""+
                    //                    ""+
                    //                    ""
                    //          ) ::    

                   new BinaryArgument[Int](   name = "digits",
                                                        arg = List("--digits"),  
                                                        valueName = "num", 
                                                        argDesc = "The number of digits after the decimal to include for counts.", 
                                                        defaultValue = Some(2)
                                                        ) ::

                   new BinaryOptionArgument[String](
                                         name = "title", 
                                         arg = List("--title"), 
                                         valueName = "title",  
                                         argDesc = "A title to be prepended to each splice junction name"
                                        ) ::
                   
                   // new FinalArgument[String](
                   //                      name = "orphanSpliceGff",
                   //                      valueName = "orphanSplices.gff.gz",
                   //                      argDesc = "The orphaned splice gff file, generated by QoRTs mergeNovelSplices. The mergeNovelSplices function is vital, as it determines orphan/ambig status and assigns each junction locus with its own unique ID." // description
                   //                     ) ::
                   new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output bed file, or '-' to write to stdout. If the filename ends with \".gz\" or \".zip\" then the file will be compressed using the appropriate method." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         makeOrphanJunctionBed.run(
             sizeFactorFile = parser.get[Option[String]]("sizeFactorFile"),
             sizeFactors = parser.get[Option[List[Double]]]("sizeFactors"),
             filenames = parser.get[Option[List[String]]]("filenames"),
             sampleList = parser.get[Option[String]]("sampleList"),
             //orphanSpliceGff = parser.get[String]("orphanSpliceGff"),
             title = parser.get[Option[String]]("title"),
             ignoreSizeFactors = parser.get[Boolean]("ignoreSizeFactors"),
             outfile = parser.get[String]("outfile"),
             infilePrefix = parser.get[String]("infilePrefix"),
             infileSuffix = parser.get[String]("infileSuffix"),
             stranded = parser.get[Boolean]("stranded"),
             digits = parser.get[Int]("digits"),
             //parser.get[Boolean]("includeFullSpliceNames"),
             calcMean = parser.get[Boolean]("calcMean"),
             //parser.get[Boolean]("nonflatgtf"),
             rgb = parser.get[Option[String]]("rgb"),
             trackTitle = parser.get[String]("trackTitle"),
             additionalTrackOptions = parser.get[String]("additionalTrackOptions"),
             junctionTypes = parser.get[List[String]]("junctionTypes")
             //parser.get[Boolean]("skipAnnotatedJunctions"),
             //parser.get[Boolean]("skipNovelJunctions")
           );
         }
     }
   }
  
  final val INDEX_CHROM = 0;
  final val INDEX_STRAND = 1;
  final val INDEX_START = 2;
  final val INDEX_END = 3;
  final val INDEX_CT = 4;
  //optional:
  final val INDEX_ID = 5;
  final val INDEX_GENESET = 6;
  
  def run(sizeFactorFile : Option[String], 
          sizeFactors : Option[List[Double]], 
          filenames : Option[List[String]], 
          sampleList : Option[String],
          //orphanSpliceGff : String,
          title : Option[String],
          ignoreSizeFactors : Boolean,
          outfile : String, infilePrefix : String, infileSuffix : String, 
          stranded : Boolean,
          digits : Int,  calcMean : Boolean,
          rgb : Option[String],
          trackTitle : String,
          additionalTrackOptions : String,
          junctionTypes : List[String]){

    val samples = getSampleList(sampleList, filenames);
    val sf = getSizeFactors(samples, ignoreSizeFactors, sizeFactorFile, sizeFactors);
    val infiles = if(filenames.isEmpty){
      samples.map(infilePrefix + _ + infileSuffix);
    } else {
      filenames.get.map(infilePrefix + _ + infileSuffix).toVector;
    }
    
    if(infiles.length != samples.length){
      error("ERROR: filenames length != # samples");
    }
    
    val rgbline = if(rgb.isEmpty) " " else " itemRgb=\"On\" ";
    
    run_helper(title, infiles, outfile , stranded , sf, rgb , digits , 
               ! calcMean ,  
               Some("track name="+trackTitle+" description="+trackTitle+" "+rgbline+additionalTrackOptions), junctionTypes);
  }
  
  def getSampleList(sampleList : Option[String], filenames : Option[List[String]]) : Vector[String] = {
    if(! sampleList.isEmpty){
      if(sampleList.get.endsWith(".txt") | sampleList.get.endsWith(".txt.gz") | sampleList.get.endsWith(".txt.zip") | sampleList.get == "-"){
        val rawlines = getLinesSmartUnzip(sampleList.get, true).toVector;
        val cells = rawlines.map(_.split("\\s+").toVector);
        if(cells.head.contains("sample.ID")){
          val sampleCol = cells.head.indexOf("sample.ID");
          return cells.tail.map(c => c(sampleCol));
        } else {
          return cells.map(c => c(0));
        }
      } else {
        return (sampleList.get.split(",").toVector);
      }
    } else if(! filenames.isEmpty) {
      return (filenames.get.toVector);
    } else {
      error("either --sampleList or --fileNames must be set!");
      return Vector();
    }
  }
  
  def getSizeFactors(samples : Vector[String], ignoreSizeFactors : Boolean, sizeFactorFile : Option[String], sizeFactors : Option[List[Double]]) : Option[Vector[Double]] = {
    if(ignoreSizeFactors){
      return None;
    } else {
      if(sizeFactors.isEmpty){
        if(sizeFactorFile.isEmpty){
          return None;
        } else {
          val rawlines = getLinesSmartUnzip(sizeFactorFile.get).toVector;
          val cells = rawlines.map(_.split("\\s+").toVector);
          val sfmap = (if(cells.head.contains("sample.ID")){
              val sampleCol = cells.head.indexOf("sample.ID");
              val sfCol = cells.head.indexOf("size.factor");
              if(cells.exists(_.length < math.max(sampleCol, sfCol)+1 )){
                error("Error: Size factor file formatting error: less than "+(math.max(sampleCol, sfCol)+1)+" columns found for line " + cells.indexWhere(_.length < 2));
              }
              cells.tail.map(c => (c(sampleCol), string2double(c(sfCol))));
            } else {
              if(cells.exists(_.length < 2)){
                error("Error: Size factor file formatting error: less than 2 columns found for line " + cells.indexWhere(_.length < 2));
              }
              cells.map(c => (c(0), string2double(c(1))));
          }).toMap;
          
          val sf = samples.map(s => {
            sfmap.get(s) match {
              case Some(f) => f;
              case None => {
                error("FATAL ERROR: Sample "+s+" not found in size factor file.");
                -1.0;
              }
            }
          })
          return Some(sf);
        }
      } else {
        if(sizeFactors.get.length != samples.length){
          error("FATAL ERROR: # of samples != # of size factors.");
        }
        return Some(sizeFactors.get.toVector);
      }
    }
  }

   def normSizeFactors(sf : Vector[(String,Double)], calcMean : Boolean) : Vector[(String,Double)]= {
     if(calcMean){
       sf.map(p => (p._1, p._2 * sf.length))
     } else {
       sf;
     }
   }
  
  /*
  def runOld(opt_title : Option[String], input : String, outfile : String, infilePrefix : String, infileSuffix : String, gff : String, 
          stranded : Boolean, sizeFactors : Option[List[Double]], opt_rgb : Option[String], digits : Int, includeFullSpliceNames : Boolean, 
          calcSum : Boolean, ignoreSizeFactors : Boolean, nonflatgtf : Boolean) {
    
    val initialpairlist : Vector[(String,Double)] = getSizeFactors(input, ignoreSizeFactors, sizeFactors);
    
    val infileList = initialpairlist.map(_._1).toList;
    val sfList = initialpairlist.map(_._2).toList;
    
    run_helper(opt_title : Option[String], infileList.map(infilePrefix + _ + infileSuffix), outfile , gff , stranded , Some(sfList), opt_rgb , digits , includeFullSpliceNames , calcSum , nonflatgtf );
    
    /*if(infileIsTable){
      val tableLinesTemp = getLinesSmartUnzip(infile.head, true).toList;
      val table : List[(String,Double)] = (if(tableLinesTemp.head.substring(0,9) == "sample.ID") {
        if(tableLinesTemp.head.substring(0,21) != "sample.ID	size.factor"){
          error("Error: first two columns of table must be sample.ID and size.factor. The table must be tab-delimited.");
        }
        tableLinesTemp.tail.map(line => { 
          val cells = line.split("	");
          (cells(0),string2double(cells(1)));
        }).toVector;
      } else {
        tableLinesTemp.map(line => { 
          val cells = line.split("	");
          (cells(0),string2double(cells(1)));
        }).toVector;
      }).toList
      
      val infileList = table.map(_._1);
      val sfList = table.map(_._2);
      
      run_helper(opt_title : Option[String], infileList.map(infilePrefix + _ + infileSuffix), outfile , gff , stranded , Some(sfList), opt_rgb , digits , includeFullSpliceNames , calcSum , nonflatgtf );
    } else {
      run_helper(opt_title : Option[String], infile.map(infilePrefix + _ + infileSuffix), outfile , gff , stranded , opt_sizeFactor, opt_rgb , digits , includeFullSpliceNames , calcSum , nonflatgtf );
    }*/
  }*/
   
  def getIdFromCells(cells : Seq[String]) : String = {
    if(cells.length > INDEX_ID){
      cells(INDEX_ID);
    } else {
      cells(INDEX_CHROM) + ":" + cells(INDEX_START)+"-"+cells(INDEX_END)+cells(INDEX_STRAND);
    }
  }
   
  def run_helper(opt_title : Option[String], infile : Vector[String], outfile : String, 
                 stranded : Boolean, 
                 opt_sizeFactor : Option[Vector[Double]], opt_rgb : Option[String], digits : Int, 
                 calcSum : Boolean, 
                 trackDefLine : Option[String],
                 junctionTypes : List[String]) {
    val rgb = if(opt_rgb.isEmpty) "0,0,0" else opt_rgb.get;

    val sizeFactorSimple = if(opt_sizeFactor.isEmpty) repToSeq(1.toDouble,infile.length) else opt_sizeFactor.get;
    val sizeFactors = if(calcSum) sizeFactorSimple else sizeFactorSimple.map(_ * infile.length.toDouble);
    if(sizeFactors.length != infile.length){
      error("Fatal error: sizeFactors must have the same length as infiles!");
    }
    val title = if(opt_title.isEmpty) "" else opt_title.get + ":";
    
    //val sjmap = makeSpliceJunctionMap(gffLines);
    
    reportln("> makeSpliceJunctionBed: Finished setup.","note");
    val lines : Vector[Iterator[String]] = infile.map(f => {
      reportln("> makeSpliceJunctionBed: opening file: \""+f+"\"", "note");
      val lineIter = getLinesSmartUnzip(f, true);
      lineIter.next; //skip the title line!
      lineIter;
    }).toVector;
    reportln("> makeSpliceJunctionBed: Initialized file iterators.","note");
    
    //val gfflines : Iterator[GtfLine] = internalUtils.GtfTool.GtfReader.getGtfReader(orphanSpliceGff, stranded = stranded, chkFmt = true, fmt_attributeBreak = "\\s+");
    
    val allCounts : Vector[(Iterator[(String,Double,Vector[String])])] = lines.zip(sizeFactors).map(pair => {
      val l = pair._1;
      val sf = pair._2;
      l.map(line => {
        val cells = line.split("	");
        (getIdFromCells(cells), string2double(cells(INDEX_CT)) / sf, cells.toVector);
      });
    });
    
    reportln("> makeSpliceJunctionBed: generated initial counts.","note");

    val counts : Iterator[(String,Double,Vector[String])] = allCounts.tail.foldLeft(allCounts.head)((soFar,curr) => {
      soFar.zip(curr).map(pair => {
        if(pair._1._1 != pair._2._1)  error("ERROR: junctionCount files do not match IDs ("+pair._1._1+") VS ("+pair._2._1+")!");
        if(pair._1._3(INDEX_CHROM)  != pair._2._3(INDEX_CHROM))  error("ERROR: junctionCount files do not match chroms ("+pair._1._3+") VS ("+pair._2._3+")!");
        if(pair._1._3(INDEX_START)  != pair._2._3(INDEX_START))  error("ERROR: junctionCount files do not match starts ("+pair._1._3+") VS ("+pair._2._3+")!");
        if(pair._1._3(INDEX_END)    != pair._2._3(INDEX_END))    error("ERROR: junctionCount files do not match ends ("+pair._1._3+") VS ("+pair._2._3+")!");
        if(pair._1._3(INDEX_STRAND) != pair._2._3(INDEX_STRAND)) error("ERROR: junctionCount files do not match strands ("+pair._1._3+") VS ("+pair._2._3+")!");
        
        (pair._1._1, pair._1._2 + pair._2._2, pair._1._3);
      });
    });
    reportln("> makeSpliceJunctionBed: calculated final counts.","note");

    writeBed(title, counts, outfile=outfile , rgb=rgb, 
             digits = digits, 
             trackDefLine = trackDefLine, 
             junctionTypes = junctionTypes.toSet);
    //convert(infile,outfile, rgb, countFilter, scoreFunction, sizeFactor, includeSpliceNames);
  }
  
  def makeSpliceJunctionMap(gffLines : Iterator[FlatGtfLine]) : scala.collection.GenMap[String,FlatGtfLine] = {
    val out = scala.collection.mutable.AnyRefMap[String,FlatGtfLine]();
    
    for(line <- gffLines){
      if(line.isSpliceJunction){
        out.update(line.getFeatureName, line);
      }
    }
    return out;
  }
  
  def writeBed(title : String, 
               counts : Iterator[(String, Double, Vector[String])], 
               outfile : String, 
               rgb : String, 
               //includeFullSpliceNames : Boolean, 
               digits : Int, 
               delim : String = "	", 
               trackDefLine : Option[String],
               junctionTypes : Set[String]
               //skipAnnotatedJunctions : Boolean,
               //skipNovelJunctions : Boolean
               ) {
    
    val filteredCounts : Iterator[(String, Double, Vector[String])] = if(! junctionTypes.contains("ALL")){
      counts.filter{ case(junctionName,ct,cells) => {
        val junctionType = junctionName.split("_")(0);
        junctionTypes.contains(junctionType);
      }}
    } else {
      counts;
    }
    
      /*
     *   final val INDEX_CHROM = 0;
  final val INDEX_STRAND = 1;
  final val INDEX_START = 2;
  final val INDEX_END = 3;
  final val INDEX_CT = 4;
  //optional:
  final val INDEX_ID = 5;
  final val INDEX_GENESET = 6;
     */
    
    val lines : Vector[((String,Int,Int),String)] = filteredCounts.map{case(junctionName, ct, cells) => {
        val chrom = cells(INDEX_CHROM);
        val strand = cells(INDEX_STRAND);
        val start =  string2int(cells(INDEX_START));
        val end = string2int(cells(INDEX_END));
        val score : Int = math.max(0,math.min(ct.round.toInt,1000));
        val len : Int = end - start;
        
        val bedStart = start - 2;
        val bedEnd = end + 1;
        
        val spliceName = getIdFromCells(cells);
        val lineTitle = title + spliceName + "(" + ("%."+digits.toString+"f").format(ct) + ")";
        
        ((chrom,bedStart,bedEnd),chrom+"	"+bedStart+"	"+bedEnd+"	"+lineTitle+"	"+score+"	"+strand+"	"+start+"	"+end+"	"+rgb+"	2	1,1	0,"+(len + 2)+"\n");
    }}.toVector.sortBy{case(triple,line) => triple};
    
    val writer : WriterUtil = openWriterSmart(outfile, true);
    
    if(! trackDefLine.isEmpty){
      writer.write(trackDefLine.get + "\n");
    }
    
    for((pair,line) <- lines){
      writer.write(line);
    }
    
    writer.close();
  }
}







