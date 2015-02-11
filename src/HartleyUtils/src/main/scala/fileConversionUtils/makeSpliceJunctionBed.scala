package fileConversionUtils

import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.optionHolder._;
import internalUtils.GtfTool._;

object makeSpliceJunctionBed {

  class converter extends CommandLineRunUtil { 
     override def priority = 25;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "makeJunctionTrack", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This utility takes the splice junction count files created by the QoRTs QC utility "+
                        "across multiple samples "+
                        "and creates a single merged splice junction 'bed' file that lists each splice junction along with the "+
                        "mean read-pair coverage counts (optionally, the mean normalized counts)."+
                        "This splice junction bed file "+
                        "can be used to visualize splice junction counts using the UCSC genome browser "+
                        "and other similar utilities."+
                        ""+
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
                    new UnaryArgument(name = "includeFullSpliceNames",
                              arg = List("--includeFullSpliceNames"), // name of value
                              argDesc = "Flag to indicate that full splice names, including gene ID, should be used." // description
                              ) :: 
                    new UnaryArgument(name = "calcMean",
                              arg = List("--calcMean"), // name of value
                              argDesc = "Flag to indicate that the splice junction counts should be averaged, rather than added up." // description
                              ) :: 
                    new UnaryArgument(name = "ignoreSizeFactors",
                              arg = List("--ignoreSizeFactors"), // name of value
                              argDesc = "Flag to indicate that this utility should ignore size factors even if they are found in the input listFile." // description
                              ) ::
                    new UnaryArgument(name = "stranded",
                              arg = List("--stranded"), // name of value
                              argDesc = "Flag to indicate that data is stranded." // description
                              ) :: 
                    new UnaryArgument(name = "nonflatgtf",
                              arg = List("--nonflatgtf"), // name of value
                              argDesc = "Flag to indicate that instead of a \"flattened\" gff file, a standard-format gtf file has been specified. If this flag is raised, it will automatically create the standard flattened gtf information, in memory. It will not be written to disk" // description
                              ) :: 
                    new UnaryArgument(name = "skipAnnotatedJunctions",
                              arg = List("--skipAnnotatedJunctions"), // name of value
                              argDesc = "If this option is used, annotated splice junctions will not be included in the output file. "+
                                        "Note: this only works if there are novel junctions in the input file."+
                                        ""+
                                        ""
                              ) :: 
                    new UnaryArgument(name = "skipNovelJunctions",
                              arg = List("--skipNovelJunctions"), // name of value
                              argDesc = "If this option is used, novel splice junctions will not be included in the output file."+
                                        ""+
                                        ""+
                                        ""
                              ) ::    
                    new BinaryArgument[Int](   name = "digits",
                                                        arg = List("--digits"),  
                                                        valueName = "num", 
                                                        argDesc = "The number of digits after the decimal to include for counts.", 
                                                        defaultValue = Some(2)
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
                                                   "If this option is set, all counts will be divided by the given normalization factors. The length must be the same as the length of infiles."+
                                                   "If sample.ID's is not specified by the --sampleList or --sampleListFile parameters, then all listed samples will be merged."
                                        ) ::
                    new BinaryOptionArgument[List[Double]](
                                         name = "sizeFactors", 
                                         arg = List("--sizeFactors"), 
                                         valueName = "val",  
                                         argDesc = "A list of double-precision floating-point values. "+
                                                   "If this or any size factor option is set,"+
                                                   " all counts will be divided by the given normalization factors."+
                                                   " The length must be the same as the number of files to merge."
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "filenames", 
                                         arg = List("--filenames"), 
                                         valueName = "file1.wig,file2.wig,file3.wig.gz,...",  
                                         argDesc = "A comma-delimited list of wiggle files to merge. "+
                                                   "This is optional, and filenames can be inferred from --infilePrefix, --infileSuffix, and the --sampleList, if those options are specified."+
                                                   ""+
                                                   ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "sampleList", 
                                         arg = List("--sampleList"), 
                                         valueName = "[sampleList.txt | - | samp1,samp2,samp3,...]",  
                                         argDesc = "Either a comma-delimited list of sample id's or a file containing a list of sample id's."+
                                                   "The file must either contain no title line, or contain a title line that includes a \"sample.ID\" column."+
                                                   ""+
                                                   ""
                                        ) ::
                   new BinaryOptionArgument[String](
                                         name = "title", 
                                         arg = List("--title"), 
                                         valueName = "title",  
                                         argDesc = "A title to be prepended to each splice junction name"
                                        ) ::
                    new FinalArgument[String](
                                         name = "gff",
                                         valueName = "flatgff.gff.gz",
                                         argDesc = "The flattened gff file." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile",
                                         argDesc = "The output bed file, or '-' to write to stdout. If the filename ends with \".gz\" or \".zip\" then the file will be compressed using the appropriate method." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         makeSpliceJunctionBed.run(
             parser.get[Option[String]]("sizeFactorFile"),
             parser.get[Option[List[Double]]]("sizeFactors"),
             parser.get[Option[List[String]]]("filenames"),
             parser.get[Option[String]]("sampleList"),
             parser.get[Option[String]]("title"),
             parser.get[Boolean]("ignoreSizeFactors"),
             parser.get[String]("outfile"),
             parser.get[String]("infilePrefix"),
             parser.get[String]("infileSuffix"),
             None,
             parser.get[String]("gff"),
             parser.get[Boolean]("stranded"),
             parser.get[Int]("digits"),
             parser.get[Boolean]("includeFullSpliceNames"),
             parser.get[Boolean]("calcMean"),
             parser.get[Boolean]("nonflatgtf"),
             parser.get[Option[String]]("rgb"),
             parser.get[String]("trackTitle"),
             parser.get[String]("additionalTrackOptions"),
             parser.get[Boolean]("skipAnnotatedJunctions"),
             parser.get[Boolean]("skipNovelJunctions")
           );
         }
     }
   }
  
  def run(sizeFactorFile : Option[String], sizeFactors : Option[List[Double]], 
          filenames : Option[List[String]], sampleList : Option[String],
          title : Option[String],
          ignoreSizeFactors : Boolean,
          outfile : String, infilePrefix : String, infileSuffix : String, 
          gffIterator : Option[Iterator[FlatGtfLine]], gff  : String, 
          stranded : Boolean,
          digits : Int, includeFullSpliceNames : Boolean, calcMean : Boolean,
          nonflatgtf : Boolean, rgb : Option[String],
          trackTitle : String,
          additionalTrackOptions : String,
          skipAnnotatedJunctions : Boolean,
          skipNovelJunctions : Boolean){

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
    
    run_helper(title, infiles, outfile , gffIterator, gff , stranded , sf, rgb , digits , includeFullSpliceNames , calcMean , nonflatgtf, Some("track name="+trackTitle+" description="+trackTitle+" "+rgbline+additionalTrackOptions), skipAnnotatedJunctions, skipNovelJunctions);
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
   
  def run_helper(opt_title : Option[String], infile : Vector[String], outfile : String, 
                 gffIterator : Option[Iterator[FlatGtfLine]], gff : String, stranded : Boolean, 
                 opt_sizeFactor : Option[Vector[Double]], opt_rgb : Option[String], digits : Int, 
                 includeFullSpliceNames : Boolean, calcSum : Boolean, nonflatgtf : Boolean, trackDefLine : Option[String],
                 skipAnnotatedJunctions : Boolean,
                 skipNovelJunctions : Boolean) {
    val rgb = if(opt_rgb.isEmpty) "0,0,0" else opt_rgb.get;

    val sizeFactorSimple = if(opt_sizeFactor.isEmpty) repToSeq(1.toDouble,infile.length) else opt_sizeFactor.get;
    val sizeFactors = if(calcSum) sizeFactorSimple else sizeFactorSimple.map(_ * infile.length.toDouble);
    if(sizeFactors.length != infile.length){
      error("Fatal error: sizeFactors must have the same length as infiles!");
    }
    val title = if(opt_title.isEmpty) "" else opt_title.get + ":";
    
    val gffLines = if(! gffIterator.isEmpty){
      gffIterator.get;
    } else if(nonflatgtf) {
      fileConversionUtils.prepFlatGtfFile.getFlatGtfLines(gff,stranded).iterator;
    } else {
      GtfReader.getFlatGtfReader(gff, stranded, true, "\\s+");
    }
    reportln("> makeSpliceJunctionBed: initialized gff reader.","note");

    val sjmap = makeSpliceJunctionMap(gffLines);
    reportln("> makeSpliceJunctionBed: build splice junction map.","note");

    reportln("> makeSpliceJunctionBed: Finished setup.","note");
    val lines : Vector[Iterator[String]] = infile.map(f => {
      reportln("> makeSpliceJunctionBed: opening file: \""+f+"\"", "note");
      getLinesSmartUnzip(f, true);
    }).toVector;
    reportln("> makeSpliceJunctionBed: Initialized file iterators.","note");
    val allCounts : Vector[(Iterator[(String,Double)])] = lines.zip(sizeFactors).map(pair => {
      val l = pair._1;
      val sf = pair._2;
      l.map(line => {
        val cells = line.split("	");
        (cells(0), string2double(cells(1)) / sf);
      });
    });
    reportln("> makeSpliceJunctionBed: generated initial counts.","note");

    val counts : Iterator[(String,Double)] = allCounts.tail.foldLeft(allCounts.head)((soFar,curr) => {
      soFar.zip(curr).map(pair => {
        if(pair._1._1 != pair._2._1) error("ERROR: junctionCount files do not match!");
        (pair._1._1, pair._1._2 + pair._2._2);
      });
    });
    reportln("> makeSpliceJunctionBed: calculated final counts.","note");

    writeBed(title, counts, outfile , rgb, sjmap , includeFullSpliceNames, digits, trackDefLine = trackDefLine, skipAnnotatedJunctions = skipAnnotatedJunctions, skipNovelJunctions = skipNovelJunctions);
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
               counts : Iterator[(String,Double)], 
               outfile : String, 
               rgb : String, 
               sjmap : scala.collection.GenMap[String,FlatGtfLine], 
               includeFullSpliceNames : Boolean, 
               digits : Int, 
               delim : String = "	", 
               trackDefLine : Option[String],
               skipAnnotatedJunctions : Boolean,
               skipNovelJunctions : Boolean) {
    
    val lines : Vector[((String,Int),String)] = counts.filter{case(junctionName, ct) =>{ 
      val featureCode = junctionName.split(":")(1).substring(0,1);
      (featureCode == "J" && (! skipAnnotatedJunctions)) || (featureCode == "N" && (! skipNovelJunctions));
    }}.map{case(junctionName, ct) =>{
        val gffLineOption : Option[FlatGtfLine] = sjmap.get(junctionName);
        if(gffLineOption.isEmpty) error("ERROR: gff does not contain splice junction named \""+junctionName+"\"! Are you using the wrong flattened gff file?");
        val gffLine : FlatGtfLine = gffLineOption.get;
        val chrom = gffLine.chromName;
        val strand = gffLine.strand;
        val start = gffLine.start;
        val end = gffLine.end;
        val score : Int = math.max(0,math.min(ct.round.toInt,1000));
        val len : Int = end - start;
        
        val bedStart = start - 2;
        val bedEnd = end + 1;
      
        val spliceName = if(includeFullSpliceNames) gffLine.getFeatureName else gffLine.getFeatureCode + gffLine.getFeaturePartNumber;
        val lineTitle = title + spliceName + "(" + ("%."+digits.toString+"f").format(ct) + ")";

        ((chrom,bedStart),chrom+"	"+bedStart+"	"+bedEnd+"	"+lineTitle+"	"+score+"	"+strand+"	"+start+"	"+end+"	"+rgb+"	2	1,1	0,"+(len + 2)+"\n");
    }}.toVector.sortBy{case(pair,line) => pair};
    
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







