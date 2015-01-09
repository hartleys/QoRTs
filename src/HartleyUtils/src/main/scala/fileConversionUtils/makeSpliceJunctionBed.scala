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
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "makeSpliceBed", 
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
                    new BinaryOptionArgument[List[Double]](
                                         name = "sizeFactors", 
                                         arg = List("--sizeFactors"), 
                                         valueName = "val",  
                                         argDesc = "A double-precision floating-point value. If this option is set, all counts will be divided by the given normalization factors. The length must be the same as the length of infiles."
                                        ) ::
                    new UnaryArgument(name = "includeFullSpliceNames",
                              arg = List("--includeFullSpliceNames"), // name of value
                              argDesc = "Flag to indicate that full splice names, including gene ID, should be used." // description
                              ) :: 
                    new UnaryArgument(name = "calcSum",
                              arg = List("--calcSum"), // name of value
                              argDesc = "Flag to indicate that the splice junction counts should be summed, rather than averaged." // description
                              ) :: 
                    new UnaryArgument(name = "infileIsTable",
                              arg = List("--infileIsTable"), // name of value
                              argDesc = "Flag to indicate that the infile isn't a junction counts file itself, but rather a tab-delimited table containing 2 columns: the junction file file path and the size factors." // description
                              ) :: 
                    new UnaryArgument(name = "stranded",
                              arg = List("--stranded"), // name of value
                              argDesc = "Flag to indicate that data is stranded." // description
                              ) :: 
                    new UnaryArgument(name = "nonflatgtf",
                              arg = List("--nonflatgtf"), // name of value
                              argDesc = "Flag to indicate that, instead of a flattened gtf file, a standard gtf file. If this flag is raised, it will automatically create the standard flattened gtf information, in memory. It will not be written to disk" // description
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
                                         name = "title", 
                                         arg = List("--title"), 
                                         valueName = "title",  
                                         argDesc = "A title to be prepended to each splice junction name"
                                        ) ::
                    new FinalArgument[List[String]](
                                         name = "infile",
                                         valueName = "infile[,infile2,infile3,...]",
                                         argDesc = "The input splice counts file, or '-' to read from stdin.  If the filename ends with \".gz\" or \".zip\" then the file will be read using the appropriate decompression method." // description
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
                                        ) :: List() );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         makeSpliceJunctionBed.run(
             parser.get[Option[String]]("title"),
             parser.get[List[String]]("infile"),
             parser.get[String]("outfile"),
             parser.get[String]("infilePrefix"),
             parser.get[String]("infileSuffix"),
             parser.get[String]("gff"),
             parser.get[Boolean]("stranded"),
             parser.get[Option[List[Double]]]("sizeFactors"),
             parser.get[Option[String]]("rgb"),
             parser.get[Int]("digits"),
             parser.get[Boolean]("includeFullSpliceNames"),
             parser.get[Boolean]("calcSum"),
             parser.get[Boolean]("infileIsTable"),
             parser.get[Boolean]("nonflatgtf")
           );
         }
     }
   }
  
  def run(opt_title : Option[String], infile : List[String], outfile : String, infilePrefix : String, infileSuffix : String, gff : String, stranded : Boolean, opt_sizeFactor : Option[List[Double]], opt_rgb : Option[String], digits : Int, includeFullSpliceNames : Boolean, calcSum : Boolean, infileIsTable : Boolean, nonflatgtf : Boolean) {
    if(infileIsTable){
      val tableLinesTemp = getLinesSmartUnzip(infile.head, true).toList;
      val table : List[(String,Double)] = if(tableLinesTemp.head.substring(0,9) == "sample.ID") {
        if(tableLinesTemp.head.substring(0,21) != "sample.ID	size.factor"){
          error("Error: first two columns of table must be sample.ID and size.factor. The table must be tab-delimited.");
        }
        tableLinesTemp.tail.map(line => { 
          val cells = line.split("	");
          (cells(0),string2double(cells(1)));
        }).toList;
      } else {
        tableLinesTemp.map(line => { 
          val cells = line.split("	");
          (cells(0),string2double(cells(1)));
        }).toList;
      }
      
      val infileList = table.map(_._1);
      val sfList = table.map(_._2);
      
      run_helper(opt_title : Option[String], infileList.map(infilePrefix + _ + infileSuffix), outfile , gff , stranded , Some(sfList), opt_rgb , digits , includeFullSpliceNames , calcSum , nonflatgtf );
    } else {
      run_helper(opt_title : Option[String], infile.map(infilePrefix + _ + infileSuffix), outfile , gff , stranded , opt_sizeFactor, opt_rgb , digits , includeFullSpliceNames , calcSum , nonflatgtf );
    }
  }

  
  def run_helper(opt_title : Option[String], infile : List[String], outfile : String, gff : String, stranded : Boolean, opt_sizeFactor : Option[List[Double]], opt_rgb : Option[String], digits : Int, includeFullSpliceNames : Boolean, calcSum : Boolean, nonflatgtf : Boolean) {
    val rgb = if(opt_rgb.isEmpty) "0,0,0" else opt_rgb.get;

    val sizeFactorSimple = if(opt_sizeFactor.isEmpty) repToSeq(1.toDouble,infile.length) else opt_sizeFactor.get;
    val sizeFactors = if(calcSum) sizeFactorSimple else sizeFactorSimple.map(_ * infile.length.toDouble);
    if(sizeFactors.length != infile.length){
      error("Fatal error: sizeFactors must have the same length as infiles!");
    }
    val title = if(opt_title.isEmpty) "" else opt_title.get + ":";
    
    val gffLines = if(nonflatgtf) {
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

    writeBed(title, counts, outfile , rgb, sjmap , includeFullSpliceNames, digits);
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
  
  def writeBed(title : String, counts : Iterator[(String,Double)], outfile : String, rgb : String, sjmap : scala.collection.GenMap[String,FlatGtfLine], includeFullSpliceNames : Boolean, digits : Int, delim : String = "	") {
    
    val lines : Vector[((String,Int),String)] = counts.filter{case(junctionName, ct) =>{ 
      val featureCode = junctionName.split(":")(1).substring(0,1);
      featureCode == "J" || featureCode == "N";
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
    
    for((pair,line) <- lines){
      writer.write(line);
    }
    
    writer.close();
  }
}







