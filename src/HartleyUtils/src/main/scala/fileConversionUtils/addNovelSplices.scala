package fileConversionUtils



import net.sf.samtools._
import internalUtils.Reporter._;
import internalUtils.stdUtils._;
import internalUtils.fileUtils._;
import internalUtils.commandLineUI._;
import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
import scala.collection.JavaConversions._;

import internalUtils.genomicUtils._;
import internalUtils.optionHolder._;

object addNovelSplices {

  class mergeNovelSplices extends CommandLineRunUtil {
     override def priority = 29;
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "mergeNovelSplices", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This utility takes the QC output from the standard QC utility run on a series of samples "+
                        "and performs two functions: first, it compiles all splice junctions across all samples and "+
                        "filters low-coverage novel splice junctions by mean coverage across all samples (optionally normalized with user-supplied size factors). "+
                        "It then assigns unique identifiers to each novel splice junction that passed this filter, and outputs a special flat gff file "+
                        "listing all exons, annotated splice junctions "+
                        "and passed-filter novel splice junctions with assigned unique identifiers for all features. "+
                        "Next, it uses these unique identifiers to create a new set of JunctionSeq-formatted count files, one for each "+
                        "input sample. This new count file will include counts for the passed-filter novel splice junctions "+
                        "in addition to the usual counts for annotated splice junctions, exons, and aggregated-genes, all listed by the assigned unique identifiers."+
                        "",   
          argList = 
                    new BinaryArgument[Double](name = "minCount",
                                                        arg = List("--minCount"),  
                                                        valueName = "num", 
                                                        argDesc = "The minimum mean normalized read coverage needed for inclusion of a novel splice junction. By default, equal to 10.0.", 
                                                        defaultValue = Some(10.0)
                                                        ) :: 
                    new UnaryArgument( name = "stranded",
                                         arg = List("--stranded","-s"), // name of value
                                         argDesc = "Flag to indicate that data is stranded." // description
                                       ) ::                 
                    new UnaryArgument( name = "stranded_fr_secondstrand",
                                         arg = List("--stranded_fr_secondstrand"), // name of value
                                         argDesc = "Nonfunctional." // description
                                       ) ::                                            
                    new UnaryArgument( name = "noGzipInput",
                                         arg = List("--noGzipInput"), // name of value
                                         argDesc = "Flag to indicate that input files are NOT be compressed into the gzip format. By default almost all input files are assumed to be compressed." // description
                                       ) ::
                    new UnaryArgument( name = "noGzipOutput",
                                         arg = List("--noGzipOutput"), // name of value
                                         argDesc = "Flag to indicate that output files should NOT be compressed into the gzip format. By default almost all output files are compressed to save space." // description
                                       ) ::
                    new FinalArgument[String](
                                         name = "infileDir",
                                         valueName = "infileDir",
                                         argDesc = "The input file directory." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "sizeFactorFile",
                                         valueName = "sizeFactorFile",
                                         argDesc = "" // description
                                        ) :: 
                    new FinalArgument[String](
                                         name = "gtfFile",
                                         valueName = "annotation.gtf.gz",
                                         argDesc = "" // description
                                        ) ::    
                    new FinalArgument[String](
                                         name = "outfileDir",
                                         valueName = "outfileDir",
                                         argDesc = "The output file directory" // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         addNovelSplices.mergeNovel(
             parser.get[String]("infileDir") + "/",
             parser.get[String]("sizeFactorFile"),
             parser.get[String]("gtfFile"),
             parser.get[String]("outfileDir"),
             parser.get[Boolean]("stranded"),
             parser.get[Double]("minCount"),
             parser.get[Boolean]("noGzipOutput"),
             parser.get[Boolean]("noGzipInput")
           );
         }
     }
   }
  
  def mergeNovel(infileDir : String, decoderFile : String, gtfFile : String, outfileDir : String, stranded : Boolean, minCount : Double, noGzipOutput : Boolean, noGzipInput : Boolean){
     //internalUtils.optionHolder.OPTION_noGzipOutput = noGzipOutput;
     
     val decoderLines = getLines(decoderFile).toVector;
     val decoderTitleLine = decoderLines.head.split("	").toVector;
     val idCol = decoderTitleLine.indexOf("sample.ID");
     val sfCol = decoderTitleLine.indexOf("size.factor");
     
     reportln("decoderTitleLine: "+decoderTitleLine,"debug");
     reportln("idCol: "+idCol,"debug");
     reportln("sfCol: "+sfCol,"debug");
     
     val sampleSF : Seq[(String,Double)] = decoderLines.tail.map((line : String) => {
       val cells = line.split("	");
       ((cells(idCol),
           string2double(cells(sfCol))));
     });
     
     val flatlines = fileConversionUtils.prepFlatGtfFile.getFlatGtfLines(gtfFile,stranded);
     
     mergeNovelHelper(infileDir , sampleSF , flatlines, outfileDir , stranded, minCount ,noGzipOutput , noGzipInput);
  }
  
  //INCOMPLETE!!!!!!!!!!!!!!!!
  
  def mergeNovelHelper(infileDir : String, sampleSF : Seq[(String,Double)], flatgff : IndexedSeq[FlatGtfLine], outfileDir : String, stranded : Boolean, minCount : Double, noGzipOutput : Boolean, noGzipInput : Boolean){
    
    val Splice_suffix = "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt" + (if(noGzipInput){""} else {".gz"});
    val NovelSplice_suffix = "/QC.spliceJunctionCounts.novelSplices.txt" + (if(noGzipInput){""} else {".gz"});
    val WithNovel_Splice_suffix = "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt" + (if(noGzipOutput){""} else {".gz"});
    val novel_gff_suffix = "/withNovel.forJunctionSeq.gff" + (if(noGzipOutput){""} else {".gz"});
    val novel_gff_file = outfileDir + novel_gff_suffix; 
    
    
    //Make a gene map out of aggregate genes.
    val geneArray : GenomicArrayOfSets[String] = buildGenomicArrayOfSets_fromGtf(stranded , (() => flatgff.iterator),  (g : FlatGtfLine) => g.isAggregateGene, (g : FlatGtfLine) => g.getFeatureAggregateGene)
    
    //read in all files and all novel splice junction counts.
    val novelSpliceMap : scala.collection.mutable.Map[GenomicInterval,Array[Int]] = scala.collection.mutable.AnyRefMap[GenomicInterval,Array[Int]]();
    
    for(i <- 0 until sampleSF.length){
      val (name, sf) = sampleSF(i);
      val lines = getLinesSmartUnzip(infileDir + name + NovelSplice_suffix);
      lines.next;
      for(line <- lines){
        val cells = line.split("	");
        val ct = string2int(cells(4));
        val iv = GenomicInterval(cells(0), cells(1).charAt(0), string2int(cells(2)), string2int(cells(3)));
        if(novelSpliceMap.contains(iv)){
          val counts = novelSpliceMap(iv);
          counts(i) = ct;
        } else {
          val counts = Array.ofDim[Int](sampleSF.length);
          counts(i) = ct;
          novelSpliceMap.update(iv, counts)
        }
      }
    }
    
    reportln("Finished reading novel splice files. Found "+novelSpliceMap.size+" novel splice junctions.","debug");
    
    //determine which splice junctions satisfy the filtering requirements.
    val pfSpliceMap : scala.collection.mutable.Map[GenomicInterval,Array[Int]] = novelSpliceMap.filter{ case (iv : GenomicInterval, counts : Array[Int]) => {
      val meanNormCounts = counts.toVector.zip(sampleSF).foldLeft(0.0)( (soFar,curr) => {
        val (ct, (name, sf)) = curr;
        soFar + (ct.toDouble * sf);
      });
      meanNormCounts > minCount;
    }}
    
    reportln("Finished filtering novel splice files. Found "+pfSpliceMap.size+" novel splice junctions that pass expression coverage filters.","debug");
    
    val passedFilter = pfSpliceMap.size;
    
    //determine which splice junctions belong to known genes.
    val geneNovelSplices : scala.collection.mutable.Map[String,Set[GenomicInterval]] = scala.collection.mutable.AnyRefMap[String,Set[GenomicInterval]]().withDefault(str => Set[GenomicInterval]());
    var addedSplices = 0;
    var dropped_orphanSplices = 0;
    var dropped_ambigSplices = 0;
    
    for((iv, counts) <- pfSpliceMap){
      val intersectingGenes : Set[String] = geneArray.findIntersectingSteps(iv).foldLeft(Set[String]())( (soFar, curr) => soFar ++ curr._2);
      if(intersectingGenes.size == 1){
        val gene = intersectingGenes.head
        geneNovelSplices.update(gene, geneNovelSplices(gene) + iv  );
        addedSplices += 1;
      } else if(intersectingGenes.size == 0){
        dropped_orphanSplices += 1;
      } else {
        dropped_ambigSplices += 1;
      }
    }
    
    reportln("Finished assigning splice junctions to genes.\n"+
             "                  Found "+addedSplices+" that fall inside one unique gene.\n"+
             "                  Found "+dropped_orphanSplices+" orphaned.\n"+
             "                  Found "+dropped_ambigSplices+" ambiguous.","debug");
    
    //Read in the existing gff:
    val gffLineMap = scala.collection.mutable.AnyRefMap[String,Vector[FlatGtfLine]]().withDefault(str => Vector[FlatGtfLine]());
    
    for(flatGtfLine <- flatgff){
      val gene = flatGtfLine.getFeatureAggregateGene;
      gffLineMap.update(gene, gffLineMap(gene) :+ flatGtfLine);
    }
    
    reportln("Finished reading flat gff. Found "+gffLineMap.size+" aggregate genes.","debug");
    
    //Add in the new splice junctions, and number them.
    for((gene, ivSet) <- geneNovelSplices){
      val ivList = ivSet.toVector.sorted;
      var elementCt = string2int(gffLineMap(gene).last.getFeaturePartNumber);
      
      for(iv <- ivList){
        elementCt += 1;
        val newGffLine = FlatOutputGtfLine.makeFlatGtfLine_feature(iv , defaultGtfCodes.JS_FEATURETYPE_NOVELSPLICE, stranded , gene, Set[String]("UNKNOWN_TX"), Set[String]("UNKNOWN_GENE_SET"), elementCt);
        gffLineMap.update(gene, gffLineMap(gene) :+ newGffLine);
      }
    }
    
    reportln("Finished adding to flat gff. Found "+gffLineMap.size+" aggregate genes.","debug");
    
    //output flat gff with all splice junctions
    
    val geneList = gffLineMap.keySet.toVector.map( g => (gffLineMap(g).head.getGenomicInterval, gffLineMap(g).head.getFeatureAggregateGene) ).sorted;
    
    val writer = openWriterSmart(novel_gff_file);
    for((geneIv, geneName) <- geneList){
      val gffLines = gffLineMap(geneName);
      for(line <- gffLines){
        writer.write(line.str + "\n");
      }
    }
    close(writer);
    
    val gffLines = GtfReader.getFlatGtfReader(novel_gff_file, stranded, true, "\\s+").toVector;
    
    reportln("Final gff feature count: "+gffLines.size+"","debug");
    
    for(i <- 0 until sampleSF.length){
      val (name,sf) = sampleSF(i);
      
      val inDir = infileDir + "/" + name + "/";
      val outDir = outfileDir + "/" + name + "/";
      val outDirFile = new java.io.File(outDir);
      if(! outDirFile.exists()){
        reportln("Creating Directory: "+ outDirFile,"note");
        outDirFile.mkdir();
        reportln("Successfully Created Directory: " + outDirFile, "note");
      }
      
      makeJunctionCountsWithNovel(inDir + Splice_suffix, inDir + NovelSplice_suffix, outDir + WithNovel_Splice_suffix, gffLines );
    }
    
    
    reportln("Done.","debug");
  }
  
  def buildGenomicArrayOfSets_fromGtf(stranded : Boolean, makeReader : (() => Iterator[FlatGtfLine]), lineFilter : ( FlatGtfLine => Boolean ) , elementExtractor : (FlatGtfLine => String)) : GenomicArrayOfSets[String] = {
    //report("reading Gtf: " + gtffile + "\n","note");
    val geneArray : GenomicArrayOfSets[String] = GenomicArrayOfSets[String](stranded);
    val gtfReader = makeReader();
    
    for(gtfLine <- gtfReader){
      if(lineFilter(gtfLine)){
        val element = elementExtractor(gtfLine);
        geneArray.addSpan(gtfLine.getGenomicInterval, element);
      }
    }
    return geneArray.finalizeStepVectors;
  }
  
  def makeJunctionCountsWithNovel(infileJS : String, infileNovel : String, outfile : String, flatgff : IndexedSeq[FlatGtfLine]){
    //val Splice_suffix = "/QC.spliceJunctionAndExonCounts.forSpliceSeq.txt.gz";
    //val NovelSplice_suffix = "/QC.spliceJunctionCounts.novelSplices.txt.gz";
    //val WithNovel_Splice_suffix = "/QC.spliceJunctionAndExonCounts.withNovel.forSpliceSeq.txt.gz";
    
    reportln("=====> Running on file:"+infileJS+"","debug");
    
    val novelMap = scala.collection.mutable.AnyRefMap[GenomicInterval,String]();
    for(gffLine <- flatgff){
      if(gffLine.isNovelSJ){
        novelMap.update(gffLine.getGenomicInterval,gffLine.getFeatureName);
      }
    }
    
    reportln("       "+novelMap.size+" novel junctions.","debug");
    
    val countMap = scala.collection.mutable.AnyRefMap[String,Int]().withDefault(str => 0);
    //Read in known junctionSeq counts.
    val knownLines = getLinesSmartUnzip(infileJS);
    //knownLines.next;
    for(line <- knownLines){
      val cells = line.split("	");
      //if(cells.length < 6)reportln("Line length < 6. Offending line is:\n\""+line+"\"","debug");
      val featureName = (cells(0));
      val ct = string2int(cells(1));
      countMap.update(featureName,ct);
    }
    
    reportln("       "+countMap.size+" known features.","debug");
    
    //Read in novel junctionSeq counts.
    val novelLines = getLinesSmartUnzip(infileNovel);
    novelLines.next;
    for(line <- novelLines){
      val cells = line.split("	");
      val ct = string2int(cells(4));
      val iv = GenomicInterval(cells(0), cells(1).charAt(0), string2int(cells(2)), string2int(cells(3)));
      
      if(novelMap.contains(iv)){
        val featureID = novelMap(iv);
        countMap.update(featureID, ct);
      }
    }
    
    reportln("       "+countMap.size+" known+novel features.","debug");
    
    //output JunctionSeq counts.
    val writer = openWriterSmart(outfile);
    for(gffLine <- flatgff){
      val id = gffLine.getFeatureName;
      val ct = countMap(id);
      writer.write(id + "	"+ct + "\n");
    }
    close(writer);
  }
}












