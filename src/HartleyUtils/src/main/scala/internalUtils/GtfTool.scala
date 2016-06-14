package internalUtils

import internalUtils.fileUtils._;
import internalUtils.Reporter._;
//import internalUtils.stdUtils._;

object GtfTool {


  
     def getFeatureID(featureName : String) : String = {
       featureName.split(":")(1);
     }
     def getAggregateGeneID(featureName : String) : String = {
       featureName.split(":")(0);
     }
     //val JS_FEATURETYPE_CODEMAP = Map(JS_FEATURETYPE_GENE -> "A", JS_FEATURETYPE_EXON -> "E", JS_FEATURETYPE_KNOWNSPLICE -> "J", JS_FEATURETYPE_NOVELSPLICE -> "N");
  
  /* object GtfCodes {
     val STD_EXON_TYPE_CODE = "exon";
     val STD_CDS_TYPE_CODE = "CDS";
   
     val GENE_ID_ATTRIBUTE_KEY = "gene_id";
     val STD_TX_ID_ATTRIBUTE_KEY = "transcript_id";
     
     val JS_TX_ID_ATTRIBUTE_KEY = "transcripts";
     val JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY = "exonic_part_number";
     val JS_GENE_SET_ATTRIBUTE_KEY = "gene_set";
     
     val JS_FEATURETYPE_GENE = "aggregate_gene";
     val JS_FEATURETYPE_EXON = "exonic_part";
     val JS_FEATURETYPE_KNOWNSPLICE = "splice_site";
     val JS_FEATURETYPE_NOVELSPLICE = "novel_splice_site";
     
     val JS_FEATURETYPE_CODEMAP = Map(JS_FEATURETYPE_GENE -> "A", JS_FEATURETYPE_EXON -> "E", JS_FEATURETYPE_KNOWNSPLICE -> "J", JS_FEATURETYPE_NOVELSPLICE -> "N");
     
     def getJSFeatureName(gtfLine : GtfLine) : String = {
       val code = JS_FEATURETYPE_CODEMAP(gtfLine.featureType);
       return gtfLine.getAttributeOrDie(GENE_ID_ATTRIBUTE_KEY) + ":" + code + gtfLine.getAttributeOrDie(JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY);
     }
     
     def getJSSubFeatureID(featureName : String) : String = {
       featureName.split(":")(1);
     }
     def getJSGeneFeatureID(featureName : String) : String = {
       featureName.split(":")(0);
     }
   }*/
   
   object defaultGtfCodes {
     val STD_EXON_TYPE_CODE = "exon";
     val STD_CDS_TYPE_CODE = "CDS";
     val BIOTYPE_ATTRIBUTE_KEY = "gene_biotype";
     
     val GENE_ID_ATTRIBUTE_KEY = "gene_id";
     val STD_TX_ID_ATTRIBUTE_KEY = "transcript_id";
     
     val JS_TX_ID_ATTRIBUTE_KEY = "tx_set";
     val JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY = "num";
     val JS_GENE_SET_ATTRIBUTE_KEY = "gene_set";
     
     val JS_FEATURETYPE_GENE = "aggregate_gene";
     val JS_FEATURETYPE_EXON = "exonic_part";
     val JS_FEATURETYPE_KNOWNSPLICE = "splice_site";
     val JS_FEATURETYPE_NOVELSPLICE = "novel_splice_site";
     val JS_AGGREGATEGENE_STRAND = "aggregateGeneStrand";
     val JS_AGGREGATEGENE_CT = "geneCt";
     val JS_AGGREGATEGENE_TXLIST = "tx_set";
     val JS_AGGREGATEGENE_TXCT = "tx_ct";
     val JS_AGGREGATEGENE_TXSTRANDS = "tx_strands";
     
     val JS_TRANSCRIPT_CDS_SPANS = "tx_cds_spans";
     
     val JS_FEATURETYPE_CODEMAP = Map(JS_FEATURETYPE_GENE -> "A", JS_FEATURETYPE_EXON -> "E", JS_FEATURETYPE_KNOWNSPLICE -> "J", JS_FEATURETYPE_NOVELSPLICE -> "N");
     
     def getJSFeatureName(gtfLine : GtfLine) : String = {
       val code = JS_FEATURETYPE_CODEMAP(gtfLine.featureType);
       return gtfLine.getAttributeOrDie(GENE_ID_ATTRIBUTE_KEY) + ":" + code + gtfLine.getAttributeOrDie(JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY);
     }
     
     def getJSSubFeatureID(featureName : String) : String = {
       featureName.split(":")(1);
     }
     def getJSGeneFeatureID(featureName : String) : String = {
       featureName.split(":")(0);
     }
   }
   
   case class GtfCodes(STD_EXON_TYPE_CODE : String = defaultGtfCodes.STD_EXON_TYPE_CODE,
                       STD_CDS_TYPE_CODE : String = defaultGtfCodes.STD_CDS_TYPE_CODE,
                       GENE_ID_ATTRIBUTE_KEY : String = defaultGtfCodes.GENE_ID_ATTRIBUTE_KEY,
                       STD_TX_ID_ATTRIBUTE_KEY : String = defaultGtfCodes.STD_TX_ID_ATTRIBUTE_KEY,
                       JS_TX_ID_ATTRIBUTE_KEY : String = defaultGtfCodes.JS_TX_ID_ATTRIBUTE_KEY,
                       JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY : String = defaultGtfCodes.JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY,
                       JS_GENE_SET_ATTRIBUTE_KEY : String = defaultGtfCodes.JS_GENE_SET_ATTRIBUTE_KEY,
                       JS_FEATURETYPE_GENE : String = defaultGtfCodes.JS_FEATURETYPE_GENE,
                       JS_FEATURETYPE_EXON : String = defaultGtfCodes.JS_FEATURETYPE_EXON,
                       JS_FEATURETYPE_KNOWNSPLICE : String = defaultGtfCodes.JS_FEATURETYPE_KNOWNSPLICE,
                       JS_FEATURETYPE_NOVELSPLICE : String = defaultGtfCodes.JS_FEATURETYPE_NOVELSPLICE,
                       JS_FEATURETYPE_CODEMAP : Map[String,String] = defaultGtfCodes.JS_FEATURETYPE_CODEMAP,
                       JS_AGGREGATEGENE_STRAND : String = defaultGtfCodes.JS_AGGREGATEGENE_STRAND,
                       JS_AGGREGATEGENE_CT : String = defaultGtfCodes.JS_AGGREGATEGENE_CT,
                       JS_AGGREGATEGENE_TXCT : String = defaultGtfCodes.JS_AGGREGATEGENE_TXCT,
                       JS_AGGREGATEGENE_TXSTRANDS : String = defaultGtfCodes.JS_AGGREGATEGENE_TXSTRANDS,
                       BIOTYPE_ATTRIBUTE_KEY : String = defaultGtfCodes.BIOTYPE_ATTRIBUTE_KEY,
                       JS_TRANSCRIPT_CDS_SPANS : String = defaultGtfCodes.JS_TRANSCRIPT_CDS_SPANS
                       ) {
     
     val KEY_SORTING : List[String] = List[String](  GENE_ID_ATTRIBUTE_KEY, 
                                                     STD_TX_ID_ATTRIBUTE_KEY,
                                                     JS_TX_ID_ATTRIBUTE_KEY,
                                                     JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY,
                                                     JS_GENE_SET_ATTRIBUTE_KEY,
                                                     JS_AGGREGATEGENE_STRAND,
                                                     JS_AGGREGATEGENE_CT,
                                                     JS_AGGREGATEGENE_TXCT,
                                                     JS_AGGREGATEGENE_TXSTRANDS,
                                                     JS_TRANSCRIPT_CDS_SPANS,
                                                     BIOTYPE_ATTRIBUTE_KEY
                                                     );
   }
   
   
   
   private var checkFmt = true;
   //private var gtfFmt_attributeBreak = "=";
   
   
   abstract class GtfLine {
     def codes : GtfCodes;
     def str : String;
     def stranded : Boolean;
     def gtfFmt_attributeBreak : String;
     def fmtError(message : String){
       error("FATAL ERROR: OutputGtfLine checkFmt: " +message+ "\n");
     }
     def cells : Array[String];
     
     def chromName : String;
     def featureSource : String;
     def featureType : String;
     def start : Int;
     def end : Int;
     def score : String;
     def strand : Char;
     def strandedStrand : Char;
     def attr : String;
     def attributeArray : Array[String];
     def attributeMap : Map[String,String];
     def getAttribute(key : String) : Option[String];
     def getAttributeOrDie(key : String) : String;
     
     
     def getGenomicInterval : internalUtils.commonSeqUtils.GenomicInterval = {
       return new internalUtils.commonSeqUtils.GenomicInterval(this.chromName, this.strand, this.start - 1, this.end);
     }
   }
   
   object FlatOutputGtfLine {
     val DEF_FEATURESOURCE = "ScalaUtils";
     val DEF_SCORE = ".";
     val DEF_ATTRBREAK = " ";
     
     val DEXSEQ_FEATURETYPE_GENE = "aggregate_gene";
     val DEXSEQ_FEATURETYPE_EXON = "exonic_part";
     val DEXSEQ_GENE_ID_ATTRIBUTE_KEY = "gene_id";
     val DEXSEQ_TX_ID_ATTRIBUTE_KEY = "transcripts";
     val DEXSEQ_EXONIC_PART_NUMBER_ATTRIBUTE_KEY = "exonic_part_number";
     
     def makeFlatGtfLine(iv : internalUtils.commonSeqUtils.GenomicInterval, featureType : String, attributeMap : Map[String,String], stranded : Boolean, codes : GtfCodes) : FlatOutputGtfLine = {
       new FlatOutputGtfLine(iv.chromName, DEF_FEATURESOURCE, featureType, iv.start + 1, iv.end, DEF_SCORE, iv.strand, attributeMap, DEF_ATTRBREAK, stranded, Some(codes.KEY_SORTING), codes);
     }
     
     def makeFlatGtfLine_aggregateGene(iv : internalUtils.commonSeqUtils.GenomicInterval, stranded : Boolean, aggregateGene : String, geneStrand : Char, geneCt : Int, txInfoMap : Map[String,Char], txCDS : Map[String,(Int,Int)], codes : GtfCodes = new GtfCodes()) : FlatOutputGtfLine = {
       val attributeMap = Map[String,String](codes.GENE_ID_ATTRIBUTE_KEY -> aggregateGene, 
                                             codes.JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY -> "000", 
                                             codes.JS_AGGREGATEGENE_STRAND -> geneStrand.toString,
                                             codes.JS_AGGREGATEGENE_CT -> geneCt.toString,
                                             codes.JS_AGGREGATEGENE_TXCT -> txInfoMap.keySet.size.toString,
                                             codes.JS_TX_ID_ATTRIBUTE_KEY -> txInfoMap.keySet.toVector.sorted.mkString("+"),
                                             codes.JS_AGGREGATEGENE_TXSTRANDS -> txInfoMap.keySet.toVector.sorted.map(txInfoMap(_).toString).mkString(","),
                                             codes.JS_TRANSCRIPT_CDS_SPANS -> txInfoMap.keySet.toVector.sorted.map((txID : String) => {
                                                                                 val (s,e) = txCDS(txID);
                                                                                 s + "_" + e;
                                                                              }).mkString(",")
                                             )
       return makeFlatGtfLine(iv,codes.JS_FEATURETYPE_GENE, attributeMap, stranded, codes);
     }
     
     def makeFlatGtfLine_feature(iv : internalUtils.commonSeqUtils.GenomicInterval, featureType : String, stranded : Boolean, aggregateGene : String, txSet : Set[String], geneSet : Set[String], featureNumber : Int, codes : GtfCodes = new GtfCodes()) : FlatOutputGtfLine = {
       var attributeMap = Map[String,String]();
       attributeMap = attributeMap.updated(codes.GENE_ID_ATTRIBUTE_KEY,aggregateGene);
       attributeMap = attributeMap.updated(codes.JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY, internalUtils.stdUtils.zeroPad(featureNumber,3));
       attributeMap = attributeMap.updated(codes.JS_TX_ID_ATTRIBUTE_KEY,txSet.mkString("+"));
       attributeMap = attributeMap.updated(codes.JS_GENE_SET_ATTRIBUTE_KEY, geneSet.mkString("+"));
       
       return makeFlatGtfLine(iv, featureType, attributeMap, stranded, codes);
     }
     
     
     
   }
   
   var ERROR_COUNT_SORTING = 0;
   class OutputGtfLine(in_chromName : String, in_featureSource : String, in_featureType : String, in_start : Int, in_end : Int, in_score : String, in_strand : Char, in_attributeMap : Map[String,String], in_gtfFmt_attributeBreak : String, in_stranded : Boolean, attribute_sorting : Option[List[String]] = None, in_codes : GtfCodes = new GtfCodes()) extends GtfLine {
     def codes = in_codes;
     def chromName = in_chromName;
     def featureSource = in_featureSource;
     def featureType = in_featureType;
     def start = in_start;
     def end = in_end;
     def score = in_score;
     def strand = if(stranded) in_strand else '.';
     def strandedStrand : Char = in_strand;
     def attributeMap = in_attributeMap;
     def gtfFmt_attributeBreak = in_gtfFmt_attributeBreak;

     def stranded = in_stranded;

     def getAttribute(key : String) : Option[String] = attributeMap.get(key);
     def getAttributeOrDie(key : String) : String = {
          attributeMap.get(key) match {
            case Some(value) => value;
            case None => fmtError("Attribute " + key + " not found!"); null;
          }
     }
     
     lazy val lz_attributeArray : Array[String] = {
       attributeMap.map((kv) => kv._1 + gtfFmt_attributeBreak + kv._2).toArray
     }
     def attributeArray : Array[String] = lz_attributeArray;
     
     lazy val lz_attr : String = {
       //This complex mechanism is designed to ensure that GTF attributes always come in a consistent and defined ordering.
       attribute_sorting match {
         case Some(sortList) => {
           if(! attributeMap.keySet.forall(attribName => sortList.contains(attribName))){
             val sp = attributeMap.keySet.span(attribName => sortList.contains(attribName));
             if(ERROR_COUNT_SORTING < 5){
               val a = sp._1.toList.mkString(",");
               val b = sp._2.toList.mkString(",");
               reportln("         Debugging notice: Not all GTF features are annotated: ("+ a +")("+ b +"). (This is ok).","debug");
               ERROR_COUNT_SORTING += 1;
             }
             lz_attributeArray.mkString("; ");
             //Sort the ones you recognize in the assigned order, the remainder go on the end in lexicographical order.
             val sortedAttrKeys   = sortList.filter(sp._1.contains(_));
             val unsortedAttrKeys = sp._2.toList.sorted;
             val keySorting = sortedAttrKeys ++ unsortedAttrKeys;
             
             keySorting.map(attribName => attribName + gtfFmt_attributeBreak + attributeMap(attribName)).mkString("; ");
           } else {
             sortList.filter(attributeMap.contains(_)).map(attribName => attribName + gtfFmt_attributeBreak + attributeMap(attribName)).mkString("; ");
           }
         }
         case None => {
           lz_attributeArray.mkString("; ");
         }
       }
     }
     def attr : String = lz_attr;
     
     lazy val lz_cells : Array[String] = {
       val out = new Array[String](9);
       out(0) = chromName;
       out(1) = featureSource;
       out(2) = featureType;
       out(3) = start.toString;
       out(4) = end.toString;
       out(5) = score.toString;
       out(6) = strand.toString;
       out(7) = ".";
       out(8) = attr;
       out;
     }
     def cells : Array[String] = lz_cells;
     
     lazy val lz_str : String = cells.mkString("	");
     def str = lz_str;
   } 
   
   class InputGtfLine(in_str : String, in_stranded : Boolean, in_gtfFmt_attributeBreak : String , in_codes : GtfCodes = new GtfCodes()) extends GtfLine {
     def codes : GtfCodes = in_codes;
     def str : String = in_str;
     def stranded : Boolean = in_stranded;
     def gtfFmt_attributeBreak : String = in_gtfFmt_attributeBreak;
     lazy val lz_cells : Array[String] = {
          val c = str.split("	");
          if(checkFmt && c.length < 9) fmtError("Line has too few columns!");
          c;
        }
     def cells : Array[String] = lz_cells;
     
     lazy val lz_chromName = internalUtils.stdUtils.cleanQuotes(cells(0));
     def chromName : String = lz_chromName;
     lazy val lz_featureSource = internalUtils.stdUtils.cleanQuotes(cells(1));
     def featureSource : String = lz_featureSource;
     lazy val lz_featureType = internalUtils.stdUtils.cleanQuotes(cells(2));
     def featureType : String = lz_featureType;
     lazy val lz_start = internalUtils.stdUtils.string2int(cells(3));
     def start : Int = lz_start;
     lazy val lz_end = internalUtils.stdUtils.string2int(cells(4));
     def end : Int = lz_end;
     lazy val lz_score = internalUtils.stdUtils.cleanQuotes(cells(5));
     def score : String = lz_score;
     lazy val lz_strand = if(stranded) cells(6).trim().charAt(0) else '.';
     def strand : Char = lz_strand;
     lazy val lz_strandedStrand = cells(6).trim().charAt(0);
     def strandedStrand : Char = lz_strandedStrand;
     lazy val lz_attr = cells(8).trim;
     def attr : String = lz_attr;
     
     //Dumb version, doesn't properly capture quotes:
     //lazy val lz_attributeArray = attr.split(";");
     //Smart version: semicolons inside quotes are allowed, but will be automatically replaced with underscores.
     lazy val lz_attributeArray = attr.split(";(?=([^\"]*\"[^\"]*\")*[^\"]*$)").map(_.replace(";","_"));
     def attributeArray : Array[String] = lz_attributeArray;
     
     lazy val lz_attributeMap : Map[String,String] = {
          var mymap = Map[String,String]();
          for(attributePair <- attributeArray){
            val pairArray = attributePair.trim.split(gtfFmt_attributeBreak,2);
            if(pairArray.length != 2) fmtError("Attribute: \""+attributePair+"\" is misformatted! Does not contain the break value: \"" + gtfFmt_attributeBreak + "\" pairString.split.length=" + pairArray.length);
          
            val (key,value) = (pairArray(0).trim,internalUtils.stdUtils.cleanQuotes(pairArray(1).trim));
            mymap = mymap + ((key,value));
          }
          mymap;
     }
     def attributeMap : Map[String,String] = lz_attributeMap;
     
     
     def getAttribute(key : String) : Option[String] = attributeMap.get(key);
     def getAttributeOrDie(key : String) : String = {
          attributeMap.get(key) match {
            case Some(value) => value;
            case None => fmtError("Attribute " + key + " not found!"); null;
          }
     }
   }
   
   trait StdGtfLine extends GtfLine {
     def isExon : Boolean = featureType == codes.STD_EXON_TYPE_CODE;
     def isCDS : Boolean  = featureType == codes.STD_CDS_TYPE_CODE;
     def getTxID : String = getAttributeOrDie(codes.STD_TX_ID_ATTRIBUTE_KEY);
     def getGeneID : String = getAttributeOrDie(codes.GENE_ID_ATTRIBUTE_KEY);
   }
   
   trait FlatGtfLine extends GtfLine {
     def isSpliceJunction : Boolean = isNovelSJ || isKnownSJ;
     def isNovelSJ : Boolean = featureType == codes.JS_FEATURETYPE_NOVELSPLICE;
     def isKnownSJ : Boolean = featureType == codes.JS_FEATURETYPE_KNOWNSPLICE;
     def isExonicPart : Boolean = featureType == codes.JS_FEATURETYPE_EXON;
     def isAggregateGene : Boolean = featureType == codes.JS_FEATURETYPE_GENE;
     
     //def getName : String = 
     def getFeatureCode : String = codes.JS_FEATURETYPE_CODEMAP(featureType);
     def getFeatureAggregateGene : String = getAttributeOrDie(codes.GENE_ID_ATTRIBUTE_KEY);
     def getFeaturePartNumber : String = getAttributeOrDie(codes.JS_EXONIC_PART_NUMBER_ATTRIBUTE_KEY);
     def getFeatureName : String = getFeatureAggregateGene + ":" + getFeatureCode + getFeaturePartNumber;
     def getTxSet : String = getAttributeOrDie(codes.JS_TX_ID_ATTRIBUTE_KEY);
     
     def getDexSeqStr : Option[String] = {
       if(isAggregateGene){
         val dex_attrMap = Map[String,String](FlatOutputGtfLine.DEXSEQ_GENE_ID_ATTRIBUTE_KEY -> getFeatureAggregateGene);
         val dex_attr = dex_attrMap.map((kv) => kv._1 + gtfFmt_attributeBreak + kv._2).mkString("; ")
         
         val out = new Array[String](9);
         out(0) = chromName;
         out(1) = "QoRT_forDEXSeq";
         out(2) =  FlatOutputGtfLine.DEXSEQ_FEATURETYPE_GENE;
         out(3) = start.toString;
         out(4) = end.toString;
         out(5) = score.toString;
         out(6) = strand.toString;
         out(7) = ".";
         out(8) = dex_attr;
         
         return(Some(out.mkString("	")));
       } else if(isExonicPart){
         val dex_attrMap = Map[String,String](FlatOutputGtfLine.DEXSEQ_GENE_ID_ATTRIBUTE_KEY -> getFeatureAggregateGene, 
                                              FlatOutputGtfLine.DEXSEQ_TX_ID_ATTRIBUTE_KEY -> getTxSet, 
                                              FlatOutputGtfLine.DEXSEQ_EXONIC_PART_NUMBER_ATTRIBUTE_KEY -> getFeaturePartNumber
                                              );
         val dex_attr = dex_attrMap.map((kv) => kv._1 + " " + kv._2).mkString("; ")
         
         val out = new Array[String](9);
         out(0) = chromName;
         out(1) = "QoRT_forDEXSeq";
         out(2) = FlatOutputGtfLine.DEXSEQ_FEATURETYPE_EXON;
         out(3) = start.toString;
         out(4) = end.toString;
         out(5) = score.toString;
         out(6) = strand.toString;
         out(7) = ".";
         out(8) = dex_attr;
         
         return(Some(out.mkString("	")));
       } else return None;
     }
   }
    
   class StdInputGtfLine(in_str : String, in_stranded : Boolean, in_gtfFmt_attributeBreak : String , in_codes : GtfCodes = new GtfCodes()) extends InputGtfLine(in_str : String, in_stranded : Boolean, in_gtfFmt_attributeBreak : String , in_codes : GtfCodes) with StdGtfLine;
   class FlatInputGtfLine(in_str : String, in_stranded : Boolean, in_gtfFmt_attributeBreak : String , in_codes : GtfCodes = new GtfCodes()) extends InputGtfLine(in_str : String, in_stranded : Boolean, in_gtfFmt_attributeBreak : String , in_codes : GtfCodes) with FlatGtfLine;
   
   class StdOutputGtfLine(in_chromName : String, in_featureSource : String, in_featureType : String, in_start : Int, in_end : Int, in_score : String, in_strand : Char, in_attributeMap : Map[String,String], in_gtfFmt_attributeBreak : String, in_stranded : Boolean, attribute_sorting : Option[List[String]] = None, in_codes : GtfCodes = new GtfCodes()) extends OutputGtfLine(in_chromName : String, in_featureSource : String, in_featureType : String, in_start : Int, in_end : Int, in_score : String, in_strand : Char, in_attributeMap : Map[String,String], in_gtfFmt_attributeBreak : String, in_stranded : Boolean, attribute_sorting : Option[List[String]], in_codes : GtfCodes) with StdGtfLine;
   class FlatOutputGtfLine(in_chromName : String, in_featureSource : String, in_featureType : String, in_start : Int, in_end : Int, in_score : String, in_strand : Char, in_attributeMap : Map[String,String], in_gtfFmt_attributeBreak : String, in_stranded : Boolean, attribute_sorting : Option[List[String]] = None, in_codes : GtfCodes = new GtfCodes()) extends OutputGtfLine(in_chromName : String, in_featureSource : String, in_featureType : String, in_start : Int, in_end : Int, in_score : String, in_strand : Char, in_attributeMap : Map[String,String], in_gtfFmt_attributeBreak : String, in_stranded : Boolean, attribute_sorting : Option[List[String]], in_codes : GtfCodes) with FlatGtfLine;
   
   /*
   case class GtfLine(str : String, stranded : Boolean, gtfFmt_attributeBreak : String){
        def fmtError(message : String){
          error("FATAL ERROR: readChromsFromFile_gtfStyle checkFmt: " +message+ "\n   Offending Line is:\n \""+str+"\"");
        } 
 
        
        lazy val cells : Array[String] = {
          val c = str.split("	");
          if(checkFmt && c.length < 9) fmtError("Line has too few columns!");
          c;
        }
        lazy val chromName : String = internalUtils.stdUtils.cleanQuotes(cells(0));
        lazy val featureType : String = internalUtils.stdUtils.cleanQuotes(cells(2));
        lazy val start : Int = internalUtils.stdUtils.string2int(cells(3));
        lazy val end : Int = internalUtils.stdUtils.string2int(cells(4));
        lazy val score : String = internalUtils.stdUtils.cleanQuotes(cells(5));
        lazy val strand : Char = if(stranded) cells(6).charAt(0) else '.';
        lazy val attr : String = cells(8).trim;

        lazy val attributeArray : Array[String] = attr.split(";");
        lazy val attributeMap : Map[String,String] = {
          var mymap = Map[String,String]();
          for(attributePair <- attributeArray){
            val pairArray = attributePair.trim.split(gtfFmt_attributeBreak,2);
            if(pairArray.length != 2) fmtError("Attribute: \""+attributePair+"\" is misformatted! Does not contain the break value: \"" + gtfFmt_attributeBreak + "\" pairString.split.length=" + pairArray.length);
          
            val (key,value) = (pairArray(0).trim,internalUtils.stdUtils.cleanQuotes(pairArray(1).trim));
            mymap = mymap + ((key,value));
          }
          mymap;
        }
        
        def getAttribute(key : String) : Option[String] = {
          attributeMap.get(key);
        }
        def getAttributeOrDie(key : String) : String = {
          attributeMap.get(key) match {
            case Some(value) => value;
            case None => fmtError("Attribute " + key + " not found!"); null;
          }
        }
    }*/

  class GtfReader(lines : Iterator[String], stranded : Boolean, chkFmt : Boolean, fmt_attributeBreak : String) extends Iterator[GtfLine] { 
    def hasNext = lines.hasNext;
    def next = new InputGtfLine(lines.next, stranded, fmt_attributeBreak);
    
    //checkFmt = chkFmt;
    //gtfFmt_attributeBreak = fmt_attributeBreak;
  }
  
  class StdGtfReader(lines : Iterator[String], stranded : Boolean, chkFmt : Boolean, fmt_attributeBreak : String, codes : GtfCodes = new GtfCodes()) extends Iterator[StdGtfLine] {
    def hasNext = lines.hasNext;
    def next = new StdInputGtfLine(lines.next, stranded, fmt_attributeBreak, codes);
  }
  class FlatGtfReader(lines : Iterator[String], stranded : Boolean, chkFmt : Boolean, fmt_attributeBreak : String, codes : GtfCodes = new GtfCodes()) extends Iterator[FlatGtfLine] {
    def hasNext = lines.hasNext;
    def next = new FlatInputGtfLine(lines.next, stranded, fmt_attributeBreak, codes);
  }
  
  object GtfReader {
    def getGtfReader(infile : String, stranded : Boolean, chkFmt : Boolean, fmt_attributeBreak : String) : GtfReader = {
      new GtfReader(skipCommentedInitialLines(internalUtils.fileUtils.getLinesSmartUnzip(infile)), stranded , chkFmt, fmt_attributeBreak);
    }
    
    def getFlatGtfReader(infile : String, stranded : Boolean, chkFmt : Boolean, fmt_attributeBreak : String, codes : GtfCodes = new GtfCodes()) : Iterator[FlatGtfLine] = {
      new FlatGtfReader(skipCommentedInitialLines(internalUtils.fileUtils.getLinesSmartUnzip(infile)), stranded, chkFmt, fmt_attributeBreak, codes);
    }
    def getStdGtfReader(infile : String, stranded : Boolean, chkFmt : Boolean, fmt_attributeBreak : String, codes : GtfCodes = new GtfCodes()) : Iterator[StdGtfLine] = {
      new StdGtfReader(skipCommentedInitialLines(internalUtils.fileUtils.getLinesSmartUnzip(infile)), stranded, chkFmt, fmt_attributeBreak, codes);
    }
  }
  
  def writeGtfLine(gtfLine : GtfLine, writer : WriterUtil){
    writer.write(gtfLine.str + "\n");
  }
  
}