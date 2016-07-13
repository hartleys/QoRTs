package internalUtils

import internalUtils.fileUtils._;
import internalUtils.Reporter._;
//import internalUtils.stdUtils._;

object VcfTool {

   abstract class VcfMetaLine {
     def key : String;
     def value : String;
     override def toString() : String = "##"+key+"="+value;
   }
   case class VcfAnnoLine(k : String, v : String) extends VcfMetaLine {
     def key : String = k;
     def value : String = v;
   }
   case class VcfInfoLine(ID : String, Number : String, Type : String, Description : String, Source:Option[String] = None, Version:Option[String] = None) extends VcfMetaLine {
     def key : String = "INFO";
     def value : String = "<" + 
                          "ID="+ ID+""+
                          ",Number="+Number+""+
                          ",Type="+Type+""+
                          ",Description="+Description+
                          (if(Source.isEmpty) "" else ",Source="+Source.get)+
                          (if(Version.isEmpty) "" else ",Version="+Version.get)+
                          ">";
   }
   case class VcfFilterLine(ID : String, Description : String) extends VcfMetaLine {
     def key : String = "FILTER";
     def value : String = "<"+
                          "ID="+ ID+""+
                          ",Description="+Description+
                          ">";
   }
   case class VcfFormatLine(ID : String, Number : String, Type : String, Description : String) extends VcfMetaLine {
     def key : String = "FORMAT";
     def value : String = "<" + 
                          "ID="+ ID+""+
                          ",Number="+Number+""+
                          ",Type="+Type+""+
                          ",Description="+Description+
                          ">";
   }
   
   private def parseVcfMetadataValue(s : String) : Seq[(String,String)] = {
     if(s.length < 2 || s.head != '<' || s.last != '>'){
       error("FATAL ERROR: Malformed metadata line in VCF file: value is too short or isn't bound by angle-brackets (errCode VcfTool:parseVcfMetadataValue:65)!");
     }
     val valCells = internalUtils.stdUtils.parseTokens(s.init.tail,',');
     val valPairs = valCells.map((v) => {
       val pair = internalUtils.stdUtils.parseTokens(v,'=');
       if(pair.length != 2) error("FATAL ERROR: Malformed metadata line in VCF file (errCode VcfTool:76)!");
       (pair(0),pair(1));
     });
     return valPairs;
   }
   
   def readVcfMetadata(lines : Iterable[String], sampleLine : String) : VcfMetadata = {
     val metaLines : Seq[VcfMetaLine] = lines.map((s : String) => readVcfMetaLine(s)).toSeq;
     val (info : Seq[VcfInfoLine], filter : Seq[VcfFilterLine], format : Seq[VcfFormatLine], anno : Seq[VcfAnnoLine]) = {
       metaLines.foldLeft((Seq[VcfInfoLine](),Seq[VcfFilterLine](),Seq[VcfFormatLine](),Seq[VcfAnnoLine]()))((soFar,curr) => {
         curr match {
           case x : VcfInfoLine => (soFar._1 :+ x, soFar._2, soFar._3, soFar._4);
           case x : VcfFilterLine => (soFar._1, soFar._2 :+ x, soFar._3, soFar._4);
           case x : VcfFormatLine => (soFar._1 , soFar._2, soFar._3 :+ x, soFar._4);
           case x : VcfAnnoLine => (soFar._1, soFar._2, soFar._3, soFar._4 :+ x);
         }
       }) 
     }
     val sampleCells = sampleLine.split("\t");
     if(sampleCells.length < 9) error("FATAL ERROR: Malformed Vcf Header line. Less than 9 columns! (errCode VcfTool:readVcfMetadata:72)\n offending line:\""+sampleLine+"\"");
     val sampleID = sampleCells.slice(9,sampleCells.length);
     
     return VcfMetadata(info, filter, format, anno, sampleID);
   }
   
   def readVcfMetaLine(line : String) : VcfMetaLine = {
     if(line.substring(0,2) != "##"){
       error("FATAL ERROR: Impossible state! readVcfMetaLine has been given a line that does not start with \"##\"! (errcode VcfTool:readVcfMetaLine:78)");
     }
     val cells = line.substring(2).split("=",2);
     if(cells(0) == "INFO"){
       readVcfInfoLine(cells(1));
     } else if(cells(0) == "FILTER"){
       readVcfFilterLine(cells(1));
     } else if(cells(0) == "FORMAT"){
       readVcfFormatLine(cells(1));
     } else {
       VcfAnnoLine(cells(0), cells(1));
     }
   }
   private def readVcfInfoLine(v : String) : VcfInfoLine = {
       val valMap = parseVcfMetadataValue(v).toMap;
       if(! valMap.contains("ID")) error("FATAL ERROR: Malformed INFO metadata line in VCF file: no ID key! (errCode VcfTool:readVcfInfoLine:84)");
       if(! valMap.contains("Number")) error("FATAL ERROR: Malformed INFO metadata line in VCF file: no Number key! (errCode VcfTool:readVcfInfoLine:85)");
       if(! valMap.contains("Type")) error("FATAL ERROR: Malformed INFO metadata line in VCF file: no Type key! (errCode VcfTool:readVcfInfoLine:86)");
       if(! valMap.contains("Description")) error("FATAL ERROR: Malformed INFO metadata line in VCF file: no Description key! (errCode VcfTool:readVcfInfoLine:87)");
       return VcfInfoLine(valMap("ID"), valMap("Number"), valMap("Type"), valMap("Description"), valMap.get("Source"), valMap.get("Version"));
   }
   private def readVcfFilterLine(v : String) : VcfFilterLine = {
       val valMap = parseVcfMetadataValue(v).toMap;
       if(! valMap.contains("ID")) error("FATAL ERROR: Malformed FILTER metadata line in VCF file: no ID key! (errCode VcfTool:readVcfFilterLine:84)");
       if(! valMap.contains("Description")) error("FATAL ERROR: Malformed FILTER metadata line in VCF file: no Description key! (errCode VcfTool:readVcfFilterLine:87)");
       return VcfFilterLine(valMap("ID"), valMap("Description"));
   }
   private def readVcfFormatLine(v : String) : VcfFormatLine = {
       val valMap = parseVcfMetadataValue(v).toMap;
       if(! valMap.contains("ID")) error("FATAL ERROR: Malformed FORMAT metadata line in VCF file: no ID key! (errCode VcfTool:readVcfFormatLine:84)");
       if(! valMap.contains("Number")) error("FATAL ERROR: Malformed FORMAT metadata line in VCF file: no Number key! (errCode VcfTool:readVcfFormatLine:85)");
       if(! valMap.contains("Type")) error("FATAL ERROR: Malformed FORMAT metadata line in VCF file: no Type key! (errCode VcfTool:readVcfFormatLine:86)");
       if(! valMap.contains("Description")) error("FATAL ERROR: Malformed FORMAT metadata line in VCF file: no Description key! (errCode VcfTool:readVcfFormatLine:87)");
       return VcfFormatLine(valMap("ID"), valMap("Number"), valMap("Type"), valMap("Description"));
   }
   
   case class VcfMetadata(info : Seq[VcfInfoLine],filter : Seq[VcfFilterLine], format : Seq[VcfFormatLine], anno : Seq[VcfAnnoLine], sampleID : Seq[String]) {
     lazy val metaLines   : Seq[VcfMetaLine] = info ++ filter ++ format ++ anno;

     lazy val infoMap : Map[String,VcfInfoLine] = info.map((v) =>{
       (v.ID,v);
     }).toMap;
     lazy val filterMap : Map[String,VcfFilterLine] = filter.map((v) =>{
       (v.ID,v);
     }).toMap;
     lazy val formatMap : Map[String,VcfFormatLine] = format.map((v) =>{
       (v.ID,v);
     }).toMap;
     
   }
   
   abstract class VcfLine {
     def CHROM : String;
     def POS : Int;
     def ID : String;
     def REF : String;
     def ALT: String;
     def QUAL : Double;
     def FILTER: String;
     def INFO: String;
     def FORMAT : String;
     def GENOTYPES: Seq[String];
     def metadata : VcfMetadata;
     
     override def toString() : String = CHROM +"\t"+POS+"\t"+ID+"\t"+REF+"\t"+ALT+"\t"+QUAL+"\t"+FILTER+"\t"+INFO+"\t"+FORMAT+"\t"+GENOTYPES.mkString("\t");
     
     def fmt : Seq[String];

     lazy val idSeq : Seq[String] = ID.split(";");
     lazy val altSeq : Seq[String] = ALT.split(",");
     
     lazy val fmtMeta : Seq[VcfFormatLine] = fmt.map((f : String) =>{
       metadata.formatMap.get(f) match {
         case Some(x) => x;
         case None => {
           error("FATAL VCF PARSING ERROR: FORMAT ID not found!\n offending line:\""+toString()+"\"");
           null;
         }
       }
     })
     lazy val fmtInfo : Seq[(String,String)] = fmtMeta.map((m : VcfFormatLine)=> (m.Number,m.Type));
   }
   
   case class InputVcfLine(line : String, meta : VcfMetadata) extends VcfLine {
     def metadata = meta;
     lazy val cells = {
       val c = line.split("\t");
       if(c.length < 9) error("FATAL ERROR: Vcf Line has fewer than 9 columns:\n  Offending line: \""+line+"\"");
       c;
     }
     def CHROM : String = cells(0);
     lazy val pos : Int = internalUtils.stdUtils.string2int(cells(1));
     def POS : Int = pos;
     def ID : String = cells(2);
     def REF : String = cells(3);
     def ALT : String = cells(4);
     lazy val qual : Double = internalUtils.stdUtils.string2double(cells(5));
     def QUAL : Double = qual;
     def FILTER : String = cells(6);
     def INFO : String = cells(7);
     def FORMAT : String = cells(8);
     def GENOTYPES : Seq[String] = cells.slice(9,cells.length);
     
     //def passFilter : Boolean = 
     
     lazy val FMT : Seq[String] = FORMAT.split(":");
     def fmt : Seq[String] = FMT;
     
     lazy val genoTableBySample : Seq[Seq[String]] = GENOTYPES.map((g : String) => {
       val out = g.split(":").toSeq;
       if(out.length != fmt.length) error("FATAL ERROR: Vcf genotype has the wrong number of elements.\n offending line: \""+line+"\"");
       out;
     });
     lazy val genoTable : Seq[Seq[String]] = internalUtils.stdUtils.transposeMatrix(genoTableBySample);
   }
   //case class OutputVcfLine(chrom : String, pos : Int, id : String, ref : String, alt : String, qual : Double, filter : String, Info : String, genotypes : Seq[String]){
     
   //}
   
   /*
   abstract class VcfMetadata {
     def getMetaLines   : Seq[VcfMetaLine];
     def getInfoLines   : Seq[VcfInfoLine];
     def getFilterLines : Seq[VcfFilterLine];
     def getFormatLines : Seq[VcfFormatLine];
     def getInfoMap : Map[String,VcfInfoLine] = getInfoLines.map((v) =>{
       (v.ID,v);
     }).toMap;
     def getFormatMap : Map[String,VcfFormatLine] = getFormatLines.map((v) =>{
       (v.ID,v);
     }).toMap;
     def getFilterMap : Map[String,VcfFilterLine] = getFilterLines.map((v) =>{
       (v.ID,v);
     }).toMap;
   }
   */
   
   abstract class VcfNumberField {
     def isInteger : Boolean;
     def isSpecial : Boolean;
     def isKnown : Boolean;
   }
   case class VcfNumberFieldDot() extends VcfNumberField {
     def isInteger : Boolean = false;
     def isSpecial : Boolean = true;
     def isKnown : Boolean = false;
   }
   abstract class VcfNumberFieldKnown extends VcfNumberField {
     def isKnown : Boolean = true;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int;
   }
   case class VcfNumberFieldInt(value : Int) extends VcfNumberFieldKnown {
     def isInteger : Boolean = true;
     def isSpecial : Boolean = false;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int = value;
   }
   case class VcfNumberFieldA() extends VcfNumberFieldKnown {
     def isInteger : Boolean = false;
     def isSpecial : Boolean = true;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int = altAlleleCt;
   }
   case class VcfNumberFieldR() extends VcfNumberFieldKnown {
     def isInteger : Boolean = false;
     def isSpecial : Boolean = true;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int = altAlleleCt + 1;
   }
   case class VcfNumberFieldG() extends VcfNumberFieldKnown {
     def isInteger : Boolean = false;
     def isSpecial : Boolean = true;
     def getFieldCount(altAlleleCt : Int, genotypeCt : Int) : Int = genotypeCt;
   }
   
   abstract class VcfVal {
     def isInt : Boolean = false;
     def isFloat : Boolean = false;
     def isChar : Boolean = false;
     def isString : Boolean = false;
   }
   case class VcfInt(v : Int) extends VcfVal {
     override def toString() = v.toString();
     override def isInt : Boolean = true;
   }
   case class VcfFloat(v : Double) extends VcfVal {
     override def toString() = v.toString();
     override def isFloat : Boolean = true;
   }
   case class VcfChar(v : Char) extends VcfVal {
     override def toString() = v.toString();
     override def isChar : Boolean = true;
   }
   case class VcfString(v : String) extends VcfVal {
     override def toString() = v.toString();
     override def isString : Boolean = true;
   }
   
   def readVcfVal(v : String, fmt : String) : VcfVal = {
     if(fmt == "Integer"){
       return VcfInt(internalUtils.stdUtils.string2int(v));
     } else if(fmt == "Float"){
       return VcfFloat(internalUtils.stdUtils.string2double(v));
     } else if(fmt == "Char"){
       if(v.length > 1) error("FATAL ERROR: Malformed VCF. Found 'char' formatted field with length > 1!");
       return VcfChar(v.charAt(0));
     } else if(fmt == "String"){
       return VcfString(v);
     } else {
       error("FATAL ERROR: Malformed VCF. Unrecognized format: \""+fmt+"\"");
       return null;
     }
   } 
   
   
   def getVcfReader(lines : Iterator[String]) : (VcfMetadata, Iterator[InputVcfLine]) = {
     val (metaLines, bodyLines) = internalUtils.stdUtils.splitIterator(lines, (ln : String) => {
       ln.startsWith("##");
     });
     if(! bodyLines.hasNext) error("FATAL ERROR: Vcf File does not have header or body lines.");
     val headerLine = bodyLines.next;
     if(! headerLine.startsWith("#")) error("FATAL ERROR: Vcf File header line not present or does not start with \"#\"");
     if(! bodyLines.hasNext) error("FATAL ERROR: Vcf File does not have body lines.");
     
     val meta : VcfMetadata = readVcfMetadata(metaLines,headerLine);
     
     val vcfLines : Iterator[InputVcfLine] = bodyLines.map( (line : String) => {
       InputVcfLine(line,meta);
     });
     
     return ((meta,vcfLines));
   }
   
  
}















