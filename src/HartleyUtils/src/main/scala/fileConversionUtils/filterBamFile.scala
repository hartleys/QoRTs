package fileConversionUtils


import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.File;
import internalUtils.stdUtils._;
import internalUtils.Reporter._;
import internalUtils.commandLineUI._;

import scala.collection.JavaConverters._;

import internalUtils.genomicUtils._;
import qcUtils.QCUtility;

import net.sf.samtools._

import internalUtils.commonSeqUtils._;
import internalUtils.genomicAnnoUtils._;
import internalUtils.GtfTool._;
//import scala.collection.JavaConversions._;
import internalUtils.optionHolder._;

import internalUtils.fileUtils._;
import java.util.regex.Pattern;

object filterBamFile {

  class cmdFilterBam extends CommandLineRunUtil {
     val parser : CommandLineArgParser = 
       new CommandLineArgParser(
          command = "filterBamFile", 
          quickSynopsis = "", 
          synopsis = "", 
          description = "This simple utility filters a BAM file."+
                        ""+
                        "WARNING: THIS SUB-UTILITY IS BETA! NOT FOR GENERAL USE!"+
                        ""+
                        ""+
                        ""+
                        "",   
          argList = 
                    new BinaryOptionArgument[String](
                                         name = "dropSequenceRegex", 
                                         arg = List("--dropSequenceRegex"), 
                                         valueName = "regex",  
                                         argDesc = ""
                                        ) ::
                    new UnaryArgument(   name = "singleEnded", 
                                         arg = List("--singleEnded"), // name of value
                                         argDesc = "Flag for single-end data. Note that many other options do not apply in this case (for example: option --countPairsTogether does nothing in single-end mode)" 
                                       ) ::
                    new UnaryArgument(   name = "dropIndels", 
                                         arg = List("--dropIndels"), // name of value
                                         argDesc = "Drops reads that contain indels." 
                                       ) ::
                    new BinaryArgument[Int](   name = "homopolymerOvrMismatch", 
                                         arg = List("--homopolymerOvrMismatch"), // name of value
                                         valueName = "num", 
                                         argDesc = "",
                                         defaultValue = Some(-1)
                                       ) ::
                    new BinaryArgument[Int](name = "minMAPQ",
                                           arg = List("--minMAPQ"),  
                                           valueName = "num", 
                                           argDesc = "Filter reads for the given minimum MAPQ. Set to 0 to turn off mapq filtering.", 
                                           defaultValue = Some(255)
                                           ) :: 

                    new BinaryOptionArgument[String](
                                         name = "keepSequenceRegex", 
                                         arg = List("--keepSequenceRegex"), 
                                         valueName = "regex",  
                                         argDesc = ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "refMatchFile", 
                                         arg = List("--refMatchFile"), 
                                         valueName = "refMatchFile.txt.gz",  
                                         argDesc = ""
                                        ) ::
                    new BinaryOptionArgument[String](
                                         name = "ovrMatchFile", 
                                         arg = List("--ovrMatchFile"), 
                                         valueName = "ovrMatchFile.txt.gz",  
                                         argDesc = ""
                                        ) ::
                    new BinaryOptionArgument[List[String]](
                                         name = "genomeFa", 
                                         arg = List("--genomeFa"), 
                                         valueName = "genomeFa.fa.gz",  
                                         argDesc = ""
                                        ) ::
                    new FinalArgument[String](
                                         name = "infile",
                                         valueName = "infile.bam",
                                         argDesc = "The input BAM file." // description
                                        ) ::
                    new FinalArgument[String](
                                         name = "outfile",
                                         valueName = "outfile.bam",
                                         argDesc = "The output BAM file." // description
                                        ) :: internalUtils.commandLineUI.CLUI_UNIVERSAL_ARGS );
      
     def run(args : Array[String]) {
       val out = parser.parseArguments(args.toList.tail);
      
       if(out){
         filterBamFile.run(
             infile = parser.get[String]("infile"),
             outfile = parser.get[String]("outfile"),
             isSingleEnd = parser.get[Boolean]("singleEnded"),
             minMAPQ = parser.get[Int]("minMAPQ"),
             keepSequenceRegex = parser.get[Option[String]]("keepSequenceRegex"),
             dropSequenceRegex = parser.get[Option[String]]("dropSequenceRegex"),
             refMatchFile = parser.get[Option[String]]("refMatchFile"),
             ovrMatchFile = parser.get[Option[String]]("ovrMatchFile"),
             genomeFa =  parser.get[Option[List[String]]]("genomeFa"),
             homopolymerOvrMismatch = parser.get[Int]("homopolymerOvrMismatch"),
             dropIndels = parser.get[Boolean]("dropIndels")
           );
         }
     }
   }
  
 // def test() {
  def run(infile : String, outfile : String, isSingleEnd : Boolean, minMAPQ : Int, keepSequenceRegex : Option[String], dropSequenceRegex : Option[String],
            refMatchFile : Option[String] = None,
            ovrMatchFile : Option[String] = None,
            genomeFa : Option[List[String]] = None,
            homopolymerOvrMismatch : Int = -1,
            dropIndels : Boolean = false){
    
     val reader : SAMFileReader = if(infile == "-"){
       new SAMFileReader(System.in);
     } else {
       new SAMFileReader(new File(infile));
     }
     val header = reader.getFileHeader()
     val qcu = new QcBamFilter(outputBamFile = outfile, isSingleEnd= isSingleEnd, keepSequenceRegex=keepSequenceRegex,dropSequenceRegex=dropSequenceRegex, header=header,
                                refMatchFile = refMatchFile, ovrMatchFile=ovrMatchFile, genomeFa=genomeFa,homopolymerOvrMismatch=homopolymerOvrMismatch,
                                dropIndels = dropIndels);
     runOnFile(
                   reader = reader, 
                   qcBTW =qcu,
                   testRun  =false, 
                   keepMultiMapped = false, 
                   readGroup  = None,
                   minMAPQ = minMAPQ,
                   isSingleEnd = isSingleEnd,resort = genomeFa.isDefined,
                   unsorted = true,checkForAlignmentBlocks = true);
    
  }
  
  
  class QcBamFilter(outputBamFile : String, isSingleEnd : Boolean, 
                    keepSequenceRegex : Option[String], dropSequenceRegex : Option[String], header : SAMFileHeader,
                    refMatchFile : Option[String] = None,
                    ovrMatchFile : Option[String] = None,
                    genomeFa : Option[List[String]] = None,
                    homopolymerOvrMismatch : Int = -1,
                    dropIndels : Boolean) extends QCUtility[Unit] {
    
    val genomeSeq = if(genomeFa.isEmpty) null else buildEfficientGenomeSeqContainer(genomeFa.get);
    val bufferSize = 10000;
    val refOverhangBuffer = 15;
    
    val keepRegex : Option[Pattern] = keepSequenceRegex match {
      case Some(rsr) => {
        Some(Pattern.compile(rsr));
      }
      case None => None;
    }
    val dropRegex : Option[Pattern] = dropSequenceRegex match {
      case Some(rsr) => {
        Some(Pattern.compile(rsr));
      }
      case None => None;
    }
    val bamWriter = openWriterSmart(outputBamFile, allowStdout = true);
    val matchWriter = refMatchFile match {
      case Some(f) => Some(openWriterSmart(f));
      case None => None;
    }
    val ovrWriter = ovrMatchFile match {
      case Some(f) => Some(openWriterSmart(f));
      case None => None;
    }
    
    bamWriter.write(header.getTextHeader());
    
    val ovrHomoPolyCountMap = scala.collection.mutable.AnyRefMap[String,Int]().withDefault{s => 0 };
    
    def runOnReadPair(r1 : SAMRecord, r2 : SAMRecord, readNum : Int){
      
      if(homopolymerOvrMismatch != -1 && readNum % 100000 == 0){
        printOvrHomoPolyCountMap
      }
      
      if(isSingleEnd){
        var keepRead = true;
        keepRegex match {
          case Some(rx) => {
            if(! (rx.matcher(r1.getReadString()).matches())){
              keepRead = false;
            }
          }
          case None => {
            //do nothing!
          }
        }
        dropRegex match {
          case Some(rx) => {
            if(rx.matcher(r1.getReadString()).matches()){
              keepRead = false;
            }
          }
          case None => {
            //do nothing!
          }
        }
        
        if(keepRead){
          bamWriter.write(r1.getSAMString());
          bamWriter.flush();
        }
      } else {
        var keepRead = true;
        keepRegex match {
          case Some(rx) => {
            if(! (rx.matcher(r1.getReadString()).matches() || rx.matcher(r2.getReadString()).matches())){
              keepRead = false;
            }
          }
          case None => {
            //do nothing!
          }
        }
        dropRegex match {
          case Some(rx) => {
            if(rx.matcher(r1.getReadString()).matches() || rx.matcher(r2.getReadString()).matches()){
              keepRead = false;
            }
          }
          case None => {
            //do nothing!
          }
        }
        
        val cb1 = getCigarBlocksFromRead(r1).toVector;
        val cb2 = getCigarBlocksFromRead(r2).toVector;
        
        if(dropIndels){
          if(cb1.exists{ cb => cb.op == CigarOperator.DELETION || cb.op == CigarOperator.INSERTION || cb.op == CigarOperator.SKIPPED_REGION}){
            keepRead = false;
          }
        }
        
        if(keepRead && refMatchFile.isDefined){

          val minStart = math.min(r1.getAlignmentStart()-1,r2.getAlignmentStart()-1);
          val chrom = r1.getReferenceName()
          if( readNum % 1000000 == 1) genomeSeq.reportBufferStatus;
          genomeSeq.shiftBufferTo(chrom,math.max(minStart-bufferSize,0));
          val cb1f = cb1.filter(c => c.op.consumesReadBases() && c.op.consumesReferenceBases())
          val cb2f = cb2.filter(c => c.op.consumesReadBases() && c.op.consumesReferenceBases())
          val cb1sc = (getLeadClipping(r1),getTailClipping(r1))
          val cb2sc = (getLeadClipping(r2),getTailClipping(r2))
          
          
          val cb1del = cb1.filter(c => (! c.op.consumesReadBases()) && c.op.consumesReferenceBases())
          val cb2del = cb2.filter(c => (! c.op.consumesReadBases()) && c.op.consumesReferenceBases())

          val titleString1raw = r1.getReadName()+" "+r1.getReferenceName()+":"+r1.getAlignmentStart()+" "+r1.getCigarString()+(if(r1.getReadNegativeStrandFlag()){"-"}else{"+"})
          val titleString2raw = repString(" ",r1.getReadName().length())+" "+r2.getReferenceName()+":"+r2.getAlignmentStart()+" "+r2.getCigarString()+(if(r2.getReadNegativeStrandFlag()){"-"}else{"+"})
          val titleLen = math.max(titleString1raw.length(),titleString2raw.length());
          val titleString1 = titleString1raw + repString(" ",1 + titleLen - titleString1raw.length())
          val titleString2 = titleString2raw + repString(" ",1 + titleLen - titleString2raw.length())
          val bufferString = repString(" ",1 + titleLen)
          
          val headBufferString = repString(" ",refOverhangBuffer)
          
          for( (r,ts,cb,cbdel,(leadClip,tailClip)) <- Seq((r1,titleString1,cb1f,cb1del,cb1sc),(r2,titleString2,cb2f,cb2del,cb2sc))){
            val matchMap = cb.flatMap{ c => {
              (c.readStart until (c.readStart + c.len) ).zip(
              genomeSeq.getSeqForInterval(chrom,c.refStart,c.refEnd).toVector.zip(
                  getSeqStringFromBlock(r,c).toVector
              ))
            }}.map{ case (pos,(refchar,readchar)) => {
              if(readchar != 'N' && refchar != readchar){
                (pos,refchar)
              } else {
                (pos,' ')
              }
            }}.toMap.withDefault((p) => '*')
            val mismatchString = Range(0,r.getReadLength()).map{ p => matchMap(p)}.mkString("")
            val offsetMM = Range(0,r.getReadLength()).flatMap{ readPos => {
              val idx = cbdel.indexWhere{ c => c.readStart == readPos }
              if(idx != -1){
                genomeSeq.getSeqForInterval(chrom,cbdel(idx).refStart,cbdel(idx).refEnd).toVector
              } else {
                if( readPos < leadClip){
                  genomeSeq.getSeqForInterval(chrom,r.getAlignmentStart() - (leadClip - readPos) - 1,r.getAlignmentStart()  - (leadClip - readPos)).toVector
                } else if(readPos >= r.getReadLength() - tailClip) {
                  genomeSeq.getSeqForInterval(chrom,r.getAlignmentEnd() + (readPos - (r.getReadLength() - tailClip)),r.getAlignmentEnd() + (readPos - (r.getReadLength() - tailClip)) + 1).toVector
                } else {
                  Some(matchMap(readPos));
                }
              }
            }}.mkString("")
            val offsetRead = Range(0,r.getReadLength()).flatMap{ readPos => {
              val idx = cbdel.indexWhere{ c => c.readStart == readPos }
              if(idx != -1){
                repString("-",cbdel(idx).len).toVector;
              } else {
                if( readPos < leadClip || readPos >= r.getReadLength() - tailClip ){
                  Some(r.getReadString().charAt(readPos).toLower);
                } else {
                  Some(r.getReadString().charAt(readPos));
                }
              }
            }}.mkString("")
            
            //matchWriter.get.write(ts+r.getReadString()+"\n");
            //matchWriter.get.write(bufferString + mismatchString+"\n");
            matchWriter.get.write(ts +headBufferString+ offsetRead+headBufferString+"\n");
            matchWriter.get.write(bufferString +genomeSeq.getSeqForInterval(chrom,r.getUnclippedStart()-refOverhangBuffer-1,r.getUnclippedStart()-1)+ 
                                                offsetMM+
                                                genomeSeq.getSeqForInterval(chrom,r.getUnclippedEnd(),r.getUnclippedEnd()+refOverhangBuffer)+"\n");
          }
        }
        if(keepRead && ovrMatchFile.isDefined){
          if(r1.getAlignmentStart() < r2.getAlignmentEnd() && r2.getAlignmentStart() < r1.getAlignmentEnd()){
            val minStart = math.min(r1.getAlignmentStart()-1,r2.getAlignmentStart()-1);
            val chrom = r1.getReferenceName()
            if( readNum % 1000000 == 1) genomeSeq.reportBufferStatus;
            genomeSeq.shiftBufferTo(chrom,math.max(minStart-bufferSize,0));
            val cb1f = cb1.filter(c => c.op.consumesReadBases() && c.op.consumesReferenceBases())
            val cb2f = cb2.filter(c => c.op.consumesReadBases() && c.op.consumesReferenceBases())
            val cb1sc = (getLeadClipping(r1),getTailClipping(r1))
            val cb2sc = (getLeadClipping(r2),getTailClipping(r2))
            
            val cb1del = cb1.filter(c => (! c.op.consumesReadBases()) && c.op.consumesReferenceBases())
            val cb2del = cb2.filter(c => (! c.op.consumesReadBases()) && c.op.consumesReferenceBases())
  
            val titleString1raw = (if(r1.getReadNegativeStrandFlag()){repString(" ",r1.getReadName().length())}else{r1.getReadName()})+"/1 "+r1.getReferenceName()+":"+r1.getAlignmentStart()+" "+r1.getCigarString()+(if(r1.getReadNegativeStrandFlag()){"-"}else{"+"})
            val titleString2raw = (if(r2.getReadNegativeStrandFlag()){repString(" ",r1.getReadName().length())}else{r1.getReadName()})+"/2 "+r2.getReferenceName()+":"+r2.getAlignmentStart()+" "+r2.getCigarString()+(if(r2.getReadNegativeStrandFlag()){"-"}else{"+"})
            val titleLen = math.max(titleString1raw.length(),titleString2raw.length());
            val titleString1 = titleString1raw + repString(" ",1 + titleLen - titleString1raw.length())
            val titleString2 = titleString2raw + repString(" ",1 + titleLen - titleString2raw.length())
            val bufferString = repString(" ",1 + titleLen)
            
            val headBufferString = repString(" ",refOverhangBuffer)
            
            val r1Offset = if(r1.getUnclippedStart() < r2.getUnclippedStart()) 0 else r1.getUnclippedStart() - r2.getUnclippedStart();
            val r2Offset = if(r2.getUnclippedStart() < r1.getUnclippedStart()) 0 else r2.getUnclippedStart() - r1.getUnclippedStart();
            val rdata = Seq((r1,titleString1,cb1f,cb1del,cb1sc,r1Offset),(r2,titleString2,cb2f,cb2del,cb2sc,r2Offset)).sortBy(x => {x._1.getReadNegativeStrandFlag()});
            //for( (r,ts,cb,cbdel,(leadClip,tailClip), offset) <- Seq((r1,titleString1,cb1f,cb1del,cb1sc,r1Offset),(r2,titleString2,cb2f,cb2del,cb2sc,r2Offset)).sortBy(x => {x._1.getReadNegativeStrandFlag()})){
            val outStrings = rdata.map{ case (r,ts,cb,cbdel,(leadClip,tailClip), offset) => {
              val matchMap = cb.flatMap{ c => {
                (c.readStart until (c.readStart + c.len) ).zip(
                genomeSeq.getSeqForInterval(chrom,c.refStart,c.refEnd).toVector.zip(
                    getSeqStringFromBlock(r,c).toVector
                ))
              }}.map{ case (pos,(refchar,readchar)) => {
                if(readchar != 'N' && refchar != readchar){
                  (pos,refchar)
                } else {
                  (pos,' ')
                }
              }}.toMap.withDefault((p) => '*')
              val mismatchString = Range(0,r.getReadLength()).map{ p => matchMap(p)}.mkString("")
              val offsetMM = Range(0,r.getReadLength()).flatMap{ readPos => {
                val idx = cbdel.indexWhere{ c => c.readStart == readPos }
                if(idx != -1){
                  genomeSeq.getSeqForInterval(chrom,cbdel(idx).refStart,cbdel(idx).refEnd).toVector
                } else {
                  if( readPos < leadClip){
                    genomeSeq.getSeqForInterval(chrom,r.getAlignmentStart() - (leadClip - readPos) - 1,r.getAlignmentStart()  - (leadClip - readPos)).toVector
                  } else if(readPos >= r.getReadLength() - tailClip) {
                    genomeSeq.getSeqForInterval(chrom,r.getAlignmentEnd() + (readPos - (r.getReadLength() - tailClip)),r.getAlignmentEnd() + (readPos - (r.getReadLength() - tailClip)) + 1).toVector
                  } else {
                    Some(matchMap(readPos));
                  }
                }
              }}.mkString("")
              val offsetRead = Range(0,r.getReadLength()).flatMap{ readPos => {
                val idx = cbdel.indexWhere{ c => c.readStart == readPos }
                if(idx != -1){
                  repString("-",cbdel(idx).len).toVector;
                } else {
                  if( readPos < leadClip || readPos >= r.getReadLength() - tailClip ){
                    Some(r.getReadString().charAt(readPos).toLower);
                  } else {
                    Some(r.getReadString().charAt(readPos));
                  }
                }
              }}.mkString("")
              
              
              
              //ovrWriter.get.write(ts+r.getReadString()+"\n");
              //ovrWriter.get.write(bufferString + mismatchString+"\n");
              /*ovrWriter.get.write(ts +repString(" ",offset)+headBufferString+ offsetRead+headBufferString+"\n");
              ovrWriter.get.write(bufferString +repString(" ",offset)+genomeSeq.getSeqForInterval(chrom,r.getUnclippedStart()-refOverhangBuffer-1,r.getUnclippedStart()-1)+ 
                                                  offsetMM+
                                                  genomeSeq.getSeqForInterval(chrom,r.getUnclippedEnd(),r.getUnclippedEnd()+refOverhangBuffer)+"\n");
              */
              (ts,repString(" ",offset)+headBufferString+ offsetRead+headBufferString,
               bufferString,repString(" ",offset)+genomeSeq.getSeqForInterval(chrom,r.getUnclippedStart()-refOverhangBuffer-1,r.getUnclippedStart()-1)+ 
                                                  offsetMM+
                                                  genomeSeq.getSeqForInterval(chrom,r.getUnclippedEnd(),r.getUnclippedEnd()+refOverhangBuffer))
            }}
            if(homopolymerOvrMismatch == -1){
              outStrings.foreach{ case (ts, readSeq, buf, refSeq) => {
                ovrWriter.get.write(ts+readSeq+"\n");
                ovrWriter.get.write(buf+refSeq+"\n");
              }}
            } else {
              val r1seqraw = outStrings(0)._2
              val r2seqraw = outStrings(1)._2
              val len = math.max(r1seqraw.length,r2seqraw.length);
              val r1seq = r1seqraw + repString(" ",len-r1seqraw.length);
              val r2seq = r2seqraw + repString(" ",len-r2seqraw.length);
              val hasOverlappedHomoPoly1 = r1seq.indices.dropRight(homopolymerOvrMismatch).indexWhere{ i => {
                (r1seq.substring(i,i+homopolymerOvrMismatch).forall( c => c == r1seq.charAt(i) && c.isUpper && c != 'N') && r2seq.substring(i,i+homopolymerOvrMismatch).forall(c => (c != ' ' && c.isUpper)))
              }}
              val hasOverlappedHomoPoly2 = r1seq.indices.dropRight(homopolymerOvrMismatch).indexWhere{ i => {
                (r2seq.substring(i,i+homopolymerOvrMismatch).forall( c => c == r2seq.charAt(i) && c.isUpper && c != 'N') && r1seq.substring(i,i+homopolymerOvrMismatch).forall(c => (c != ' ' && c.isUpper)))
              }}
              
              if(hasOverlappedHomoPoly1 != -1 || hasOverlappedHomoPoly2 != -1){
                val hasOverlappedHomoPolyMis1 = r1seq.indices.dropRight(homopolymerOvrMismatch).indexWhere{ i => {
                  (r1seq.substring(i,i+homopolymerOvrMismatch).forall( c => c == r1seq.charAt(i) && c.isUpper && c != 'N') && r2seq.substring(i,i+homopolymerOvrMismatch).forall(c => (c != ' ' && c.isUpper)) && r1seq.substring(i,i+homopolymerOvrMismatch) != r2seq.substring(i,i+homopolymerOvrMismatch))
                }}
                val hasOverlappedHomoPolyMis2 = r1seq.indices.dropRight(homopolymerOvrMismatch).indexWhere{ i => {
                  (r2seq.substring(i,i+homopolymerOvrMismatch).forall( c => c == r2seq.charAt(i) && c.isUpper && c != 'N') && r1seq.substring(i,i+homopolymerOvrMismatch).forall(c => (c != ' ' && c.isUpper)) && r1seq.substring(i,i+homopolymerOvrMismatch) != r2seq.substring(i,i+homopolymerOvrMismatch))
                }}
                outStrings.foreach{ case (ts, readSeq, buf, refSeq) => {
                  ovrWriter.get.write(ts+readSeq+"\n");
                  ovrWriter.get.write(buf+refSeq+"\n");
                }}
                
                val summaryString = if(hasOverlappedHomoPolyMis1 != -1 && hasOverlappedHomoPolyMis2 != -1){
                  "MISMATCH:BOTH/r"+(if(r1.getReadNegativeStrandFlag()) "2" else "1")+"+"+"."+r1seq.charAt(hasOverlappedHomoPolyMis1)+"/"+(if(r1.getReadNegativeStrandFlag()) "1" else "2")+"-"+"."+r2seq.charAt(hasOverlappedHomoPolyMis2)
                } else if(hasOverlappedHomoPolyMis1 != -1){
                  "MISMATCH:r"+(if(r1.getReadNegativeStrandFlag()) "2" else "1")+"+"+"."+r1seq.charAt(hasOverlappedHomoPolyMis1)
                  //ovrWriter.get.write("r1"+(if(r1.getReadNegativeStrandFlag()) "-" else "+")+":("+hasOverlappedHomoPolyMis1+")("+r1seq.charAt(hasOverlappedHomoPolyMis1)+")"+"\n");
                  //ovrHomoPolyCountMap((1,r1seq.charAt(hasOverlappedHomoPoly1))) = ovrHomoPolyCountMap((1,r1seq.charAt(hasOverlappedHomoPoly1))) + 1;
                } else if(hasOverlappedHomoPolyMis2 != -1){
                  "MISMATCH:r"+(if(r1.getReadNegativeStrandFlag()) "1" else "2")+"-"+"."+r2seq.charAt(hasOverlappedHomoPolyMis2)
                  //ovrWriter.get.write("r2"+(if(r2.getReadNegativeStrandFlag()) "-" else "+")+":("+hasOverlappedHomoPoly2+")("+r1seq.charAt(hasOverlappedHomoPolyMis1)+")"+"\n");
                  //ovrHomoPolyCountMap((1,r1seq.charAt(hasOverlappedHomoPoly1))) = ovrHomoPolyCountMap((1,r1seq.charAt(hasOverlappedHomoPoly1))) + 1;
                } else if(hasOverlappedHomoPolyMis1 != -1){
                  "MATCH:r"+(if(r1.getReadNegativeStrandFlag()) "2" else "1")+"+"+"."+r1seq.charAt(hasOverlappedHomoPoly1)
                } else {
                  "MATCH:r"+(if(r1.getReadNegativeStrandFlag()) "1" else "2")+"-"+"."+r2seq.charAt(hasOverlappedHomoPoly2)
                }
                ovrHomoPolyCountMap(summaryString) = ovrHomoPolyCountMap(summaryString) + 1;
                
                ovrWriter.get.write(summaryString+"\n");
                 
                
              } else {
                keepRead = false;
              }
            }
            
          } else {
            keepRead = false;
          }
        }
        
        if(keepRead){
          if(r1.getAlignmentStart() < r2.getAlignmentStart()){
            bamWriter.write(r1.getSAMString());
            bamWriter.write(r2.getSAMString());
          } else {
            bamWriter.write(r2.getSAMString());
            bamWriter.write(r1.getSAMString());
          }
          bamWriter.flush();
        }
      }
    }
    def getSeqStringFromBlock(r : SAMRecord, cb : CigarBlock) : String = {
      r.getReadString().substring(cb.readStart,cb.readEnd);
    }
    
    def writeOutput(outfile : String, summaryWriter : WriterUtil, docWriter : DocWriterUtil = null){
      bamWriter.close();
      if(refMatchFile.isDefined){
        matchWriter.get.close();
      }
      if(ovrWriter.isDefined){
        ovrWriter.get.close();
        printOvrHomoPolyCountMap;
      }
    }
    
    
    val baseSeq = Seq(Seq("A+","C+","G+","T+"),Seq("T-","G-","C-","A-"));
    
    def printOvrHomoPolyCountMap {
      
      
      
      reportln("MISMATCH:"+ovrHomoPolyCountMap.keys.toVector.sorted.filter{_.startsWith("MIS")}.map{ s => {
                             s.drop(9)+":"+ovrHomoPolyCountMap(s)
                           }}.mkString(",") + "\n"+
              "MATCH:    "+ovrHomoPolyCountMap.keys.toVector.sorted.filter{! _.startsWith("MIS")}.map{ s => {
                             s.drop(6)+":"+ovrHomoPolyCountMap(s)
                           }}.mkString(","),"note");
    }
    def getUtilityName : String = "bamToWig";

  }
  

  def runOnFile(  
                   reader : SAMFileReader, 
                   qcBTW : QCUtility[Unit],
                   testRun  : Boolean, 
                   keepMultiMapped : Boolean = false, 
                   readGroup : Option[String] = None,
                   minMAPQ : Int,
                   isSingleEnd : Boolean,resort : Boolean = false,
                   unsorted : Boolean = true,checkForAlignmentBlocks : Boolean = true){
    
    val peekCt = 2000;

    val (samFileAttributes, recordIter) = initSamRecordIterator(reader, peekCt);
    
    val pairedIter : Iterator[(SAMRecord,SAMRecord)] = 
      if(isSingleEnd){
        if(testRun) samRecordPairIterator_withMulti_singleEnd(recordIter, true, 200000) else samRecordPairIterator_withMulti_singleEnd(recordIter);
      } else {
        if(unsorted){
          if(resort){
            if(testRun) samRecordPairIterator_resorted(recordIter, true, 200000) else samRecordPairIterator_resorted(recordIter)
          } else if(testRun) samRecordPairIterator_unsorted(recordIter, true, 200000) else samRecordPairIterator_unsorted(recordIter)
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
      if( isSortedByPosition & (! unsorted )){ reportln("Warning: Based on the first "+peekCt+" reads, SAM/BAM file appears to be sorted by read position. If this is so, you should probably omit the '--nameSorted' option, as errors may follow.","warn"); }
      if( ((! isSortedByPosition) & ( unsorted ))){ reportln("Note: Sam lines do not appear to be sorted by read position (based on the first "+peekCt+" reads).","note"); }
      //if( ((! isSortedByNameLexicographically) & (! unsorted ))) error("FATAL ERROR: SAM/BAM file is not sorted by name (based on the first "+peekCt+" reads)! Either sort the file by name, or sort by read position and use the \"--coordSorted\" option.");
    }
    
    //reportln("SAMRecord Reader Generated. Read length: "+readLength+".","note");

    val coda : Array[Int] = internalUtils.commonSeqUtils.getNewCauseOfDropArray;
    val coda_options : Array[Boolean] = internalUtils.commonSeqUtils.CODA_DEFAULT_OPTIONS.toArray;
    if(isSingleEnd) CODA_SINGLE_END_OFF_OPTIONS.foreach( coda_options(_) = false );
    if(keepMultiMapped) coda_options(internalUtils.commonSeqUtils.CODA_NOT_UNIQUE_ALIGNMENT) = false;
    if(! readGroup.isEmpty) coda_options(internalUtils.commonSeqUtils.CODA_NOT_MARKED_RG) = true;
    if(checkForAlignmentBlocks) coda_options(internalUtils.commonSeqUtils.CODA_NO_ALN_BLOCKS) = true;

    
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
  
  
}










