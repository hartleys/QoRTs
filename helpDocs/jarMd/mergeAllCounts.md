# QoRTs: Quality Of Rna-seq Tool Set
Version 0.1.12 ([back to index](index.html))

([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "mergeAllCounts"

## USAGE:

    java [Java Options] -jar QoRTs.jar mergeAllCounts [options] infileDir decoderFile outfile


## DESCRIPTION:

This tool uses a replicate decoder to merge count/wiggle data of all techical replicates in a dataset, producing sample\-wise counts\. You must supply a replicate decoder which indicates which replicates are technical replicates of which samples\. This tool will then merges each sample's technical replicates using the "mergeCounts" function\.

## REQUIRED ARGUMENTS:
### infileDir:

> The top-level directory in which all the QC output can be found. This concatenated with the qc.data.dir column must equal the path to the raw QC output directories (String)


### decoderFile:

> The decoder file, which must conform to the requirements of the QoRT decoder specification. In particular it MUST have two specific columns: 
"sample.ID": This utility will merge the count data output from all bamfiles that have the same sample.ID
and
"qc.data.dir" (OR "unique.ID"): This must be the file path to the output data directory, from the infileDir file location. (String)


### outfile:

> The output file directory. (String)



## OPTIONAL ARGUMENTS:
### --sampleID sampid:

> Optional: the id of the specific sample that is to be merged. By default, this utility will merge all samples found in the decoder. With this option selected, it will ONLY merge the one sample named here. (String)

### --additionalTrackOptions "track options":

> More options for the wiggle tracks. For more information refer to the wiggle track definition on the UCSC genome browser website. (String)

### --wiggleWindow val:

> The window size of the alternate-size wiggle track, if applicable. (Int)

### --mergeFiles file1[,file2,...]:

> A comma-delimited list of strings, indicating which file types to attempt to merge. By default, this utility autodetects the presence of all mergable qc files and merges all standard files. Valid codes are:DESeq,DEXSeq,JunctionSeq,NovelSplice,KnownSplice,WiggleTrack,WiggleTrackAltWin (CommaDelimitedListOfStrings)

### --verbose:

> Flag to indicate that debugging information and extra progress information should be sent to stderr. (flag)

### --quiet:

> Flag to indicate that only errors and warnings should be sent to stderr. (flag)

## AUTHORS:

Stephen W\. Hartley, Ph\.D\. <stephen\.hartley@nih\.gov>

## LEGAL:

 This software is "United States Government Work" under the terms of the United States Copyright  Act\.  It was written as part of the authors' official duties for the United States Government and  thus cannot be copyrighted\.  This software is freely available to the public for use without a  copyright notice\.  Restrictions cannot be placed on its present or future use\.  Although all reasonable efforts have been taken to ensure the accuracy and reliability of the  software and data, the National Human Genome Research Institute \(NHGRI\) and the U\.S\. Government  does not and cannot warrant the performance or results that may be obtained by using this software  or data\.  NHGRI and the U\.S\. Government disclaims all warranties as to performance, merchantability  or fitness for any particular purpose\.  In any work or product derived from this material, proper attribution of the authors as the source  of the software or data should be made, using "NHGRI Genome Technology Branch" as the citation\.  NOTE: This package includes \(internally\) the sam\-1\.113\.jar library from picard tools\. That package uses the MIT license, which can be accessed using the command:  java \-jar thisjarfile\.jar help samjdkinfo

