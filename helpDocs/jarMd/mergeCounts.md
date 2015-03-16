# QoRTs: Quality Of Rna-seq Tool Set
> Version 0.2.7 (Updated Mon Mar 16 17:02:31 EDT 2015)

> ([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "mergeCounts"

## USAGE:

    java [Java Options] -jar QoRTs.jar mergeCounts [options] infileDirs outfilePrefix


## DESCRIPTION:

This utility merges count, wiggle, and similar data from multiple QoRTs QC runs\. This is intended for use in merging the data from multiple technical replicates of the same sample/library\.This tool will then merge all count data \(including gene\-level, exon\-level, and known/novel splice\-junction\) counts, as well as wiggle files, assuming all files use the standard naming conventions \(for example, the fwd\-strand wiggle files must be named: "QC\.wiggle\.fwd\.wig\.gz", etc\)\.If any files are missing, they will be skipped\.

## REQUIRED ARGUMENTS:
### infileDirs:

> The replicates' QC output directories (the output directory used with the initial 'QC' step), as a comma-delimited list (no whitespace). (String)


### outfilePrefix:

> The output file prefix (or directory) (String)



## OPTIONAL ARGUMENTS:
### --wiggleWindow val:

> The window size of the alternate-size wiggle track, if applicable. (Int)

### --additionalTrackOptions "track options":

> More options for the wiggle tracks. For more information refer to the wiggle track definition on the UCSC genome browser website. (String)

### --mergeFiles filetype1[,filetype2,...]:

> A comma-delimited list of strings, indicating which file types to attempt to merge. By default, this utility autodetects the presence of all mergable qc files and merges all standard files. Valid codes are: DESeq, DEXSeq, JunctionSeq, NovelSplice, KnownSplice, WiggleTrack, WiggleTrackAltWin (CommaDelimitedListOfStrings)

### --trackTitle options:

> The prefix of the title of the merged wiggle tracks. (String)

### --verbose:

> Flag to indicate that debugging information and extra progress information should be sent to stderr. (flag)

### --quiet:

> Flag to indicate that only errors and warnings should be sent to stderr. (flag)

## AUTHORS:

Stephen W\. Hartley, Ph\.D\. stephen\.hartley \(at nih dot gov\)

## LEGAL:

 This software is "United States Government Work" under the terms of the United States Copyright  Act\.  It was written as part of the authors' official duties for the United States Government and  thus cannot be copyrighted\.  This software is freely available to the public for use without a  copyright notice\.  Restrictions cannot be placed on its present or future use\.  Although all reasonable efforts have been taken to ensure the accuracy and reliability of the  software and data, the National Human Genome Research Institute \(NHGRI\) and the U\.S\. Government  does not and cannot warrant the performance or results that may be obtained by using this software  or data\.  NHGRI and the U\.S\. Government disclaims all warranties as to performance, merchantability  or fitness for any particular purpose\.  In any work or product derived from this material, proper attribution of the authors as the source  of the software or data should be made, using "NHGRI Genome Technology Branch" as the citation\.  NOTE: This package includes \(internally\) the sam\-1\.113\.jar library from picard tools\. That package uses the MIT license, which can be accessed using the command:  java \-jar thisjarfile\.jar help samjdkinfo

