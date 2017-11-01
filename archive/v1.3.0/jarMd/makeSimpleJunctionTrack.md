# QoRTs: Quality Of Rna-seq Tool Set
> Version 1.3.0 (Updated Fri Oct 20 11:56:37 EDT 2017)

> ([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "makeSimpleJunctionTrack"

## USAGE:

    java [Java Options] -jar QoRTs.jar makeSimpleJunctionTrack [options] indir fileType outdir


## DESCRIPTION:

This utility converts QoRTs splice junction count files into bed format\. Unlike makeJunctionTrack and makeOrphanJunctionTrack, this utility is not designed to compile multiple samples or replicates together\. It is a simple converter from a QoRTs junction count file to a bed file\. The count files come in 3 types: known, novel, and orphan\.

## REQUIRED ARGUMENTS:
### indir:

> The data directory in which to find the splice junction count files. (String)


### fileType:

> The type of splice junction counts to compile. Must be one of: [known,novel,orphan]. (String)


### outdir:

> The location to output the bed file. Traditionally the same as indir. (String)



## OPTIONAL ARGUMENTS:
### --rgb r,g,b:

> The rgb color for all the bed file lines. (String)

### --sizeFactor val:

> A double-precision floating-point value. If this option is set, all counts will be divided by the given normalization factor. (Double)

### --includeSpliceNames:

> Flag to indicate that splice names should be used as well as splice counts. (flag)

### --digits num:

> The number of digits after the decimal to include in counts. CURRENTLY NOT IMPLEMENTED! (Int)

### --filterMin val:

> If this option is set, then all bed lines with a count LESS THAN the given value will be dropped. (Double)

### --filterMax val:

> If this option is set, then all bed lines with a count GREATER THAN the given value will be dropped. (Double)

### --noGzip:

> Flag to indicate whether whether input and output data is/will be gzip-compressed. (flag)

### --maxIdentifierLength num:

> The max number of characters a junction ID can have. If longer than this, the middle will be truncated. This is to prevent browser errors when the length of a bed line name is greater than 255 characters. (Int)

### --outfileSuffix file.bed.gz:

> The name of the output file. (String)

### --verbose:

> Flag to indicate that debugging information and extra progress information should be sent to stderr. (flag)

### --quiet:

> Flag to indicate that only errors and warnings should be sent to stderr. (flag)

## AUTHORS:

Stephen W\. Hartley, Ph\.D\. stephen\.hartley \(at nih dot gov\)

## LEGAL:

 This software is "United States Government Work" under the terms of the United States Copyright  Act\.  It was written as part of the authors' official duties for the United States Government and  thus cannot be copyrighted\.  This software is freely available to the public for use without a  copyright notice\.  Restrictions cannot be placed on its present or future use\.  Although all reasonable efforts have been taken to ensure the accuracy and reliability of the  software and data, the National Human Genome Research Institute \(NHGRI\) and the U\.S\. Government  does not and cannot warrant the performance or results that may be obtained by using this software  or data\.  NHGRI and the U\.S\. Government disclaims all warranties as to performance, merchantability  or fitness for any particular purpose\.  In any work or product derived from this material, proper attribution of the authors as the source  of the software or data should be made, using "NHGRI Genome Technology Branch" as the citation\.  NOTE: This package includes \(internally\) the sam\-1\.113\.jar library from picard tools\. That package uses the MIT license, which can be accessed using the command:  java \-jar thisjarfile\.jar help samjdkinfo

