# QoRTs: Quality Of Rna-seq Tool Set
> Version 0.3.3 (Updated Thu May 21 12:39:22 EDT 2015)

> ([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "mergeWig"

## USAGE:

    java [Java Options] -jar QoRTs.jar mergeWig [options] outfile


## DESCRIPTION:

This utility merges multiple '\.wig' wiggle files into a single summary '\.wig' wiggle file\. Optionally it can be used to display the mean read\-pair coverage of each window across all input wiggle files rather than the sum\. Also optionally, the mean/sum can be weighted by a set of user\-supplied normalization factors\.
Note: Either the '\-\-filenames' or the '\-\-sampleList' option must be set\! The sampleList option is generally used with the \-\-infilePrefix and \-\-infileSuffix options\.

## REQUIRED ARGUMENTS:
### outfile:

> The name of the output wiggle file, or '-' to write to stdout. 
If this ends with ".gz" or ".zip", then the file will automatically be compressed using the appropriate method. (String)



## OPTIONAL ARGUMENTS:
### --filenames file1.wig,file2.wig,file3.wig.gz,...:

> A comma-delimited list of wiggle files to merge. This is optional, and filenames can be inferred from --infilePrefix, --infileSuffix, and the --sampleList, if those options are specified. Either this option OR --sampleList MUST BE SPECIFIED. (CommaDelimitedListOfStrings)

### --sampleList [sampleList.txt | - | samp1,samp2,samp3,...]:

> Either a comma-delimited list of sample id's or a '.txt' file containing a list of sample id's. The file must either contain no title line, or contain a title line that includes a "sample.ID" column. Either this option OR --filenames MUST BE SPECIFIED. Note that the sample list file must end with the extension '.txt'. Similarly, the sample names CANNOT end with '.txt'. (String)

### --infilePrefix infilePrefix:

> A file prefix for all input wiggle files. Used with the --sampleList parameter. (String)

### --infileSuffix infileSuffix:

> A file suffix for all input wiggle files. Used with the --sampleList parameter. (String)

### --sizeFactorFile val:

> A file containing (at least) two columns: a list of sample ID's and their double-precision floating-point size factors. The first line must include at least two columns: "sample.ID" and "size.factor"If this option is set, all counts will be divided by the given normalization factors. The length must be the same as the length of infiles.If sample.ID's is not specified by the --sampleList or --sampleListFile parameters, then all listed samples will be merged. (String)

### --sizeFactors val,val,val,...:

> A list of double-precision floating-point values. If this or any size factor option is set, all counts will be divided by the given normalization factors. The length must be the same as the number of files to merge. (CommaDelimitedListOfDoubles)

### --makeNegative:

> Flag to indicate that every counting bin value should be multiplied by -1 (flag)

### --calcMean:

> Flag to indicate that the mean average should be calculated, rather than the sum. (flag)

### --trackTitle options:

> The title of the new merged track. (String)

### --additionalTrackOptions options:

> Additional track definition options, added to the track definition line. See the UCSC documentation for more information. (String)

### --ignoreSizeFactors:

> Flag to indicate that this utility should ignore size factors even if they are found in the input listFile. (flag)

### --verbose:

> Flag to indicate that debugging information and extra progress information should be sent to stderr. (flag)

### --quiet:

> Flag to indicate that only errors and warnings should be sent to stderr. (flag)

## AUTHORS:

Stephen W\. Hartley, Ph\.D\. stephen\.hartley \(at nih dot gov\)

## LEGAL:

 This software is "United States Government Work" under the terms of the United States Copyright  Act\.  It was written as part of the authors' official duties for the United States Government and  thus cannot be copyrighted\.  This software is freely available to the public for use without a  copyright notice\.  Restrictions cannot be placed on its present or future use\.  Although all reasonable efforts have been taken to ensure the accuracy and reliability of the  software and data, the National Human Genome Research Institute \(NHGRI\) and the U\.S\. Government  does not and cannot warrant the performance or results that may be obtained by using this software  or data\.  NHGRI and the U\.S\. Government disclaims all warranties as to performance, merchantability  or fitness for any particular purpose\.  In any work or product derived from this material, proper attribution of the authors as the source  of the software or data should be made, using "NHGRI Genome Technology Branch" as the citation\.  NOTE: This package includes \(internally\) the sam\-1\.113\.jar library from picard tools\. That package uses the MIT license, which can be accessed using the command:  java \-jar thisjarfile\.jar help samjdkinfo

