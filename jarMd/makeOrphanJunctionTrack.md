# QoRTs: Quality Of Rna-seq Tool Set
> Version 1.3.6 (Updated Tue Sep 25 11:21:46 EDT 2018)

> ([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "makeOrphanJunctionTrack"

## USAGE:

    java [Java Options] -jar QoRTs.jar makeOrphanJunctionTrack [options] outfile


## DESCRIPTION:

This utility takes the 'orphan' splice junction count files created by the QoRTs QC utility \(optionally across multiple samples\) and creates a single merged splice junction 'bed' file that lists each splice junction along with the read\-pair coverage counts\. It can optionally calculate the mean counts, and/or normalize the counts using the supplied normalization size factors\.The output splice junction bed file can be used to visualize splice junction counts using the UCSC genome browser, IGV, or other similar utilities\.Note: Either the '\-\-filenames' or the '\-\-sampleList' option MUST be set\! The sampleList option is generally used with the \-\-infilePrefix and \-\-infileSuffix options to determine the input filenames\.

## REQUIRED ARGUMENTS:
### outfile:

> The output bed file, or '-' to write to stdout. If the filename ends with ".gz" or ".zip" then the file will be compressed using the appropriate method. (String)



## OPTIONAL ARGUMENTS:
### --filenames file1,file2,file3,...:

> A comma-delimited list of wiggle files to merge. Allows input files to be specified manually. Alternatively, filenames can be inferred from --infilePrefix, --infileSuffix, and the --sampleList, if those options are specified. Either this parameter OR --sampleList MUST BE SPECIFIED. (CommaDelimitedListOfStrings)

### --sampleList [sampleList.txt | - | samp1,samp2,samp3,...]:

> Either a comma-delimited list of sample id's or a file containing a list of sample id's. The file must either contain no title line, or contain a title line that includes a "sample.ID" column. Either this option OR --filenames MUST BE SPECIFIED. Note if the sample list is a file then it must end with the extension '.txt' (String)

### --infilePrefix infilePrefix:

> A file prefix for all input junction count files. By default the full file path should be specified by the infile parameter. (String)

### --infileSuffix infileSuffix:

> A file suffix for all input junction count files. By default the full file path should be specified by the infile parameter. (String)

### --sizeFactorFile val:

> A file containing (at least) two columns: a list of sample ID's and their double-precision floating-point size factors. The first line must include at least two columns: "sample.ID" and "size.factor"If this option is set, all counts will be divided by the given normalization factors. Sample ID's will be matched from the --sampleList parameter. (String)

### --sizeFactors val:

> A list of double-precision floating-point values. If this or any size factor option is set, all counts will be divided by the given normalization factors. The length must be the same as the number of files to merge. Must have the same length and ordering as the --sampleList or --filenames parameter. (CommaDelimitedListOfDoubles)

### --calcMean:

> Flag to indicate that the splice junction counts should be averaged, rather than added up. (flag)

### --stranded:

> Flag to indicate that data is stranded. (flag)

### --junctionTypes junctionTypes:

> Whether to include ambiguous junctions, orphaned junctions, or both. Comma-delimited list (no spaces!). (CommaDelimitedListOfStrings)

### --rgb r,g,b:

> The rgb color for all the bed file lines. (String)

### --trackTitle title:

> The title of the track. By default this will be the same as the title parameter. (String)

### --additionalTrackOptions options:

> Additional track definition options, added to the track definition line. See the UCSC documentation for more information. (String)

### --ignoreSizeFactors:

> Flag to indicate that this utility should ignore size factors even if they are found in the input listFile. (flag)

### --digits num:

> The number of digits after the decimal to include for counts. (Int)

### --title title:

> A title to be prepended to each splice junction name (String)

### --verbose:

> Flag to indicate that debugging information and extra progress information should be sent to stderr. (flag)

### --quiet:

> Flag to indicate that only errors and warnings should be sent to stderr. (flag)

## AUTHORS:

Stephen W\. Hartley, Ph\.D\. stephen\.hartley \(at nih dot gov\)

## LEGAL:

 This software is "United States Government Work" under the terms of the United States Copyright  Act\.  It was written as part of the authors' official duties for the United States Government and  thus cannot be copyrighted\.  This software is freely available to the public for use without a  copyright notice\.  Restrictions cannot be placed on its present or future use\.  Although all reasonable efforts have been taken to ensure the accuracy and reliability of the  software and data, the National Human Genome Research Institute \(NHGRI\) and the U\.S\. Government  does not and cannot warrant the performance or results that may be obtained by using this software  or data\.  NHGRI and the U\.S\. Government disclaims all warranties as to performance, merchantability  or fitness for any particular purpose\.  In any work or product derived from this material, proper attribution of the authors as the source  of the software or data should be made, using "NHGRI Genome Technology Branch" as the citation\.  NOTE: This package includes \(internally\) the sam\-1\.113\.jar library from picard tools\. That package uses the MIT license, which can be accessed using the command:  java \-jar thisjarfile\.jar help samjdkinfo

