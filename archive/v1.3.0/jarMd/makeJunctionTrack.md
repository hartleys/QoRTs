# QoRTs: Quality Of Rna-seq Tool Set
> Version 1.3.0 (Updated Fri Oct 20 11:56:37 EDT 2017)

> ([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "makeJunctionTrack"

## USAGE:

    java [Java Options] -jar QoRTs.jar makeJunctionTrack [options] flatgff.gff.gz outfile


## DESCRIPTION:

This utility takes the splice junction count files created by the QoRTs QC utility across multiple samples and creates a single merged splice junction 'bed' file that lists each splice junction along with the mean read\-pair coverage counts \(optionally, the mean normalized counts\)\.This splice junction bed file can be used to visualize splice junction counts using the UCSC genome browser and other similar utilities\. Note: Either the '\-\-filenames' or the '\-\-sampleList' option must be set\! The sampleList option is generally used with the \-\-infilePrefix and \-\-infileSuffix options\. Also note: This command only compiles the named splice junctions\. For other unnamed splice junctions such as novel splice junctions with low coverage, or novel splice junctions that bridge multiple genes, use the makeAltJunctionTrack command instead\.

## REQUIRED ARGUMENTS:
### flatgff.gff.gz:

> The flattened gff file. (String)


### outfile:

> The output bed file, or '-' to write to stdout. If the filename ends with ".gz" or ".zip" then the file will be compressed using the appropriate method. (String)



## OPTIONAL ARGUMENTS:
### --filenames file1,file2,file3,...:

> A comma-delimited list of wiggle files to merge. This is optional, and filenames can be inferred from --infilePrefix, --infileSuffix, and the --sampleList, if those options are specified. Either this option OR --sampleList MUST BE SPECIFIED. (CommaDelimitedListOfStrings)

### --sampleList [sampleList.txt | - | samp1,samp2,samp3,...]:

> Either a comma-delimited list of sample id's or a file containing a list of sample id's. The file must either contain no title line, or contain a title line that includes a "sample.ID" column. Either this option OR --filenames MUST BE SPECIFIED. Note if the sample list is a file then it must end with the extension '.txt' (String)

### --infilePrefix infilePrefix:

> A file prefix for all input junction count files. By default the full file path should be specified by the infile parameter. (String)

### --infileSuffix infileSuffix:

> A file suffix for all input junction count files. By default the full file path should be specified by the infile parameter. (String)

### --sizeFactorFile val:

> A file containing (at least) two columns: a list of sample ID's and their double-precision floating-point size factors. The first line must include at least two columns: "sample.ID" and "size.factor"If this option is set, all counts will be divided by the given normalization factors. The length must be the same as the length of infiles.If sample.ID's is not specified by the --sampleList or --sampleListFile parameters, then all listed samples will be merged. (String)

### --sizeFactors val:

> A list of double-precision floating-point values. If this or any size factor option is set, all counts will be divided by the given normalization factors. The length must be the same as the number of files to merge. (CommaDelimitedListOfDoubles)

### --calcMean:

> Flag to indicate that the splice junction counts should be averaged, rather than added up. (flag)

### --stranded:

> Flag to indicate that data is stranded. (flag)

### --rgb r,g,b:

> The rgb color for all the bed file lines. (String)

### --trackTitle title:

> The title of the track. By default this will be the same as the title parameter. (String)

### --additionalTrackOptions options:

> Additional track definition options, added to the track definition line. See the UCSC documentation for more information. (String)

### --includeFullSpliceNames:

> Flag to indicate that full splice names, including gene ID, should be used. (flag)

### --ignoreSizeFactors:

> Flag to indicate that this utility should ignore size factors even if they are found in the input listFile. (flag)

### --nonflatgtf:

> Flag to indicate that instead of a "flattened" gff file, a standard-format gtf file has been specified. If this flag is raised, it will automatically create the standard flattened gtf information, in memory. It will not be written to disk (flag)

### --skipAnnotatedJunctions:

> If this option is used, annotated splice junctions will not be included in the output file. Note: this only works if there are novel junctions in the input file. (flag)

### --skipNovelJunctions:

> If this option is used, novel splice junctions will not be included in the output file. (flag)

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

