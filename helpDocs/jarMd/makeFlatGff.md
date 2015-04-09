# QoRTs: Quality Of Rna-seq Tool Set
> Version 0.2.13 (Updated Thu Apr  9 11:56:46 EDT 2015)

> ([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "makeFlatGff"

## USAGE:

    java [Java Options] -jar QoRTs.jar makeFlatGff [options] gtffile flatgfffile


## DESCRIPTION:

When running the QC command, QoRT first generates a set of non\-overlapping exonic fragments out of all the exons in the genome annotation gtf file\. It then assigns each exonic fragment a unique identifier\. Similarly, it assigns every splice junction its own unique identifier\. This command can be used to write that data to file\.
It can also be used to produce a flattened gff file that adheres to the specifications used by DEXSeq\.

## REQUIRED ARGUMENTS:
### gtffile:

> The gtf annotation file. This tool was designed to use the standard gtf annotations provided by Ensembl, but other annotations can be used as well. Note: if the file ends in .zip or .gz the compression method will be auto-detected and read accordingly. (String)


### flatgfffile:

> The output destination for the "flattened" gff annotation file to be created, or '-' to write to stdout. Note: if the filename ends in ".zip" or ".gz" the corresponding compression method will be applied. (String)



## OPTIONAL ARGUMENTS:
### --stranded:

> The strandedness mode. Note that to precisely replicate DEXSeq behavior, always use the --stranded mode regardless of the strandedness of your dataset. However: for most purposes it is usually safer to use the same strandedness mode as your dataset. Otherwise, genes that overlap on different strands will not be identified as such, and instead these reads will simply be ignored as "ambiguous" during the counting step. This may lead to misleading read counts. (flag)

### --DEXSeqFmt:

> Flag to indicate that the output gff file should be formatted for use with DEXSeq. (flag)

### --verbose:

> Flag to indicate that debugging information and extra progress information should be sent to stderr. (flag)

### --quiet:

> Flag to indicate that only errors and warnings should be sent to stderr. (flag)

## AUTHORS:

Stephen W\. Hartley, Ph\.D\. stephen\.hartley \(at nih dot gov\)

## LEGAL:

 This software is "United States Government Work" under the terms of the United States Copyright  Act\.  It was written as part of the authors' official duties for the United States Government and  thus cannot be copyrighted\.  This software is freely available to the public for use without a  copyright notice\.  Restrictions cannot be placed on its present or future use\.  Although all reasonable efforts have been taken to ensure the accuracy and reliability of the  software and data, the National Human Genome Research Institute \(NHGRI\) and the U\.S\. Government  does not and cannot warrant the performance or results that may be obtained by using this software  or data\.  NHGRI and the U\.S\. Government disclaims all warranties as to performance, merchantability  or fitness for any particular purpose\.  In any work or product derived from this material, proper attribution of the authors as the source  of the software or data should be made, using "NHGRI Genome Technology Branch" as the citation\.  NOTE: This package includes \(internally\) the sam\-1\.113\.jar library from picard tools\. That package uses the MIT license, which can be accessed using the command:  java \-jar thisjarfile\.jar help samjdkinfo

