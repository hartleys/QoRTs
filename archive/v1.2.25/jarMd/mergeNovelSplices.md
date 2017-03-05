# QoRTs: Quality Of Rna-seq Tool Set
> Version 1.2.25 (Updated Sat Mar  4 13:22:59 EST 2017)

> ([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "mergeNovelSplices"

## USAGE:

    java [Java Options] -jar QoRTs.jar mergeNovelSplices [options] infileDir decoderFile annotation.gtf.gz outfileDir


## DESCRIPTION:

This utility takes the QC output from the standard QC utility run on a series of samples and performs two functions: first, it compiles all splice junctions across all samples and filters low\-coverage novel splice junctions by mean coverage across all samples \(optionally normalized with user\-supplied size factors\)\. It then assigns unique identifiers to each novel splice junction that passed this filter, and outputs a special flat gff file listing all exons, annotated splice junctions and passed\-filter novel splice junctions with assigned unique identifiers for all features\. Next, it uses these unique identifiers to create a new set of JunctionSeq\-formatted count files, one for each input sample\. This new count file will include counts for the passed\-filter novel splice junctions in addition to the usual counts for annotated splice junctions, exons, and aggregated\-genes, all listed by the assigned unique identifiers\.


## REQUIRED ARGUMENTS:
### infileDir:

> The input file directory. All samples should be contained inside this directory, in a subdirectory with the same name as the sample's sample.ID. (String)


### decoderFile:

> This tab-delimited file must contain a column titled 'sample.ID'. Optionally, it may also contain a column titled 'size.factor', allowing you to specify size factors generated using DESeq, EdgeR, DEXSeq, CuffLinks, or similar utilities.  If the size.factor column is missing, then the size factors will automatically be calculated using the gene-level read counts for all non-aggregate (ie. non-overlapping) genes, using the standard 'Geometric' normalization method used by DESeq, DESeq2, and CuffDiff, as described in Anders et.al. (PMC3218662). Any columns other than 'sample.ID' and 'size.factor' will be ignored. (String)


### annotation.gtf.gz:

> An input gtf file, containing the reference transcript annotation. A number of transcript annotations are available from ensembl, UCSC, or RefSeq. (String)


### outfileDir:

> The output file directory. This can be the same as the input file directory, in which case this utility will simply place the merged novel/known count files in each sample's subdirectory. (String)



## OPTIONAL ARGUMENTS:
### --minCount num:

> The minimum mean normalized read coverage needed for inclusion of a novel splice junction. By default, equal to 10.0. (Double)

### --minSpan len:

> The minimum (genomic) distance threshold for novel splice junctions. 'Novel splice junctions' that span a distance smaller than this value will be IGNORED. This can be useful because many aligners do not distinguish between deletions and splice junctions. The default is 10 bp. (Int)

### --stranded:

> Flag to indicate that data is stranded. This MUST be the same as the strandedness of the original QoRTs QC run. (flag)

### --stranded\_fr\_secondstrand:

> Nonfunctional, as the strandedness rule will have already been applied. (flag)

### --noGzipInput:

> Flag to indicate that input files are NOT be compressed into the gzip format. By default almost all input files are assumed to be compressed. (flag)

### --noGzipOutput:

> Flag to indicate that output files should NOT be compressed into the gzip format. By default almost all output files are compressed to save space. (flag)

### --verbose:

> Flag to indicate that debugging information and extra progress information should be sent to stderr. (flag)

### --quiet:

> Flag to indicate that only errors and warnings should be sent to stderr. (flag)

## AUTHORS:

Stephen W\. Hartley, Ph\.D\. stephen\.hartley \(at nih dot gov\)

## LEGAL:

 This software is "United States Government Work" under the terms of the United States Copyright  Act\.  It was written as part of the authors' official duties for the United States Government and  thus cannot be copyrighted\.  This software is freely available to the public for use without a  copyright notice\.  Restrictions cannot be placed on its present or future use\.  Although all reasonable efforts have been taken to ensure the accuracy and reliability of the  software and data, the National Human Genome Research Institute \(NHGRI\) and the U\.S\. Government  does not and cannot warrant the performance or results that may be obtained by using this software  or data\.  NHGRI and the U\.S\. Government disclaims all warranties as to performance, merchantability  or fitness for any particular purpose\.  In any work or product derived from this material, proper attribution of the authors as the source  of the software or data should be made, using "NHGRI Genome Technology Branch" as the citation\.  NOTE: This package includes \(internally\) the sam\-1\.113\.jar library from picard tools\. That package uses the MIT license, which can be accessed using the command:  java \-jar thisjarfile\.jar help samjdkinfo

