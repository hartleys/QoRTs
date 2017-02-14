# QoRTs: Quality Of Rna-seq Tool Set
> Version1.1.8 (Updated Wed Jul 13 13:35:56 EDT 2016)

> ([back to help base](../index.html))

## General Help

## DESCRIPTION:

This jar-file contains the data processing module of the software package QoRTs, which is intended for use with Paired-End or Single-End High-Throughput RNA-Seq data. This tool can perform a number of different functions to assist in assessing the data quality, detecting errors or biases, performing analyses, data cleaning, data visualization, and data formatting.

NOTE: if you run into OutOfMemoryExceptions, try adding the java options: "-Xmx8G"

## GENERAL SYNTAX:

    java [java_options] -jar QoRTs.jar COMMAND [options]

## COMMANDS:
### [QC](QC.html)

> This utility runs a large battery of QC / data processing tools on a single given sam or bam file. This is the primary function of the QoRTs utility. All analyses are run via a single pass through the sam/bam file.

### [mergeCounts](mergeCounts.html)

> This utility merges count, wiggle, and similar data from multiple QoRTs QC runs. This is intended for use in merging the data from multiple technical replicates of the same sample/library.This tool will then merge all count data (including gene-level, exon-level, and known/novel splice-junction) counts, as well as wiggle files, assuming all files use the standard naming conventions (for example, the fwd-strand wiggle files must be named: "QC.wiggle.fwd.wig.gz", etc).If any files are missing, they will be skipped.

### [mergeAllCounts](mergeAllCounts.html)

> This tool uses a replicate decoder to merge count/wiggle data of all techical replicates in a dataset, producing sample-wise counts. You must supply a replicate decoder which indicates which replicates are technical replicates of which samples. This tool will then merges each sample's technical replicates using the "mergeCounts" function.

### [bamToWiggle](bamToWiggle.html)

> Generates '.wig' wiggle files, suitable for use with the UCSC genome browser or similar tools. Wiggle files contain depth-of-coverage counts for all equally-sized windows across the entire genome. These depth-of-coverage counts can be either be by read (the default) or by read-pair.

### [makeJunctionTrack](makeJunctionTrack.html)

> This utility takes the splice junction count files created by the QoRTs QC utility across multiple samples and creates a single merged splice junction 'bed' file that lists each splice junction along with the mean read-pair coverage counts (optionally, the mean normalized counts).This splice junction bed file can be used to visualize splice junction counts using the UCSC genome browser and other similar utilities. Note: Either the '--filenames' or the '--sampleList' option must be set! The sampleList option is generally used with the --infilePrefix and --infileSuffix options. Also note: This command only compiles the named splice junctions. For other unnamed splice junctions such as novel splice junctions with low coverage, or novel splice junctions that bridge multiple genes, use the makeAltJunctionTrack command instead.

### [makeOrphanJunctionTrack](makeOrphanJunctionTrack.html)

> This utility takes the 'orphan' splice junction count files created by the QoRTs QC utility (optionally across multiple samples) and creates a single merged splice junction 'bed' file that lists each splice junction along with the read-pair coverage counts. It can optionally calculate the mean counts, and/or normalize the counts using the supplied normalization size factors.The output splice junction bed file can be used to visualize splice junction counts using the UCSC genome browser, IGV, or other similar utilities.Note: Either the '--filenames' or the '--sampleList' option MUST be set! The sampleList option is generally used with the --infilePrefix and --infileSuffix options to determine the input filenames.

### [mergeNovelSplices](mergeNovelSplices.html)

> This utility takes the QC output from the standard QC utility run on a series of samples and performs two functions: first, it compiles all splice junctions across all samples and filters low-coverage novel splice junctions by mean coverage across all samples (optionally normalized with user-supplied size factors). It then assigns unique identifiers to each novel splice junction that passed this filter, and outputs a special flat gff file listing all exons, annotated splice junctions and passed-filter novel splice junctions with assigned unique identifiers for all features. Next, it uses these unique identifiers to create a new set of JunctionSeq-formatted count files, one for each input sample. This new count file will include counts for the passed-filter novel splice junctions in addition to the usual counts for annotated splice junctions, exons, and aggregated-genes, all listed by the assigned unique identifiers.


### [mergeWig](mergeWig.html)

> This utility merges multiple '.wig' wiggle files into a single summary '.wig' wiggle file. Optionally it can be used to display the mean read-pair coverage of each window across all input wiggle files rather than the sum. Also optionally, the mean/sum can be weighted by a set of user-supplied normalization factors.
Note: Either the '--filenames' or the '--sampleList' option must be set! The sampleList option is generally used with the --infilePrefix and --infileSuffix options.

### [makeFlatGff](makeFlatGff.html)

> When running the QC command, QoRT first generates a set of non-overlapping exonic fragments out of all the exons in the genome annotation gtf file. It then assigns each exonic fragment a unique identifier. Similarly, it assigns every splice junction its own unique identifier. This command can be used to write that data to file.
It can also be used to produce a flattened gff file that adheres to the specifications used by DEXSeq.

### [generateSamplePlots](generateSamplePlots.html)

> This simple function invokes R and generates a simple, single-replicate plots (or a similar pdf report) given a single replicate's QoRTs QC output.

### [longReadClassifier](longReadClassifier.html)

> This function Classifies long reads (such as those from pacbio or oxford nanopore sequencers) It determines which reads are matches to the known isoforms supplied in the gene annotation.Note that many of the options listed below are currently nonfunctional. They are placeholders for future improvements and/or holdovers from the QC function, upon which this function is built.WARNING! THIS FUNCTION IS BETA!

### [makeSimpleJunctionTrack](makeSimpleJunctionTrack.html)

> This utility converts QoRTs splice junction count files into bed format. Unlike makeJunctionTrack and makeOrphanJunctionTrack, this utility is not designed to compile multiple samples or replicates together. It is a simple converter from a QoRTs junction count file to a bed file. The count files come in 3 types: known, novel, and orphan.

## AUTHORS:

Stephen W. Hartley, Ph.D. stephen.hartley (at nih dot gov)
## LEGAL:

    This software is "United States Government Work" under the terms 
    of the United States Copyright Act. It was written as part of 
    the authors' official duties for the United States Government 
    and thus cannot be copyrighted. This software is freely 
    available to the public for use without a copyright notice. 
    Restrictions cannot be placed on its present or future use.
    Although all reasonable efforts have been taken to ensure the 
    accuracy and reliability of the software and data, the National 
    Human Genome Research Institute (NHGRI) and the U.S. Government 
    does not and cannot warrant the performance or results that may 
    be obtained by using this software or data. NHGRI and the U.S. 
    Government disclaims all warranties as to performance, 
    merchantability  or fitness for any particular purpose.
    In any work or product derived from this material, proper 
    attribution of the authors as the source of the software or data 
    should be made, using "NHGRI Genome Technology Branch" as the 
    citation.
    NOTE: This package includes (internally) the sam-1.113.jar 
    library from picard tools. That package uses the MIT license, 
    which can be accessed using the command:
     java -jar thisjarfile.jar help samjdkinfo
