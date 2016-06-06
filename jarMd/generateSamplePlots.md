# QoRTs: Quality Of Rna-seq Tool Set
> Version 1.1.2 (Updated Mon Jun  6 15:18:44 EDT 2016)

> ([back to main](../index.html)) ([back to java-utility help](index.html))

## Help for java command "generateSamplePlots"

## USAGE:

    java [Java Options] -jar QoRTs.jar generateSamplePlots [options] qcDataDir


## DESCRIPTION:

This simple function invokes R and generates a simple, single\-replicate plots \(or a similar pdf report\) given a single replicate's QoRTs QC output\.

## REQUIRED ARGUMENTS:
### qcDataDir:

> The qc directory in which all the QC files are contained. (String)



## OPTIONAL ARGUMENTS:
### --makePdf:

> Flag to indicate that you want a pdf multi-plot report to be generated. (flag)

### --noPng:

> Flag to indicate that you do NOT want the primary single-png multi-plot to be generated. (flag)

### --makeSeparatePngs:

> Flag to indicate that you want a battery of separate pngs to be generated. (flag)

### --uniqueID id:

> The ID of the replicate. This will be only used for the plot labels. (String)

### --verbose:

> Flag to indicate that debugging information and extra progress information should be sent to stderr. (flag)

### --quiet:

> Flag to indicate that only errors and warnings should be sent to stderr. (flag)

## AUTHORS:

Stephen W\. Hartley, Ph\.D\. stephen\.hartley \(at nih dot gov\)

## LEGAL:

 This software is "United States Government Work" under the terms of the United States Copyright  Act\.  It was written as part of the authors' official duties for the United States Government and  thus cannot be copyrighted\.  This software is freely available to the public for use without a  copyright notice\.  Restrictions cannot be placed on its present or future use\.  Although all reasonable efforts have been taken to ensure the accuracy and reliability of the  software and data, the National Human Genome Research Institute \(NHGRI\) and the U\.S\. Government  does not and cannot warrant the performance or results that may be obtained by using this software  or data\.  NHGRI and the U\.S\. Government disclaims all warranties as to performance, merchantability  or fitness for any particular purpose\.  In any work or product derived from this material, proper attribution of the authors as the source  of the software or data should be made, using "NHGRI Genome Technology Branch" as the citation\.  NOTE: This package includes \(internally\) the sam\-1\.113\.jar library from picard tools\. That package uses the MIT license, which can be accessed using the command:  java \-jar thisjarfile\.jar help samjdkinfo

