# QoRTs: Quality of RNA-seq Tool-Set
v1.3.6
Revised Tue Sep 25 11:21:46 EDT 2018 (Release)

The QoRTs software package is a fast, efficient, and portable multifunction toolkit designed to assist in
the analysis, quality control, and data management of RNA-Seq datasets. Its primary function is to aid
in the detection and identification of errors, biases, and artifacts produced by paired-end high-throughput
RNA-Seq technology. In addition, it can produce count data designed for use with differential expression
and differential exon usage tools, as well as individual-sample and/or group-summary genome track
files suitable for use with the UCSC genome browser (or any compatible browser).

The entire QoRTs toolkit can be used in almost any operating system that supports java and R.

The most recent release of QoRTs is available on the [QoRTs github page](http://github.com/hartleys/QoRTs). Additional help and documentation is available online [here](http://hartleys.github.io/QoRTs/index.html). The creator and maintainer of QoRTs can be reached by emailing "QoRTs-contact" (at) list.nih.gov.

##Help Index:

* QoRTs vignette. User manual and general how-to. Available as a [pdf](doc/QoRTs-vignette.pdf) or [on the web](QoRTs-vignette/index.html).
* [Java jar utility help index.](jarHtml/index.html) The Java utility is responsible for the bulk of the data processing, and contains a number of useful tools. It must be run on all replicates.
* [R package help.](Rhtml/index.html) The R package is responsible for producing visualizations, graphs, plots, and summary tables. (also available as a [pdf](doc/QoRTs-reference.pdf))
* The [Example QC data](QoRTsExampleData_LATEST.tar.gz) used in the vignette.
* [Documentation for Older versions of QoRTs](archive.html)
* [Frequently Asked Questions](FAQ.html): Check here if you run into problems!

For advanced users:

* A [comprehensive walkthrough](doc/example-walkthrough.pdf) of an entire analysis pipeline. This walkthrough primarily covers the analysis and data processing. For more information on the quality control steps see the Vignette, above.
* The [example dataset and results](ftp://nhgriftp.nhgri.nih.gov/pub/outgoing/mullikin/QoRTsExample/QoRTsPipelineWalkthrough.zip) for the walkthrough.
* The [raw data files](ftp://nhgriftp.nhgri.nih.gov/pub/outgoing/mullikin/QoRTsExample/QoRTsPipelineWalkthroughData.zip) for the example walkthough (File is ~2gb).

##Citing QoRTs:
If you use QoRTs and find it helpful, you can cite it in your publications as:

Hartley SW, Mullikin JC. [**QoRTs: a comprehensive toolset for quality control and data processing of RNA-Seq experiments.**](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4506620/) [*BMC Bioinformatics.*](http://www.biomedcentral.com/bmcbioinformatics) doi: 10.1186/s12859-015-0670-5


##INSTALLATION:
The latest version of the R package can be installed in R using the command:

    > install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz",
                       repos=NULL, 
                       type="source");

Alternatively, the latest STABLE version can be installed using the command:

    > install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",
                       repos=NULL, 
                       type="source");

The stable version may lack new features available in the latest version, but has undergone much more testing and is less likely to have problems.

The java jar file does not need to be installed, it can simply be downloaded [here](http://hartleys.github.io/QoRTs/QoRTs.jar) and used directly. The latest stable version is available [here](http://hartleys.github.io/QoRTs/QoRTs-STABLE.jar).
Just execute it using the java command:
    
    java -jar /path/to/jarfile/QoRTs.jar QC input.bam anno.gtf.gz /output/dir/

See the [FAQ](FAQ.html) if you have trouble installing QoRTs.

##EXAMPLE DATA:
The example QC output used in the vignette can be found on the github main page, and installed
with the command:
    
    install.packages("http://hartleys.github.io/QoRTs/QoRTsExampleData_LATEST.tar.gz", 
                     repos = NULL, 
                     type="source")

The original bam files are too large to upload to github. See the  
[example dataset](ftp://nhgriftp.nhgri.nih.gov/pub/outgoing/mullikin/QoRTsExample/QoRTsPipelineWalkthrough.zip) with 
[example raw bam and fastq files](ftp://nhgriftp.nhgri.nih.gov/pub/outgoing/mullikin/QoRTsExample/QoRTsPipelineWalkthroughData.zip).

##COMMAND-LINE HELP:
Additional options and syntax information for the main QC java utility 
can be found using the command:

    java -jar /path/to/jarfile/QoRTs.jar QC --man

Options and information about other sub-utilities within the java package
can be found using the command:

    java -jar /path/to/jarfile/QoRTs.jar --man

And for each sub-utility:

    java -jar /path/to/jarfile/QoRTs.jar utilname --man
    
For help with individual R functions in the R utility, use the R command:

    > help(functionname);

For a full listing of all help topics for the R utility, use the R command: 

    > help(package="QoRTs");



##LEGAL:
This software is "United States Government Work" under the terms of the United 
States Copyright Act. It was written as part of the authors' official duties 
for the United States Government and thus QoRT Package User Manual cannot be 
copyrighted. This software is freely available to the public for use without a 
copyright notice. Restrictions cannot be placed on its present or future use.

Although all reasonable efforts have been taken to ensure the accuracy and 
reliability of the software and data, the National Human Genome Research 
Institute (NHGRI) and the U.S. Government does not and cannot warrant the 
performance or results that may be obtained by using this software or data. 
NHGRI and the U.S. Government disclaims all warranties as to performance, 
merchantability or fitness for any particular purpose.

In any work or product derived from this material, proper attribution of the 
authors as the source of the software or data should be made, using "NHGRI 
Genome Technology Branch" as the citation.

NOTE: The scala package includes (internally) the sam-JDK library 
(sam-1.113.jar), from picard tools, which is covered under the MIT license. 
The MIT license and copyright information can be accessed using the command:
java -jar /path/to/jarfile/QoRTs.jar --man samjdkinfo
