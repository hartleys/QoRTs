# Frequently Asked Questions
v0.3.17
Revised Thu Aug 27 13:53:48 EDT 2015

## What if none of these answers solve my problem?

No problem. Just email me at qorts.contact at gmail dot com. If possible attach the output from QoRTs, and run QoRTs with the "--verbose" option to maximize the debugging information.

## Using QoRTs with different aligners:

Different aligners often use different conventions, in particular with regards to the MAPQ field which is often used to differentiate multi-mapped and uniquely-mapped reads.

* **RNA-STAR**: QoRTs is designed to work with RNA-STAR without any special options.
* **TopHat1**: Prior to version 2.0.0, TopHat used the same MAPQ convention as RNA-STAR now uses. Thus, no special options are needed.
* **TopHat2**: TopHat2 uses a MAPQ of 50 to mark reads that are uniquely aligned. Thus, when using TopHat2 aligned data, you must set the --minMAPQ parameter to 50.
* **MapSplice**: Under many test conditions MapSplice produces malformed SAM files that violate the [current SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf). As such, QoRTs may not be able to properly assess MapSplice output, particularly for paired-ended data.
* **GSNAP**: I have been unable to find documentation describing exactly what the MAPQ field means in GSNAP output. However, setting the --minMAPQ parameter to 20 seems to work.
* **BowTie, BowTie2, BWA, NOVOALIGN, BLAT, BLAST, or ELEND**: The use of older aligners or aligners not specifically designed for use with RNA-Seq data is NOT recommended or supported. You MIGHT be able to get these alignments to work using --minMAPQ 0.

##How do I install QoRTs?
The R package can be installed in R using the command:

    > install.packages("QoRTs_0.3.14.tar.gz", repos = NULL, type="source")

or using the command-line:
    
    R CMD INSTALL QoRTs_0.3.14.tar.gz

Note: If you run into an error that looks like this:

    ERROR: cannot extract package from '0.3.9.tar.gz'

Then you may be attempting to install the entire github repository archive, rather than the actual R package archive. Download the QoRTs\_0.3.9.zip file from <a href="https://github.com/hartleys/QoRTs/releases/latest">the most recent release</a>, unzip it, and use the file "QoRTs\_0.3.9.tar.gz" to install the R package. Do NOT attempt to install the full repository source in R!

The java jar file does not need to be installed. 
Just execute it using the java command:
    
    java -jar /path/to/jarfile/QoRTs.jar QC input.bam anno.gtf.gz /output/dir/

##Non-standard Phred scores

If you see an warning that looks like this:

    WARNING WARNING WARNING: 
    SAM format check:
        Phred Qual > 41!
        You will need to set either --adjustPhredScores or --maxPhredScores
        in order to compute Phred quality metrics! QoRTs WILL throw an error
        if quality metrics are attempted!

Followed by an error that looks like this:

    ERROR! ArrayIndexOutOfBoundsException caught while attempting to store quality score metrics.
           Most likely cause: quality score found above 41!

Then you may have a problem with your Phred score encoding. This can be caused
by one of two possibilities:

(1) Certain newer Illumina platforms (reportedly) produce Phred scores in excess of 41. 
While not explicitly forbidden by the Phred encoding specification (such as it is), 
this is very uncommon at the current time. However, if a Phred score in excess of
41 is encountered then QoRTs will throw a descriptive error like this:

In this case, you will need to use the "--maxPhredScore" parameter, set to the maximum expected
Phred score. According to one (secondhand) report, Illumina stated that the absolute maximum 
(as of August 2015) was 45.

(2) Certain older Illumina platforms (~1.3-1.7) used Phred+64 encoding for their phred scores. 
While technically SAM files should always be converted to Phred+33, QoRTs does offer 
optional support for reading such malformed input files. Use the "--adjustPhredScore 31"
option to automatically subtract 31 from all phred scores, thus converting Phred+64 to Phred+33.
Note that negative Phred scores are NOT supported by QoRTs, and they WILL throw an error.
See the [related wikipedia article](https://en.wikipedia.org/wiki/FASTQ_format#Encoding) for more information. 

##How do I cite QoRTs?
If you use QoRTs and find it helpful, you can cite it in your publications as:

Hartley SW, Mullikin JC. [**QoRTs: a comprehensive toolset for quality control and data processing of RNA-Seq experiments.**](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4506620/) *BMC Bioinformatics.*](http://www.biomedcentral.com/bmcbioinformatics) doi: 10.1186/s12859-015-0670-5
