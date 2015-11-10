# QoRTs v1.0.1
(Compiled Tue Nov 10 10:05:13 EST 2015)

The [QoRTs software package](http://hartleys.github.io/QoRTs/) is a fast, efficient, and portable 
multifunction toolkit designed to assist in
the analysis, quality control, and data management of RNA-Seq datasets. Its primary function is to aid
in the detection and identification of errors, biases, and artifacts produced by paired-end high-throughput
RNA-Seq technology. In addition, it can produce count data designed for use with differential expression
and differential exon usage tools, as well as individual-sample and/or group-summary genome track
files suitable for use with the UCSC genome browser (or any compatible browser). 
The entire QoRTs toolkit can be used in almost any operating system that supports java and R.
For more information visit the [QoRTs project page](http://hartleys.github.io/QoRTs/index.html). 

##Citing QoRTs:
If you use QoRTs and find it helpful, you can cite it in your publications as:

Hartley SW, Mullikin JC. [**QoRTs: a comprehensive toolset for quality control and data processing of RNA-Seq experiments.**](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4506620/) [*BMC Bioinformatics.*](http://www.biomedcentral.com/bmcbioinformatics) doi: 10.1186/s12859-015-0670-5

##HELP AND DOCUMENTATION:
Additional help and documentation is available online 
[here](http://hartleys.github.io/QoRTs/index.html). The creator and maintainer of QoRTs can be reached by emailing "QoRTs-contact" (at) list.nih.gov.

Issues, bug reports, or feature requests can be posted to the [github issues page](https://github.com/hartleys/QoRTs/issues).

##INSTALLATION:
The R package can be installed in R using the command:

    > install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz",
                       repos=NULL, 
                       type="source");

The java jar file does not need to be installed, it can simply be downloaded and used directly.
Just execute it using the java command:
    
    java -jar /path/to/jarfile/QoRTs.jar QC input.bam anno.gtf.gz /output/dir/

See the [FAQ](http://hartleys.github.io/QoRTs/FAQ.html) if you have trouble installing QoRTs.

##EXAMPLE DATA:
The example QC output can be found on the github main page, and installed
with the command:
    
    install.packages("http://hartleys.github.io/QoRTs/QoRTsExampleData_LATEST.tar.gz", 
                     repos = NULL, 
                     type="source")

The original bam files are too large to upload to github. 
These files can be found [here](
https://dl.dropboxusercontent.com/u/103621176/qorts/exData/QoRTsFullExampleData.zip). 
(File is ~1.5gb)

There is also a [comprehensive walkthrough](http://hartleys.github.io/JunctionSeq/doc/example-walkthrough.pdf) of 
the entire analysis pipeline, run on a different 
[example dataset](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/QoRTsPipelineWalkthrough.zip) with 
[example bam files](https://dl.dropboxusercontent.com/u/103621176/pipelineWalkthrough/bamfiles.zip).

##MORE INFORMATION:
For more information, see the [QoRTs vignette](http://hartleys.github.io/QoRTs/doc/QoRTs-vignette.pdf) or the [QoRTs help and documentation page](http://hartleys.github.io/QoRTs/index.html).

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
States Copyright Act. It was written as part of the authors’ official duties 
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

The QoRTs Scala package includes (internally) the sam-JDK library (sam-1.113.jar) from picard tools, which is licensed under the MIT license:
    The MIT License
    Copyright (c) 2009 The Broad Institute
    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use,
    copy, modify, merge, publish, distribute, sublicense, and/or
    sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following
    conditions:
    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.

The MIT license and copyright information can also be accessed using the command:
\begin{framed} \begin{verbatim}
java -jar /path/to/jarfile/QoRTs.jar "?" samjdkinfo
\end{verbatim}  \end{framed}

