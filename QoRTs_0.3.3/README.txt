# QoRTs v0.3.3
(Compiled Thu May 21 12:39:22 EDT 2015)

The QoRTs software package is a fast, efficient, and portable multifunction toolkit designed to assist in
the analysis, quality control, and data management of RNA-Seq datasets. Its primary function is to aid
in the detection and identification of errors, biases, and artifacts produced by paired-end high-throughput
RNA-Seq technology. In addition, it can produce count data designed for use with differential expression
and differential exon usage tools, as well as individual-sample and/or group-summary genome track
files suitable for use with the UCSC genome browser (or any compatible browser).

The entire QoRTs toolkit can be used in almost any operating system that supports java and R.

The most recent release of QoRTs is available on the [QoRTs github page](http://github.com/hartleys/QoRTs). 

##HELP AND DOCUMENTATION:
Additional help and documentation is available online 
[here](http://hartleys.github.io/QoRTs/index.html). 

Issues, bug reports, or feature requests can be posted to the [github issues page](https://github.com/hartleys/QoRTs/issues).

##INSTALLATION:
The R package can be installed in R using the command:

    > install.packages("QoRTs_0.3.3.tar.gz", repos = NULL, type="source")

or using the command-line:
    
    R CMD INSTALL QoRTs_0.3.3.tar.gz

The java jar file does not need to be installed. 
Just execute it using the java command:
    
    java -jar /path/to/jarfile/QoRTs.jar QC input.bam anno.gtf.gz /output/dir/

##EXAMPLE DATA:
The example QC output can be found on the github main page, and installed
with the command:
    
    install.packages("QoRTsExampleData_0.3.3.tar.gz", repos = NULL, type="source")

The original bam files are too large to upload to github. 
These files, along with a comprehensive walkthrough demonstrating 
how all analysis can be run on them, can be found [here](
https://dl.dropboxusercontent.com/u/103621176/qorts/exData/QoRTsFullExampleData.zip). 
(File is ~1.5gb)

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

