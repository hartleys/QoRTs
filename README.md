# QoRTs v1.2.26
(Compiled Tue Mar  7 14:22:03 EST 2017)

The [QoRTs software package](http://hartleys.github.io/QoRTs/) is a fast, efficient, and portable 
multifunction toolkit designed to assist in
the analysis, quality control, and data management of RNA-Seq and DNA-Seq datasets. Its primary function is to aid
in the detection and identification of errors, biases, and artifacts produced by high-throughput
sequencing technology. In addition, it can produce count data designed for use with RNA-Seq differential gene expression
and differential exon usage tools, as well as individual-sample and/or group-summary genome track
files suitable for use with the UCSC genome browser (or similar browsers). 
The entire QoRTs toolkit can be used in almost any operating system that supports java and R.
For more information visit the [QoRTs project page](http://hartleys.github.io/QoRTs/index.html). 

##Citing QoRTs:
If you use QoRTs and find it helpful, you can cite it in your publications as:

Hartley SW, Mullikin JC. [**QoRTs: a comprehensive toolset for quality control and data processing of RNA-Seq experiments.**](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4506620/) [*BMC Bioinformatics.*](http://www.biomedcentral.com/bmcbioinformatics) doi: 10.1186/s12859-015-0670-5

##HELP AND DOCUMENTATION:
**Additional help and documentation is available online 
[here](http://hartleys.github.io/QoRTs/index.html).** 

If you have a question, encounter an error, or have a feature request, we recommend that you create an issue on the [github issues page](https://github.com/hartleys/QoRTs/issues).

Alternatively, the creator and maintainer of QoRTs can be reached by emailing "QoRTs-contact" (at) list.nih.gov.

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

