

library("knitr");
library("Cairo");
library("QoRTs");
library("MASS");
library("DESeq2");
library("edgeR");
library("BiocStyle");

knit('QoRTs-advanced.Rnw');

tools::texi2dvi("QoRTs-advanced.tex", pdf = TRUE);

