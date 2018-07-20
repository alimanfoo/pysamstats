FROM continuumio/miniconda:4.5.4
LABEL pysamstats - A fast utility for extracting statistics from a SAM or BAM file.

RUN conda install -c bioconda pysamstats
