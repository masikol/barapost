# Barapost toolkit

**"Barapost"** command line toolkit is designed for FASTA, FASTQ and FAST5 files sorting (i.e. separation into different files) according to taxonomic classification of nucleotide sequences stored in them.

## Motivation

Find yourself constantly having a large amount of FASTA, FASTQ or FAST5 files, where sequences that belong to different organisms are **mixed up** together, but you want them to be **separated**?

It is awful to sit in front of the computer for hours sending all these sequences to NCBI BLAST server and to rewrite large FAST(A/Q) files by hand, isn't it?

"Barapost" toolkit is the thing that can do it for you and **save your time**.

## Getting started

Since current version is too raw to make a release, you can get "barapost" in the following ways:

Way 1: go to terminal and run `git clone https://github.com/masikol/barapost.git`

Way 2: download ZIP archive (green button at the top right of this page "Clone or downlaod" -> "Download ZIP").

## The workflow and what Barapost does

1. **prober.py** -- this script sends several sequences (i.e. only a part of your data set) to NCBI BLAST server in order to determine what taxonomic units are "present" in data set. "prober.py" saves information about the best hit(s) of each sequence from probing batch.
Processing all sequences in this way takes too much time, what leads us to "barapost.py".

2. **barapost.py** -- this script firstly downloads best hits "discovered" by "prober.py" from Genbank, then uses these downloaded sequences to build a database on your local machine and finally aligns the rest of data set against builded database in order to perform texonomic annotation. Database building and "BLASTing" is performed by "BLAST+" toolkit.

3. **fastQA5-sorter.py** -- this script performs sorting (dividing into separate files) of your data set according to results of "prober.py" and "barapost.py"

![](imgs/Barapost-wokflow.png)

## More information:
  
  You can find detailed information about Barapost at [Barapost Wiki](https://github.com/masikol/barapost/wiki).