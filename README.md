# Barapost toolkit

**"Barapost"** command line toolkit is designed for determination the taxonomic position of nucleotide sequences and subsequent sorting. In other words, it performs taxonomic annotation of sets of nucleotide sequences and then separates these sets into different files.

- [Motivation](#motivation)
- [The default workflow](#the-default-workflow-looks-like)
- [Getting barapost](#getting-barapost)
- [Pre-requirements and where to get them](#pre-requirements)
- [1. prober](#prober)
- [2. barapost](#barapost)
- [3. fastQA sorter](#fastQA-sorter)
- [FAST5 sorting](#FAST5-sorting)
- [FAST5 untwisting](#FAST5-untwisting)
- [Examples of usage in combination](#Examples-of-usage-in-combination)

## Motivation

>*Rhodococcus* separately, *Pseudomonas* separately!

Find yourself constantly having a large amount of FASTA, FASTQ or FAST5 files, where sequences that belong to different organisms are **mixed up** together, but you want them to be **separated**?

It is awful to sit in front of the computer for hours sending all these sequences to NCBI BLAST server and to rewrite large FAST(A/Q) files by hand, isn't it?

"Barapost" toolkit is the thing that can do it for you and **save your time**.

## The default workflow looks like

**!** - these scripts cannot be executed by Python interpreter version < 3.0.

1. **prober.py** -- this script sends several sequences (i.e. only a part of your data set) to NCBI BLAST server in order to determine what taxonomic units are present in data set. "prober.py" saves information about the best hit of each sequence from probing batch.
Processing all sequences in this way takes too much time, what leads us to "barapost.py".

2. **barapost.py** -- this script firstly downloads best hits "discovered" by "prober.py" from Genbank, then uses these downloaded sequences to build a database on your local machine and finally aligns the rest of data set against builded database. Database building and "BLASTing" is performed by using "BLAST+" toolkit.
Results of taxonomic annotation are written in TSV file named according to name of input file(s).

3. **fastQA5-sorter.py** -- this script performs sorting (dividing into separate files) of your data set according to results of "prober.py" and "barapost.py"

![](imgs/Barapost-wokflow.png)

## Getting barapost

Since current version is too raw to make a release, you can get "barapost" in the following ways:

Way 1: go to terminal and run `git clone https://github.com/masikol/barapost.git`

Way 2: download ZIP archive (green button at the top right of this page "Clone or downlaod" -> "Download ZIP").

## Pre-requirements

1. **Python 3** (https://www.python.org/). Barapost is tested on Python interpreter version 3.8.0.

2. **BLAST+** toolkit is used by "barapost.py" for building a database and aligning.

   It can be downloaded [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/);

   **Installation:**

   - Linux: download tarball, unpack it and add `bin/` directory from unpacked tree to `PATH` variable.
     (this should also work for macOS, but I have not checked it)

   - Windows: download `.exe` executable, run it and follow installation. Do not forget to add "BLAST+" to `PATH` during installation (you will see corresponding checkbox at one of installation steps).

   "barapost.py" has been tested on Linux and Windows with BLAST+ version 2.9.0.

3. [**h5py**](https://www.h5py.org/) Python package is used by "fastQA5-sorter.py" for sorting FAST5 files.

   This package is not necessary if you do not intend to sort FAST5 files.

   To build `h5py`, `libhdf5-dev` is in turn required.

   **Installation**:

   - Linux (assumming that your packet manager is `apt` and you use `pip3` to install Python side packages):
     ```
     apt install libhdf5-dev
     pip3 install h5py
     ```

   Sorting FAST5 files have been tested only on Linux (`h5py` version 2.10.0 was used).

## prober

Version 1.12.i; 2019.11.02 edition;

### DESCRIPTION:

**prober.py** -- this script is designed for determination the taxonomic position
of nucleotide sequences by sending each of them to NCBI BLAST server and regarding the best hit.

The main goal of this script is to send a probing batch (see `-b` option) of sequences to NCBI BLAST server and discover, what Genbank records can be downloaded and used for building a database on your local machine by "barapost.py".

It means that you should not process all your data by "prober.py' -- it would take long time. "barapost.py" will process the rest of you sequences in the same way like "prober.py", but on your local computer.

Obviously, a probing batch cannot cover all variety of a data set, so some sequences can be recognized as "unknown" while processing by "barapost.py". But you always can run "prober.py" again on "unknown" sequences.

This script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.

Results of taxonomic annotation are written in TSV file named according to name of input file(s), that can be found in result directory.

"prober.py" also generates a file named `...probe_acc_list.tsv`. It contains accessions and names of Genbank records (these "best hits" mentioned above) that
  can be used for building a database on your local machine by "barapost.py".

If you have your own FASTA files that can be used as database to blast against, you can omit "prober.py" step and go to "barapost.py" (see `-l` option in "barapost.py" description).


### Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- packet size (see `-p` option): 100 sequences;
- probing batch size (see `-b` option): 200 sequences;
- algorithm (see `-a` option): `megaBlast`;
- organisms (see `-g` option): full `nt` database, i.e. no slices;
- output directory (`-o` option): directory named `barapost_result`
  nested in current directory;
- no email information (`-e` option) is send to NCBI;

### OPTIONS:

- Files that you want "prober.py" to process should be specified as positional arguments (see EXAMPLE #2 below).
  Wildcards do work: `./prober.py my_directory/*` will process all files in `'my_directory'`.

```
    -h (--help) --- show help message;

    -v (--version) --- show version;

    -d (--indir) --- directory which contains FASTQ or FASTA files meant to be processed.
          I.e. all FASTQ and FASTA files in this direcory will be processed; Files can be gzipped.

    -o (--outdir) --- output directory.
          Default value: `barapost_result`;

    -p (--packet-size) --- size of the packet, i.e. number of sequence to blast in one request.
          Value: integer number [1, 500]. Default value is 100;

    -a (--algorithm) --- BLASTn algorithm to use for aligning.
          Available values: 'megaBlast', 'discoMegablast', 'blastn'.
          Default is megaBlast;

    -g (--organisms) --- TaxIDs of organisms to align your sequences against. I.e. 'nt' database slices.
          You can find your Taxonomy IDs here: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
          Format of value (TaxIDs separated by comma): 
            <organism1_name>,<organism2_taxid>...
          See EXAMPLES #4 and #5 below.
          Spaces are not allowed.
          Default is: full 'nt' database, i.e. no slices.

    -b (--probing-batch-size) --- total number of sequences that will be sent to BLAST server
          during 'prober.py' run.
          You can specify '-b all' to process all your sequeces by 'prober.py'.
          Value: positive integer number.
          Default value is 200;

    -e (--email) --- your email. Please, specify your email when you run "prober.py",
        so that the NCBI can contact you if there is a problem. See EXAMPLE #2 below.
```

- More clearly, functionality of `-g` option is totally equal to "Organism" text boxes on this BLASTn page:
    https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome.


### EXAMPLES:

Note for Windows users: `./prober.py` won't work on Windows -- type `py -3 prober.py` instead.

Sure, you can do the same thing on Unix-like systems, but you might face problems with path completions if you call Python interpreter explicitly. Therefore I recommend to make .py-file executable (by running `chmod +x prober.py`) and run it as it is shown in examples below.

  1. Process all FASTA and FASTQ files in working directory with default settings:

`./prober.py`

  2. Process all files in the working directory that start with "some_my_fasta".
Provide NCBI with your email. Use default settings:

`./prober.py some_my_fasta* -e my.email@smth.com`

  3. Process one file with default settings:

`./prober.py reads.fastq`

  4. Process a FASTQ file and a FASTA file with discoMegablast, packet size of 100 sequences.
Search only among Erwinia sequences (551 is Erwinia taxid):

`./prober.py reads_1.fastq.gz some_sequences.fasta -a discoMegablast -p 100 -g 551`

  5. Process all FASTQ and FASTA files in directory named `some_dir`. Process 300 sequences, packet size is 100 sequnces (3 packets will be sent).
Search only among Escherichia (taxid 561) and viral (taxid 10239) sequences:

`./prober.py -d some_dir -g 561,10239 -o outdir -b 300 -p 100`


## barapost

Version 3.5.i; 2019.11.08 edition;

### DESCRIPTION:

**barapost.py** -- this script is designed for taxonomic annotation of nucleotide sequences by "BLASTing" each of them with 'blastn' script from "BLAST+" toolkit
and regarding the best hit.

"barapost.py" is meant to be used just after 'prober.py'.

"barapost.py" downloads records-hits from Genbank according to results (`...probe_acc_list.tsv`) generated by "prober.py", builds a database (on your local machine) which consists of downloaded sequences, and continues aligning the rest of your with "BLAST+" toolkit.

Script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.

"barapost.py" writes (actually appends) it's results in the same TSV file as "prober.py" does.

Files processed by this script are meant to be sorted afterwards by 'fastQA5_sorter.py'.

If you have your own FASTA files that can be used as database to blast against, you can omit "prober.py" step and go to "barapost.py" (see `-l` option).

### Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- packet size (see `-p` option): 100 sequences;
- algorithm (see `-a` option): 'megaBlast';
- numbers of threads to launch (`-t` option): 1 thread.


### OPTIONS:

- Files that you want "barapost.py" to process should be specified as positional arguments (see EXAMPLE #2 below).
  Wildcards do work: `./barapost.py my_directory/*` will process all files in `'my_directory'`.

```
    -h (--help) --- show help message;

    -v (--version) --- show version;

    -r (--taxannot-resdir) --- directory that contains results of taxonomic annotation.
        This is prober's output directory (see prober's `-o` option).
        If you omit 'prober.py' and use only your own FASTA files to create a database,
        this option will behave just like simple output-directory-specifying option.
        Default value is `barapost_result`.

    -d (--indir) --- directory which contains FASTQ or FASTA files meant to be processed.
        I.e. all FASTQ and FASTA files in this direcory will be processed; Files might be gzipped.

    -p (--packet-size) --- size of the packet, i.e. number of sequence to blast in one request.
        Value: integer number [1, 500]. Default value is 100;

    -a (--algorithm) --- BLASTn algorithm to use for aligning.
        Available values: 'megaBlast', 'discoMegablast', 'blastn'.
        Default is megaBlast;

    -l (--local-fasta-to-db) --- your own FASTA file that will be added to downloaded database
        (or used instead of it if you omit 'prober.py' step);

    -t (--threads) --- number of threads to launch;
```

### Notes about using your own FASTA files as database:

1. Besides using `-l` option, you can specify your own FASTA files using accession TSV file generated by "prober.py". To do this, just write your FASTA file's path to this TSV file in new line.

2. "makeblastdb" utility from "BLAST+" toolkit considers first word (it separates words by spaces) of sequence ID in FASTA file as sequence accession. Naturally, duplicated accessions are not allowed. Therefore, in order to avoid this duplication, "barapost.py" uses modified sequence IDs of your own sequences in FASTA files while database creating. It adds custom accession number in the beginning of sequence IDs. This custom accessions have following format: OWN_SEQ_<N>, where <N> is an integer number. Actually, it is order number of this sequence (I mean order of adding to database). These modified sequence IDs are used only in database -- your own FASTA files will be kept intact.

3. If you include SPAdes or a5 assembly FASTA file in the database with "barapost.py", sequence IDs will be modified in a specific (i.e. *ad hoc*) way. If there are **more than one** assembly file generated by **one** assembler (e.g. two files named "contigs.fasta" generated by SPAdes), paths to these "contigs.fasta" files will be added to sequence IDs while database creation. So, sequence IDs will look like, e.g. for SPAdes:

    `OWN_SEQ_4 /some/happy/path/contigs.fasta_NODE_3_length_546787_cov_102.642226`

    Sequences from assembly files affect on sorting process in their own specific way (see "Notes about sorting" section in "fastQA5-sorter.py" decsription below, note number 4).

### EXAMPLES:

Note for Windows users: `./barapost.py` won't work on Windows -- type `py -3 barapost.py` instead.

Sure, you can do the same thing on Unix-like systems, but you might face problems with path completions if you call Python interpreter explicitly. Therefore I recommend to make .py-file executable (by running `chmod +x barapost.py`) and run it as it is shown in examples below.

  1. Process all FASTA and FASTQ files in working directory with default settings:

`./barapost.py`

  2. Process all files in the working directory that start with "some_my_fasta". Use default settings:

`./barapost.py some_my_fasta*`

  3. Process one FASTQ file with default settings.
     File `reads.fastq` has been already processed by "prober.py".
     Results of "prober.py" work are in directory `prober_outdir`:

`./barapost.py reads.fastq -r prober_outdir`

  4. Process FASTQ file and FASTA file with discoMegablast, packet size of 100 sequences.
     Files `reads.fastq.gz` and `another_sequences.fasta` have been already processed by "prober.py".
     Results of "prober.py" work are in directory `prober_outdir`:

`./barapost.py reads.fastq.gz another_sequences.fasta -a discoMegablast -p 100 -r prober_outdir`

  5. Process all FASTQ and FASTA files in directory named `some_dir`.
    All these files have been already processed by "prober.py".
    Results of "prober.py" work are in directory `prober_outdir`:

`/barapost.py -d some_dir -r prober_outdir`

  6. Process file named `some_reads.fastq` **omitting prober step**. Sequence(s) from file `my_own_sequence.fasta` will be included to the database.
     Packet size is 50 sequences. Launch 4 threads.

`./barapost.py some_reads.fastq -p 50 -l my_own_sequence.fasta -t 4`


## fastQA5 sorter
(fast**Q**, fast**A** and fast**5** sorter)

Version 3.3.a; 2019.11.08 edition;

### DESCRIPTION:

**fastQA5-sorter.py** -- this script is designed for sorting (dividing into separate files) FASTQ and FASTA files processed by "barapost.py".

Moreover, it can sort FAST5 files according to taxonomic annotation of FASTQ files, that are in turn results of basecalling of these FAST5 files. See "FAST5 sorting" and "FAST5 untwisting" sections below.

"fastQA5-sorter.py" is meant to be used just after "barapost.py".

### Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- sorting sensitivity (see `-s` option): `"genus"`;
- output directory (`-o` option): directory named `fastQA_sorter_result_<date_and_time_of_run>`
  nested in current directory;
- minimum mean quality of a read to keep (`-q` option): 20 (Phred33);
- length filtering (`-m` option) is disabled by default;
- "FAST5 untwisting" is disaled by default;

### OPTIONS:

- Files that you want "fastQA5-sorter.py" to process should be specified as positional arguments (see EXAMPLE #2 below).
  Wildcards do work: `./fastQA5-sorter.py my_directory/*` will process all files in `'my_directory'`.

```
    -h (--help) --- show help message;

    -v (--version) --- show version;

    -r (--taxannot-resdir) --- directory that contains results of taxonomic annotation.
        This is prober's output directory (see prober's `-o` option).
        Or, if you omit prober step, this is barapost's output directory (see barapost's `-r` option).
        Default value is "barapost_result", since it is the default name of
        output directory generated by "prober.py".

    -d (--indir) --- directory which contains FASTQ and/or FASTA files
        (files can be gzipped) meant to be sorted;

    -o (--outdir) --- output directory;

    -s (--sorting-sensitivity) --- the lowest taxonomy rank that will be used
        in names of result files;
        Available values: "genus", "species", "strain".
        Default value is "genus";

    -q (--min-ph33-qual) --- minimum mean Phred33 quality of a read to keep.
        Reads of lower quality will be written to separate "trash" file.
        Default value: 20;

    -m (--min_seq_len) --- minimum length of a sequence to keep.
        Shorter sequences will be written to separate "trash" file.
        Length filtering is disabled by default;

    -u (--untwist-fast5) --- flag option. If specified, FAST5 files will be
        sorted considering that they and "corresponding" FASTQ files contain
        different reads (like after basecalling performed by Guppy).
        For details, see "FAST5 untwisting" section below.
        Disabled by default;
```

### Notes about sorting:

1. Sorting sensitivity by taxonomic rank is quite symbolic. Common Genbank record name usually looks like this: `'Erwinia amylovora strain E-2 chromosome, complete genome'`. Sorting sensitivity "genus" in this case provides you with sorted file named `"Erwinia.fastq.gz"`. Analogically, "species" leads to `"Erwinia_amylovora.fastq.gz"`. But if you scpecify "strain", you'll have file named `"Erwinia_amylovora_strain_E-2_chromosome.fastq.gz."`, since information after strain ID is quite unpredictable.

2. If Genbank record name looks like `'Shigella sp. PAMC 28760 chromosome, complete genome'`. If you scpecify sorting sensitivity "genus", you will get file named `"Shigella.fastq.gz"`. But "species" and "strain" ones both will provide you with file named `"Shigella_sp._PAMC_28760_chromosome.fastq.gz"`.

3. Hits corresponding to phage sequences are processed in a bit different way. For example, we have Genbank record `'Escherichia phage vB_EcoM_JS09, complete genome'`. If you'll sort by genus or species, you'll get file named `"Escherichia_phage.fastq.gz"`. And "strain" sorting sensitivity will provide you with `"Escherichia_phage_vB_EcoM_JS09.fastq.gz"` one.

4. If you include SPAdes or a5 assembly FASTA files (files that contain contigs) to your database, sequences that hit them will be sorted in a specific way. There can be two situations:

    a) there is **one** FASTA file with assembly (or there are two files, but one of them was generated by SPAdes, and another one -- by a5). In this case if you sort your sequnces by genus, you'll get files named `"SPAdes_assembly_NODE.fastq.gz"` and/or `"a5_assembly_scaffold.fastq.gz"` depending on the assembler. It means that all sequences that hit contigs in your assembly will be placed in one file. If sorting sensitivity is "species" or "strain", you'll get separate files for each NODE (or scaffold, if you've used a5). For example: `"SPAdes_assembly_NODE_6.fastq.gz"` or `"a5_assembly_scaffold_8.fastq.gz"`.

    b) there **more than one** FASTA file with assembly generated by one assembler. For example, you have two SPAdes outputs: `"outdir_1/contigs.fasta"` and `"outdir_2/contigs.fasta"`. Paths to these files will be used to name sorted files. In this case if you sort your sequnces by genus, you'll get files named `"SPAdes_assembly__outdir_1_contigs.fasta_NODE.fastq.gz"` and `"SPAdes_assembly__outdir_2_contigs.fasta_NODE.fastq.gz"`. Situation with a5 is the same -- instead of "SPAdes" "a5" will be written (and "scaffold" instead of "NODE"). As you see, path separators are replaced by underscores in order not to held a bacchanalia in file system. If sorting sensitivity is "species" or "strain", you'll get separate files for each NODE (or scaffold, if you've used a5), just as in 'a)' section above, but paths will be included (e.g. `"SPAdes_assembly__outdir_1_contigs.fasta_NODE_3.fastq.gz"`)

    "fastQA5-sorter.py" decides whether it is a SPAdes assembly file or a a5 one by looking at sequence IDs in these files. If ID is like `"NODE_1_length_245432_cov_23.5412"` -- probably it is SPAdes work. If ID is like `"scaffold_1"` -- it looks just like a5 output. "fastQA5-sorter.py" regards only the first sequnce ID in file.

5. All sequences that do not pass quality and/or length controle will be written to one "trash"-file named in the following way:

    `qual_less_Q<min_quality>_len_less_<min_length>.fastq.gz`.

### EXAMPLES:

Note for Windows users: `./fastQA5-sorter.py` won't work on Windows -- type `py -3 fastQA5-sorter.py` instead.

Sure, you can do the same thing on Unix-like systems, but you might face problems with path completions if you call Python interpreter explicitly. Therefore I recommend to make .py-file executable (by running `chmod +x fastQA5-sorter.py`) and run it as it is shown in examples below.

  1. Process all FASTA, FASTQ and FAST5 files in working directory with default settings:

`./fastQA5-sorter.py`

  2. Process all files in the working directory that start with "some_my_fastq". Place reads with mean Phred33 quality < 15 into separate "trash" file:

`./fastQA5-sorter.py some_my_fastq* -q 15`

  2. Process one FASTQ file with default settings.
     File `reads.fastq` has been already processed by "barapost.py".
     Results of "barapost.py" work are in directory `prober_outdir`:

`./fastQA5-sorter.py reads.fastq.gz -r prober_outdir/`

  3. Process a FASTQ file and a FASTA file, place results in `outdir` directory.
     Files `reads.fastq.gz` and `another_sequences.fasta` have been already processed by "barapost.py".
     Results of "barapost.py" work are in directory `prober_outdir`:

`./fastQA5-sorter.py reads.fastq.gz another_sequences.fasta -o outdir -r prober_outdir/`

  4. Process all FASTQ, FASTA and FAST5 files in directory named `dir_with_seqs`. Sort by species.
     All these files have been already processed by "barapost.py". Perform "FAST5 untwisting".
     Results of "barapost.py" work are in directory `prober_outdir`:

`./fastQA5-sorter.py -d dir_with_seqs -o outdir -r prober_outdir/ -s species -u`


## FAST5 sorting

FAST5 files can be sorted by fastQA5-sorter.py.

**!** - Barapost toolkit does **not** perform basecalling of nanopore data.

Since it is recommended to keep your FAST5 files in order to re-basecall them later, with more accurate (e.g. more sensible for base modifications) basecall algorithms, it worth following the pipeline below:

1. Basecall FAST5 files and get FASTQ.
2. Perform taxonomic annotation of obtained FASTQ files.
3. Sort source FAST5 files according to this taxonomic annotation.
4. Keep sorted FAST5 files in order to re-basecall them later.

Therefore, you can pass FAST5 files to "fastQA5-sorter.py" just as FASTQ or FASTA files and they will be sorted as well.

### Example:

Assuming you have already performed basecalling of file "some_reads.fast5" and have file "some_reads.fastq". You may act as follows to sort source FAST5 file:

Taxonomic annotation:
```
./prober.py some_reads.fastq
./barapost.py some_reads.fastq
```
Sorting:
```
./fastQA5-sorter.py some_reads.fast5
```

## FAST5 untwisting

The problem is following: basecallers (popular Guppy, in particular) often missasign names of input FAST5 and output FASTQ files. In result, source **FAST5** and basecalled **FASTQ** files **contain different reads**. Therefore, straitforward sorting of FAST5 files, that relies on names of "corresponding" FASTQ files (that have ondergone taxonomic annotation) is, in general, impossible.

In fastQA-sorter.py, this issue is solved by developing a "FAST5 untwisting" procedure (it can be enabled by specifying `-u` flag).

"Untwisting" is performed by creating a DBM index file that maps FAST5 files and each read in it to TSV file containing taxonomic annotation information about this read. Subsequent sorting goes on according to this index file.

"Untwisting" procedure also determines, if all reads in input FAST5 files have undergone taxonomic annotation and gives you IDs of missing ones if there are any.

One obvious disadvantage: you may need to perform taxonomic annotation of all your FASTQ files to sort some (maybe not all) FAST5 files from the same data set.

Here another problem arises: how to find out, in which FASTQ file(s) are your reads from given FAST5 file placed? You can find this information in "sequencing_summary" file which is often generated by basecaller (at least, Guppy behaves so). But these files are rather bulky and not very enjoyable to use (and often lack essential information, like names of FASTQ files).

Therefore I will develop an auxiliary tool which will make much more handy summary about how reads are distributed between FAST5 and FASTQ files.

## Examples of usage in combination:

1. You can place all .py-files provided with by this toolkit in a directory that contains some FASTA and FASTQ files and run whole "pipeline" with default settings:

`./prober.py && ./barapost.py && ./fastQA5-sorter.py`

2. You can try these scripts on toy test dataset named `some_reads.fastq` (there are 4 reads):

`./prober.py some_reads.fastq -o some_outdir -g 561,10239 -p 2 -b 2`

`./barapost.py some_reads.fastq -r some_outdir`

`./fastQA5-sorter.py some_reads.fastq -r some_outdir -o some_sorted_reads`
