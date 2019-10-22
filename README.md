# Barapost toolkit

**"Barapost"** command line toolkit is designed for determinating the taxonomic position of nucleotide sequences and subsequent sorting. In other words, it performs taxonomic annotation of sets of nucleotide sequences.

- [Motivation](#motivation)
- [The default workflow](#the-default-workflow-looks-like)
- [Getting barapost](#getting-barapost)
- [1. prober](#prober)
- [2. barapost](#barapost)
- [3. fastQA sorter](#fastQA-sorter)
- [FAST5 sorting](#FAST5-sorting)
- [Examples of usage in combination](#Examples-of-usage-in-combination)

## Motivation

>*Rhodococcus* separately, *Pseudomonas* separately!

Everyone can encounter a situation, when there is a large amount of FASTA or FASTQ files, where sequences that belong to different organisms are **mixed up** together, but you want them to be **separated**.

It is awful to to sit in front of the computer for hours sending all these sequences to NCBI BLAST server and to rewrite large FAST(A/Q) files by hand, isn't it?

"Barapost" toolkit is the thing that can do it for you and **save your time**.

## The default workflow looks like

**!** - these scripts cannot be executed by Python interpreter version < 3.0. They have been tested on Python interpreter version 3.8.0.

1. **prober-v1-12d.py** -- this script sends several sequences (aka probing batch) to NCBI BLAST server in order to determine what taxonomic units are present in data set. "prober-v1-12d.py" saves information about the best hit of each sequence from probing batch.
Process all sequences in this way takes too much time, what leads us to "barapost-v3-5c.py".

2. **barapost-v3-5c.py** -- this script firstly downloads best hits discovered by "prober-v1-12d.py" from Genbank, then uses these downloaded sequences to build a database on your local machine and finally aligns the rest of data set against builded database. Database building and "BLASTing" is performed by using "BLAST+" toolkit.
Results are written in TSV file named `...results.tsv`.

3. **fastQA5-sorter-v3-1a.py** -- this script performs sorting (dividing into separate files) of your data set according to results of "prober-v1-12d.py" and "barapost-v3-5c.py"

## Getting barapost

Since current version is too raw to make a release, you can get "barapost" in the following ways:

Way 1: go to terminal and run `git clone https://github.com/masikol/barapost.git`

Way 2: download ZIP archive (green button at the top right of this page "Clone or downlaod" -> "Download ZIP").

## prober

Version 1.12.d; 2019.10.22 edition;

### DESCRIPTION:

**prober-v1-12d.py** -- this script is designed for determinating the taxonomic position
of nucleotide sequences by sending each of them to NCBI BLAST server and regarding the best hit.

The main goal of this script is to send a probing batch of sequences to NCBI BLAST server
and discover, what Genbank records can be downloaded and used for building a database
on your local machine by "barapost-v3-5c.py".

This script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.

Results of the work of this script are written to TSV files, that can be found in result directory:

1) There is a file named `...acc_list.tsv`. It contains accessions and names of Genbank records that
    can be used for building a database on your local machine by "barapost-v3-5c.py".

2) There is a file named `...result.tsv`. It contains full result of "BLASTing".
    Results of barapost-v3-5c.py's work will be appended to this file.

Files processed by this script are meant to be processed afterwards by "barapost-v3-5c.py".

If you have your own FASTA files that can be used as database to blast against, you can omit "prober-v1-12d.py" step and go to "barapost-v3-5c.py" (see `-l` option in "barapost-v3-5c.py" description).


### Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- packet size (see `-p` option): 100 sequences;
- probing batch size (see `-b` option): 200 sequences;
- algorithm (see `-a` option): `megaBlast`;
- organisms (see `-g` option): full `nt` database, i.e. no slices;
- output directory (`-o` option): directory named `prober_result`
  nested in current directory;
- no email information (`-e` option) is send to NCBI;

Dedication of this script is to send small batch (see `-b` option) of sequences to NCBI BLAST server.
It means that you should not process all your data by "prober-v1-12d.py' -- it would take long time.

Instead of this you should process some sequences by "prober-v1-12d.py" -- it will determine,
what Genbank records (genomes, if you want) are present in your data and then go to "barapost-v3-5c.py".

"barapost-v3-5c.py" will process the rest of you sequences in the same way like "prober-v1-12d.py", but on your local computer.
"barapost-v3-5c.py" uses 'BLAST+' toolkit for this purpose. It will be much faster.

Obviously, a probing batch cannot cover all variety of a data set,
so some sequences can be recognized as "unknown" while processing by "barapost-v3-5c.py".
But you always can run "prober-v1-12d.py" again on "unknown" sequences.

### OPTIONS:

- Files that you want "prober-v1-12d.py" to process should be specified as positional arguments (see EXAMPLE #2 below).
  Wildcards do work: `./prober-v1-12d.py my_directory/*` will process all files in `'my_directory'`.

```
    -h (--help) --- show help message;

    -d (--indir) --- directory which contains FASTQ or FASTA files meant to be processed.
          I.e. all FASTQ and FASTA files in this direcory will be processed; Files can be gzipped.

    -o (--outdir) --- output directory;

    -p (--packet-size) --- size of the packet, i.e. number of sequence to blast in one request.
          Value: integer number [1, 500]. Default value is 100;

    -a (--algorithm) --- BLASTn algorithm to use for aligning.
          Available values: 'megaBlast', 'discoMegablast', 'blastn'.
          Default is megaBlast;

    -g (--organisms) --- 'nt' database slices, i.e. organisms that you expect to see in result files.
          Format of value (TaxIDs separated by comma): 
            <organism1_name>,<organism2_taxid>...
          See EXAMPLES #4 and #5 below.
          Spaces are not allowed.
          Default is: full 'nt' database, i.e. no slices.

    -b (--probing-batch-size) --- number of sequences that will be aligned on BLAST server
          during 'prober-v1-12d.py' work.
          You can specify '-b all' to process all your sequeces by 'prober-v1-12d.py'.
          Value: positive integer number.
          Default value is 200;

    -e (--email) --- your email. Please, specify your email when you run "prober-v1-12d.py",
        so that the NCBI can contact you if there is a problem. See EXAMPLE #2 below.
```

- More clearly, functionality of `-g` option is totally equal to "Organism" text boxes on this BLASTn page:
    https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome.
- You can find your Taxonomy IDs here: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi.


### EXAMPLES:

Note for Windows users: `./prober-v1-12d.py` won't work on Windows -- type `python prober-v1-12d.py` instead.

Sure, you can do the same thing on Unix-like systems, but you might face problems with path completions if you call Python interpreter explicitly. Therefore I recommend to make .py-file executable (by running `chmod +x prober-v1-12d.py`) and run it as it is shown in examples below.

  1. Process all FASTA and FASTQ files in working directory with default settings:

`./prober-v1-12d.py`

  2. Process all files in the working directory that start with "some_my_fasta".
Provide NCBI with your email. Use default settings:

`./prober-v1-12d.py some_my_fasta* -e my.email@smth.com`

  3. Process one file with default settings:

`./prober-v1-12d.py reads.fastq`

  4. Process a FASTQ file and a FASTA file with discoMegablast, packet size of 100 sequences.
Search only among Erwinia sequences (551 is Erwinia taxid):

`./prober-v1-12d.py reads_1.fastq.gz some_sequences.fasta -a discoMegablast -p 100 -g 551`

  5. Process all FASTQ and FASTA files in directory named `some_dir`. Process 300 sequences, packet size is 100 sequnces (3 packets will be sent).
Search only among Escherichia (taxid 561) and viral (taxid 10239) sequences:

`./prober-v1-12d.py -d some_dir -g 561,10239 -o outdir -b 300 -p 100`


## barapost

Version 3.5.c; 2019.10.22 edition;

### DESCRIPTION:

**barapost-v3-5c.py** -- this script is designed for determinating the taxonomic position
of nucleotide sequences by "BLASTing" each of them with 'blastn' script from "BLAST+" toolkit
and regarding the best hit.

"barapost-v3-5c.py" is meant to be used just after 'prober-v1-12d.py'.

"barapost-v3-5c.py" downloads records-hits from Genbank according to results (`...acc_list.tsv`)
of work of "prober-v1-12d.py", builds an indexed local database which consists of
downloaded sequences, and continues aligning with "BLAST+" toolkit in order to save time.

script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.

"barapost-v3-5c.py" writes it's results in the same TSV file as "prober-v1-12d.py" does.

FASTQ files processed by this script are meant to be sorted afterwards by 'fastQA_sorted.py'.

`BLAST+` toolkit (including blastn, makeblastdb and makembindex) can
be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download);


### Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- packet size (see '-p' option): 100 sequences;
- algorithm (see '-a' option): 'megaBlast';


### OPTIONS:

- Files that you want "barapost-v3-5c.py" to process should be specified as positional arguments (see EXAMPLE #2 below).
  Wildcards do work: `./barapost-v3-5c.py my_directory/*` will process all files in `'my_directory'`.

```
    -h (--help) --- show help message;

    -r (--prober-result-dir) --- result directory generated by script "prober-v1-12d.py".
        This is directory specified to 'prober-v1-12d.py' by '-o' option.
        If you omit 'prober-v1-12d.py' and use your own FASTA files
        to create a database, this directory may not exist before start of 'barapost-v3-5c.py'
        (i.e. it will be a simple output directory).
        Default value is "prober_result", since it is the default name of
        output directory generated by "prober-v1-12d.py"

    -d (--indir) --- directory which contains FASTQ or FASTA files meant to be processed.
        I.e. all FASTQ and FASTA files in this direcory will be processed; Files might be gzipped.

    -p (--packet-size) --- size of the packet, i.e. number of sequence to blast in one request.
        Value: integer number [1, 500]. Default value is 100;

    -a (--algorithm) --- BLASTn algorithm to use for aligning.
        Available values: 'megaBlast', 'discoMegablast', 'blastn'.
        Default is megaBlast;

    -l (--local-fasta-to-db) --- your own FASTA file that will be added to downloaded database
        or used instead of it if you omit 'prober-v1-12d.py' step;

    -t (--threads) --- number of threads to launch;
```

### Notes about using your own FASTA files as database:

1. Besides using `-l` option, you can specify your own FASTA files using accession TSV file generated by "prober-v1-12d.py". To do this, just write your FASTA file's path to this TSV file in new line.

2. "makeblastdb" utility from "BLAST+" toolkit considers first word (it separates words by spaces) of sequence ID in FASTA file as sequence accession. Naturally, duplicated accessions are not allowed. Therefore, in order to avoid this duplication, "barapost-v3-5c.py" uses modified sequence IDs of your own sequences in FASTA files while database creating. It adds custom accession number in the beginning of sequence IDs. This custom accessions have following format: OWN_SEQ_<N>, where <N> is an integer number. Actually, it is order number of this sequence (I mean order of adding to database). Do not worry: these modified sequence IDs are used only in database -- your own FASTA files will be kept intact.

3. If you include SPAdes or a5 assembly FASTA file in the database with "barapost-v3-5c.py", sequence IDs will be modified in a specific (i.e. *ad hoc*) way. If there are **more than one** assembly file generated by **one** assembler (e.g. two files named "contigs.fasta" generated by SPAdes), paths to these "contigs.fasta" files will be added to sequence IDs while database creation. So, sequence IDs will look like, e.g. for SPAdes:

    `OWN_SEQ_4 /some/happy/path/contigs.fasta_NODE_3_length_546787_cov_102.642226`

    Sequences from assembly files affect on sorting process in their own specific way (see "Notes about sorting" section below, note number 4).

### EXAMPLES:

Note for Windows users: `./barapost-v3-5c.py` won't work on Windows -- type `python barapost-v3-5c.py` instead.

Sure, you can do the same thing on Unix-like systems, but you might face problems with path completions if you call Python interpreter explicitly. Therefore I recommend to make .py-file executable (by running `chmod +x barapost-v3-5c.py`) and run it as it is shown in examples below.

  1. Process all FASTA and FASTQ files in working directory with default settings:

`./barapost-v3-5c.py`

  2. Process all files in the working directory that start with "some_my_fasta". Use default settings:

`./barapost-v3-5c.py some_my_fasta*`

  3. Process one FASTQ file with default settings.
     File `reads.fastq` has been already processed by "prober-v1-12d.py".
     Results of "prober-v1-12d.py" work are in directory `prober_outdir`:

`./barapost-v3-5c.py reads.fastq -r prober_outdir`

  4. Process FASTQ file and FASTA file with discoMegablast, packet size of 100 sequences.
     Files `reads.fastq.gz` and `another_sequences.fasta` have been already processed by "prober-v1-12d.py".
     Results of "prober-v1-12d.py" work are in directory `prober_outdir`:

`./barapost-v3-5c.py reads.fastq.gz another_sequences.fasta -a discoMegablast -p 100 -r prober_outdir`

  5. Process all FASTQ and FASTA files in directory named `some_dir`.
    All these files have been already processed by "prober-v1-12d.py".
    Results of "prober-v1-12d.py" work are in directory `prober_outdir`:

`/barapost-v3-5c.py -d some_dir -r prober_outdir`

  6. Process file named `some_reads.fastq`. This file has been already processed by "prober-v1-12d.py".
     Results of "prober-v1-12d.py" work are in directory `prober_outdir`. Sequence from file `my_own_sequence.fasta` will be included to the database.
     Packet size is 50 sequences. Launch 4 threads.

`./barapost-v3-5c.py some_reads.fastq -p 50 -l my_own_sequence.fasta -t 4 -r prober_outdir`


## fastQA5 sorter
(fast**Q**, fast**A** and fast**5** sorter)

Version 3.1.a; 2019.10.22 edition;

### DESCRIPTION:

**fastQA5-sorter-v3-1a.py** -- this script is designed for sorting (dividing into separate files) FASTQ and FASTA files processed by "barapost-v3-5c.py".

Moreover, it can sort FAST5 files according to taxonomical annotation of FASTQ files, that are result of basecalling these FAST5 files.

"fastQA5-sorter-v3-1a.py" is meant to be used just after "barapost-v3-5c.py".

### Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- sorting sensitivity (see `-s` option): `"genus"`;
- output directory (`-o` option): directory named `fastQA_sorter_result_<date_and_time_of_run>`
  nested in current directory;

### OPTIONS:

- Files that you want "fastQA5-sorter-v3-1a.py" to process should be specified as positional arguments (see EXAMPLE #2 below).
  Wildcards do work: `./fastQA5-sorter-v3-1a.py my_directory/*` will process all files in `'my_directory'`.

```
    -h (--help) --- show help message;

    -r (--prober-result-dir) --- result directory genearted by script "prober-v1-12d.py"
        This is directory specified to "prober-v1-12d.py" with '-o' option.
        Default value is "prober_result", since it is the default name of
        output directory generated by "prober-v1-12d.py".

    -d (--indir) --- directory which contains FASTQ and/or FASTA files
        (files can be gzipped) meant to be sorted;

    -o (--outdir) --- output directory;

    -s (--sorting-sensitivity) --- the lowest taxonomy rank that will be used
        in names of result files;
        Available values: "genus", "species", "strain".
        Default value is "genus";

    -q (--min-ph33-qual) --- minimum mean Phred33 quality of a read to keep.
        Reads of lower quality will be ignored.
        Default value: 20;
```

### Notes about sorting:

1. Sorting sensitivity by taxonomic rank is quite symbolic. Common Genbank record name usually looks like this: `'Erwinia amylovora strain E-2 chromosome, complete genome'`. Sorting sensitivity "genus" in this case provides you with sorted file named `"Erwinia.fastq.gz"`. Analogically, "species" leads to `"Erwinia_amylovora.fastq.gz"`. But if you scpecify "strain", you'll have file named `"Erwinia_amylovora_strain_E-2_chromosome.fastq.gz."`, since information after strain ID is quite unpredictable.

2. If Genbank record name looks like `'Shigella sp. PAMC 28760 chromosome, complete genome'`. If you scpecify sorting sensitivity "genus", you will get file named `"Shigella.fastq.gz"`. But "species" and "strain" ones both will provide you with file named `"Shigella_sp._PAMC_28760_chromosome.fastq.gz"`.

3. Hits corresponding to phage sequences are processed in a bit different way. For example, we have Genbank record `'Escherichia phage vB_EcoM_JS09, complete genome'`. If you'll sort by genus or species, you'll get file named `"Escherichia_phage.fastq.gz"`. And "strain" sorting sensitivity will provide you with `"Escherichia_phage_vB_EcoM_JS09.fastq.gz"` one.

4. If you include SPAdes or a5 assembly FASTA files to your database, sequences that hit them will be sorted in a specific way. There are two situations:

    a) there is **one** FASTA file with assembly (or there are two files, but one of them was generated by SPAdes, and another -- by a5). In this case if you sort your sequnces by genus, you'll get files named `"SPAdes_assembly_NODE.fastq.gz"` and/or `"a5_assembly_scaffold.fastq.gz"` depending on the assembler. It means that all sequences that hit contigs in your assembly will be placed in one file. If sorting sensitivity is "species" or "strain", you'll get separate files for each NODE (or scaffold, if you've used a5). For example: `"SPAdes_assembly_NODE_6.fastq.gz"` or `"a5_assembly_scaffold_8.fastq.gz"`.

    b) there **more than one** FASTA file with assembly generated by one assembler. For example, you have two SPAdes outputs: `"outdir_1/contigs.fasta"` and `"outdir_2/contigs.fasta"`. Paths to these files will be used to name sorted files. In this case if you sort your sequnces by genus, you'll get files named `"SPAdes_assembly__outdir_1_contigs.fasta_NODE.fastq.gz"` and `"SPAdes_assembly__outdir_2_contigs.fasta_NODE.fastq.gz"`. Situation with a5 is the same -- instead of "SPAdes" "a5" will be written (and "scaffold" instead of "NODE"). As you see, path separators are replaced by underscores in order not to held a bacchanalia in file system. If sorting sensitivity is "species" or "strain", you'll get separate files for each NODE (or scaffold, if you've used a5), just as in 'a)' section above, but paths will be included (e.g. `"SPAdes_assembly__outdir_1_contigs.fasta_NODE_3.fastq.gz"`)

"fastQA_sorted.py" decides whether it is a SPAdes assembly file or a a5 one by looking at sequence IDs in these files. If ID is like `"NODE_1_length_245432_cov_23.5412"` -- probably it is SPAdes work. If ID is like `"scaffold_1"` -- it looks just like a5 output. "fastQA5-sorter-v3-1a.py" regards only the first sequnce ID in file.

### EXAMPLES:

Note for Windows users: `./fastQA5-sorter-v3-1a.py` won't work on Windows -- type `python fastQA5-sorter-v3-1a.py` instead.

Sure, you can do the same thing on Unix-like systems, but you might face problems with path completions if you call Python interpreter explicitly. Therefore I recommend to make .py-file executable (by running `chmod +x fastQA5-sorter-v3-1a.py`) and run it as it is shown in examples below.

  1. Process all FASTA and FASTQ files in working directory with default settings:

`./fastQA5-sorter-v3-1a.py`

  2. Process all files in the working directory that start with "some_my_fastq". Ignore reads with mean Phred33 quality < 15. The rest of settings are default:

`./fastQA5-sorter-v3-1a.py some_my_fastq* -q 15`

  2. Process one FASTQ file with default settings.
     File `reads.fastq` has been already processed by "barapost-v3-5c.py".
     Results of "barapost-v3-5c.py" work are in directory `prober_outdir`:

`./fastQA5-sorter-v3-1a.py reads.fastq.gz -r prober_outdir/`

  3. Process a FASTQ file and a FASTA file, place results in `outdir` directory.
     Files `reads.fastq.gz` and `another_sequences.fasta` have been already processed by "barapost-v3-5c.py".
     Results of "barapost-v3-5c.py" work are in directory `prober_outdir`:

`./fastQA5-sorter-v3-1a.py reads_1.fastq.gz some_sequences_2.fasta -o outdir -r prober_outdir/`

  4. Process all FASTQ and FASTA files in directory named `dir_with_seqs`. Sort by genus.
     All these files have been already processed by "barapost-v3-5c.py".
     Results of "barapost-v3-5c.py" work are in directory `prober_outdir`:

`./fastQA5-sorter-v3-1a.py -d dir_with_seqs -o outdir -r prober_outdir/ -s genus`


## FAST5 sorting

FAST5 files can be sorted by fastQA5-sorter-v3-1a.py.

**!** - Barapost toolkit does **not** perform basecalling of nanopore data.

Since it is recommended to keep your FAST5 files in order to re-basecall them later, with more accurate (e.g. more sensible for base modifications) basecall algorithms, it worth following the pipeline below:

1. Basecall FAST5 files and get FASTQ.
2. Perform taxonomical annotation of obtained FASTQ files.
3. Sort source FAST5 files according to this taxonomical annotation.
4. Keep sorted FAST5 files in order to re-basecall them later.

Therefore, you can pass FAST5 files to "fastQA5-sorter-v3-1a.py" just as FASTQ or FASTA files and they will be sorted as well.

### Pre-requirements:

Sorting of FAST5 files have been tested only on Linux.

[h5py](https://www.h5py.org/) Python package is used to perform this feature. To build it, `libhdf5-dev` is in turn required. To install these soft, run following commands (assumming that your packet manager is `apt` and you use `pip3` to install Python side packages):

```
apt install libhdf5-dev
pip3 install h5py
```

### The constraint on names of FAST5 and FASTQ files

The short and reliable rule souds so: **do not rename your FASTQ and FAST5 files after basecalling** if you intend to sort these FAST5 files with "fastQA5-sorter-v3-1a.py".

Details:

- Software that drives nanopore sequencing devices (like MinKNOW) yields FAST5 files named as follows:
  `FAK94973_e6f2851ddd414655574208c18f2f51e590bf4b27_1.fast5`
  (at least, MinKNOW does so, I am not sure about other programs).

- Basecallers (Guppy, for instance), can drop `FAK94973` part of file name or replace it with some other string (e.g. `fastq_runid`), but `e6f2851ddd414655574208c18f2f51e590bf4b27_1` part is kept intact. Therefore the second one (let's call it "check string") can be used to identify basecalled FASTQ file and associate it with source FAST5 file unambiguously.

- Experience tells that check string (more precisely, it's hash-like part before the underscore) has length of 41. Assuming that this situation can vary from case to case, "fastQA5-sorter-v3-1a.py" will match check string that contains hash-like part with length at least 30 characters (or more). Also I have not seen uppercase letters in such hash-like character sequences, but it is better to assume and forsee their existance.

Algorithm of finding taxonomical annotation infornation for a particular FAST5 file sorting:

1. "fastQA5-sorter-v3-1a.py" tries to match the pattern described above with regex. Exact pattern:

    `[a-zA-Z0-9]{30,}_[0-9]+`
  
2. If mathing string is substring of the name of directory with results of taxonomic annotation, FAST5 file of our interest will be associated with this directory and will be sorted according to TSV file in it.

3. If there is no match (for example, the user renames his/her FAST5 files just after sequensing) whole FAST5 file name will be used to find corresponding directory with results of taxonomic annotation. For example, if FAST5 file is named `my_favorite_reads_324egf.fast5`, then `my_favorite_reads_324egf` will be used as check string. Again, a directory which name contains check string is considered as corresponding to this FAST5 file.

### Example:

Assuming you have already performed basecalling of file "some_reads.fast5" and have file "some_reads.fastq". You may act as follows to sort source FAST5 file:

Taxonomical annotation:
```
./prober-v1-12d.py some_reads.fastq
./barapost-v3-5c.py some_reads.fastq
```
Sorting:
```
./fastQA5-sorter-v3-1a.py some_reads.fast5
```

## Examples of usage in combination:

1. You can place all .py-files provided with by this toolkit in a directory that contains some FASTA and FASTQ files and run whole "pipeline" with default settings:

`./prober-v1-12d.py && ./barapost-v3-5c.py && ./fastQA5-sorter-v3-1a.py`

2. You can try these scripts on test dataset named `some_reads.fastq` (there are 4 reads):

`./prober-v1-12d.py some_reads.fastq -o some_outdir -g Escherichia,561+viruses,10239 -p 2 -b 2`

`./barapost-v3-5c.py some_reads.fastq -r some_outdir`

`./fastQA5-sorter-v3-1a.py some_reads.fastq -r some_outdir -o some_sorted_reads`
