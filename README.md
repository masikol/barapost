# barapost toolkit

## 1. prober.py

Version 1.2; 31.08.2019 edition

### DESCRIPTION:

"prober.py" -- script is designed for determinating the taxonomic position
of nucleotide sequences by blasting each of them and regarding the best hit.

The main goal of this script is to send a probing batch of sequences to BLAST server
and discover, what Genbank records can be downloaded and used for further processing
on your local machine by "barapost.py".

This script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.

Results of the work of this script are written to TSV files, that can be found in result directory:

1) There is a file named `"...acc_list.tsv"`. It contains accessions and names of Genbank records that
    can be used for further processing on your local machine by "barapost.py".

2) There is a file named `"...result.tsv"`. It contains full result of blasting.
    Results of barapost.py's work will be appended to this file.

FASTQ files processed by this script are meant to be processed afterwards by "barapost.py".


### Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- packet size (see `-p` option): 100 sequences;
- probing batch size (see `-b` option): 200 sequences;
- algorithm (see `-a` option): `megaBlast`;
- organisms (see `-g` option): full `nt` database, e.i. no slices;
- output directory (`-o` option): directory named `"barapost_result"`
  nested in current directory;

Default behavior of this script is to send certain batch (see `-b` option) of sequences to BLAST server.
It means that you should not process all your data by "prober.py' -- it would take long time.

Instead of this you should process some sequences by "prober.py" -- it will determine,
what Genbank records (genomes, if you want) are present in your data and then go to "barapost.py".

"barapost.py" will process the rest of you sequences in the same way like "prober.py", but on your local computer.
"barapost.py" uses 'blast+' toolkit for this purpose. It would be much faster.

Obviously, a probing batch cannot cover all variety of a data set,
so some sequences can be recognized as "unknown" while processing by "barapost.py".
But you always can run "prober.py" again on "unknown" sequences.

### OPTIONS:

```
    -h (--help) --- show help message;

    -f (--infile) --- input FASTQ or FASTA file. File can be gzipped.
            You can specify multiple input files with this option (see EXAMPLES #2);

    -d (--indir) --- directory which contains FASTQ or FASTA files meant to be processed.
            E.i. all FASTQ and FASTA files in this direcory will be processed; Files can be gzipped.

    -o (--outdir) --- output directory;

    -p (--packet-size) --- size of the packet, e.i. number of sequence to blast in one request.
            Value: integer number [1, 500]. Default value is 100;

    -a (--algorithm) --- BLASTn algorithm to use for aligning.
            Available values: 'megaBlast', 'discoMegablast', 'blastn'.
            Default is megaBlast;

    -g (--organisms) --- 'nt' database slices, e.i. organisms that you expect to see in result files.
            Format of value: 
              <organism1_name>,<organism1_taxid>+<organism2_name>,<organism2_taxid>+...
            See EXAMPLES #2 and #3 below.
            Spaces are not allowed. Number of organisms can be from 1 to 5 (5 is maximum).
            Default is: full 'nt' database, e.i. no slices.

    -b (--probing-batch-size) --- number of sequences that will be aligned on BLAST server.
            After that a local database will be builded according to results of probing blasting.
            More clearly: records-hits will be downloaded from Genbank and will be used
            as local database. Further blasting against this database will be performed
            on local machine with 'blast+' toolkit.
            Value: positive integer number. Default value is 200;
```

- More clearly, functionality of `-g` option is totally equal to "Organism" text boxes on this BLASTn page:
    'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome'.
- You can find your Taxonomy IDs here: 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi'.


### EXAMPLES:

  1) Process one file with default settings:

`./prober.py -f reads.fastq`

  2) Process a FASTQ file and a FASTA file with discoMegablast, packet size of 5 sequences.
Search only among Erwinia sequences:

`./prober.py -f reads_1.fastq.gz -f some_sequences.fasta -a discoMegablast -p 5 -g Erwinia,551`

  3) Process all FASTQ and FASTA files in directory named `some_dir`.
Search only among Escherichia and viral sequences:

`./prober.py -d some_dir -g Escherichia,561+viruses,10239 -o outdir`

## 2. barapost.py

Version 2.4; 31.08.2019 edition

### DESCRIPTION:

barapost.py -- script is designed for determinating the taxonomic position
of nucleotide sequences by blasting each of them with 'blastn' from 'blast+' toolkit
and regarding the best hit.

"barapost.py" is meant to be used just after 'prober.py'.

"barapost.py" downloads records-hits from Genbank according to results
of work of 'prober.py' script, builds an indexed local database which consists of
downloaded sequences, and continues aligning with 'blast+' toolkit in order to save time.

Script processes FASTQ and FASTA (as well as '.fastq.gz' and '.fasta.gz') files.

"barapost.py" writes it's results in the same TSV file as 'prober.py' does.

FASTQ files processed by this script are meant to be sorted afterwards by 'fastQA_sorted.py'.

`blast+` toolkit (including blastn, makeblastdb and makembindex) can
be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download);


### Default parameters:

Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- packet size (see '-p' option): 100 sequences;
- algorithm (see '-a' option): 'megaBlast';


### OPTIONS (* means mandatory option):

```
    -h (--help) --- show help message;

  * -r (--prober-result-dir) --- result directory genearted by script "prober.py"
        This is directory specified to "prober.py" by '-o' option.
        Or named 'prober_result' if you've ran 'prober.py' with default options.

    -f (--infile) --- input FASTQ or FASTA file. File can be gzipped.
        You can specify multiple input files with this option (see EXAMPLES #2);

    -d (--indir) --- directory which contains FASTQ or FASTA files meant to be processed.
        E.i. all FASTQ and FASTA files in this direcory will be processed; Files can be gzipped.

    -o (--outdir) --- output directory;

    -p (--packet-size) --- size of the packet, e.i. number of sequence to blast in one request.
        Value: integer number [1, 500]. Default value is 100;

    -a (--algorithm) --- BLASTn algorithm to use for aligning.
        Available values: 'megaBlast', 'discoMegablast', 'blastn'.
        Default is megaBlast;
```


### EXAMPLES:

1) Process one FASTQ file with default settings.
     File `reads.fastq` has been already processed by "prober.py".
     Results of "prober.py" work are in directory `prober_outdir`:

`./barapost.py -f reads.fastq -r prober_outdir`

2) Process FASTQ file and FASTA file with discoMegablast, packet size of 5 sequences.
     Files `reads.fastq.gz` and `another_sequences.fasta` have been already processed by "prober.py".
     Results of "prober.py" work are in directory `prober_outdir`:

`./barapost.py -f reads.fastq.gz -f another_sequences.fasta -a discoMegablast -r prober_outdir`

3) Process all FASTQ and FASTA files in directory named `some_dir`.
    All these files have been already processed by "prober.py".
    Results of "prober.py" work are in directory `prober_outdir`:

`/barapost.py -d some_dir -r prober_outdir`


## 3. fastQA_sorter.py
(fast**Q** and fast**A** sorter)

Version 2.1; 31.08.2019 edition

### DESCRIPTION:

fastQA_forter.py -- script designed for sorting FASTQ and FASTAfiles processed by barapost.py.

'fastQA_forter.py' is meant to be used just after "barapost.py".

- Separate FASTQ or FASTA files should be specified with -f(--infile) option (see EXAMPLES below).
- If -d (--indir) option is specified, all FASTQ and FASTA files in directory specified by this option will be sorted.
- If no separate FASTQ or FASTA files and not input directory is specified,
fastQA_sorter will process all FASTQ and FASTA files in current directory.


### Default parameters:

- all FASTQ and FASTA files in current directory will be processed;
- sorting sensitivity (see `-s` option): `"species"`;
- output directory (`-o` option): directory named `"fastQA_sorter_result_<date_and_time_of_run>"`
  nested in current directory;

### OPTIONS (* means mandatory option):

```
    -h (--help) --- show help message;

  * -r (--prober-result-dir) --- result directory genearted by script "prober.py"
        This is directory specified to "prober.py" by '-o' option.
        Or named 'prober_result' if you've ran 'prober.py' with default options.

    -f (--infile) --- input FASTQ or FASTA file (can be gzipped);

    -d (--indir) --- directory which contains FASTQ and/or FASTA files
        (files can be gzipped) meant to be sorted;

    -o (--outdir) --- output directory;

    -s (--sorting-sensitivity) --- sorting sensitivity;
        E.i. the lowest taxonomy rank that will be used in names of resut files;
        Available values: 'genus', 'species', 'strain';
```


### EXAMPLES:

1) Process one FASTQ file with default settings.
     File `reads.fastq` has been already processed by "barapost.py".
     Results of "barapost.py" work are in directory `prober_outdir`:

`./fastq_sorter.py -f reads.fastq.gz -r prober_outdir/`

2) Process a FASTQ file and a FASTA file, place results in `outdir` directory.
     Files `reads.fastq.gz` and `another_sequences.fasta` have been already processed by "barapost.py".
     Results of "barapost.py" work are in directory `prober_outdir`:

`./fastq_sorter.py -f reads_1.fastq.gz -f some_sequences_2.fasta -o outdir -r prober_outdir/`

3) Process all FASTQ and FASTA files in directory named `dir_with_seqs`. Sort by genus.
     All these files have been already processed by "barapost.py".
     Results of "barapost.py" work are in directory `prober_outdir`:

`./fastq_sorter.py -d dir_with_seqs -o outdir -r prober_outdir/ -s genus`


## Example of usage in combination:

You can try this scripts on test dataset named `some_reads.fastq` (there are 4 reads):

`./prober.py -f some_reads.fastq -o some_outdir -g Escherichia,561+viruses,10239`

`./barapost.py -f some_reads.fastq -r some_outdir`

`./fastQA_sorter.py -f some_reads.fastq -r some_outdir -o some_sorted_reads`
