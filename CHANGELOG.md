# "Barapost" changelog

## 2020-02-27 edition.

- prober, barapost: network errors handling (sleep 30 sec);
- prober, barapost: genus name parsing fixed;

### Version changes:

1. prober: `1.17.b --> 1.17.c`
2. barapost: `3.11.c --> 3.11.d`;

## 2020-02-26 evening edition.

- prober, barapost: "syntetic construct" fix;
- prober, barapost: performance fix;

### Version changes:

1. prober: `1.17.a --> 1.17.b`
2. barapost: `3.11.b --> 3.11.c`;

## 2020-02-26 morinig edition.

- barapost: -d option bug fixed;

### Version changes:

1. barapost: `3.11.a --> 3.11.b`;

## 2020-02-26 night edition.

- all: new taxonomy system embedded;
- sorter: sorting sensitivity can now vary from domain to species;

### Version changes:

1. prober: `1.16.b --> 1.17.a`;
2. barapost: `3.10.b --> 3.11.a`;
3. sorter: `4.4.b --> 4.5.a`

## 2020-02-23 edition.

- sorter: separate filter for alignment identity and coverage added;

### Version changes:

1. sorter: `4.4.a --> 4.4.b`

## 2020-02-22 edition.

- sorter: alignment identity and coverage filters added;
- prober, barapost: fasta reading bug fixed;

### Version changes:

1. prober: `1.16.a --> 1.16.b`;
2. barapost: `3.10.a --> 3.10.b`;
3. sorter: `4.3.a --> 4.4.a`

## 2020-02-21 edition.

- all: code restructirized, adding directory with executable scripts to PATH variable is necessary;
- prober, barapost: `-a` option syntax changed;
- sorter: `-s` and `-z` options syntax changed;
- several bug fixes;

### Version changes:

1. prober: `1.15.b --> 1.16.a`;
2. barapost: `3.9.a --> 3.10.a`;
3. sorter: `4.2.b --> 4.3.a`

## 2020-01-29 evening edition.

- barapost.py now downloads all replicons related to records "discovered" by prober.py (other chromosomes, plasmids) and adds them to database;
- prober, barapost: file renaming procedure fixed;

### Version changes:

1. prober: `1.15.a --> 1.15.b`;
2. barapost: `3.8.a --> 3.9.a`; 

## 2020-01-29 edition.

- sorter: default `-q` value switched to 10;

### Version changes:

1. sorter: `4.2.a --> 4.2.b`;

## 2020-01-28 edition.

- prober, barapost: average read quality is now calculated in a correct manner (by mean error propability);
- prober, barapost: taxonomy parsing modified;
- sorter: default `-q` value switched to 15;
- sorter: parallel progress bar fixed;
- sorter: correct handling of custom sequence IDs in database;

### Version changes:

1. prober: `1.14.e --> 1.15.a`;
2. barapost: `3.7.g --> 3.8.a`;
3. sorter: `4.1.b --> 4.2.a`;

## 2020-01-26 edition.

- barapost: shelve's "Service temporarily unavailable" bug fixed;

### Version changes:

1. barapost: `3.7.g --> 3.7.e`;

## 2020-01-15 edition.

- barapost: accessions with '.1' terminus returned by blastn hit taxonomy properly;

### Version changes:

1. barapost: `3.7.f --> 3.7.g`;

## 2020-01-13 edition.

- prober: 'blastsrv4.REAL'-error handling on resumption fixed;

### Version changes:

1. prober: `1.14.d --> 1.14.e`;

## 2020-01-10 edition.

- barapost, prober: fmt_seq_id function fixed;
- barapost: more secure database index handling;

### Version changes:

1. prober: `1.14.c --> 1.14.d`;
2. barapost: `3.7.e --> 3.7.f`;

## 2019-12-29 edition.

- barapost, sorter: task distributing among processes improved;

### Version changes:

1. barapost: `3.7.d --> 3.7.e`;
2. sorter: `4.1.a --> 4.1.b`

## 2019-12-26 evening edition.

- prober, barapost: shelve open-keys bugs fixed;

### Version changes:

1. prober: `1.14.b --> 1.14.c`;
2. barapost: `3.7.b --> 3.7.c`;

## 2019-12-26 evening edition.

- barapost: shelve open mode bug fixed;

### Version changes:

1. barapost: `3.7.b --> 3.7.c`;

## 2019-12-26 edition.

- prober, barapost: lineage downloading optimized;

### Version changes:

1. prober: `1.14.a --> 1.14.b`;
2. barapost: `3.7.a --> 3.7.b`;


## 2019-12-22 edition. Classification algorithm modified

- Full lineages are now used for classification instead of just hit definitions;
- If several best hits have equal Bit scores, lowest common ancestor (LCA) of these hits will be determined and used for classification;

### Version changes:

1. prober: `1.13.c --> 1.14.a`;
2. barapost: `3.6.d --> 3.7.a`;
3. sorter: `4.0.a --> 4.1.a`;

## 2019-12-19 edition

### Noteworthy changes:

1. prober: '-x' sequence pruning fixed;
2. prober: resuming bug fixed;

### Version changes:

1. prober: `1.13.b --> 1.13.c`;

## 2019-12-12 edition

### Noteworthy changes:

1. prober: seq IDs formatting fixed;

### Version changes:

1. prober: `1.13.a --> 1.13.b`;

## 2019-12-10 evening edition (prober's '-x' option added)

### Noteworthy changes:

1. prober: '-x' ('--max-seq-len') option enabled. It is now possible to prune sequences that prober sends to NCBI in order to spare NCBI BLAST servers;

### Version changes:

1. prober: `1.12.m --> 1.13.a`;

## 2019-12-10 edition (bug fix)

### Noteworthy changes:

1. prober, barapost: error while processing .fasta and .fasta.gz files fixed;

### Version changes:

1. prober: `1.12.k --> 1.12.m`;
2. barapost: `3.6.c --> 3.6.d`;

## 2019-12-07 edition (parallel sorting enabled)

### Noteworthy changes:

1. sorter: parallel FASTA and FASTQ sorting enabled;
2. sorter: parallel `-u` index creating enabled;
3. sorter: `-z (--gzip)` option added;

Parallel FAST5 sorting is not embedded and perhaps won't be -- it gives no performance profit. The point is that writing to FAST5 files takes much more time than 'calculating'. Thus threads mostly just stay in a queue for writing rather than doinig their work.

### Version changes:

1. prober: `1.12.l --> 1.12.k`;
2. barapost: `3.6.b --> 3.6.c`;
3. sorter: `3.4.b. --> 4.0.a`;
4. summarizer: `1.1.b --> 1.1.c`;

## 2019-11-19 edition (bug-fixing and performance update):

### Noteworthy changes:

1. barapost: "inc_lock 1-thr" bug fixed;
2. barapost: rebuild-db "no -l no acc" bug fixed;
3. prober: waiting for NCBI response no more makes a mess in terminal;

### Version changes:

1. prober: `1.12.k --> 1.12.l`;
2. barapost: `3.6.a --> 3.6.b`;
3. sorter: `3.4.a --> 3.4.b`;
4. summarizer: `1.1.a --> 1.1.b`;

## 2019-11-17 edition:

### Noteworthy changes:

- sorter: singleFAST5 processing enabled;
- barapost, sorter, summarizer: concise status bar embedded;

### Version changes:

1. prober: `1.12.j --> 1.12.k`;
2. barapsot: `3.5.j --> 3.6.a`;
3. sorter: `3.3.d --> 3.4.a`;
4. summarizer: `1.0.d --> 1.1.a`;