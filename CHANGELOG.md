# "Barapost" changelog

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