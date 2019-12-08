# "Barapost" changelog

## 2019-12-07 edition (parallel sorting enabled)

### Noteworthy changes:

1. sorter: parallel FASTA and FASTQ sorting enabled;
2. sorter: parallel `-u` index creating enabled;
3. sorter: `-z (--gzip)` option added;

Parallel FAST5 sorting is not embedded and perhaps won't be -- it gives no performance profit. The point is that writing to FAST5 files takes much more time than 'calculating'. Thus threads mostly just stay in a queue for writing rather than doinig their work.

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