# Barapost changelog

## 2025-02-05 edition.

### barapost-prober

A patch update.

Now, barapost-prober won't crush if it cannot parse a taxonomy ID from a GenBank page. It will arrogantly ignore such a record and give it empty taxonomy.

## Vesrion changes:

- barapost-prober: `1.24.a --> 1.24.b`

## 2023-05-26 edition.

### barapost-binning

- Now the script won't just terminate if it does not find some reads in the index. The script will just omit the reads.

### barapost-local

- Fix a bug in search for related replicons. The thing is, since Python 3.9 the `xml.etree.ElementTree.Element.getchildren` method is removed, and barapost-local used it. Now it does not, and the bug is gone.

- A note. Previously, Barapost was tested on `h5py` package version 2.10.0. However, its dependency, `numpy`, does not have attribute `typeDict` since version 1.21.0, and `h5py` 2.10.0 used this attribute. This used to cause Barapost failure (error message `AttributeError: module 'numpy' has no attribute 'typeDict'`). Now, Barapost is tested on `h5py` 3.8.0, and the inconsistency is gone.

## Vesrion changes:

- barapost-local: `3.18.e --> 3.18.f`
- barapost-binning: `4.9.b --> 4.9.c`

## 2022-08-02 edition.

- barapost-binning: fix a bug that didn't allow to recover taxonomy.

## Vesrion changes:

- barapost-binning: `4.9.a --> 4.9.b`

## 2021-11-19 edition.

- barapost-binning: this script now does not terminate with error if it cannot find the classification for an input sequence. Instead, it now places all these sequences into a `classification_not_found` file. This file is different from the `unknown` file: an `unknowm` file contains sequences, for which similarity was not found during the classification. And a `classification_not_found` contains sequences, which have not been classified at all.

## Vesrion changes:

- barapost-binning: `4.8.d --> 4.9.a`

## 2021-06-24 edition.

- barapost-prober: removed option `-e`, since it is the **program author's** email that [should](https://www.ncbi.nlm.nih.gov/books/NBK25497/) be sent to BLAST server within requests, not user's one.

## Vesrion changes:

- barapost-prober: `1.23.b --> 1.24.a`

## 2021-06-24 edition.

- barapost-local: fixed bug that would cause barapost-local to misinterpret user-defined path in `hits_to_download.tsv` as GenBank accession number;
- barapost-binnning: improved status messages: it now reports info e.g. "File #2 out of 23 is binned".

## Vesrion changes:

- barapost-local: `3.18.d --> 3.18.e`
- barapost-binning: `4.8.c --> 4.8.d`

## 2021-02-01 edition.

- barapost-local: yet another database managing bug fixed;

## Vesrion changes:

- barapost-local: `3.18.c --> 3.18.d`

## 2021-01-29 edition.

- addition to 2021-01-27 edition: non-breaking spaces are now replaced too;
- improved behavior of `barapost-local` and `barapost-binning` when they use assembly files generated by SPAdes and a5 (+ bug fix);

## Vesrion changes:

- barapost-local: `3.18.b --> 3.18.c`
- barapost-binning: `4.8.b --> 4.8.c`

## 2021-01-27 edition.

- all scripts: barapost now does not allow inappropriate characters (such as semisolons, parentheses or spaces) to appear in names of ouput files;

## Vesrion changes:

- barapost-prober: `1.23.a --> 1.23.b`
- barapost-local: `3.18.a --> 3.18.b`
- barapost-binning: `4.8.a --> 4.8.b`

## 2021-01-18 edition.

- barapost-binnig: `-z` option removed: now barapost-binnig does not compress output files;
- all: new format of log files;
- all: fixed error that would cause barapost-prober and barapost-local to fail if an input file is broken and utf-8 codec cannot decode it. Now the programs don't terminate execution and just continue working after this error;

## Vesrion changes:

- barapost-prober: `1.22.a --> 1.23.a`
- barapost-local: `3.17.a --> 3.18.a`
- barapost-binning: `4.7.b --> 4.8.a`

## 2021-01-11 edition.

- barapost-binning. Taxonomy recovering added: now if barapost-binning fails to find taxonomy for a sequence in `taxonomy.tsv`, it is able to recover it and continue binning;

## Vesrion changes:

- barapost-binning: `4.7.a --> 4.7.b`

## 2021-01-09 edition.

- Replaced DBM taxanomy files with plain TSV keping backward capability;

## Vesrion changes:

- barapost-prober: `1.21.j --> 1.22.a`
- barapost-local: `3.16.g --> 3.17.a`
- barapost-binning: `4.6.g --> 4.7.a`

## 2020-12-28 edition.

- Added handling of error "db type is dbm.gnu, but the module is not available" and suggestion for fixing;
- A little bit of refactoring;

## Vesrion changes:

- barapost-prober: `1.21.i --> 1.21.j`
- barapost-local: `3.16.f --> 3.16.g`
- barapost-binning: `4.6.f --> 4.6.g`

## 2020-11-16 edition.

- barapost-prober: fixed error that would cause prober to compose packets improperly with `-c 1` specified if fastq file is being processed;
- barapost-local: fixed a bug that would cause barapost-local to improperly spread sequences among threads;

## Vesrion changes:

- barapost-prober: `1.21.h --> 1.21.i`
- barapost-local: `3.16.e --> 3.16.f`

## 2020-11-13 edition.

- barapost-local, barapost-prober: made parsing file `hits_to_download` less prone to errors;
- barapost-binning: fixed a bug that would cause gzip to fail to compress some files dure to naming errors, such as brackets of plus-characters (+) in file names;

## Vesrion changes:

- barapost-prober: `1.21.g --> 1.21.h`
- barapost-local: `3.16.d --> 3.16.e`
- barapost-binning: `4.6.e --> 4.6.f`


## 2020-10-30 edition.

- barapost-local: fixed bug that would cause barapost-binning to exit with `taxonomy error 994` on lambda phage sequence.
- barapost-binning: naming of output files named after SPAdes or a5 assembly file was changed: assembler name now goes after contig name;
- barapost-local: fixed a bug that would cause BLAST+ to be unable to create a database due to non-ASCII characters in fasta titles of reference sequences;
- barapost-binnging: fixed messed stdout on parallel binning;

## Vesrion changes:

- barapost-barapost: `3.16.b --> 3.16.c --> 3.16.d`
- barapost-binning: `4.6.c --> 4.6.d --> 4.6.e`

## 2020-10-29 edition.

- increased durability of prober and barapost-local to instability of NCBI servers.

## Vesrion changes:

- barapost-prober: `1.21.f --> 1.21.g`
- barapost-barapost: `3.16.a --> 3.16.b`

## 2020-10-23 edition.

- barapost-prober: fixed a bug that would cause prober to stop after the first iteration if `-c 1` is specified.

## Vesrion changes:

- barapost-prober: `1.21.e --> 1.21.f`

## 2020-09-21 edition.

- barapost-local now in advance adds sequence of nanopore lambda phage control sequence (DNA-CS) to reference database;
- barapost-binning: fixed bug that would cause barapost-binning to crash while renaming extant result directory;
- barapost-local: minor cosmetical improvements;

## Vesrion changes:

- barapost-local: `3.15.b--> 3.15.c --> 3.16.a`;
- barapost-binning: `4.6.b--> 4.6.c`;

## 2020-09-03 edition.

- barapost-local: fixed `-s`-assitiated bug that would cause barapost-local to crash in absence of pre-created `-r` directory with `hits_to_download.tsv` file.

## Version changes:

- barapost-local: `3.15.a --> 3.15.b`;

## 2020-09-02 edition.

- barapost-local: algorithm of searching for related replicons improved: it no longer tries to insanely download over 7000 transcripts (like for this [fungus](https://www.ncbi.nlm.nih.gov/nuccore?LinkName=biosample_nuccore&from_uid=7457167)). Instead, barapost now requests for GenBank records given ID of corresponding record in the Assembly database and downloads only complete genomes/chromosomes or at least scaffolds.
- barapost-local: entertaining conda-like spinning thing added to indicate that the script is actually working while searching for related replicons.
- barapost-local, barapost-binning: option renamed: `--taxannot-resdir --> --annot-resdir`. Much more convenient.

## Version changes:

- barapost-local: `3.14.l --> 3.15.a`;
- barapost-binning: `4.6.a --> 4.6.b`;

## 2020-09-01 edition.

Small fix: removed code left from debug.

## Version changes:

- barapost-prober: `1.21.d --> 1.21.e`;
- barapost-local `3.14.k --> 3.14.l`;

## 2020-06-29 edition
- barapost-local: `-n` flag added, which cause the script not to output "trash" sequences;
- all: scripts print Python interpreter version at startup;
- prober: fixed bug that would cause prober incorrectly split probing batch into packets;

## Version changes:

- barapost-prober: `1.21.c --> 1.21.d`;
- barapost-local `3.14.j --> 3.14.k`;
- barapost-binning: `4.5.g --> 4.6.a`;

## 2020-06-26 edition

- barapost-local: (`i`) fixed bug that would cause barapost-local not to found assembly files specified with `-l` option;
- barapost-local: (`j`) fixed bug that would cause blastn not to accept `-l` fasta records with spaces in headers;

## Version changes:

- barapost-local `3.14.h --> 3.14.i --> 3.14.j`;

## 2020-06-12 edition

- README updated;
- barapost-local: `-l` bug fixed;
- minor renaming changes in source files;

## Version changes:

- barapost-prober: `1.21.b --> 1.21.c`;
- barapost-local `3.14.f --> 3.14.g --> 3.14.h`;
- barapost-binning: `4.5.f --> 4.5.g`;

## 2020-06-11

All scripts are renamed.

- prober.py --> barapost-prober.py
- barapost.py --> barapost-local.py
- fastQA5-sorter.py --> barapost-binning.py

## 2020-06-08 edition

- barapost-local: fixed bug that did not allow to download many (200+) genomes due to too long URL;

## Version changes:

- barapost-local `3.14.e --> 3.14.f`;

## 2020-06-05 edition

- barapost-binning. Bug fix: binning fast5 files according to classification made on fasta files no longer fails;
- barapost-prober, barapost-local: minor modifications;

## Version changes:

- barapost-prober: `1.21.a --> 1.21.b`;
- barapost-local `3.14.d --> 3.14.e`;
- barapost-binning: `4.5.e --> 4.5.f`;

## 2020-05-26 edition

- barapost-prober: barapost-prober now will split packet into two if BLAST error encounters, leaving sequences intact. barapost-prober now prunes sequences twofold only if packet consists of the only sequence;
- barapost-local: one more minor multiprocessing bug fixed;

## Version changes:

- barapost-prober: `1.20.b -> 1.21.a`;
- barapost-local: `3.14.c --> 3.14.d`;

## 2020-05-23

- barapost-prober: bug emerged due to embedding `-c` option fixed;
- barapost-local: multiprocessing bug fixed;

## Version changes:

1. barapost-prober: `1.20.a -> 1.20.b`;
2. barapost-local: `3.14.b --> 3.14.c`;

## 2020-05-21

- barapost-local: minor performance fix;

## Version changes:

1. barapost-local: `3.14.a --> 3.14.b`;

## 2020-05-21

- barapost-prober: `-c` option added, enabling better performance if input sequences are binned by length within input file;
- barapost-prober: "no-hits" error fixed;

## Version changes:

1. barapost-prober: `1.19.b -> 1.20.a`;

## 2020-05-13

- barapost-local: searching for related replicons fixed. barapost-local no longer downloads identical sequences from GenBank and RefSeq;
- barapost-local now prints amount of processed files whilst wokring;

## Version changes:

1. barapost-local: `3.13.b -> 3.14.a`;

## 2020-05-10

- barapost-local now firstly tries to download reference sequences with GNU wget, and if there are no wget, downloads them with standart Python tools;
- barapost-prober: minor durability fix;

## Version changes:

1. barapost-prober: `1.19.a -> 1.19.b`;
2. barapost-local: `3.13.a -> 3.13.b`;

## 2020-05-06 edition.

- barapost-prober: behaviour `-x` is modified: barapost-prober now prunes sequences from both ends instead of cutting off 3'-end and leaving 5'-end;
- barapost-prober: minor durability fix;

### Version changes:

1. barapost-prober: `1.18.e -> 1.19.a`;

## 2020-05-01 edition.

- barapost-prober, barapost-local: got rid of GI numbers in `hits_to_download.tsv`. Backward compatibility with old version of such files is keeped.
- barapost-local: a possibility to directly specify accessions of GenBank records meant to be downloaded and included in database added (`-s` option).
- all: minor performance fixes;

### Version changes:

1. barapost-prober: `1.18.d -> 1.18.e`;
2. barapost-local: `3.12.a --> 3.13.a`;
3. barapost-binning: `4.5.d --> 4.5.e`;

## 2020-04-24 edition.

- barapost-prober: minor fixes. Network error handling and saving BLAST results in .txt files are modified;

### Version changes:

1. barapost-prober: `1.18.c -> 1.18.d`;

## 2020-04-07 edition.

- barapost-local: searching for related replicons modified. Restriction to 20 first links on page removed;

### Version changes:

1. barpost: `3.11.j -> 3.12.a`;

## 2020-04-06 edition.

- barapost-prober, barapost-local: value of HIT_NAME' field in `classification.tsv` is changed: full definition of GenBank record is now specified in the file instead of `<Genus> <species>`;

### Version changes:

1. barapost-prober: `1.18.b --> 1.18.c`
2. barpost: `3.11.i -> 3.11.j`;

## 2020-03-18 edition.

- barapost-prober: taxonomy parsing fixed;
- barapost-prober, barapost-local: bug fixed -- resuming processing fastq file(s) works fine now;

### Version changes:

1. barapost-prober: `1.18.a --> 1.18.b`
2. barpost: `3.11.h -> 3.11.i`;

## 2020-03-07 edition.

- barapost-prober: missing sequences on resumption and changing packet size fixed;

### Version changes:

1. barapost-prober: `1.17.g --> 1.18.a`

## 2020-03-03 edition.

- all: restriction on duplicated basenames of input files added;

### Version changes:

1. barapost-prober: `1.17.f --> 1.17.g`
2. barapost-local: `3.11.g --> 3.11.h`;
3. barapost-binning: `4.5.c --> 4.5.d`;

## 2020-03-02 edition.

- barapost-local, barapost-prober: GET requests durability "sleep"-bug fixed;

### Version changes:

1. barapost-prober: `1.17.e --> 1.17.f`
2. barapost-local: `3.11.f --> 3.11.g`;

## 2020-02-28 afternoon edition.

- all: taxonomy parsing improved, absence of genus name allowed;
- barapost-local, barapost-prober: GET requests durability bug fixed;
- barapost-binning excludes characters casuing gzip errors from names of result files;

### Version changes:

1. barapost-prober: `1.17.d --> 1.17.e`
2. barapost-local: `3.11.e --> 3.11.f`;
3. barapost-binning: `4.5.b --> 4.5.c`;

## 2020-02-28 night edition.

- barapost-local: blastn optimization options added;
- all: bug leading to improper naming of binned files fixed;

### Version changes:

1. barapost-prober: `1.17.c --> 1.17.d`
2. barapost-local: `3.11.d --> 3.11.e`;
3. barapost-binning: `4.5.a --> 4.5.b`;

## 2020-02-27 edition.

- barapost-prober, barapost-local: network errors handling (sleep 30 sec);
- barapost-prober, barapost-local: genus name parsing fixed;

### Version changes:

1. barapost-prober: `1.17.b --> 1.17.c`
2. barapost-local: `3.11.c --> 3.11.d`;

## 2020-02-26 evening edition.

- barapost-prober, barapost-local: "syntetic construct" fix;
- barapost-prober, barapost-local: performance fix;

### Version changes:

1. barapost-prober: `1.17.a --> 1.17.b`
2. barapost-local: `3.11.b --> 3.11.c`;

## 2020-02-26 morinig edition.

- barapost-local: -d option bug fixed;

### Version changes:

1. barapost-local: `3.11.a --> 3.11.b`;

## 2020-02-26 night edition.

- all: new taxonomy system embedded;
- barapost-binning: binning sensitivity can now vary from domain to species;

### Version changes:

1. barapost-prober: `1.16.b --> 1.17.a`;
2. barapost-local: `3.10.b --> 3.11.a`;
3. barapost-binning: `4.4.b --> 4.5.a`

## 2020-02-23 edition.

- barapost-binning: separate filter for alignment identity and coverage added;

### Version changes:

1. barapost-binning: `4.4.a --> 4.4.b`

## 2020-02-22 edition.

- barapost-binning: alignment identity and coverage filters added;
- barapost-prober, barapost-local: fasta reading bug fixed;

### Version changes:

1. barapost-prober: `1.16.a --> 1.16.b`;
2. barapost-local: `3.10.a --> 3.10.b`;
3. barapost-binning: `4.3.a --> 4.4.a`

## 2020-02-21 edition.

- all: code restructirized, adding directory with executable scripts to PATH variable is necessary;
- barapost-prober, barapost-local: `-a` option syntax changed;
- barapost-binning: `-s` and `-z` options syntax changed;
- several bug fixes;

### Version changes:

1. barapost-prober: `1.15.b --> 1.16.a`;
2. barapost-local: `3.9.a --> 3.10.a`;
3. barapost-binning: `4.2.b --> 4.3.a`

## 2020-01-29 evening edition.

- barapost-local.py now downloads all replicons related to records "discovered" by barapost-prober.py (other chromosomes, plasmids) and adds them to database;
- barapost-prober, barapost-local: file renaming procedure fixed;

### Version changes:

1. barapost-prober: `1.15.a --> 1.15.b`;
2. barapost-local: `3.8.a --> 3.9.a`; 

## 2020-01-29 edition.

- barapost-binning: default `-q` value switched to 10;

### Version changes:

1. barapost-binning: `4.2.a --> 4.2.b`;

## 2020-01-28 edition.

- barapost-prober, barapost-local: average read quality is now calculated in a correct manner (by mean error propability);
- barapost-prober, barapost-local: taxonomy parsing modified;
- barapost-binning: default `-q` value switched to 15;
- barapost-binning: parallel progress bar fixed;
- barapost-binning: correct handling of custom sequence IDs in database;

### Version changes:

1. barapost-prober: `1.14.e --> 1.15.a`;
2. barapost-local: `3.7.g --> 3.8.a`;
3. barapost-binning: `4.1.b --> 4.2.a`;

## 2020-01-26 edition.

- barapost-local: shelve's "Service temporarily unavailable" bug fixed;

### Version changes:

1. barapost-local: `3.7.g --> 3.7.e`;

## 2020-01-15 edition.

- barapost-local: accessions with '.1' terminus returned by blastn hit taxonomy properly;

### Version changes:

1. barapost-local: `3.7.f --> 3.7.g`;

## 2020-01-13 edition.

- barapost-prober: 'blastsrv4.REAL'-error handling on resumption fixed;

### Version changes:

1. barapost-prober: `1.14.d --> 1.14.e`;

## 2020-01-10 edition.

- barapost-local, barapost-prober: fmt_seq_id function fixed;
- barapost-local: more secure database index handling;

### Version changes:

1. barapost-prober: `1.14.c --> 1.14.d`;
2. barapost-local: `3.7.e --> 3.7.f`;

## 2019-12-29 edition.

- barapost-local, barapost-binning: task distributing among processes improved;

### Version changes:

1. barapost-local: `3.7.d --> 3.7.e`;
2. barapost-binning: `4.1.a --> 4.1.b`

## 2019-12-26 evening edition.

- barapost-prober, barapost-local: shelve open-keys bugs fixed;

### Version changes:

1. barapost-prober: `1.14.b --> 1.14.c`;
2. barapost-local: `3.7.b --> 3.7.c`;

## 2019-12-26 evening edition.

- barapost-local: shelve open mode bug fixed;

### Version changes:

1. barapost-local: `3.7.b --> 3.7.c`;

## 2019-12-26 edition.

- barapost-prober, barapost-local: lineage downloading optimized;

### Version changes:

1. barapost-prober: `1.14.a --> 1.14.b`;
2. barapost-local: `3.7.a --> 3.7.b`;


## 2019-12-22 edition. Classification algorithm modified

- Full lineages are now used for classification instead of just hit definitions;
- If several best hits have equal Bit scores, lowest common ancestor (LCA) of these hits will be determined and used for classification;

### Version changes:

1. barapost-prober: `1.13.c --> 1.14.a`;
2. barapost-local: `3.6.d --> 3.7.a`;
3. barapost-binning: `4.0.a --> 4.1.a`;

## 2019-12-19 edition

### Noteworthy changes:

1. barapost-prober: '-x' sequence pruning fixed;
2. barapost-prober: resuming bug fixed;

### Version changes:

1. barapost-prober: `1.13.b --> 1.13.c`;

## 2019-12-12 edition

### Noteworthy changes:

1. barapost-prober: seq IDs formatting fixed;

### Version changes:

1. barapost-prober: `1.13.a --> 1.13.b`;

## 2019-12-10 evening edition (barapost-prober's '-x' option added)

### Noteworthy changes:

1. barapost-prober: '-x' ('--max-seq-len') option enabled. It is now possible to prune sequences that barapost-prober sends to NCBI in order to spare NCBI BLAST servers;

### Version changes:

1. barapost-prober: `1.12.m --> 1.13.a`;

## 2019-12-10 edition (bug fix)

### Noteworthy changes:

1. barapost-prober, barapost-local: error while processing .fasta and .fasta.gz files fixed;

### Version changes:

1. barapost-prober: `1.12.k --> 1.12.m`;
2. barapost-local: `3.6.c --> 3.6.d`;

## 2019-12-07 edition (parallel binning enabled)

### Noteworthy changes:

1. barapost-binning: parallel FASTA and FASTQ binning enabled;
2. barapost-binning: parallel `-u` index creating enabled;
3. barapost-binning: `-z (--gzip)` option added;

Parallel FAST5 binning is not embedded and perhaps won't be -- it gives no performance profit. The point is that writing to FAST5 files takes much more time than 'calculating'. Thus threads mostly just stay in a queue for writing rather than doinig their work.

### Version changes:

1. barapost-prober: `1.12.l --> 1.12.k`;
2. barapost-local: `3.6.b --> 3.6.c`;
3. barapost-binning: `3.4.b. --> 4.0.a`;
4. summarizer: `1.1.b --> 1.1.c`;

## 2019-11-19 edition (bug-fixing and performance update):

### Noteworthy changes:

1. barapost-local: "inc_lock 1-thr" bug fixed;
2. barapost-local: rebuild-db "no -l no acc" bug fixed;
3. barapost-prober: waiting for NCBI response no more makes a mess in terminal;

### Version changes:

1. barapost-prober: `1.12.k --> 1.12.l`;
2. barapost-local: `3.6.a --> 3.6.b`;
3. barapost-binning: `3.4.a --> 3.4.b`;
4. summarizer: `1.1.a --> 1.1.b`;

## 2019-11-17 edition:

### Noteworthy changes:

- barapost-binning: singleFAST5 processing enabled;
- barapost-local, barapost-binning, summarizer: concise status bar embedded;

### Version changes:

1. barapost-prober: `1.12.j --> 1.12.k`;
2. barapsot: `3.5.j --> 3.6.a`;
3. barapost-binning: `3.3.d --> 3.4.a`;
4. summarizer: `1.0.d --> 1.1.a`;

## Unfortunately, I hadn't been recording the changelog before 2019-11-17 :(
