# GAPFILLER

A de novo local assembler for paired reads.


## Overview

GapFiller is a seed-and-extend local assembler to fill the gap within paired reads.
A set of coupled reads (seeds) is required, as well as a set of reads to fill the gap between the
two mates.
The only restriction on the reads library is that seeds are inner oriented (paired-reads) or 
outer oriented (mate-pairs).
Insert size refers to the largest distance between the seeds. 
GapFiller works better on short inserts (i.e., less than 1 kbp), however in principle it can be 
used to assemble larger ones (i.e., mate-pairs inserts).
GapFiller can be used whenever a sequence is to be assembled starting from reads lying on its 
ends, provided a loose estimate of sequence length.

## Compiling

GapFiller is written in C++, to compile the package you need the g++ compiler.
The LIBBOOST library (version >= 1.40) and the ZLIB library are required. 

To configure and build the package, type the following command lines:

```
$ unzip /path/to/gapfiller.zip
$ cd /path/to/gapfiller
$ ./configure; make
```

## Running

Move to the source folder:

```
$ cd /path/to/gapfiller/src
```

and enter the following command to get a list of parameters with a brief description:

```
$ ./GapFiller --help
```


### Basic command

GapFiller requires the input reads to be provided (in FASTA/FASTQ or SAM/BAM format) as well as
the insert size.
For example, the following command line

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION
```

produces three files: 

- "GapFiller_output.fasta"        : contigs obtained by filling the gap between pairs of sequences
                                    belonging to "seeds_1.fasta" and "seeds_2.fasta".
- "GapFiller_output_trash.fasta"  : contigs for which the seed mate was not found.
- "GapFiller_output.stats"        : statistics on the extensions.

Default behavior is to consider paired reads as inner-oriented (paired-ends). If outer-oriented 
sequences are used instead (mate-pairs), it is possible to specify it with option --mate-pairs:

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --mate-pairs
```
	
The above command will reverse-complement reads of "seeds_1.fasta" instead of "seeds_2.fasta".

### Fasta/fastq input specification

Seed files provided with --seed1 and --seed2 must contain pairs of sequences in the same order.
That is, the n-th paired read is constituted by the n-th read of seed 1 file and the n-th 
read of seed 2 file, respectively.

To provide a different set of reads to be used for assembly, type the command

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --query reads.fasta
```
	
The above command will produce a set of contigs produced by filling the sequence between seed 
reads from files "seeds_1.fasta" and "seeds_2.fasta" using reads from "reads.fasta".

Seed reads are NOT used for assembly if option --query is provided. 
To use reads from seed files also for assembly, type the following command

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --query reads.fasta --query seeds_1.fasta \
	--query seeds_2.fasta
```

Options --seed1, --seed2, and --query can be repeated several times.
Supported file formats are fasta and fastq, possibly compressed with gzip or bzip2.

### Sam/bam input specification

Alignment file can be provided as input by typing the following 

```
$ ./GapFiller --seed-sam seeds.sam --seed-ins AVG_INSERT_SIZE --seed-var INSERT_SIZE_VARIATION 
```

The only request is that "seeds.sam" is SORTED BY NAME (e.g. with samtools sort -n).
To use a different set of reads for assembly, type

```
$ ./GapFiller --seed-sam seeds.sam --seed-ins AVG_INSERT_SIZE --seed-var INSERT_SIZE_VARIATION \
	--query-sam reads.sam
```

Reads within "reads.sam" can be provided in any order.

It is possible to combine different file formats for query reads. Namely, both the following commands will work:

```
$ ./GapFiller --seed1 seeds_1.fasta --seed1 seeds_2.fasta --seed-ins AVG_INSERT_SIZE 
	--seed-var INSERT_SIZE_VARIATION --query-sam reads.sam --query reads.fasta
```

```
$ ./GapFiller --seed-sam seeds.sam --seed-ins AVG_INSERT_SIZE --seed-var INSERT_SIZE_VARIATION \
	--query-sam reads.sam --query reads.fasta
```

It is NOT possible to specify different seed file formats (e.g., combine fasta and sam).

### Output files

It is possible to specify a different output file name with option --output-prefix as follows:

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --output-prefix contigs
```
	
The above command will produce "contigs.fasta", "contigs_trash.fasta", and "contigs.stats".

Each contig in "contigs.fasta" reports the flag "pairFound", meaning that the paired read has been 
successfully filled. 
Contigs printed to "contigs_trash.fasta" do not correspond to insert sequences because the seed mate 
was not found. There are various reasons why this happened and contigs are flagged accordingly with
one of the following:

- noMoreExt 	: extension is halted because no more reads are available (low coverage)
- repSeq 		: extension is halted because a bifurcation occurs (exiting from a repeat)
- lengthExceed 	: maximum insert size is reached without finding seed 2
- readCycle		: a read has been re-used within the same contig

To possibly re-use reads within the same contig, add option --no-read-cycle.
In this case, no contig will be flagged as "readCycle".

"contigs.fasta" and "contigs_trash.fasta" be compressed with gzip (bzip2) by setting the option 
--gz (--bz2).
 
File "contigs.stats" contains statistics on the extensions:

- read_extended          : number of seeds considered for extension, in slots of 10000
- mean_extension_length  : average number of extension steps
- StopbyNoMoreExt        : number of contigs flagged "noMoreExt"
- StopByRepSeq           : number of contigs flagged "repSeq"
- StopByLengthExceed     : number of contigs flagged "lengthExceed" 
- StopByPairFound        : number of contigs flagged "pairFound"
- StopbyReadCycle        : number of contigs flagged "readCycle"

### Parameters for hashing

The hash table used for reads fingerprint depends on parameter k.
Set option --k to change hash size.
The fingerprint is computed on a subsequence of each read, called block. To change block length 
set option --block-length.
For example:

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --k 10 --block-length 20
```
	
Setting k to a lower value (default: 12) decreases the hash table size but increases the search 
time. 

Block length is the size of the subsequence on which the fingerprint is computed. Overlaps are 
searched among reads with the same fingerprint. Because the hash function has low false positive 
rate, with high probability reads with the same fingerprint share a block-length subsequence.
For this reason, a larger value of block length will increase specificity but can reduce 
sensitivity (reads that overlap for a large fraction and with few mismatches may not be found). 

### Parameters for extension

Depending on overlap length and similarity among reads, GapFiller can produce different results.
If reads are short or high error-affected, it may be preferable to let the extension phase to be
less strict.
For example, the following command

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --overlap 25 --mismatch-rate 10
```

will consider two reads as overlapping if the overlap length is at least 25 bp long (default: 30)
and if the number of mismatches among the two sequences is at most 10% of overlap length (default:
5%).

If either coverage is high or overlap length has been set to a low value, it can be convenient 
to increase the number of overlapping reads required to safely extend a contig. 

With the following command

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --overlap 25 --mismatch-rate 10 --extThreshold 5
```
	
at least 5 reads (default: 2) reads are required to extend the current contig. The consensus 
sequence is computed only in the positions covered by at least 5 overlapping reads.

### Layout computation

For further analyses, it is possible to specify the creation of a layout file that stores 
information on reads used to assemble each "pairFound" contig.
Given a contig C, its layout contains information on ID, position, length, and strand of each read
used to assemble C.

The following command:

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --store-layout
```
	
will produce two additional files: 

- "GapFiller_output.layout"       : binary file containing non-redundant information on contigs
                                    layout.
- "GapFiller_output.temp_layout"  : temporary file used to compute "GapFiller_output.layout".

Option --store-layout is incompatible with --query and --query-sam.
This means that the seed reads are exactly the reads used for contig assembly.

### Other parameters

For small test the following options may be useful: --limit and --verbose.
The former is the number of pairs considered (contigs assembly is halted as soon as --limit is
reached) whereas the latter prints detailed information on STDOUT.
--verbose option should be used in combination with --limit and STDOUT should be redirected to 
file:

```
$ ./GapFiller --seed1 seeds_1.fasta --seed2 seeds_2.fasta --seed-ins AVG_INSERT_SIZE \
	--seed-var INSERT_SIZE_VARIATION --limit 100 --verbose > GapFiller_output.STDOUT
```

## Contact

For feedback or questions about this repository, please contact [Francesca Nadalin](mailto:francesca@ebi.ac.uk).
