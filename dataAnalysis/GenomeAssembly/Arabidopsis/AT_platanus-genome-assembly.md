---
title: "Short reads assembly using Platanus"
layout: single
author: Arun Seetharam
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Platanus

This tutorial covers the short read genome assembly for _Arabidopsis thaliana_ using `platanus` genome assembler. This assembly program can use both paired-end and mate-pair data, progressively, for contigging and for scaffolding to generate an assembly. This was initially designed for highly heterozygous genomes, but performs equally well for inbred lines as well. For more information about this assembly program, checkout the its [publication](https://genome.cshlp.org/content/24/8/1384.long).

## Dataset

For this tutorial, we will use _Arabidopsis thaliana_, Ler strain, short read dataset, from the BioProject [PRJNA311266](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA311266). This assembly was published in 2016 by [_Zapata_ et. al](https://www.pnas.org/content/113/28/E4052). The dataset has 4 libraries, with both PE and MP reads (2 each). Although, it does have long reads (PacBio), we will not be using it for this tutorial. The SRA info can be found [here](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN04457953).  For a 117Mb genome, this is quite a bit of coverage.


**Table 1: Dataset used for the Assembly**

| Run        | Instrument | Layout        | Insert (bp) | ReadLength | TotalReads   | Bases (Mbp) |
|------------|------------|---------------|------------:|------------|-------------:|------------:|
| SRR3157034 | HiSeq 2000 | paired-end    | 0           | 100x2      | 93,446,768   | 17,823      |
| SRR3166543 | HiSeq 2000 | paired-end    | 0           | 100x2      | 162,362,560  | 30,968      |
| SRR3156163 | HiSeq 2000 | mate-pair     | 8,000       | 100x2      | 51,332,776   | 9,790       |
| SRR3156596 | HiSeq 2000 | mate-pair     | 20,000      | 100x2      | 61,030,552   | 11,640      |

## Download

We will use the [sra-toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) to get the fastq files needed. After you copy the SRR ids for the files, we will save it as `srr.ids`

`srr.ids`

```
SRR3156163
SRR3156596
SRR3157034
SRR3166543
```

```bash
module load sra-toolkit
while read line; do
  fastq-dump --split-files --origfmt ${line};
done<srr.ids
```
<details>
  <summary>sra-toolkit stdout</summary>

```
Read 51332776 spots for /work/GIF/arnstrm/ncbi/public/sra/SRR3156163.sra
Written 51332776 spots for /work/GIF/arnstrm/ncbi/public/sra/SRR3156163.sra
Read 61030552 spots for /work/GIF/arnstrm/ncbi/public/sra/SRR3156596.sra
Written 61030552 spots for /work/GIF/arnstrm/ncbi/public/sra/SRR3156596.sra
Read 93446768 spots for /work/GIF/arnstrm/ncbi/public/sra/SRR3157034.sra
Written 93446768 spots for /work/GIF/arnstrm/ncbi/public/sra/SRR3157034.sra
Read 162362560 spots for /work/GIF/arnstrm/ncbi/public/sra/SRR3166543.sra
Written 162362560 spots for /work/GIF/arnstrm/ncbi/public/sra/SRR3166543.sra
```
</details>

Typically, we will start by checking the quality of reads after you obtain the data (using Fastqc), but for this tutorial, since we already know the quality, and to keep it focused on `platanus` assembler, we will skip it.

The compiled versions of the `platanus` assembler is available on the assembly [homepage](http://platanus.bio.titech.ac.jp). You can download them directly and use the binaries for this tutorial (no compiling or installation required).

```bash
wget -O platanus http://platanus.bio.titech.ac.jp/?ddownload=145
wget -O platanus_trim http://platanus.bio.titech.ac.jp/?ddownload=153
wget -O platanus_internal_trim http://platanus.bio.titech.ac.jp/?ddownload=154
# make them executables
chmod +x platanus*
```
## Overview

Steps for running the assembly are as follows:

1. Trim the reads using `platanus` trimming programs. The paired-end and mate-pair reads should be trimmed using different programs as they are normally have different orientation (paired-end=innie; mate-pair=outie).
2. Run contig generation step using only trimmed paired-end reads.
3. The contigs are then scaffolded using trimmed mate-pairs.
4. Gap-closing to fill the gaps using both types of reads and finish the assembly process.

## Organization

The assembly process is simple and it does not create a lot of files. So to keep the organization simple as well, we will run everything in a single directory.

```
usda/
├── platanus
├── platanus_internal_trim
├── platanus_trim
├── SRR3156163_1.fastq
├── SRR3156163_2.fastq
├── SRR3156596_1.fastq
├── SRR3156596_2.fastq
├── SRR3157034_1.fastq
├── SRR3157034_2.fastq
├── SRR3166543_1.fastq
├── SRR3166543_2.fastq
└── srr.ids
```


## Assembly process

As mentioned before, there are 4 stages. We will walk you through each step.

### 1. Trimming

Since the trimming is different for each type of library (paired-end and mate-pair), we will create a file of filenames (fofn) to process them separately.

`pe.fofn`

```
SRR3157034_1.fastq
SRR3157034_2.fastq
SRR3166543_1.fastq
SRR3166543_2.fastq
```

`mp.fofn`
```
SRR3156163_1.fastq
SRR3156163_2.fastq
SRR3156596_1.fastq
SRR3156596_2.fastq
```

Now run the trimming as follows (here `-i` is for input file and `-t` for number of threads):

```
platanus_trim -i pe.fofn -t 16
platanus_internal_trim -i mp.fofn -t 16
```

<details>
  <summary>pe trimming stdout</summary>

```
Running with trim adapter mode
Checking files:
  SRR3157034_1.fastq SRR3157034_2.fastq  (100%)
Checking files:
  SRR3166543_1.fastq SRR3166543_2.fastq  (100%)

Number of trimmed read with adapter:
NUM_OF_TRIMMED_READ(FORWARD) = 136950407
NUM_OF_TRIMMED_BASE(FORWARD) = 369118390
NUM_OF_TRIMMED_READ(REVERSE) = 136950019
NUM_OF_TRIMMED_BASE(REVERSE) = 402829853
NUM_OF_TRIMMED_PAIR(OR) = 136951669
NUM_OF_TRIMMED_PAIR(AND) = 136948757


Number of trimmed read because of low quality or too short (< 11bp):
NUM_OF_TRIMMED_READ(FORWARD) = 122955524
NUM_OF_TRIMMED_BASE(FORWARD) = 6257718860
NUM_OF_TRIMMED_READ(REVERSE) = 109237679
NUM_OF_TRIMMED_BASE(REVERSE) = 7507282392
NUM_OF_TRIMMED_PAIR(OR) = 146417481
NUM_OF_TRIMMED_PAIR(AND) = 85775722


#### PROCESS INFORMATION ####
User Time:         266.78 min
System Time:         1.68 min
VmPeak:           0.833 GByte
VmHWM:            0.144 GByte
Execution time:     59.32 min
```
</details>

<details>
  <summary>mp trimming stdout</summary>

```
Running with trim internal adapter mode
Checking files:
  SRR3156163_1.fastq SRR3156163_2.fastq  (100%)
Checking files:
  SRR3156596_1.fastq SRR3156596_2.fastq  (100%)

Number of trimmed read with internal adapter:
NUM_OF_TRIMMED_READ(FORWARD) = 0
NUM_OF_TRIMMED_BASE(FORWARD) = 0
NUM_OF_TRIMMED_READ(REVERSE) = 0
NUM_OF_TRIMMED_BASE(REVERSE) = 0
NUM_OF_TRIMMED_PAIR(OR) = 0
NUM_OF_TRIMMED_PAIR(AND) = 0


Number of trimmed read with adapter:
NUM_OF_TRIMMED_READ(FORWARD) = 28379276
NUM_OF_TRIMMED_BASE(FORWARD) = 621458117
NUM_OF_TRIMMED_READ(REVERSE) = 28378747
NUM_OF_TRIMMED_BASE(REVERSE) = 205825175
NUM_OF_TRIMMED_PAIR(OR) = 28379567
NUM_OF_TRIMMED_PAIR(AND) = 28378456


Number of trimmed read because of low quality or too short (< 11bp):
NUM_OF_TRIMMED_READ(FORWARD) = 62990345
NUM_OF_TRIMMED_BASE(FORWARD) = 2047046061
NUM_OF_TRIMMED_READ(REVERSE) = 60871687
NUM_OF_TRIMMED_BASE(REVERSE) = 3391343381
NUM_OF_TRIMMED_PAIR(OR) = 79275459
NUM_OF_TRIMMED_PAIR(AND) = 44586573


#### PROCESS INFORMATION ####
User Time:         234.61 min
System Time:         0.79 min
VmPeak:           0.872 GByte
VmHWM:            0.132 GByte
Execution time:     36.22 min
```
</details>


After the run is complete, you will see the files with `.trimmed` and `.int_trimmed` extension. It will also display progress while trimming is running.

### 2. Contig generation:

The main executable for this assembler is `platanus`. This comes with 3 options. `assemble`, `scaffold` and `gap_close`. For this step, we need the `assemble` option.

```bash
platanus assemble \
   -o platanus \
   -f {SRR3157034,SRR3166543}_?.fastq.trimmed \
   -t 12 \
   -m 115 \
   -tmp $TMPDIR
```

We used most of the default options, except for `-t` number of threads (set to 12), `-m` max memory un GB to use for assembly step (set to 115GB), and the `-tmp` temp directory to use (set to scratch space). Note that we used a regular expression (pattern) to provide all the files as input. We could have also provided each filename individually, but using patterns will reduce errors and makes it easier to read.

<details>
  <summary>assembly stdout</summary>

```
Platanus version: 1.2.4
./platanus assemble -o platanus -f SRR3157034_1.fastq.trimmed SRR3157034_2.fastq.trimmed SRR3166543_1.fastq.trimmed SRR3166543_2.fastq.trimmed -t 12 -m 115 -tmp /scratch/arnstrm/504932

K = 32, saving kmers from reads...
AVE_READ_LEN=95.7062

KMER_EXTENSION:
K=32, KMER_COVERAGE=167.299 (>= 44), COVERAGE_CUTOFF=44
K=42, KMER_COVERAGE=141.444, COVERAGE_CUTOFF=44, PROB_SPLIT=10e-inf
K=52, KMER_COVERAGE=115.588, COVERAGE_CUTOFF=44, PROB_SPLIT=10e-13.035
K=62, KMER_COVERAGE=89.7333, COVERAGE_CUTOFF=34, PROB_SPLIT=10e-10.2256
K=72, KMER_COVERAGE=63.8782, COVERAGE_CUTOFF=18, PROB_SPLIT=10e-10.4294
K=82, KMER_COVERAGE=38.0231, COVERAGE_CUTOFF=5, PROB_SPLIT=10e-10.485
K=85, KMER_COVERAGE=30.2666, COVERAGE_CUTOFF=2, PROB_SPLIT=10e-11.0475
K=86, KMER_COVERAGE=27.6811, COVERAGE_CUTOFF=2, PROB_SPLIT=10e-10.2631
loading kmers...
connecting kmers...
removing branches...
BRANCH_DELETE_THRESHOLD=0.5
NUM_CUT=52771
NUM_CUT=325
NUM_CUT=0
TOTAL_NUM_CUT=53096
mapping reads on de Bruijn Graph nodes...
TOTAL_MAPPED_READS=301160221
TOTAL_UNMAPPED_READS=49230885
TOTAL_SHORT_READS(<32)=2449016
NUM_DELETE_NODE(reads are unmapped)=127838
NUM_CUT_NODE=14250
extracting reads (containing kmer used in contig assemble)...
K = 42, loading kmers from contigs...
K = 42, saving additional kmers(not found in contigs) from reads...
COVERAGE_CUTOFF = 44
loading kmers...
connecting kmers...
removing branches...
BRANCH_DELETE_THRESHOLD=0.5
NUM_CUT=32977
NUM_CUT=262
NUM_CUT=0
TOTAL_NUM_CUT=33239
mapping reads on de Bruijn Graph nodes...
TOTAL_MAPPED_READS=297870532
TOTAL_UNMAPPED_READS=48839229
TOTAL_SHORT_READS(<42)=6130361
NUM_DELETE_NODE(reads are unmapped)=34658
NUM_CUT_NODE=9411
extracting reads (containing kmer used in contig assemble)...
K = 52, loading kmers from contigs...
K = 52, saving additional kmers(not found in contigs) from reads...
COVERAGE_CUTOFF = 44
loading kmers...
connecting kmers...
removing branches...
BRANCH_DELETE_THRESHOLD=0.5
NUM_CUT=11612
NUM_CUT=58
NUM_CUT=1
NUM_CUT=0
TOTAL_NUM_CUT=11671
mapping reads on de Bruijn Graph nodes...
TOTAL_MAPPED_READS=294187616
TOTAL_UNMAPPED_READS=48525575
TOTAL_SHORT_READS(<52)=10126931
NUM_DELETE_NODE(reads are unmapped)=14651
NUM_CUT_NODE=6611
extracting reads (containing kmer used in contig assemble)...
K = 62, loading kmers from contigs...
K = 62, saving additional kmers(not found in contigs) from reads...
COVERAGE_CUTOFF = 34
loading kmers...
connecting kmers...
removing branches...
BRANCH_DELETE_THRESHOLD=0.5
NUM_CUT=8760
NUM_CUT=85
NUM_CUT=0
TOTAL_NUM_CUT=8845
mapping reads on de Bruijn Graph nodes...
TOTAL_MAPPED_READS=292269452
TOTAL_UNMAPPED_READS=46285943
TOTAL_SHORT_READS(<62)=14284727
NUM_DELETE_NODE(reads are unmapped)=11077
NUM_CUT_NODE=6253
extracting reads (containing kmer used in contig assemble)...
K = 72, loading kmers from contigs...
K = 72, saving additional kmers(not found in contigs) from reads...
COVERAGE_CUTOFF = 18
loading kmers...
connecting kmers...
removing branches...
BRANCH_DELETE_THRESHOLD=0.5
NUM_CUT=11089
NUM_CUT=73
NUM_CUT=1
NUM_CUT=0
TOTAL_NUM_CUT=11163
mapping reads on de Bruijn Graph nodes...
TOTAL_MAPPED_READS=290273837
TOTAL_UNMAPPED_READS=42592196
TOTAL_SHORT_READS(<72)=19974089
NUM_DELETE_NODE(reads are unmapped)=16011
NUM_CUT_NODE=5055
extracting reads (containing kmer used in contig assemble)...
K = 82, loading kmers from contigs...
K = 82, saving additional kmers(not found in contigs) from reads...
COVERAGE_CUTOFF = 5
loading kmers...
connecting kmers...
removing branches...
BRANCH_DELETE_THRESHOLD=0.5
NUM_CUT=35099
NUM_CUT=316
NUM_CUT=2
NUM_CUT=0
TOTAL_NUM_CUT=35417
mapping reads on de Bruijn Graph nodes...
TOTAL_MAPPED_READS=287485676
TOTAL_UNMAPPED_READS=38455436
TOTAL_SHORT_READS(<82)=26899010
NUM_DELETE_NODE(reads are unmapped)=23565
NUM_CUT_NODE=5757
extracting reads (containing kmer used in contig assemble)...
K = 85, loading kmers from contigs...
K = 85, saving additional kmers(not found in contigs) from reads...
COVERAGE_CUTOFF = 2
loading kmers...
connecting kmers...
removing branches...
BRANCH_DELETE_THRESHOLD=0.5
NUM_CUT=92929
NUM_CUT=705
NUM_CUT=5
NUM_CUT=0
TOTAL_NUM_CUT=93639
mapping reads on de Bruijn Graph nodes...
TOTAL_MAPPED_READS=286286607
TOTAL_UNMAPPED_READS=36890605
TOTAL_SHORT_READS(<85)=29662910
NUM_DELETE_NODE(reads are unmapped)=18033
NUM_CUT_NODE=4914
extracting reads (containing kmer used in contig assemble)...
K = 86, loading kmers from contigs...
K = 86, saving additional kmers(not found in contigs) from reads...
COVERAGE_CUTOFF = 2
loading kmers...
connecting kmers...
removing branches...
BRANCH_DELETE_THRESHOLD=0.5
NUM_CUT=76051
NUM_CUT=455
NUM_CUT=5
NUM_CUT=0
TOTAL_NUM_CUT=76511
LENGTH_CUTOFF = 172
COVERAGE_CUTOFF = 4
removing erroneous nodes...
NUM_REMOVED_NODES=62182
NUM_REMOVED_NODES=1787
NUM_REMOVED_NODES=6
NUM_REMOVED_NODES=0
TOTAL_NUM_REMOVED_NODES=63975
AVE_KMER_COV_REMOVING_BUBBLE=27.7394
removing bubbles...
BUBBLE_IDENTITY_THRESHOLD=0.1
NUM_REMOVED_BUBBLES=2332
NUM_REMOVED_BUBBLES=1
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=2333
mapping reads on de Bruijn Graph nodes...
TOTAL_MAPPED_READS=284641512
TOTAL_UNMAPPED_READS=37462911
TOTAL_SHORT_READS(<86)=30735699
NUM_DELETE_NODE(reads are unmapped)=613
NUM_CUT_NODE=2839
assemble completed!

#### PROCESS INFORMATION ####
VmPeak:         101.557 GByte
VmHWM:            5.544 GByte
```
</details>

There are only 3 output files from this step:
```
platanus_contig.fa : assembled contiguous sequences
platanus_contigBubble.fa : merged and removed bubble sequences
platanus_32merFrq.tsv : kmer frequency distribution
```


### 3. Scaffolding

In this step, we will use the contigs create in the previous step, along with trimmed reads (PE and MP) to scaffold the initial set of contigs. We will run the command as follows:

```bash
platanus scaffold \
    -o platanus \
    -c platanus_contig.fa \
    -b platanus_contigBubble.fa \
    -IP1 {SRR3157034,SRR3166543}_?.fastq.trimmed \
    -OP2 SRR3156163_?.fastq.int_trimmed \
    -n2 7000 \
    -a2 8000 \
    -d2 1000 \
    -OP3 SRR3156596_?.fastq.int_trimmed \
    -n3 19000 \
    -a3 20000 \
    -d3 2000 \
    -t 12 \
    -tmp $TMPDIR
```

Again, most of the default options were used. It requires file from previous steps (contigs and bubbles file, using options `-c` and `-b`, respectively) and reads (both paired-end and mate-pair) as input. For paired-end reads, the insert size was 0, so no other information is provided. However, the 2 set of files we have for mate-pair have different insert sizes and should be included while running scaffolding. We do this by using `-OP2` and `-OP3` options (you can provide as many as you want, depending on number of libraries you use, here we only need 2). For each library, we will also provide `-n` minimum insert size, `-a` average insert size and `-d` standard deviation. The trailing number is used for identifying library with its insert size tags. _Eg.,_ `-n2`, `-a2`, and `-d2` is associated with library `-OP2` and, `-n3`, `-a3`, and `-d3` with library `-OP3`. Just like in previous step, we will change `-t` number of threads (set to 12), and the `-tmp` temp directory to use (set to scratch space).

There will 3 output files again:

```
platanus_scaffold.fa: assembled sequences with gaps
platanus_scaffoldBubble.fa: removed bubble sequences
platanus_scaffoldComponent.tsv: table describing contig joins
```

<details>
  <summary>scaffolding stdout</summary>

```
Platanus version: 1.2.4
./platanus scaffold -o platanus -c platanus_contig.fa -b platanus_contigBubble.fa -IP1 SRR3157034_1.fastq.trimmed SRR3157034_2.fastq.trimmed SRR3166543_1.fastq.trimmed SRR3166543_2.fastq.trimmed -OP2 SRR3156163_1.fastq.int_trimmed SRR3156163_2.fastq.int_trimmed -n2 7000 -a2 8000 -d2 1000 -OP3 SRR3156596_1.fastq.int_trimmed SRR3156596_2.fastq.int_trimmed -n3 19000 -a3 20000 -d3 2000 -t 12 -tmp /scratch/arnstrm/505110

K=32, making hash table...
K=32, making hash table...
CONTIG_AVERAGE_COVERAGE = 214.328
mapping bubbles on contigs...
[LIBRARY 1]
mapping reads...
TOTAL_PAIR = 176420061
MAPPED_PAIR = 120493353 (0.682991)
MAPPED_IN_DIFFERENT_CONTIGS = 9222646 (0.0522766)
MAPPED_IN_SAME_CONTIG = 111270707 (0.630715)
AVERAGE_COVERAGE = 161.193
[LIBRARY 2]
mapping reads...
TOTAL_PAIR = 43437115
MAPPED_PAIR = 13187994 (0.303611)
MAPPED_IN_DIFFERENT_CONTIGS = 10891832 (0.250749)
MAPPED_IN_SAME_CONTIG = 2296162 (0.0528618)
AVERAGE_COVERAGE = 17.3053
Average insert size specified: AVE = 8000
[LIBRARY 3]
mapping reads...
TOTAL_PAIR = 35013593
MAPPED_PAIR = 3295364 (0.0941167)
MAPPED_IN_DIFFERENT_CONTIGS = 3245212 (0.0926843)
MAPPED_IN_SAME_CONTIG = 50152 (0.00143236)
AVERAGE_COVERAGE = 4.30912
Average insert size specified: AVE = 20000
estimating insert-size...
PEAK = 171
LOWER_LIMIT (permissible range to estimate AVE_INS)= 43
UPPER_LIMIT (permissible range to estimate AVE_INS)= 299
AVE_INS = 174
SD_INS = 23
[LIBRARY 1]
AVE_INS = 174, SD_INS = 23
saving overlaps... (LEN_CUTOFF=46)
destructing mapper objects...
[LIBRARY 1]
AVE_INS = 174, SD_INS = 23
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=17 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 46
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=23 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=312
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=312
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=1931
NUM_SPLIT_LINK (not enough mapped pairs)=1
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=1932
scaffolding...
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=1 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=3
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=3
NUM_SPLIT_LINK (not originate from heterozygosity)=717 (COVERAGE_THRESHOLD)
scaffolding...
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=1 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 69
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=1
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=1
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=1892
NUM_SPLIT_LINK (not enough mapped pairs)=1
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=1893
scaffolding...
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=1
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=1
NUM_SPLIT_LINK (not originate from heterozygosity)=9 (COVERAGE_THRESHOLD)
scaffolding...
[LIBRARY 2]
AVE_INS = 8000, SD_INS = 1000
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=24 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 2000
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=298
NUM_SPLIT_LINK (not enough mapped pairs)=1
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=299
deleting edges from repeat contigs...
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=294
NUM_SPLIT_LINK (not enough mapped pairs)=1
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=295
scaffolding...
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=5 (COVERAGE_THRESHOLD)
scaffolding...
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=10 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 3000
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=8
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=8
deleting edges from repeat contigs...
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=8
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=8
scaffolding...
linking scaffolds (MIN_LINK = 7)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=0 (COVERAGE_THRESHOLD)
scaffolding...
[LIBRARY 3]
AVE_INS = 20000, SD_INS = 2000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=42 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 4000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=40
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=40
deleting edges from repeat contigs...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=40
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=40
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=0 (COVERAGE_THRESHOLD)
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=4 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 6000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=5
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=5
deleting edges from repeat contigs...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=5
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=5
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=0 (COVERAGE_THRESHOLD)
scaffolding...
[LIBRARY 1]
AVE_INS = 174, SD_INS = 23
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=3 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 46
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=12
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=12
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=3026
NUM_SPLIT_LINK (not enough mapped pairs)=1
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=3027
deleting edges from repeat contigs...
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=1
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=1
NUM_SPLIT_LINK (not originate from heterozygosity)=100 (COVERAGE_THRESHOLD)
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=0 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 69
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=4180
NUM_SPLIT_LINK (not enough mapped pairs)=2
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=4182
deleting edges from repeat contigs...
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=5 (COVERAGE_THRESHOLD)
scaffolding...
Library1 PHYSICAL_COVERAGE=146
SUM_SHORT_LIBRARY_PHYSICAL_COVERAGE=146
Library2 PHYSICAL_COVERAGE=736
Library3 PHYSICAL_COVERAGE=463
SUM_LONG_LIBRARY_PHYSICAL_COVERAGE=1199
checking erroneous scaffold using long libraries...
spliting low coverage links...
NUM_SPLIT_LINK(low coverage)= 143
[LIBRARY 2]
AVE_INS = 8000, SD_INS = 1000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=3 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 2000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=20
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=20
deleting edges from repeat contigs...
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=0 (COVERAGE_THRESHOLD)
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=1 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 3000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=16
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=16
deleting edges from repeat contigs...
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=0 (COVERAGE_THRESHOLD)
scaffolding...
Library1 PHYSICAL_COVERAGE=146
Library2 PHYSICAL_COVERAGE=736
SUM_SHORT_LIBRARY_PHYSICAL_COVERAGE=882
Library3 PHYSICAL_COVERAGE=463
SUM_LONG_LIBRARY_PHYSICAL_COVERAGE=463
[LIBRARY 3]
AVE_INS = 20000, SD_INS = 2000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=3 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 4000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=20
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=20
deleting edges from repeat contigs...
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=0 (COVERAGE_THRESHOLD)
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_OVERLAP_CONTIGS=1 (CONTAINED_HETERO)
TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = 6000
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
removing erroneous edges...
NUM_SPLIT_LINK (not enough mapped pairs)=1
NUM_SPLIT_LINK (not enough mapped pairs)=0
TOTAL_SPLIT_LINK (not enough mapped pairs)=1
deleting edges from repeat contigs...
scaffolding...
linking scaffolds (MIN_LINK = 3)
sorting links in contigID order...
estimating contig distances...
constructing scaffold graph
NUM_REMOVED_BUBBLES=0 (COVERAGE_THRESHOLD)
removing bubbles... (MAX_BUBBLE_IDENTITY = 0.1)
NUM_REMOVED_BUBBLES=0
TOTAL_NUM_REMOVED_BUBBLES=0
NUM_SPLIT_LINK (not originate from heterozygosity)=0 (COVERAGE_THRESHOLD)
scaffolding...
Library1 PHYSICAL_COVERAGE=146
Library2 PHYSICAL_COVERAGE=736
Library3 PHYSICAL_COVERAGE=463
SUM_SHORT_LIBRARY_PHYSICAL_COVERAGE=1345
SUM_LONG_LIBRARY_PHYSICAL_COVERAGE=0
writing scaffold files...
scaffold completed!

#### PROCESS INFORMATION ####
VmPeak:          12.303 GByte
VmHWM:            2.066 GByte
```
</details>

### 4. Gap closing
The final step is to close the gaps (reduce the N content in the genome introduced during scaffolding). We will need the scaffolds generated in the previous step. For this we will use the `gap_close` module. We will run it as follows:

```bash
platanus gap_close \
    -o platanus \
    -c platanus_scaffold.fa \
    -IP1 {SRR3157034,SRR3166543}_?.fastq.trimmed \
    -OP2 SRR3156163_?.fastq.int_trimmed \
    -OP3 SRR3156596_?.fastq.int_trimmed \
    -t 16 \
    -tmp $TMPDIR
```
Most of the default options were used as usual, except for `-t` number of threads (set to 12), and the `-tmp` temp directory to use (set to scratch space). The 2 sets of PE libraries were provided as input using `IP1` option and the 2 mate-pair libraries separately using `-OP2` and `-OP3` options.  

<details>
  <summary>gap-closing stdout</summary>

```
Platanus version: 1.2.4
./platanus gap_close -o platanus -c platanus_scaffold.fa -IP1 SRR3157034_1.fastq.trimmed SRR3157034_2.fastq.trimmed SRR3166543_1.fastq.trimmed SRR3166543_2.fastq.trimmed -OP2 SRR3156163_1.fastq.int_trimmed SRR3156163_2.fastq.int_trimmed -OP3 SRR3156596_1.fastq.int_trimmed SRR3156596_2.fastq.int_trimmed -t 16 -tmp /scratch/arnstrm/505111

K=32, making hash table...
[PAIR_LIBRARY 1]
mapping reads...
TOTAL_PAIR = 176420061
MAPPED_IN_SAME_CONTIG = 114172231 (0.647161)
estimating insert-size...
PEAK = 170
LOWER_LIMIT (permissible range to estimate AVE_INS)= 43
UPPER_LIMIT (permissible range to estimate AVE_INS)= 298
AVE_INS = 173
SD_INS = 24
mapping reads that cover small gaps...
[PAIR_LIBRARY 2]
mapping reads...
TOTAL_PAIR = 43437115
MAPPED_IN_SAME_CONTIG = 9243937 (0.212812)
estimating insert-size...
PEAK = 8227
LOWER_LIMIT (permissible range to estimate AVE_INS)= 2057
UPPER_LIMIT (permissible range to estimate AVE_INS)= 14397
AVE_INS = 8967
SD_INS = 1236
mapping reads that cover small gaps...
[PAIR_LIBRARY 3]
mapping reads...
TOTAL_PAIR = 35013593
MAPPED_IN_SAME_CONTIG = 1663835 (0.0475197)
estimating insert-size...
PEAK = 21776
LOWER_LIMIT (permissible range to estimate AVE_INS)= 5444
UPPER_LIMIT (permissible range to estimate AVE_INS)= 38108
AVE_INS = 23799
SD_INS = 4754
mapping reads that cover small gaps...
making hash table of gaps...
making consensus sequences to close small gaps...
NUM_GAP=18477
NUM_CLOSED_GAP=5061
[PAIR_LIBRARY 1]
saving reads covering gaps...
loading reads covering gaps...
assembling localized reads...
NUM_GAPS = 13416
NUM_NOT_CLOSED_GAPS (too many reads are mapped comapering to coverage)= 91
NUM_CLOSED_GAPS_USING_DE_BRUIJN = 4229
NUM_CLOSED_GAPS_USING_OVERLAP_LAYOUT_CONSENSUS = 785
[PAIR_LIBRARY 2]
saving reads covering gaps...
loading reads covering gaps...
assembling localized reads...
NUM_GAPS = 8402
NUM_NOT_CLOSED_GAPS (too many reads are mapped comapering to coverage)= 20
NUM_CLOSED_GAPS_USING_DE_BRUIJN = 26
NUM_CLOSED_GAPS_USING_OVERLAP_LAYOUT_CONSENSUS = 94
[PAIR_LIBRARY 3]
saving reads covering gaps...
loading reads covering gaps...
assembling localized reads...
NUM_GAPS = 8282
NUM_NOT_CLOSED_GAPS (too many reads are mapped comapering to coverage)= 15
NUM_CLOSED_GAPS_USING_DE_BRUIJN = 3
NUM_CLOSED_GAPS_USING_OVERLAP_LAYOUT_CONSENSUS = 0
[ALL LIBRARY]
assembling localized reads...
NUM_GAPS = 8279
NUM_NOT_CLOSED_GAPS (too many reads are mapped comapering to coverage)= 122
NUM_CLOSED_GAPS_USING_DE_BRUIJN = 17
NUM_CLOSED_GAPS_USING_OVERLAP_LAYOUT_CONSENSUS = 0
TOTAL_NUM_CLOSED_GAPS = 10215
gap_close completed!!

#### PROCESS INFORMATION ####
VmPeak:          12.710 GByte
VmHWM:            2.065 GByte
```
</details>

## Benchmark

The entire assembly was run on Intel(R) Xeon(R) CPU E5-2650 0 @ 2.00GHz machine with 12 processors and with 128gb RAM (HPC Condo cluster, free nodes). The table shows the time taken for each step (`[h]:mm:ss` format)

**Table 2: Time (real, user and sys) used for each step of the assembly**

| step                   | real     | user     | sys      |
|:-----------------------|----------|----------|----------|
| fastq_dump             | 0:00:00  | 0:00:00  | 0:00:00  |
| platanus_trim          | 0:59:19  | 4:26:47  | 0:01:41  |
| platanus_internal_trim | 0:36:13  | 3:54:37  | 0:00:48  |
| assembly               | 3:53:35  | 22:14:12 | 0:56:56  |
| scaffold               | 0:25:35  | 1:59:41  | 0:04:53  |
| gapclosing             | 5:58:32	| 55:58:05 | 0:14:05  |

The stats for the assemblies obtained after each step (`assemble`, `scaffold` and `gap_close`) are as follows:

**Table 3: Assemblathon stats for contigs, scaffolds and final assembly.**

| Metrics                       | Contig      | Scaffold    | Final       |
|:------------------------------|------------:|------------:|------------:|
| Number of scaffolds           | 249,090     | 18,915      | 18,915      |
| Total size of scaffolds       | 146,125,822 | 128,269,400 | 127,986,775 |
| Percentage genome represented | 124.9%      | 109.6%      | 109.4%      |
| Longest scaffold              | 153,225     | 15,741,082  | 15,721,187  |
| Shortest scaffold             | 86          | 100         | 86          |
| Number of scaffolds > 1K nt   | 19,518      | 2,799       | 2,784       |
| Number of scaffolds > 10K nt  | 3,108       | 230         | 230         |
| Number of scaffolds > 100K nt | 1           | 77          | 77          |
| Number of scaffolds > 1M nt   | 0           | 25          | 25          |
| Number of scaffolds > 10M nt  | 0           | 1           | 1           |
| Mean scaffold size            | 587         | 6,781       | 6,766       |
| Median scaffold size          | 101         | 204         | 203         |
| N50 scaffold length           | 6,442       | 3,885,840   | 3,861,540   |
| L50 scaffold count            | 5,705       | 9           | 9           |
| NG50 scaffold length          | 8,795       | 4,991,028   | 4,978,789   |
| LG50 scaffold count           | 3,773       | 8           | 8           |
| N50- NG50 difference          | 2,353       | 1,105,188   | 1,117,249   |
| scaffold %A                   | 31.6%       | 29.9%       | 30.0%       |
| scaffold %C                   | 18.7%       | 17.5%       | 17.6%       |
| scaffold %G                   | 18.3%       | 17.5%       | 17.5%       |
| scaffold %T                   | 31.4%       | 29.9%       | 30.0%       |
| scaffold %N                   | 0.0%        | 5.3%        | 5.0%        |
| scaffold %non-ACGTN           | 0.0%        | 0.0%        | 0.0%        |


[Arabidopsis data set Info](Arabidopsis_background.md)
[Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
