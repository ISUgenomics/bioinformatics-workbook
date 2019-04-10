# Short read assembly using Platanus

This tutorial covers genome assembly for a sample dataset using `platanus` genome assembler. This assembly program uses short-reads (paired-end and mate-pair) as input data and works well for highly heterozygous genomes.

## Dataset

For this tutorial, we will use _Arabidopsis thaliana_, Ler strain, short read dataset, from the BioProject [PRJNA311266](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA311266). This dataset as 4 libraries, with both PE and MP reads. Although, it does have long reads (PacBio), we will not be using it for this tutorial. The SRA info can be found [here](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN04457953).

Table 1: Dataset used for the Assembly

| Run        | Instrument | Layout        | Insert (bp) | ReadLength | TotalReads   | Bases (Mbp) |
|------------|------------|---------------|------------:|------------|-------------:|------------:|
| SRR3157034 | HiSeq 2000 | paired-end    | 0           | 100x2      | 93,446,768   | 17,823      |
| SRR3166543 | HiSeq 2000 | paired-end    | 0           | 100x2      | 162,362,560  | 30,968      |
| SRR3156163 | HiSeq 2000 | mate-pair     | 8,000       | 100x2      | 51,332,776   | 9,790       |
| SRR3156596 | HiSeq 2000 | mate-pair     | 20,000      | 100x2      | 61,030,552   | 11,640      |

## Download

We will use the [sra-toolkit]() to get the fastq files needed. After you copy the SRR ids for the files, we will save it as `srr.ids`

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

For a typical project, once you obtain the data, you will run the quality check (using Fastqc) and then proceed with the assembly steps. For this tutorial, to keep it focused on `platanus` assembler, we will skip it.

The compiled versions of the `platanus` assembler are available at the assembly [homepage](http://platanus.bio.titech.ac.jp). You can download them directly and use the binaries for this tutorial.

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
├── mp.fofn
├── pe.fofn
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

We used most of the default options, except for `-t` number of threads (set to 12), `-m` max memory un GB to use for assembly step (set to 115GB), and the `-tmp` temp directory to use (set to scratch space).

<details>
  <summary>assembly stdout</summary>

```
```
</details>

### 3. Scaffolding


<details>
  <summary>scaffolding stdout</summary>

```
```
</details>
### 4. Gap closing


<details>
  <summary>gap-closing stdout</summary>

```
```
</details>

### Benchmark

The entire assembly was run on Intel(R) Xeon(R) CPU E5-2650 0 @ 2.00GHz machine with 12 processors and with 128gb RAM (HPC Condo cluster, free nodes). The table shows the time taken for each step (`[h]:mm:ss` format)


| step                   | real    | user    | sys     |
|------------------------|---------|---------|---------|
| fastq_dump             | 0:00:00 | 0:00:00 | 0:00:00 |
| platanus_trim          | 0:59:19 | 4:26:47 | 0:01:41 |
| platanus_internal_trim | 0:36:13 | 3:54:37 | 0:00:48 |
| assembly               | 0:00:00 | 0:00:00 | 0:00:00 |
| scaffold               | 0:00:00 | 0:00:00 | 0:00:00 |
| gapclosing             | 0:00:00 | 0:00:00 | 0:00:00 |
