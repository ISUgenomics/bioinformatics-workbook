---
title: Arabidopsis thaliana dataset
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Arabidopsis thaliana dataset

## Background
*Arabidopsis thaliana* which has a genome size of approximately 135 Mb.  

Here is a link to the article written about these data and experiment: [Chromosome-level assembly of Arabidopsis thaliana Ler reveals the extent of translocation and inversion polymorphisms](https://www.pnas.org/content/113/28/E4052.long).  PacBio reads were taken from [Pacific Biosciences Model Organism Genome Sequencing-Arabidopsis thaliana P4C2](https://www.ncbi.nlm.nih.gov/sra/?term=SRX533608).

## Files


| Run        | Instrument | Layout        | Insert (bp) | ReadLength | TotalReads   | Bases (Mbp) |
|------------|------------|---------------|------------:|------------|-------------:|------------:|
| SRR3157034 | HiSeq 2000 | paired-end    | 0           | 100x2      | 93,446,768   | 17,823      |
| SRR3166543 | HiSeq 2000 | paired-end    | 0           | 100x2      | 162,362,560  | 30,968      |
| SRR3156163 | HiSeq 2000 | mate-pair     | 8,000       | 100x2      | 51,332,776   | 9,790       |
| SRR3156596 | HiSeq 2000 | mate-pair     | 20,000      | 100x2      | 61,030,552   | 11,640      |
| SRR1284771 | PacBio RSII | Single       | NA          | unknown    | 163,482      | 2,3       |


## How to download the data from SRA


Downloading from SRA will be performed using the [sra-toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc).
* create a file with the SRA ids and name it `srr.ids`

  ```
  SRR3156163
  SRR3156596
  SRR3157034
  SRR3166543
  ```

* Assuming the sra-toolkit is installed then load the module and run the following bash script on the command line.

  ```bash
  module load sra-toolkit
  while read line; do
    fastq-dump --split-files --origfmt ${line};
  done<srr.ids
  ```

* The pacbio data has to be downloaded separately

  ```
  fastq-dump --table SEQUENCE --origfmt SRR3156160
  ```

  **Note:** sra-toolkit will create a folder named ```ncbi``` in your home directory ```/home/userid/ncbi```  If you have a disk storage limit on your home directory (most supercomputers do), you will want to move that folder to a different location and then create a softlink in your home folder.

  **Error example:** 2019-04-16T19:43:49 fastq-dump.2.8.1 err: unknown while writing file within file system module - unknown system error errno=Disk quota exceeded(122)


## Assembly statistics for Arabidopsis

| Statistics | value|
| :-- | :-- |
|Total sequence length	|118,890,721|
|Total ungapped length	|117,113,196|
|Gaps between scaffolds|	0|
|Number of scaffolds|	30|
|Scaffold N50	|22,588,203|
|Scaffold L50|	3|
|Number of contigs|	525|
|Contig N50	|1,193,183|
|Contig L50	|27|
|Total number of chromosomes and plasmids|	7|
|Number of component sequences (WGS or clone)	|30|


[Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
