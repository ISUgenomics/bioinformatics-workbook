---
title: Spades Genome assembly of a Bacillus thuringiensis
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


## Organizing the project folder

Make a main directory where all the work will be performed.

```
mkdir BTspades
cd BTspades
```


## Downloading the dataset

Make a subdirectory for the raw data files.

```
mkdir 00_rawdata
```

Download the raw data from EBI.
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/001/SRR2093871/SRR2093871_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/001/SRR2093871/SRR2093871_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/002/SRR2093872/SRR2093872_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/002/SRR2093872/SRR2093872_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/006/SRR2093876/SRR2093876_subreads.fastq.gz
```

## QC on the data

### Illumina data
If you only have a few files to check, it is easier to start an interactive node and run them quickly in the terminal.

Start interactive node.
```
salloc -N 1 -t 4:00:00
```

Load the fastqc module and execute fastqc for each file.
```
module load fastqc

fastqc --version
FastQC v0.11.3

fastqc SRR2093871_1.fastq.gz  
fastqc SRR2093871_2.fastq.gz  
fastqc SRR2093872_1.fastq.gz  
fastqc SRR2093872_2.fastq.gz
```

Results

* [Fastqc results for SRR2093871_1](https://isugenomics.github.io/bioinformatics-workbook/dataAnalysis/GenomeAssembly/BT/assets/SRR2093871_1_fastqc.html)
* [Fastqc results for SRR2093871_2](https://isugenomics.github.io/bioinformatics-workbook/dataAnalysis/GenomeAssembly/BT/assets/SRR2093871_2_fastqc.html)
* [Fastqc results for SRR2093872_1](https://isugenomics.github.io/bioinformatics-workbook/dataAnalysis/GenomeAssembly/BT/assets/SRR2093872_1_fastqc.html)
* [Fastqc results for SRR2093872_2](https://isugenomics.github.io/bioinformatics-workbook/dataAnalysis/GenomeAssembly/BT/assets/SRR2093872_2_fastqc.html)

FASTQC indicates there are no major issues with these data.  

For datasets with many files see our tutorial on QC [TODO Add tutorial on QC]


### PacBio or long read data

```
SRR2093876_subreads.fastq.gz
```

## Genome assembly using Spades

Change back into the main folder BTspades and make a subdirectory for the Spades assembly files

```
mkdir 01_spades
cd 01_spades
```

spades.sub header used on XSEDE bridges

```
#!/bin/bash
#SBATCH -J spades4.e%j
#SBATCH -o spades4.o%j
#SBATCH --mem 1000GB
#SBATCH -p LM
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 48:00:00

source /usr/share/Modules/init/bash
module load spades/3.11.1
```

Running Spades

```
module load spades
spades.py \
-k 21,33,55,77,99,127 \
-t 16 -m 1000 \
--pacbio /pylon5/mc48o5p/severin/spadesTutorial/00_rawdata/SRR2093876_subreads.fastq.gz \
--pe1-1 /pylon5/mc48o5p/severin/spadesTutorial/00_rawdata/SRR2093871_1.fastq.gz  \
--pe1-2 /pylon5/mc48o5p/severin/spadesTutorial/00_rawdata/SRR2093871_2.fastq.gz  \
--mp1-1 /pylon5/mc48o5p/severin/spadesTutorial/00_rawdata/SRR2093872_1.fastq.gz  \
--mp1-2 /pylon5/mc48o5p/severin/spadesTutorial/00_rawdata/SRR2093872_2.fastq.gz  \
--careful \
-o BT
```


|SPAdes| files output| from assembly|
|--|--|--|
|K21|K33|K55|
|K77|K99|assembly_graph.fastg|
|assembly_graph_with_scaffolds.gfa|before_rr.fasta|contigs.fasta|
|contigs.paths|corrected|dataset.info|
|input_dataset.yaml|misc|mismatch_corrector|
|params.txt|**scaffolds.fasta**|scaffolds.paths|
|spades.log|tmp|warnings.log|

## Assembly statistics




```
new_Assemblathon scaffolds.fasta

---------------- Information for assembly './scaffolds.fasta' ----------------


                                         Number of scaffolds       1334
                                     Total size of scaffolds    6768928
                                            Longest scaffold     395916
                                           Shortest scaffold        100
                                 Number of scaffolds > 1K nt        104   7.8%
                                Number of scaffolds > 10K nt         58   4.3%
                               Number of scaffolds > 100K nt         23   1.7%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size       5074
                                        Median scaffold size        306
                                         N50 scaffold length     218044
                                          L50 scaffold count         12
                                                 scaffold %A      31.95
                                                 scaffold %C      17.95
                                                 scaffold %G      18.09
                                                 scaffold %T      32.01
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0

                Percentage of assembly in scaffolded contigs       0.0%
              Percentage of assembly in unscaffolded contigs     100.0%
                      Average number of contigs per scaffold        1.0
Average length of break (>25 Ns) between contigs in scaffold          0
```

---

* [Bacillus thuringiensis data set Info](BT_background.md)
* [Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
