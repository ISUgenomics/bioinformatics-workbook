---
title: "Introduction to Data Acquisition"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


# What dataset will we use throughout this text?

As much as possible, we will be using Arabidopsis data from multiple NCBI BioProject that contains datasets for many of the most common data analyses.  The following BioProjects were chosen.  

| SeqType              | Platform | ReadType   | BioProject                                                                                | Experiment                                                                                                                                                                 |
|----------------------|----------|------------|-------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ChIP-seq             | Illumina | single     | [PRJNA316877](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA316877)                  | [Requirement for flap endonuclease 1 (FEN1) to maintain genomic stability and transcriptional gene silencing in Arabidopsis](https://www.ncbi.nlm.nih.gov/pubmed/27231839) |
| ChIP-seq             | Illumina | paired     | [PRJNA349052](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA349052)                  | [Centromere location in Arabidopsis is unaltered by extreme divergence in CENH3 protein sequence]()                                                                        |
| ATAC-seq             | Illumina | paired     | [PRJNA394532](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA394532)                  | [ATAC-seq profiling of open chromatin in the root tips]()                                                                                                                  |
| RNA-seq              | Illumina | single     | [PRJNA312637](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA312637)                  | [RNA-seq analysis of transcriptomes in cae2-1, CA1-1 and cae2-1 CA1-1 Arabidopsis genotypes](https://www.ncbi.nlm.nih.gov/pubmed/27911772)                                 |
| RNA-seq              | Illumina | paired     | [PRJNA348194](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA348194)                  | [Analysis of gene expression in a ATRX loss-of-function line](https://www.ncbi.nlm.nih.gov/pubmed/28684426)                                                                |
| ncRNA                | SOLiD    | single     | [PRJNA169627](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA169627)                  | [Deep sequencing of small RNAs]()                                                                                                                                          |
| microRNA             | Illumina | single     | [PRJNA355875](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA355875)                  | [Differential expression of microRNAs in wildtype versus DCL1 mutants in Arabidopsis thaliana](https://www.ncbi.nlm.nih.gov/pubmed/28407097)                               |
| Long Reads           | PacBio   | long-reads | [PRJNA314706](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA314706)                  | [Diploid Arabidopsis thaliana genome sequencing and assembly](https://www.ncbi.nlm.nih.gov/pubmed/27749838)                                                                |
| DNAseq               | Illumina | paired-end | [PRJEB13889](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB13889)                    | [Genome stability under UV-B in Arabidopsis   thaliana](https://www.nature.com/articles/ncomms13522/)                                                                                                                    |
| DNAseq               | Illumina | mate-pair  | [SRX1434948](https://www.ncbi.nlm.nih.gov/sra/SRX1434948/)                                | [Arabidopsis thaliana Genome sequencing and assembly](https://www.ncbi.nlm.nih.gov/pubmed/27711162)                                                                        |
| 16s-rRNA             | Illumina | paired-end | [MG-RAST:4457768.3-4459735.3](https://docs.qiime2.org/2017.10/tutorials/moving-pictures/) | [Moving pictures of the human microbiome (“Moving Pictures” tutorial)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r50)                          |
| Shotgun metagenomics | Illumina | paired-end | [ERX2017035](https://www.ncbi.nlm.nih.gov/sra/ERX2017035/)                                | [A case of hepatic brucelloma studied by next generation sequencing]()                                                                                                     |

#### Why Arabidopsis?
Arabidopsis is one of several model organisms where significant amounts of data has been collected on a wide variety of bioinformatic data analysis problems.  Additional example datasets from a variety of organism will also be provided as problem sets to explore.



---


[Previous](fileTransfer/sra.md){: .btn  .btn--primary}
[Table of contents](../index-bk.md){: .btn  .btn--primary}
