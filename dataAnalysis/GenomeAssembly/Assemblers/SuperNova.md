---
title: Genome Assembly
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Introduction to SuperNova genome assembler


## Background

* Main software website links
  * [SuperNova home page](https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/welcome)
  * [SuperNOva Paper](http://genome.cshlp.org/content/27/5/757.short)

* SuperNova is intended to run only Illumina linked-read libraries that were produced using the 10x Chromium instrument.
* SuperNova should be run using 38-56x coverage of the genome. If the coverage is far from the recommended range SuperNova will exit.
* At most 2.14 billion reads are allowed.


## Installing SuperNova
-------------------


## How to run SuperNova
-------------------

* **SuperNova parameter considerations**
```
--id      A unique run ID string
--fastqs  Path of the FASTQ folder
--sample  (optional) Can be used to select only a single sample of those specified in the sample sheet supplied to `fastq`. Default is all the samples
--description (optional) Description of the data set.
--maxreads  Target number of reads to be used. You should estimate this number to achieve 56x coverage.
--localcores  (optional) Limits concurrent sections of SuperNova to use the specified number of cores.
--localmem (optional) Limits memory use
```
