---
title: Genome Assembly
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Introduction to SPAdes genome assembler


SPAdes is not intended for larger genomes (e.g. mammalian size genomes)

### Background reading about the Assembling program SPAdes

* [SPAdes 3.13.0 Manual](http://cab.spbu.ru/files/release3.13.0/manual.html)
* [Center for Algorithmic Biotechnology](http://cab.spbu.ru/software/spades/#benchmark)
* [SPAdes Genome paper 2012](http://cab.spbu.ru/software/spades/#benchmark)


SPAdes is a very memory intensive program.  In multithreaded mode (-t 16), you will want at least 500 Gigabytes if not 1000 Gigabytes of RAM.  


## Errors


#### Not enough memory
These two errors are extremely unhelpful for diagnosing the problem.  If you get them, try a computer with much higher memory or running spades single threaded (-t 1).  Also, I have seen multiple err code message numbers while attempting multithreaded spades on a low memory machine (<256Gb RAM)

```
libgomp: Thread creation failed: Resource temporarily unavailable
```
```
finished abnormally, err code: 1
```

#### Remove spaces in fastq files

If you still are running into errors, you may also consider removing spaces in fastq sequence headers as follows.

```
gunzip SRR2093871_1.fastq.gz &
gunzip SRR2093871_2.fastq.gz &
gunzip SRR2093872_1.fastq.gz &
gunzip SRR2093872_2.fastq.gz &
perl -i -pe 's/ /_/g' SRR2093871_1.fastq &
perl -i -pe 's/ /_/g' SRR2093871_2.fastq &
perl -i -pe 's/ /_/g' SRR2093872_2.fastq &
perl -i -pe 's/ /_/g' SRR2093872_1.fastq &
gzip SRR2093871_1.fastq &
gzip SRR2093871_2.fastq &
gzip SRR2093872_1.fastq &
gzip SRR2093872_2.fastq &
```
