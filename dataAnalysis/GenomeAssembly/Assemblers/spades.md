---
title: Introduction to SPAdes genome assembler
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


## Background

* Main software website links
  * [SPAdes 3.13.0 Manual](http://cab.spbu.ru/files/release3.13.0/manual.html)
  * [Center for Algorithmic Biotechnology](http://cab.spbu.ru/software/spades/#benchmark)
  * [SPAdes Genome paper 2012](http://cab.spbu.ru/software/spades/#benchmark)

  SPAdes is not intended for larger genomes (e.g. mammalian size genomes). SPAdes is a very memory intensive program.  In multithreaded mode (-t 16), you will want at least 500 Gigabytes if not 1000 Gigabytes of RAM.  

## How to run SPAdes

* #### SPAdes parameter considerations

  ```bash
  --careful		tries to reduce number of mismatches and short indels
  --tmp-dir	<dirname>	directory for temporary files
  -t/--threads	<int>		number of threads
  -m/--memory	<int>		RAM limit for SPAdes in Gb (terminates if exceeded) [default 250 (Gb)]
  --sanger	<filename>	file with Sanger reads
  --pacbio	<filename>	file with PacBio reads
  --nanopore	<filename>	file with Nanopore reads
  ```

* #### Example SLURM Job

  ```bash
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


## Expected files generated during assembly

|SPAdes| files output| from assembly|
|--|--|--|
|K21|K33|K55|
|K77|K99|assembly_graph.fastg|
|assembly_graph_with_scaffolds.gfa|before_rr.fasta|contigs.fasta|
|contigs.paths|corrected|dataset.info|
|input_dataset.yaml|misc|mismatch_corrector|
|params.txt|**scaffolds.fasta**|scaffolds.paths|
|spades.log|tmp|warnings.log|

## SLURM standard output

* An example of the output log can be found here: [spades.log](dataAnalysis/GenomeAssembly/Assemblers/logs/spades.log)

## Errors

* #### Not enough memory
  These two errors are extremely unhelpful for diagnosing the problem.  If you get them, try a computer with much higher memory or running spades single threaded (-t 1).  Also, I have seen multiple err code message numbers while attempting multithreaded spades on a low memory machine (<256Gb RAM)

  ```bash
  libgomp: Thread creation failed: Resource temporarily unavailable
  ```
  ```bash
  finished abnormally, err code: 1
  ```
* #### Remove spaces in fastq files

  If you still are running into errors, you may also consider removing spaces in fastq sequence headers as follows.

  ```bash
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

## Further Reading

* Main software website links
  * [SPAdes 3.13.0 Manual](http://cab.spbu.ru/files/release3.13.0/manual.html)
  * [Center for Algorithmic Biotechnology](http://cab.spbu.ru/software/spades/#benchmark)
  * [SPAdes Genome paper 2012](http://cab.spbu.ru/software/spades/#benchmark)


[Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
