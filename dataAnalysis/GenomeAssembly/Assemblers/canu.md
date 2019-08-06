---
title:  Introduction to Canu genome assembler
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---



## Background

* Main software website links
  * [Canu Manual](http://canu.readthedocs.io/en/latest/history.html)
  * [Canu Quick Start Guide](http://canu.readthedocs.io/en/latest/quick-start.html)
  * [Canu Genome Paper 2017](https://genome.cshlp.org/content/27/5/722.full)

  Canu is based off of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page) and is designed for noisy long-read data from [PacBio](http://www.pacb.com/) and [NanoPore](https://www.nanoporetech.com/).  More on the history of Canu can be found .

## Canu Steps

* The Canu Assembler has three steps
  * Correct   (-correct)
  * Trim      (-trim)
  * Assemble  (-assemble)

  These steps can be run individual if you specifiy the (-step) parameter as defined in parenthesis above or these steps will be all run if no step parameter is specified.

## How to run Canu

* #### Canu parameter considerations

  ```
  Trio binning
    * If you have sequence data from the parents and the F1 offspring
  Raw Data type requires a different parameter to read each datatype
    * PacBio
    * Nanopore
  Coverage
    * Low Coverage parameters can be set to improve assembly output depending on the sequencing technology See [Canu Quick Start Guide](http://canu.readthedocs.io/en/latest/quick-start.html) for more details.
  Canu Basics
    * -p is the assembly prefix and this is the name that will be prefixed to all output Files
    * -d is the directory that it will make and write all the files to.
  input file types (multiple files can be listed after this parameter but should be of the same type)
      * -pacbio-raw
      * -pacbio-corrected
      * -nanopore-raw
      * -nanopore-corrected
  ```

* #### Running Canu

  The name of the module will vary here and you should check to see what version you are using.  

* #### Example SLURM Job

  ```bash
  #!/bin/bash
  #SBATCH --nodes=1
  #SBATCH --time=24:00:00
  #SBATCH --mem=64G
  #SBATCH --mail-user=YOUREMAILADDRESS
  #SBATCH --mail-type=begin
  #SBATCH --mail-type=end
  #SBATCH --error=JobName.%J.err
  #SBATCH --output=JobName.%J.out
  module load canu

  canu -p Bt2 -d Bt2_assembly genomeSize=6.2m  -pacbio-raw SRR2093876_subreads.fastq.gz
  ```

## Expected files generated during assembly

|Files|output|from assembly|
|--|--|--|
|Bt2.contigs.fasta|Bt2.contigs.gfa|Bt2.contigs.layout|
|Bt2.contigs.layout.readToTig|Bt2.contigs.layout.tigInfo|Bt2.correctedReads.fasta.gz|
|Bt2.gkpStore|Bt2.gkpStore.err|Bt2.gkpStore.gkp|
|Bt2.report|Bt2.trimmedReads.fasta.gz|Bt2.unassembled.fasta|
|Bt2.unitigs.bed|Bt2.unitigs.fasta|Bt2.unitigs.gfa|
|Bt2.unitigs.layout|Bt2.unitigs.layout.readToTig|Bt2.unitigs.layout.tigInfo|
|canu.out|canu-logs|canu-scripts|
|correction|trimming|unitigging|

## SLURM standard output
* An example of the output log can be found here: [canu.out](logs/canu.out)

## Errors

## Further Reading

* Main software website links
  * [Canu Manual](http://canu.readthedocs.io/en/latest/history.html)
  * [Canu Quick Start Guide](http://canu.readthedocs.io/en/latest/quick-start.html)
  * [Canu Genome Paper 2017](https://genome.cshlp.org/content/27/5/722.full)


[Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
