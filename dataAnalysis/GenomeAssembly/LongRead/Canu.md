---
title: Genome Assembly
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Genome assembly with Canu

Canu is based off of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page) and is designed for noisy long-read data from [PacBio](http://www.pacb.com/) and [NanoPore](https://www.nanoporetech.com/).  More on the history of Canu can be found [here](http://canu.readthedocs.io/en/latest/history.html).


## Learning Objectives

Upon completion of this section on experimental design the learner will be able to

* Download data from SRA
* Evaluate the correct parameters to use when running canu for genome assembly
* Execute Canu and obtain a the assembled genome of


## Canu Steps

The Canu Assembler has three steps
* Correct   (-correct)
* Trim      (-trim)
* Assemble  (-assemble)

These steps can be run individual if you specifiy the (-step) parameter as defined in parenthesis above or these steps will be all run if no step parameter is specified.

## Canu parameter considerations

* Trio binning
  * If you have sequence data from the parents and the F1 offspring
* Raw Data type requires a different parameter to read each datatype
  * PacBio
  * Nanopore
* Coverage
  * Low Coverage parameters can be set to improve assembly output depending on the sequencing technology See [Canu Quick Start Guide](http://canu.readthedocs.io/en/latest/quick-start.html) for more details.


## An Example with canu

Let us assemble a bacterial genome [Bacillus thuringiensis strain: HS18-1](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2093876). This particular genome is interesting because it shows insecticidal activity against Diptera.  The raw data contains 1.5G of data and the genome size is 6.4 mb.  The assembly can be found on [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_001182785.1).  Therefore, we will have over 200x coverage, plenty for a genome assembly with canu.  Canu recommends 30-60x minimum coverage.

#### Download the data from ENA

The European Nucleotide Archive has the files already in fastq format and it is easy to download.
The main page with this project data is [PRJNA288953](https://www.ebi.ac.uk/ena/data/view/PRJNA288953)
We are only going to take the pacbio data which is the 3rd file link under the '''Read Files''' tab.


```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/006/SRR2093876/SRR2093876_subreads.fastq.gz
```

#### Running Canu

The name of the module will vary here and you should check to see what version you are using.  In this tutorial we used version 1.7 with java/1.8.0_121

```
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

#### parameters explained

You can read all of this from the quick start and documentation for Canu but here are the basics.

* -p is the assembly prefix and this is the name that will be prefixed to all output Files
* -d is the directory that it will make and write all the files to.
* input file types (multiple files can be listed after this parameter but should be of the same type)
  * -pacbio-raw
  * -pacbio-corrected
  * -nanopore-raw
  * -nanopore-corrected


#### OUTPUT
Canu took about half an hour to run this dataset.  The output looked like this.

```
Bt2_assembly]$ ls
Bt2.contigs.fasta   Bt2.contigs.layout.readToTig  Bt2.gkpStore      Bt2.report                 Bt2.unitigs.bed    Bt2.unitigs.layout            canu-logs     correction
Bt2.contigs.gfa     Bt2.contigs.layout.tigInfo    Bt2.gkpStore.err  Bt2.trimmedReads.fasta.gz  Bt2.unitigs.fasta  Bt2.unitigs.layout.readToTig  canu.out      trimming
Bt2.contigs.layout  Bt2.correctedReads.fasta.gz   Bt2.gkpStore.gkp  Bt2.unassembled.fasta      Bt2.unitigs.gfa    Bt2.unitigs.layout.tigInfo    canu-scripts  unitigging
```

You can get a sense of progress on your larger datasets by looking at the canu-scripts folder.  The final script is canu.12.sh.  

```
Bt2_assembly/canu-scripts/
canu.01.out  canu.02.out  canu.03.out  canu.04.out  canu.05.out  canu.06.out  canu.07.out  canu.08.out  canu.09.out  canu.10.out  canu.11.out  canu.12.out
canu.01.sh   canu.02.sh   canu.03.sh   canu.04.sh   canu.05.sh   canu.06.sh   canu.07.sh   canu.08.sh   canu.09.sh   canu.10.sh   canu.11.sh   canu.12.sh
```

The assembled contigs are in Bt2.contigs.fasta.  Recall we set the -p prefix parameter to Bt2.  Other useful files include:

* Bt2.unassembled.fasta   These are reads that weren't assembled.  
* Bt2.report              This file contains information about each step.


Below I grabbed the contigs that appear to be Circular representing plasmids.
```
grep ">" Bt2_assembly/Bt2.contigs.fasta  | grep Circular=yes
>tig00000137 len=515061 reads=2993 covStat=742.81 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
>tig00000138 len=344327 reads=1984 covStat=496.27 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
>tig00000162 len=99169 reads=378 covStat=251.39 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
>tig00000163 len=102856 reads=497 covStat=185.89 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
>tig00000167 len=52597 reads=520 covStat=-94.48 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
>tig00000168 len=14459 reads=116 covStat=-11.79 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
>tig00000169 len=24282 reads=281 covStat=-75.26 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
>tig00000175 len=8489 reads=267 covStat=-146.76 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
>tig00002887 len=15142 reads=37 covStat=32.30 gappedBases=no class=contig suggestRepeat=no suggestCircular=yes
```


## Additional Reading

Below are links to canu documentation, github repo and other relevant websites related to canu.

* [Canu Tutorial](http://canu.readthedocs.io/en/latest/tutorial.html) Includes most of the links below and is fairly comprehensive
* [Canu Quick Start](http://canu.readthedocs.io/en/latest/quick-start.html)  When you are too impatient to read the documentation
* [Canu Github](https://github.com/marbl/canu)
* [Canu Faq Page](https://canu.readthedocs.io/en/latest/faq.html#) Really useful information about Canu
* [Read the Canu Paper](http://biorxiv.org/content/early/2016/08/24/071282)
