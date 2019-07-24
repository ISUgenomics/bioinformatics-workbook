---
title: Canu genome assembly of Bacillus thuringiensis
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


## Learning Objectives

Upon completion of this section on Genome Assembly you will be able to

* Download data from SRA  
* Evaluate the correct parameters to use when running canu for genome assembly
* Execute Canu and obtain a the assembled genome of [Bacillus thuringiensis](https://en.wikipedia.org/wiki/Bacillus_thuringiensis)

Let us assemble a bacterial genome [Bacillus thuringiensis strain: HS18-1](https://www.ncbi.nlm.nih.gov/sra/?term=SRR2093876). This particular genome is interesting because it shows insecticidal activity against Diptera.  The raw data contains 1.5G of data and the genome size is 6.4 mb.  The assembly can be found on [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_001182785.1).  Therefore, we will have over 200x coverage, plenty for a genome assembly with canu.  Canu recommends 30-60x minimum coverage.  The paper for this assembly can be read [here](https://www.sciencedirect.com/science/article/pii/S0168165615300961)

#### Download the data from ENA

The European Nucleotide Archive has the files already in fastq format and it is easy to download.
The main page with this project data is [PRJNA288953](https://www.ebi.ac.uk/ena/data/view/PRJNA288953)
We are only going to take the pacbio data which is the 3rd file link under the '''Read Files''' tab.


```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/006/SRR2093876/SRR2093876_subreads.fastq.gz
```

#### Running Canu

The name of the module will vary here and you should check to see what version you are using.  In this tutorial we used version 1.8 with java/1.8.0_121

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

You can get a sense of progress on your larger datasets by looking at the canu-scripts folder.  The final script is canu.12.sh.  

```
Bt2_assembly/canu-scripts/
canu.01.out  canu.02.out  canu.03.out  canu.04.out  canu.05.out  canu.06.out  canu.07.out  canu.08.out  canu.09.out  canu.10.out  canu.11.out  canu.12.out
canu.01.sh   canu.02.sh   canu.03.sh   canu.04.sh   canu.05.sh   canu.06.sh   canu.07.sh   canu.08.sh   canu.09.sh   canu.10.sh   canu.11.sh   canu.12.sh
```

The assembled contigs are in Bt2.contigs.fasta.  Recall we set the -p prefix parameter to Bt2.  Other useful files include:

* Bt2.unassembled.fasta   These are reads that weren't assembled.  
* Bt2.report              This file contains information about each step.


Below I grabbed the contigs that appear to be Circular representing plasmids.  This is the same number of plasmids reported by the [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_001182785.1) assembly.

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

Below are links to canu documentation, GitHub repo and other relevant websites related to canu.

* [Canu Tutorial](http://canu.readthedocs.io/en/latest/tutorial.html) Includes most of the links below and is fairly comprehensive
* [Canu Quick Start](http://canu.readthedocs.io/en/latest/quick-start.html)  When you are too impatient to read the documentation
* [Canu Github](https://github.com/marbl/canu)
* [Canu Faq Page](https://canu.readthedocs.io/en/latest/faq.html#) Really useful information about Canu
* [Read the Canu Paper](http://biorxiv.org/content/early/2016/08/24/071282)

---

* [Bacillus thuringiensis data set Info](BT_background.md)
* [Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
