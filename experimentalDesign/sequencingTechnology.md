---
title: "Sequencing Technology"
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Sequencing Technology
There are now three main sequencing technologies that are available and commonly used: Illumina, PacBio and Oxford Nanopore.  Understanding the assumptions and limitations of each of these technologies can aid in planning the experimental design.

---
## **[Illumina](https://www.illumina.com/systems/sequencing-platforms/novaseq/specifications.html)**  


Illumina raw data are short (100-300bp) in size and of high quality for reads shorter than 200 bps.  Quality scores for bases on reads between 250-300bp usually are of significant lower quality.  The quality of the read diminishes as the length of the read increases.  This trend of of quality does not change with the length of the run.

### Number of fragments to expect to pass the filter

The NovaSeq 6000 is Illumina's latest machine and has significantly higher output than previous generations of sequencing machines.  The amount of output is tied to the `flow cell type`.  The numbers we need to know are below represented as fragments (single or paired).  For the paired case, if you want to know the number of reads then just multiply by 2.

|NovaSeq 6000 System| | | M =| Millions of Fragments|
| -- | -- | -- |-- | -- |  
|Flow Cell Type	|SP	|S1	|S2	|S4|
|Number of fragments	|650–800 M|	1300–1600 M|	3300 - 4100 M|	8000 - 10000 M|


### Approximate number of samples you could run with each type of flow cell by application

This assumes 60-80X coverage per genome run.

| NovaSeq 6000 System|  |  |  | |
| -- | -- | -- |-- | -- |  
| Flow Cell Type| SP 	| S1| 	S2| 	S4|
| 3Gb Genomes per Run	| 4	| 8	| 20	| 48|
| 1Gb Genomes per Run	| 12	| 24	| 60	| 144|
| Exomes per Run| 	40	| 80	| 200| 	500|
|Transcriptomes per Run	|32|	64	|164	|400|

### Read lengths and output at that read length

| Flow Cell Type|	SP	|S1	|S2	|S4|
| -- | -- | -- |-- | -- |  
|1 × 35 bp	|No	|No	|No	|280-350 Gb|
|2 × 50 bp	|65–80 Gb	|134–167 Gb	|333–417 Gb|	No|
|2 × 100 bp	|134–167 Gb	|266–333 Gb	|667–833 Gb|	1600–2000 Gb|
|2 × 150 bp	|200–250 Gb	|400–500 Gb	|1000–1250 Gb|	2400–3000 Gb|
|2 x 250 bp	|325-400 Gb	|No	|No|No|

#### Some useful links to Illumina related information

* [Sequencing rates at a service provider](http://www.biotech.iastate.edu/biotechnology-service-facilities/dna-facility/#rates)

#### Video explanation

* [![Illlumina Video](https://img.youtube.com/vi/fCd6B5HRaZ8/0.jpg)](https://www.youtube.com/watch?v=fCd6B5HRaZ8)


---
## **[PacBio](https://www.pacb.com/wp-content/uploads/Sequel-II-System-v8.0-and-SMRT-Link-v8.0-Technical-Overview-Customer-Training.pdf)**


PacBio raw data are long (~13,000-20,000bp) with max read lengths around 300,000 bp.

* HiFi = High Fidelity reads which are produced from a high fidelity polymerase making for longer more accurate reads by reading the same molecule multiple times.  
* CLR = Continuous Long Reads, read once but capable of reading much longer reads.

| System | Gb | Millions of Reads|
| -- | -- | --|
|Sequel II | ~100  | ~400 | HiFi|
|Sequel II | ~50  | ~40 | CLR |
|Sequel I  | ~15 | ~0.5  | HiFi |


#### Some useful links to Pacbio related information

*  [Table of Application-Options-and-Sequencing-Recommendations](https://www.pacb.com/wp-content/uploads/Overview-Sequel-Systems-Application-Options-and-Sequencing-Recommendations.pdf)
*  [Preparing samples](https://www.pacb.com/wp-content/uploads/Technical-note-Preparing-Samples-for-PacBio-Whole-Genome-Sequencing-for-de-novo-Assembly-collection-and-storage.pdf)
*  [DNA extraction for Pacbio](https://www.pacb.com/wp-content/uploads/Technical-Note-Preparing-DNA-for-PacBio-HiFi-Sequencing-Extraction-and-Quality-Control.pdf)
*  [Intro to PacBio from UCDavis](https://dnatech.genomecenter.ucdavis.edu/pacbio-library-prep-sequencing/)
* [Sequencing rates at a service provider](https://dnatech.genomecenter.ucdavis.edu/industry-rates/)
* [Reference for estimated output](https://i2.wp.com/www.dnalinkseqlab.com/wp-content/uploads/2019/07/Sequel-ii-Sequel-i-chart-2-2.jpg?w=1182&ssl=1)

#### Video explanation

* [![Pacbio Video](https://img.youtube.com/vi/v8p4ph2MAvI/0.jpg)](https://www.youtube.com/watch?v=v8p4ph2MAvI&feature=youtu.be)




#### Notes
- Multiplex up to 48 microbial samples per SMRT Cell 8M

---
## **[Nanopore](https://nanoporetech.com/products/comparison)**

Nanopore raw data are long (10,000 - 30,000 bp) with the longest confirmed read of 2.3 million bases.  Nanopore is the fastest evolving of the three sequencing technologies and therefore this data is continuously becoming outdated.  In December of 2020, a huge jump in base calling quality was [announced](https://twitter.com/nanopore/status/1334238633008754688) with the mean above Q20 (99.13%) using the base caller [Bonito](https://github.com/nanoporetech/bonito)


| System | Gb | Millions of Reads|
| -- | -- | --|
|Minion| ~40 | ~2.5 |
|Promethion| ~180| 11.5 |


#### Some useful links to Nanopore related information

 [Nanopore Information](https://nanoporetech.com/products/comparison)
 - This paper provides a nice overview of MinIon sequencing technologies and uses [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1103-0)
 - [Opportunities and challenges in long-read sequencing data analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1935-5#:~:text=Nanopore%20sequencing%20provides%20the%20longest,kb%20genomic%20libraries%20being%20common.)

#### Video Explanation

* * [![Nanopore Video](https://img.youtube.com/vi/E9-Rm5AoZGw/0.jpg)](https://www.youtube.com/watch?v=E9-Rm5AoZGw)


---

[Sequencing rates at a service provider](http://www.biotech.iastate.edu/biotechnology-service-facilities/dna-facility/#rates)

## Funding and Cost
Most research has a strict allowance for how much sequencing and bioinformatics can be performed to answer the biological question of interest. An understanding of the following terminology can aid in determining the type and amount of sequencing that is best suited for your biological purpose.

 - **Read length:**```Short reads (50bp) are difficult to align to unique locations in a genome, so unless the experiment is for smRNA it is uncommon to use very short reads.```


 - **Paired-end** ```Both ends of the DNA fragment are sequenced.  This type of sequencing is useful for obtaining more unique alignments to a genome  For RNA-Seq experiments with a known genome, it is recommended to use at least 100bp paired-end Illumina data.  For RNA-Seq experiments without a genome or a genome of questionable quality, it recommended to use 150bp Illumina paired-end  data. ```


 - **Single-end** ```Used when the experiment has DNA fragments shorter than the length of the read.  For example, smRNA experiments are typically done with 50bp single-end data. ```


 - **Biological Replicates**  ```It is extremely important to have at least 3 replicates and preferably 5 to 10 replicates for RNA-Seq experiments to determine differential expression```

---

## Examples

In the next sections we will go over several example experimental design problems from real world examples.

[Next](costs.md){: .btn  .btn--primary}
[Previous](exp_design_index.md){: .btn  .btn--primary}
[Table of contents](exp_design_index){: .btn  .btn--primary}
