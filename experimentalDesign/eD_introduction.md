# Introduction to Experimental Design

Bioinformatics and data analysis starts with a good experimental design.  In this section we will explore the aspects of experimental design that are important to a successful bioinformatic data analysis.  


## Learning Objective
Upon completion of this section on experimental design students will be able to:

* assess the assumptions and limitations of the biological domain and the sequencing technology as it relates to a biology question of interest
* evaluate which sequencing technology will provide the best results for a given project
* determine the limitations that result from reducing the amount of data collected in order to fit within a given budget.


---

## What is your biological question?

First and foremost you should ask yourself what is it that you are trying to answer with high throughput sequencing data.  A hypothesis driven experiment will always be more insightful and easier to design than an experiment without direction.  

**Fishing Expedition:**  ```
Common feedback given by grant review panelists that feel a researcher is sequencing samples without a clear direction or hypothesis.```

There are four main aspects that should be considered to answer your biological question during experimental design: biological system, sequencing technology and cost.

---

## Biological System

- **Domain:**  ```The organism or group of organisms that a researcher or group of researchers are studying.```

- **Polyploidy** ```Multiple copies of the same genome contained in a single nucleus due to a recent whole genome duplication event ```

- **Paleopolyploidy**  ```An organism that has had a whole genome duplication in the ancient past (millions of years) but has rediploidized```

- **Repeat Content**   ```The amount of highly repetitive sequences contained in the genome.```

- **GC Content** ```The percent GC or AT in your genome can affect the quality of the sequencing data that you acquire.  Some sequencing technologies and assembly programs have trouble with AT rich genomes.```

---

## Sequencing Technology
There are now three main sequencing technologies that are available and commonly used: Illumina, PacBio and MinIon.  Understanding the assumptions and limitations of each of these technologies can aid in planning the experimental design.

- **Illumina**  
Illumina raw data are short (100-300bp) in size and of high quality for reads shorter than 200 bps.  Quality scores for bases on reads between 250-300bp usually are of lower quality.
 - **Size:** ```short```
 - **Output:** ```90,000 Gb (HiSeq 3000-4000)```
 - **Base Quality:** ```99.9% High```


- **PacBio**
PacBio raw data (subreads) are long (3,000-15,000bp) with max read lengths around 40,000-80,000bp.
 - **Size:** ```long```
 - **Output**: ```1 Gb (RSII)```
 - **Output**: ```1-5 Gb (Sequel)```
 - **Base Quality:** ```85% Medium```


- **MinIon**
MinIon raw data are long (3,000-70,000bp) with max read lengths as high as 250,000bp.
 - **Size:** ```long```
 - **Output**: ```1 Gb (RSII)```
 - **Output**: ```1-5 Gb (Sequel)```
 - **Base Quality:** ```85% Medium```
 - This paper provides a nice overview of MinIon sequencing technologies and uses [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1103-0)

---

## Funding and Cost
Most research has a strict allowance for how much sequencing and bioinformatics can be performed to answer the biological question of interest.

 - **Read length:**```Short reads (50bp) are difficult to align to unique locations in a genome, so unless the experiment is for smRNA it is uncommon to use very short reads.```


 - **Paired-end** ```Both ends of the DNA fragment are sequenced.  This type of sequencing is useful for obtaining more unique alignments to a genome  For RNA-Seq experiments with a known genome, it is recommended to use at least 100bp paired-end Illumina data.  For RNA-Seq experiments without a genome or a genome of questionable quality, it recommended to use 150bp Illumina paired-end  data. ```


 - **Single-end** ```Used when the experiment has DNA fragments shorter than the length of the read.  For example, smRNA experiments are typically done with 50bp single-end data. ```


 - **Biological Replicates**  ```It is extremely important to have at least 3 replicates and preferably 5 to 10 replicates for RNA-Seq experiments to determine differential expression```

---

## Examples

In the next sections we will go over several example experimental design problems from real world examples.

[Table of contents](README.md)
