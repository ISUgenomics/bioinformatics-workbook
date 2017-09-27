# Introduction to Experimental Design

Bioinformatics and data analysis starts with a good experimental design.  In this section we will explore the aspects of experimental design that are important to a successful bioinformatic data analysis.

## What is your biological question?

First and foremost you should ask yourself what is it that I am trying to answer with high throughput sequencing data.  A hypothesis driven experiment will always be more insightful and easier to design than an experiment without direction.  

**Fishing Expedition:**  ```
Common feedback given by grant review panelists that feel a researcher is sequencing samples without a clear direction or hypothesis.```

There are three main aspects that should be considered to answer your biological question during experimental design.

## Biological System

- **Domain:**  ```The organism or group of organisms that a researcher or group of researchers are studying.```

- **Polyploidy** ```Multiple copies of the same genome contained in a single nucleus due to a recent whole genome duplication event ```

- **Paleopolyploidy**  ```An organism that has had a whole genome duplication in the ancient past (millions of years) but has rediploidized```

- **Repeat Content**   ```The amount of highly repetitive sequences contained in the genome.```


## Sequencing Technology
There are now three main sequencing technologies that are available and commonly used: Illumina, PacBio and MinIon.  Understanding the assumptions and limitations of each of these technologies can aid in planning the experimental design.

- **Illumina**  
Illumina raw data are short (100-300bp) in size and of high quality for reads shorter than 200 bps.  Quality scores for bases on reads between 250-300bp usually are of lower quality.
 - **Size:** short
 - **Output:** 90,000 Gb (HiSeq 3000-4000)
 - **Base Quality:** 99.9% High


- **PacBio**
PacBio raw data (subreads) are long (3,000-15,000bp) with max read lengths around 40,000-80,000bp.
 - **Size:** long
 - **Output**: 1 Gb (RSII)
 - **Output**: 1-5 Gb (Sequel)
 - **Base Quality:** 85% Medium


- **MinIon**
MinIon raw data are long (3,000-70,000bp) with max read lengths as high as 250,000bp.
 - **Size:** long
 - **Output**: 1 Gb (RSII)
 - **Output**: 1-5 Gb (Sequel)
 - **Base Quality:** 85% Medium
 - This paper provides a nice overview of MinIon sequencing technologies and uses [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1103-0)

## Funding and Cost


## Examples

In the next sections we will go over several example experimental design problems from real world examples.
