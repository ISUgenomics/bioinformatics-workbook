---
title: "Experimental Design Examples"
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---



# Sex Determination in a Fish

#### Background
Aquaculture production has become increasingly important to satisfy seafood and fishery product demands. Researcher wants to develop genomic resources to improve sustainable domestic aquaculture.  They recently assembled and annotated a genome and now want to collect wild fish to explore genetic diversity in local populations and identify the sex determining region.

#### Ideal Experimental Design Elements
* ```Funding for project``` $X,XXX
* ```Sequencing Technology```: Illumina HiSeq 3000
* ```Number of lanes```: 20 lanes
* ```Length of read```: 150bp
* ```Number of Samples```: 90
* ```Coverage Depth```: 29x
* ```Year of Sequencing```: hypothetical 2017

#### Real Experimental Design Elements
* ```Funding for project``` $X,XXX
* ```Sequencing Technology```: Illumina HiSeq 2500
* ```Number of lanes```: 4 lanes
* ```Length of read```: 100bp
* ```Number of Samples```: 90
* ```Coverage Depth```: 1.9x
* ```Year of Sequencing```: 2015


# Questions  

##### How is Coverage Depth Determined?

The quick and dirty calculation is to take the read length (150) multiply by 2 if it is paired end multiply by the number of fragments the Illumina machine produces (300 Million) = 150bp * 2 * 300M = 90,000 million or 90 billion bases output from a single lane.  Take that number and divide it by your genome size.  If your genome size is 1 gigabase then you have approximately 90x coverage.

While this can give an approximation for coverage depth the reality is that the coverage depth at any given locus in the genome will be a distribution around this number with some regions having more and some regions having less coverage.

[Table of contents](https://isugenomics.github.io/bioinformatics-workbook/)
