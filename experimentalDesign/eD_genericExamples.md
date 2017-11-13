# **RNA-Seq Analysis Generic Example:** experiments with 192 samples or less

#### Background
 This example assumes genomic resources already exist for the organism under study.  For example, the genome is assembled, annotated and available.  

 Consider the following experiment where an organism with 4 strains/lines/individuals is grown under 2 conditions (control and treated) with 3 time points and the experiment is biologically replicated 8 times for a total of ```4 strains x 2 conditions x 3 time points x 8 replicates = 192 samples```.  The number of strains, conditions and time points can be changed by factors or combinations of factors of 24 (ie 2,3,4,6,8,12) depending on experimental design.

 Currently, the max number of samples that can be indexed in a same lane for RNA-Seq is 24 and the max number of lanes in a flow cell is 8 (8x24=192). Therefore, each lane can act as a replicate of the entire experiment.  This has the added advantage of avoiding lane effects and since all samples fit on a single flow cell, chip effects are also avoided.

An average estimate for the size of the genic space in a genome  assumes a genome with ~30,000 genes with an average gene size of 1000 bases.  

#### Ideal Experimental Design Elements for $136,000
* ```Genome Assembled``` Yes
* ```Cost for Sequencing Lanes```: $22,720  (8 lanes * $2,840)
* ```Cost for Library Prep```: $31,680 (192 libraries * $165)
* ```Total Sequencing Cost```: $54,400
* ```Cost for Bioinformatics```: $81,600 (1.5 x sequencing cost)
* ```Total Project Cost```: $136,000
* ```Sequencing Technology```: Illumina HiSeq 3000
* ```Assumed sequencing output```: 300 million fragments/lane
* ```Number of lanes```: 8 lanes
* ```Length of read```: 150bp
* ```Number of Samples```: 192 (4 strains x 2 conditions x 3 time points x 8 replicates)
* ```Coverage Depth per sample```: 12.5 million fragments per lane and on average ~416 fragments/gene
* ```Number of Replicates```: 8
* ```Cost and technology as of```:  2017

#### Ideal Experimental Design Elements for less than $20,000
* ```Genome Assembled``` Yes
* ```Cost for Sequencing Lanes```: $2,840  (1 lanes * $2,840)
* ```Cost for Library Prep```: $3960 (24 libraries * $165)
* ```Total Sequencing Cost```: $6,800
* ```Cost for Bioinformatics```: $10,200 (1.5 x sequencing cost)
* ```Total Project Cost```: $17,000
* ```Sequencing Technology```: Illumina HiSeq 3000
* ```Assumed sequencing output```: 300 million fragments/lane
* ```Number of lanes```: 1 lanes
* ```Length of read```: 150bp
* ```Number of Samples```: 24 (1 strains x 3 conditions  8 replicates)
* ```Coverage Depth per sample```: 12.5 million fragments per lane and on average ~416 fragments/gene
* ```Number of Replicates```: 8
* ```Cost and technology as of```:  2017

# Genome Assembly Generic Example


[Table of contents](index.md)
