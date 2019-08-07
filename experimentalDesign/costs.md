---
title: "Experimental Design"
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

{% include toc %}

## Funding and Cost
Most research has a strict allowance for how much sequencing and bioinformatics can be performed to answer the biological question of interest. An understanding of the following terminology can aid in determining the type and amount of sequencing that is best suited for your biological purpose.
{: style="text-align: justify"}

 - **Read length:**```Short reads (50bp) are difficult to align to unique locations in a genome, so unless the experiment is for smRNA it is uncommon to use very short reads.```
 {: style="text-align: justify"}


 - **Paired-end** ```Both ends of the DNA fragment are sequenced.  This type of sequencing is useful for obtaining more unique alignments to a genome  For RNA-Seq experiments with a known genome, it is recommended to use at least 100bp paired-end Illumina data.  For RNA-Seq experiments without a genome or a genome of questionable quality, it recommended to use 150bp Illumina paired-end  data. ```
 {: style="text-align: justify"}


 - **Single-end** ```Used when the experiment has DNA fragments shorter than the length of the read.  For example, smRNA experiments are typically done with 50bp single-end data. ```
 {: style="text-align: justify"}


 - **Biological Replicates**  ```It is extremely important to have at least 3 replicates and preferably 5 to 10 replicates for RNA-Seq experiments to determine differential expression```
 {: style="text-align: justify"}

---

## Examples

In the next sections we will go over several example experimental design problems from real world examples.
{: style="text-align: justify"}



[Next](eD_genericExamples.md){: .btn  .btn--primary}
[Previous](sequencingTechnology.md){: .btn  .btn--primary}
[Table of contents](exp_design_index){: .btn  .btn--primary}
