---
title: "Glossary"

permalink: /glossary.html
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

{% include toc %}

<!--template
<span style="color: #518cc2; font-weight:bold">TERM</span> -
<span style="color: #24376b;">
definition
</span><br>
<span style="color: #8997c1;"> _#hashtag1_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#hashtag2_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#hashtag3_ </span> &nbsp;|&nbsp;
<br>
<br>
-->


## A


<span style="color: #518cc2; font-weight:bold">Alignment</span> [of sequencing reads to reference genome] -
<span style="color: #24376b;">
in RNA-seq analysis, means detecting the presence and location (position) of individual reads (*fragments of a sequence obtained from sequencing experiment*) within the sequence space of a reference genome. In particular, sequenced RNA can be aligned to the reference genome to identify genes and get information about gene expression.
</span><br>
<span style="color: #8997c1;"> _#RNA-seq_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#reads mapping_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#gene expression_ </span> &nbsp;|&nbsp;
<br>
<br>

## B



## C



## D



## E



## F



## G

<span style="color: #518cc2; font-weight:bold">Genome Index</span> -
<span style="color: #24376b;">
is a data structure of a compressed full-text file containing the reference genome (e.g., .fna file). Using the genome index makes it efficient when searching a substring (e.g., matching reads) in a large text. Programs such as <a href="https://daehwankimlab.github.io/hisat2/manual/" target="_blank">HISAT2  ⤴</a> (<i>hisat2-build</i> indexer) build the reference genome index using the <a href="https://en.wikipedia.org/wiki/FM-index" target="_blank">FM-index  ⤴</a> approach, where the data is both compressed and indexed to reasonably fit within a computer's memory.
</span><br>
<span style="color: #8997c1;"> _#RNA-seq_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#mapping_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#alignment_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#sequencing_ </span> &nbsp;|&nbsp;
<br>
<br>

## H



## I



## J



## K



## L



## M

<span style="color: #518cc2; font-weight:bold">Mapping</span> [sequencing reads to reference genome] -
<span style="color: #24376b;">
in RNA-seq analysis, means detecting the presence and location (position) of individual reads (*fragments of a sequence obtained from sequencing experiment*) within the sequence space of a reference genome. In particular, sequenced RNA can be mapped to the reference genome to identify genes and get information about gene expression.
</span><br>
<span style="color: #8997c1;"> _#RNA-seq_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#reads alignment_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#gene expression_ </span> &nbsp;|&nbsp;
<br>
<br>

## N



## O



## P



## Q



## R



## S

<span style="color: #518cc2; font-weight:bold">Soft-Clip</span> [of reads] -
<span style="color: #24376b;">
means ignoring the terminal fragments (ends) of the reads that do not match perfectly to the reference genome alignment. This procedure enables higher mapping efficiency and facilitates detecting structural variants. However, this also bears the danger of incorrectly trimming the reads, leading to the misassignment of reads to repetitive regions. [<a href="https://sequencing.qcfail.com/articles/soft-clipping-of-reads-may-add-potentially-unwanted-alignments-to-repetitive-regions/" target="_blank">Learn more  ⤴</a>]
</span><br>
<span style="color: #8997c1;"> _#RNA-seq_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#reads mapping_ </span> &nbsp;|&nbsp;
<span style="color: #8997c1;"> _#reads trimming_ </span> &nbsp;|&nbsp;
<br>
<br>

## T



## U



## V



## W



## X



## Y



## Z
