---
title: FastTree
layout: single
author: Jennifer Chang
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# FastTree

A phylogenetic tree is a hypothesis of the evolutionary inheritance of genes across individual taxa. Trees have been used to summarize an organism's pedigree, infer viral host-spillover events ([Zhou et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32015507/)), and determine if the red panda is a bear or a raccoon ([Slattery & O'Brien, 1995](https://pubmed.ncbi.nlm.nih.gov/8568209/))

There are a number of phylogenetic tree building programs including BEAST, MrBayes, PAUP, PhyML, RAxML, IQ-Tree, and FastTree. We are focusing on FastTree because it tends to be a faster tree building program to provide a quick draft tree. We recommend running FastTree to get a general sense of individuals in the tree, subsample down or add references taxa to create a well formed tree, and then move on to the other tree-building programs for more detailed diversity and time-scaled analysis.

## Software required

* [MAFFT](https://mafft.cbrc.jp/alignment/software/) - for aligning your sequences
* [FastTree](http://tree.bio.ed.ac.uk/software/figtree/) - for building a quick phylogenetic tree from aligned sequences
  ```
  curl -O http://www.microbesonline.org/fasttree/FastTree.c
  gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
  ```
* [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) - for viewing your phylogenetic tree

## Overview of the pipeline:

```
  sequences -> MAFFT -> FastTree -> FigTree -> pdf figure that can go into your manuscript
```

## Step 1: MAFFT

```
mafft --auto input.fasta > input_aln.fasta
```

## Step 2: FastTree

```
# For a nucleotide alignment
FastTree -nt -gtr -gamma input_aln.fasta > input.tre

# For a protein alignment
FastTree input_aln.fasta > input.tre
```

## Step 3: FigTree

Open FigTree application and use `File/Open Tree` to select your `input.tre`. Usually the tree will be hard to read until you midpoint using the command keys `Ctrl-M` and  `Ctrl-U`.

## Citations

* Katoh, K., Misawa, K., Kuma, K.I. and Miyata, T., 2002. [MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform](https://pubmed.ncbi.nlm.nih.gov/12136088). Nucleic acids research, 30(14), pp.3059-3066.
* Price, M.N., Dehal, P.S. and Arkin, A.P., 2009. [FastTree: computing large minimum evolution trees with profiles instead of a distance matrix](https://pubmed.ncbi.nlm.nih.gov/19377059). Molecular biology and evolution, 26(7), pp.1641-1650.
* Rambaut, A., 2007. [FigTree, a graphical viewer of phylogenetic trees.](http://tree.bio.ed.ac.uk/software/figtree/)
