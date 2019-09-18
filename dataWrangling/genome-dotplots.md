# DotPlots for comparing genomes

One of the primary comparative analyses that can be done once you have the genome is by visualizing the synteny with closely related species. Many characteristics for the genome can be easily highlighted with a good dotplot. Structural variations like inversion, deletion, duplication and insertions can be identified from these dotplots.

You will need two genomes for generating dotplots. A higher quality, preferably at chromosome level "reference" genome (also referred as target genome) and your genome (scaffold or contigs are okay, but chromosomes would be ideal), which is referred as query genome. We will run this tutorial using maize genomes, but can be easily applied on any other genome as well.


## Data

Download the data from [Grameme](http://ensembl.gramene.org/Zea_mays/Info/Index) and [MaizeGDB](https://www.maizegdb.org)

```bash
# B73 version 4.0 (Target)
wget ftp://ftp.gramene.org/pub/gramene/release-61/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
# Genome of interest (Query)
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML247-REFERENCE-PANZEA-1.1/Zm-CML247-REFERENCE-PANZEA-1.1.fa.gz
# maize custom TE libraries
wget https://github.com/mcstitzer/dawe_ab10_kindr/blob/master/identify_haplotypes/te_alignments/TE_12-Feb-2015_15-35.fa
gunzip *.gz
```

## Overview

![DotPlots](R/assets/dotplots.png)

Figure 1: Overview of the DotPlot construction method


## Organization

```
SyntenyTutorial
├── 1_data
│   ├── Zea_mays.B73_RefGen_v4.dna.toplevel.fa
│   └── Zm-CML247-REFERENCE-PANZEA-1.1.fa
├── 2_repeatmasking
│   ├── runRepeatMasker.sh
│   └── TE_12-Feb-2015_15-35.fa
├── 3_minimap
│   └── runMinimap.sh
└── 4_paf-processing
    └── pafCoordsDotPlotly.R
```

## Programs needed

The o

## Repeatmask the genomes

Setup a script as follows:

```
