---
title: Index page for Genome Assembly and Annotation
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Genome Assembly

---

**Genome Assembly** is the process of using DNA sequencing data to generate a representation of bases contained in chromosomes or the full genomes of an organism in the proper order and orientation.   

### [Introduction to Genome Assembly](../GenomeAssembly/Intro_GenomeAssembly.md)
  * [Introduction to Canu](../GenomeAssembly/Assemblers/canu.md)
    * [Canu on XSEDE Bridges Super Computer](../GenomeAssembly/BT/BT_Canu_bridges.md)
  * [Introduction to SPAdes](../GenomeAssembly/Assemblers/spades.md)
  * [Introduction to MaSuRCA](../GenomeAssembly/Assemblers/MaSuRCA.md)

### Genome Assembly Examples

  * [Bacillus thuringiensis data set](../GenomeAssembly/BT/BT_background.md)
    * [Canu Assembly of Bacillus thuringiensis](../GenomeAssembly/BT/BT_Canu.md)
    * [SPAdes Assembly of Bacillus thuringiensis](../GenomeAssembly/BT/BT_spades.md)
  * [Arabidopsis thaliana data set](../GenomeAssembly/Arabidopsis/Arabidopsis_background.md)
    * [MaSuRCA assembly of Arabidopsis thaliana](../GenomeAssembly/Arabidopsis/AT_MaSuRCA.md)
    * [Platanus assembly of Arabidopsis thaliana](../GenomeAssembly/Arabidopsis/AT_platanus-genome-assembly.md)
  * [Iteration of Pilon polishing and gap-filling on a genome](../GenomeAssembly/Iterate_Pilon_Genome_Polishing.md)
  * [Iterating a long read assembly to get higher contiguity by eliminating contaminant reads](../GenomeAssembly/IteratingGenomeAssemblyWithReadFiltration.md)
  * [Extracting and Assembling a Mitochondrial and Chloroplast Genome from Nuclear-targeted Nanopore](../GenomeAssembly/Organellar_Genome_Assembly_From_Nuclear-Targeted_Reads.md) 

### Tools for assessing the quality of a Genome Assembly

  * [GenomeScope to Estimate Genome Size](../GenomeAssembly/genomescope.md)
  * [Checking a genome for contamination from vectors using UniVec](../GenomeAssembly/univecContaminationCheck.md)
  * [Check a genome for PhiX contamination](../GenomeAssembly/PhiXContaminationCheck.md)

### Tools for Scaffolding assemblies
  * [Hi-C scaffolding](../GenomeAssembly/Hybrid/Scaffolding_with_HiC_Juicer.md)
  * [Hi-C scaffolding Update](../GenomeAssembly/Hybrid/Juicer_Juicebox_3dDNA_pipeline.md)
  * [Genome Scaffolding Using Synteny with Pyscaf](../GenomeAssembly/Pyscaf_Synteny_Scaffolding.md)

### Genetic Map Construction

  * [Generating Genetic Maps from GBS data](../GenomeAssembly/GeneticMaps/creating-genetic-maps.md)
  * [Using Genetic Map to create AGPs](../GenomeAssembly/GeneticMaps/scaffolding-using-genetic-maps.md)

# Genome Annotation

---

**Genome Annotation** has two separate but related definitions but is often used to mean both:

1. The process of identifying the location of genes by predicting the coding regions in a genome and generating gene models that represent the structure of a gene (start, stop, intron-exon boundaries, regulatory sequences, repeats).

2. The process of assigning a function to the gene models (gene names, protein products, domain structure)


### Introduction to Genome Annotation
  * [Calling Genome Methylation with Nanopolish and Comparing Promoter Methylation among Samples](Calling_Genome_Methylation_with_Nanopore.md)
  * [Introduction to Maker Gene Prediction](Intro_To_Maker.md)
  * [Introduction to Braker2 Gene Prediction](Intro_to_Braker2.md)
  * [Tutorial for NCBI PGAP](PGAP_tutorial.md)
  * [Motif Identification and Finding with MEME and FIMO](MEME_Motif_Finding_In_Genomes.md)  
  * [Assign GO Terms to Proteins using Deep Learning with DeepGoPlus](DeepGoPlus_AI_Functional_Prediction_of_Proteins.md)
  * [How to Identify the Secreted Protein from an Annotation and Predict Subcellular Localization](Secreted_Protein_Prediction_with_SignalP_and_TMHMM.md)
  * [How to Functionally Classify Proteins using ProtTrans, a Kinase Example ](Protein_Classification_with_ProtTrans)
