---
title: "Table of Contents"

permalink: /list.html
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

### [Command Line Basics and Useful Programs](https://isugenomics.github.io/bioinformatics-workbook/Appendix/programs)
* <span style="color:lightblue">Introduction to Unix</span>
  * [Unix Basics 1](../Appendix/Unix/unix-basics-1.md)  
  * [Introduction to Job Scheduling: SLURM](../Appendix/Unix/01_slurm-basics.md)
  * [Introduction to GNU parallel](../Appendix/GNUparallel/GNU_parallel_examples.md)
  * [Introduction to HPC](../Appendix/Unix/unix-basics-6HPC.md)
  * [Unix Admin commands](../Appendix/Unix/unix-basics-2admin.md)
  * [Grep](../Appendix/Unix/unix-basics-3grep.md)
  * [Sed](../Appendix/Unix/unix-basics-4sed.md)
  * [Introduction to regular expression](../introduction/introduction-to-regular-expressions.md)
* <span style="color:lightblue">HPC</span>
  * [How to Install a Program](../Appendix/HPC/guide-for-installing-various-types-of-programs-in-linux.md)
  * <span style="color:lightblue">SSH</span>
    * [SSH Shortcuts](../Appendix/HPC/ssh-shortcuts.md)
    * [Password-less SSH](../Appendix/HPC/password-less-ssh-login.md)
  * <span style="color:lightblue">SLURM</span>
    * [SLURM Cheatsheet](../Appendix/HPC/SLURM/slurm-cheatsheat.md)
    * [Creating SLURM job submission scripts](../Appendix/HPC/SLURM/creating-slurm-job-submission-scripts-for-condo.md)
    * [Submitting dependency jobs using SLURM](../Appendix/HPC/SLURM/submitting-dependency-jobs-using-slurm.md)
  * <span style="color:lightblue">XSEDE</span>
      * [XSEDE supercomputer](../Appendix/HPC/xsede/xsede.md)
      * [XSEDE supercell storage](../Appendix/HPC/xsede/using-psc-supercell-storage-for-bridges-and-greenfield.md)
    * <span style="color:lightblue">Containers</span>
      * [Introduction to Containers](../Appendix/HPC/Containers/Intro_Singularity.md)
      * [Creating Containers Using Singularity](../Appendix/HPC/Containers/creatingContainers.md)
      * [Modifying Existing Containers](../Appendix/HPC/Containers/modifyingExistingContainers.md)
* <span style="color:lightblue">Bioawk</span>
  * [Bioawk Basics](../Appendix/Unix/bioawk-basics.md)
* <span style="color:lightblue">[Viewing Files In Remote Machine Without Downloading](../Appendix/HPC/viewing-files-in-remote-machine-without-downloading-locally.md)</span>

### [Project Management](https://isugenomics.github.io/bioinformatics-workbook/projectManagement/projectManagement_index)
  * [Introduction to Project Management](https://isugenomics.github.io/bioinformatics-workbook/projectManagement/Intro_projectManagement)
  * <span style="color:lightblue">Project Management Tools</span>
    * [Introduction to Slack](../Appendix/slack.md)
    * [Introduction to GitHub](../Appendix/github/introgithub.md)
    * [Introduction to Markdown](../Appendix/Markdown.md)
    * [Zenhub Download](https://www.zenhub.com/extension)


### [Introduction to BLAST](../dataAnalysis/blast/blast_index.md)
* [BLAST example: finding a gene in a genome](../dataAnalysis/blast/blastExample.md)
* [Running BLAST in parallel](../dataAnalysis/blast/running-blast-jobs-in-parallel.md)

### [Experimental Design](../experimentalDesign/exp_design_index.md)
* [Biological system](../experimentalDesign/bio_sys.md)
* [Sequencing Technology](../experimentalDesign/sequencingTechnology.md)
* [Costs](../experimentalDesign/costs.md)
* [Generic Examples of Experimental Design](../experimentalDesign/eD_genericExamples.md)
* [List of Biology exceptions and irregularities](../Appendix/biology_tidbits.md)

### [Data Acquisition and Wrangling](../dataAcquisition/dAc_introduction.md)
* <span style="color:lightblue">File Transfer</span>
  * [File Transfer using wget](../dataAcquisition/fileTransfer/downloading-files-via-wget.md)
  * [File Transfer using Globus](../dataAcquisition/fileTransfer/file-transfer-using-globus-connect-personal-gcp.md)
  * [File Transfer using irods](../dataAcquisition/fileTransfer/getting-data-from-iplant-via-irods.md)
  * [File Transfer using SRA toolkit](../dataAcquisition/fileTransfer/sra.md)
* [Data Sets Used in the tutorials](../dataAcquisition/dataSets.md)
* <span style="color:lightblue">FASTA manipulation</span>
  * [Determining Sequence length](../dataWrangling/fastaq-manipulations/calculate-sequence-lengths-in-a-fasta-file.md)
  * [Converting FastQ to FASTA](../dataWrangling/fastaq-manipulations/converting-fastq-format-to-fasta.md)
  * [FASTQ Quality trimming](../dataWrangling/fastaq-manipulations/fastq-quality-trimming.md)
  * [Retrieving FASTA file using sequence ID](../dataWrangling/fastaq-manipulations/retrieve-fasta-sequences-using-sequence-ids.md)
* <span style="color:lightblue">Manipulating Excel data sheets</span>
  * [Create Workbook from Multiple Text Files](../dataWrangling/microsoftExcel/import-multiple-text-files-as-separate-worksheets-in-excel.md)
  * [Export multiple worksheets as separate text files ](../dataWrangling/microsoftExcel/export-multiple-worksheets-as-separate-text-files-in-excel.md)
  * [Create Index for All Worksheets](../dataWrangling/microsoftExcel/generate-index-sheet-linking-all-spreadsheets-in-excel.md)
  * [Merge two spreadsheets using a common column](../dataWrangling/microsoftExcel/Merge_two_spreadsheets_using_a_common_column_in_Excel.md)
* <span style="color:lightblue">Data Management</span>
  * [Deposition of Data to NCBI SRA](../dataWrangling/NCBI_Data_Submission.md)

### [Bioinformatics terminology](../introduction/terminology_index.md)
* [Reads, Contigs, Scaffolds and Chromosome](../introduction/dataTerminology.md)
* [File Formats](../introduction/fileFormats.md)
* [Fasta Quality Score](../introduction/fastqquality-score-encoding.md)

### [RNA Sequencing](../dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-intro.md)
* [RNA-Seq Example with a Genome Assembly](../dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-using-a-genome.md)
* [RNA-Seq Example without a Genome Assembly](../dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-without-a-genome.md)
* [Different Expression Analysis:DESeq2](../dataAnalysis/RNA-Seq/RNA-SeqIntro/Differential-Expression-Analysis.md)
* [10x genomics single-cell RNAseq analysis from SRA data using Cell Ranger and Seurat](../dataAnalysis/RNA-Seq/Single_Cell_RNAseq/Chromium_Cell_Ranger.md)


### [Genome Assembly and Annotation](../dataAnalysis/GenomeAnnotation/annotation_and_assembly_index.md)
#### <span style="color:lightblue">[Introduction to Genome Assembly](../dataAnalysis/GenomeAssembly/Intro_GenomeAssembly.md)</span>
* [Introduction to Canu](../dataAnalysis/GenomeAssembly/Assemblers/canu.md)
  * [Canu on XSEDE Bridges Super Computer](../dataAnalysis/GenomeAssembly/BT/BT_Canu_bridges.md)
* [Introduction to SPAdes](../dataAnalysis/GenomeAssembly/Assemblers/spades.md)
* [Introduction to MaSuRCA](../dataAnalysis/GenomeAssembly/Assemblers/MaSuRCA.md)  

#### <span style="color:lightblue">Genome Assembly Examples</span>
* [Bacillus thuringiensis data set](../dataAnalysis/GenomeAssembly/BT/BT_background.md)
  * [Canu Assembly of Bacillus thuringiensis](../dataAnalysis/GenomeAssembly/BT/BT_Canu.md)
  * [SPAdes Assembly of Bacillus thuringiensis](../dataAnalysis/GenomeAssembly/BT/BT_spades.md)
* [Arabidopsis thaliana data set](../dataAnalysis/GenomeAssembly/Arabidopsis/Arabidopsis_background.md)
  * [MaSuRCA assembly of Arabidopsis thanliana](../dataAnalysis/GenomeAssembly/Arabidopsis/AT_MaSuRCA.md)
  * [Short read assembly using Platanus](../dataAnalysis/GenomeAssembly/Arabidopsis/AT_platanus-genome-assembly.md)

#### <span style="color:lightblue">Tools for assessing the quality of a Genome Assembly</span>
* [GenomeScope to Estimate Genome Size](../dataAnalysis/GenomeAssembly/genomescope.md)
* [Checking a genome for contamination from vectors using UniVec](../dataAnalysis/GenomeAssembly/univecContaminationCheck.md)
* [Check a genome for PhiX contamination](../dataAnalysis/GenomeAssembly/PhiXContaminationCheck.md)

#### <span style="color:lightblue">Tools for Scaffolding assemblies</span>  
* [Genome Scaffolding Using Synteny with Pyscaf](../dataAnalysis/GenomeAssembly/Pyscaf_Synteny_Scaffolding.md)
* [Hi-C scaffolding](../dataAnalysis/GenomeAssembly/Hybrid/Scaffolding_with_HiC_Juicer.md)

#### <span style="color:lightblue">Genetic Map Construction</span>

* [Generating Genetic Maps from GBS data](../dataAnalysis/GenomeAssembly/GeneticMaps/creating-genetic-maps.md)
* [Using Genetic Map to create AGPs](../dataAnalysis/GenomeAssembly/GeneticMaps/scaffolding-using-genetic-maps.md)

#### <span style="color:lightblue">Introduction to Genome Annotation</span>
* [Introduction to Maker Gene Prediction](../dataAnalysis/GenomeAnnotation/Intro_To_Maker.md)
* [Introduction to Braker2 Gene Prediction](../dataAnalysis/GenomeAnnotation/Intro_to_Braker2.md)
* [Motif Identification and Finding with MEME and FIMO](../dataAnalysis/GenomeAnnotation/MEME_Motif_Finding_In_Genomes.md)  


### [Comparative Genomics](../dataAnalysis/ComparativeGenomics/compGenomics_index.md)
* [Gene Orthology, Synteny, and Visualzation with Opscan, Iadhore and Circos](../dataAnalysis/ComparativeGenomics/Gene_Orthology_And_Synteny.md)
* [Gene Orthology, Synteny, and Visualzation with Orthofinder, Iadhore and Circos](../dataAnalysis/ComparativeGenomics/OrthofinderSynteny.md)
* [Gene Overlap Significance Testing with R Gene_overlap Package](../dataAnalysis/ComparativeGenomics/Gene_Category_Overlap_Fishers_exact_testing.md)  
* [Phylostratiophraphy:Determining the LCA of all Genes in a Genome](../dataAnalysis/ComparativeGenomics/phylostratr.md)


### [Variant Discovery](../dataAnalysis/VariantCalling/variant-calling-index.md)
* [FreeBayes variant calling workflow for DNA-Seq](../dataAnalysis/VariantCalling/freebayes-dnaseq-workflow.md)
* [GATK Best Practices Workflow for DNA-Seq](../dataAnalysis/VariantCalling/gatk-dnaseq-best-practices-workflow.md)
* [SNP calling for GBS data using Stacks pipeline](../dataAnalysis/VariantCalling/gbs-data-snp-calling-using-stacks.md)
* [SNP calling for GBS data using Tassel pipeline (GBS.v2)](../dataAnalysis/VariantCalling/gbs-data-snp-calling-using-tassel.md)

### [Metagenomics](../dataAnalysis/Metagenomics/metagenomics_index.md)
* [Quality control, assembly and mapping](../dataAnalysis/Metagenomics/MetagenomicsP1.md)
* [Qiime2](../dataAnalysis/Metagenomics/Qiime2.md)
* [DADA2](../dataAnalysis/Metagenomics/Dada2.md)

### [Genome Repeat Identification](../dataAnalysis/ComparativeGenomics/Repeats_index.md)
* [Helitron Identification in a Genome Sequence](../dataAnalysis/ComparativeGenomics/Helitron_Scanner.md)
* [DNA Transposon Annotation with Inverted-Repeats Finder](../dataAnalysis/ComparativeGenomics/InvertedRepeatsFinderForDNATransposonAnnotation.md)
* [LTR Retrotransposon Annotation with LTR-Finder](../dataAnalysis/ComparativeGenomics/LTRFinder.md)  
* [Repeat Annotation from Next-gen Sequencing Reads Using RepeatExplorer](../dataAnalysis/ComparativeGenomics/RepeatExplorer.md)
* [De-Novo Repeat Identification and Annotation from Genome Assemblies using RepeatModeler and RepeatMasker](../dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.md)
* [Tandem Duplication Annotation in a Genome Assembly Using Mummer and RedTandem](../dataAnalysis/ComparativeGenomics/Tandem_Duplication_Detection.md)

### [ATAC-Sequencing](../dataAnalysis/ATAC-seq/ATAC-index.md)
* [ATAC-seq](../dataAnalysis/ATAC-seq/ATAC_tutorial.md)

### [Data Visuallization](../Appendix/dataVisualization_index.md)
* [Viewing Files Remotely without Transferring](../Appendix/HPC/viewing-files-in-remote-machine-without-downloading-locally.md)
* [Setting up an R and RStudio Environment](../dataWrangling/R/r-setup.md)
* [Creating Boxplots in R](../dataWrangling/R/generate-boxplots.md)
* [Creating Heatmaps in R](../dataWrangling/R/generate_heatmaps.md)
* [Visulaize Gaps in the Genome assemblies](../dataWrangling/R/visualize-gaps-in-genomes.md)
