---
title: "Table of Contents"

permalink: /list.html
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

### Command Line Basics and Useful Programs
* <span style="color:lightblue">Introduction to Unix</span>
  * [Unix Basics 1](../Appendix/unix-basics-1.md)  
  * [Unix Basics 2](../Appendix/unix-basics-2.md)
  * [Unix Basics 3 (grep)](../Appendix/unix-basics-3.md)
  * [Unix Basics 4 (sed)](../Appendix/unix-basics-4.md)
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

* <span style="color:lightblue">GitHub</span>
  * [Introduction to GitHub](../Appendix/github/introgithub.md)
  * [Some Helpful Commands For Your New Repository](../Appendix/github/github2.md)
  * [Best Practices on Github](../Appendix/github/githubBasics.md)
* <span style="color:lightblue">Bioawk</span>
  * [Bioawk Basics](../Appendix/bioawk-basics.md)
* [Slack](../Appendix/slack.md)
* [Markdown](../Appendix/Markdown.md)
* [Viewing Files In Remote Machine Without Downloading](../Appendix/HPC/viewing-files-in-remote-machine-without-downloading-locally.md)

### Experimental Design
* [Biological system](../experimentalDesign/bio_sys.md)
* [Sequencing Technology](../experimentalDesign/sequencing.md)
* [Costs](../experimentalDesign/costs.md)
* [Generic Examples of Experimental Design](../experimentalDesign/eD_genericExamples.md)
* [List of Biology exceptions and irregularities](../Appendix/biology_tidbits.md)

### Data Acquisition and Wrangling
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
* <span style="color:lightblue">Manuplating Excel data sheets</span>
  * [Create Workbook from Multiple Text Files](../dataWrangling/microsoftExcel/import-multiple-text-files-as-separate-worksheets-in-excel.md)
  * [Export multiple worksheets as separate text files ](../dataWrangling/microsoftExcel/export-multiple-worksheets-as-separate-text-files-in-excel.md)
  * [Create Index for All Worksheets](../dataWrangling/microsoftExcel/generate-index-sheet-linking-all-spreadsheets-in-excel.md)
* <span style="color:lightblue">Data Management</span>
  * [Deposition of Data to NCBI SRA](../dataWrangling/NCBI_Data_Submission.md)

### Bioinformatics terminology
* [Read, Contigs, Scafolds and Choromosome](../introduction/dataTerminology.md)
* [File Formats](../introduction/fileFormats.md)
* [Fasta Quality Score](../introduction/fastqquality-score-encoding.md)

### RNA Sequencing analysis
* [RNA-Seq Example with a Genome Assembly](../dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-using-a-genome.md)
* [RNA-Seq Example without a Genome Assembly](../dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-without-a-genome.md)
* [Different Expression Analysis:DESeq2](../dataAnalysis/RNA-Seq/RNA-SeqIntro/Differential-Expression-Analysis.md)

### Genome Assembly and Annotation
* <span style="color:lightblue">Genome Assembly</span>
  * [GenomeScope to Estimate Genome Size](../dataAnalysis/GenomeAssembly/genomescope.md)
  * [Canu for Long Read Assembly](../dataAnalysis/GenomeAssembly/LongRead/Canu.md)
    * [Canu on XSEDE Bridges Machine](../dataAnalysis/GenomeAssembly/LongRead/Canu_bridges.md)
  * [Mascurca with Pacbio and Illumina](../dataAnalysis/GenomeAssembly/Hybrid/MaSuRCA.md)
  * [Genome Scaffolding Using Synteny eith Pyscaf](../dataAnalysis/GenomeAssembly/Pyscaf_Synteny_Scaffolding.md)
* <span style="color:lightblue">Genome Annotation</span>
  * [Introduction to Maker Gene Prediction](../dataAnalysis/GenomeAnnotation/Intro_To_Maker.md)
  * [Introduction to Braker2 Gene Prediction](../dataAnalysis/GenomeAnnotation/Intro_to_Braker2.md)
  * [Motif Identification and Finding with MEME and FIMO](../dataAnalysis/GenomeAnnotation/MEME_Motif_Finding_In_Genomes.md)  

### Comparative Genomics
* [Gene Orthology, Synteny, and Visualzation with Opscan, Iadhore and Circos](../dataAnalysis/ComparativeGenomics/Gene_Orthology_And_Synteny.md)
  * [Gene Orthology, Synteny, and Visualzation with Orthofinder, Iadhore and Circos](../dataAnalysis/ComparativeGenomics/OrthofinderSynteny.md)
  * [Gene Overlap Significance Testing with R Gene_overlap Package](../dataAnalysis/ComparativeGenomics/Gene_Category_overlap_Fisher_exact_testing.md)  
  * [Phylostratiophraphy:Determining the LCA of all Genes in a Genome](../dataAnalysis/ComparativeGenomics/phylostratr.md)

### Genome Repeat Identification
* [Helitron Identification in a Genome Sequence](../dataAnalysis/ComparativeGenomics/Helitron_Scanner.md)
* [DNA Transposon Annotation with Inverted-Repeats Finder](../dataAnalysis/ComparativeGenomics/InvertedRepeatsFinderForDNATransposonAnnotation.md)
* [LTR Retrotransposon Annotation with LTR-Finder](../dataAnalysis/ComparativeGenomics/LTRFinder.md)  
* [Repeat Annotation from Next-gen Sequencing Reads Using RepeatExplorer](../dataAnalysis/ComparativeGenomics/RepeatExplorer.md)
* [De-Novo Repeat Identification and Annotation from Genome Assemblies using RepeatModeler and RepeatMasker](../dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.md)
* [Tandem Duplication Annotation in a Genome Assembly Using Mummer and RedTandem](../dataAnalysis/ComparativeGenomics/Tandem_Duplication_Detection.md)

### Data Visuallization
* [Viewing Files Remotely without Transferring](../Appendix/HPC/viewing-files-in-remote-machine-without-downloading-locally.md)
* [Creating Boxplots in R](../dataWrangling/R/generate-boxplots.md)
* [Creating Heatmaps in R](../dataWrangling/R/generate_heatmaps.md)
