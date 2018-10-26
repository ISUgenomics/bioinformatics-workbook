

# Bioinformatics Workbook


## Table of Contents

### Introduction
  * [Introduction to Bioinformatics](introduction/introduction.md)
  * Introduction to bioinformatics Terminalogy
    * [Reads, contigs and scaffolds](/introduction/dataTerminology.md)
  * [Introduction to File Formats](introduction/fileFormats.md)
      * [What is a Quality Score?](introduction/fastqquality-score-encoding.md)

### Experimental Design
  * [Introduction to Experimental Design](experimentalDesign/eD_introduction.md)
  * [Generic examples of Experimental Design](/experimentalDesign/eD_genericExamples.md)


### Data Acquisition
  * [Introduction to Data Acquisition](dataAcquisition/dAc_introduction.md)
  * Transferring data
    * [Downloading with wget](dataAcquisition/FileTransfer/downloading-files-via-wget.md)

### Data Wrangling
  * FASTA(Q) manipulations
    * [Determine Sequence Lengths](dataWrangling/fastaq-manipulations/calculate-sequence-lengths-in-a-fasta-file.md)
    * [Converting FASTQ to FASTA](dataWrangling/fastaq-manipulations/converting-fastq-format-to-fasta.md)
    * [Extract sequences based on sequence ID](dataWrangling/fastaq-manipulations/retrieve-fasta-sequences-using-sequence-ids.md)

* Microsoft Excel Tips and Tricks
    * [Create workbook from multiple text files](dataWrangling/microsoftExcel/export-multiple-worksheets-as-separate-text-files-in-excel.md)
    * [Create text files from workbook with multiple worksheets](dataWrangling/microsoftExcel/export-multiple-worksheets-as-separate-text-files-in-excel.md)
    * [Create Index for all worksheets](dataWrangling/microsoftExcel/generate-index-sheet-linking-all-spreadsheets-in-excel.md)
 * Data Management
    * [Deposition of Data to NCBI SRA](dataWrangling/NCBI_Data_Submission.md)
### Data analysis
  * Alignment
    * [Blast Command Line Basic Example](dataAnalysis/blast/blastExample.md)
  * Comparative Genomics
    * [Gene orthology, synteny, and visualization with opscan, iadhore, and circos](dataAnalysis/ComparativeGenomics/Gene_Orthology_And_Synteny.md)
    * [Gene orthology, synteny, and visualization with orthofinder, iadhore, and circos](dataAnalysis/ComparativeGenomics/OrthofinderSynteny.md)
    * [Gene overlap significance testing with R Gene_overlap package](dataAnalysis/ComparativeGenomics/Gene_Category_Overlap_Fishers_exact_testing.md)
  * MetaGenomics
    * [Intro to Metagenomics](dataAnalysis/Metagenomics/MetagenomicsP1.md)
    * [Qiime2](dataAnalysis/Metagenomics/Qiime2.md)
    * [DADA2](dataAnalysis/Metagenomics/Dada2.md)
  * Genomic Repeat Identification
    * [Helitron identification in a genome sequence](dataAnalysis/ComparativeGenomics/Helitron_Scanner.md)    
    * [DNA transposon annotation with Inverted-Repeats Finder](dataAnalysis/ComparativeGenomics/InvertedRepeatsFinderForDNATransposonAnnotation.md)
    * [LTR retrotransposon annoatation with LTR-Finder](dataAnalysis/ComparativeGenomics/LTRFinder.md)
    * [Repeat annotation from next-gen sequencing reads using RepeatExplorer](dataAnalysis/ComparativeGenomics/RepeatExplorer.md)
    * [De-novo repeat identification and annotation from genome assemblies using RepeatModeler and RepeatMasker ](dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.md)
    * [Tandem duplication annotation in a genome assembly using Mummer and RedTandem](dataAnalysis/ComparativeGenomics/Tandem_Duplication_Detection.md)
  * Genome Assembly
    * [GenomeScope to estimate genome size](dataAnalysis/GenomeAssembly/genomescope.md)
    * [Canu for long read assembly](dataAnalysis/GenomeAssembly/LongRead/Canu.md)
      * [Canu on XSEDE bridges machine](dataAnalysis/GenomeAssembly/LongRead/Canu_bridges.md)
    * [Genome scaffolding using synteny with pyscaf](dataAnalysis/GenomeAssembly/Pyscaf_Synteny_Scaffolding.md): E. coli and Arabidopsis genome scaffolding
  * Gene Annotation
    * [Intro to Maker gene prediction](dataAnalysis/GenomeAnnotation/Intro_To_Maker.md)
    * [Intro to Braker2 gene prediction](dataAnalysis/GenomeAnnotation/Intro_to_Braker2.md)
    * [Motif Identification and finding with MEME and FIMO](dataAnalysis/GenomeAnnotation/MEME_Motif_Finding_In_Genomes.md)
  * RNA-Seq analysis
    * [RNA-Seq example with a genome assembly](dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-using-a-genome.md): Arabidopsis with a genome
    * [RNA-Seq example without a genome assembly](dataAnalysis/RNA-Seq/RNA-SeqIntro/RNAseq-without-a-genome.md): Arabidopsis without a genome
    * [Differential Expression Analysis](dataAnalysis/RNA-Seq/RNA-SeqIntro/Differential-Expression-Analysis.md): DESeq2
  * Specialized sequencing analysis
    * [ATAC-seq](https://github.com/ISUgenomics/bioinformatics-workbook/blob/master/dataAnalysis/ATAC-seq/ATAC_tutorial.md) in Arabidopsis

### Data visualization
* [Viewing files remotely without transferring](Appendix/HPC/viewing-files-in-remote-machine-without-downloading-locally.md)
* [Creating Boxplots in R](dataWrangling/R/generate-boxplots.md)
* [Creating Heatmaps in R](dataWrangling/R/generate_heatmaps.md)

### Appendix
  * Useful programs
    * [Introduction to GitHub](Appendix/github/introgithub.md)
    * [Introduction to Slack](Appendix/slack.md)
    * [How to use Markdown](Appendix/Markdown.md)
  * Scripting and command line
    * [Bioawk Basics](Appendix/bioawk-basics.md)
    * [Introduction to Version Control and Github](Appendix/github/githubBasics.md)
  * Installation
    * [How do I install a program?](Appendix/HPC/guide-for-installing-various-types-of-programs-in-linux.md)
  * HPC
    * Containers
        * [Introduction to Containers](Appendix/HPC/Containers/Intro_Singularity.md)
        * [Creating Containers using Singularity](Appendix/HPC/Containers/creatingContainers.md)
        * [Modifying exiting containers](Appendix/HPC/Containers/modifyingExistingContainers.md)
    * SLURM
        * [cheatsheat SLURM](/Appendix/HPC/SLURM/slurm-cheatsheat.md)
        * [SLURM job submission dependencies](/Appendix/HPC/SLURM/submitting-dependency-jobs-using-slurm.md)
    * PBS-Torque
        * [cheatsheat Torque](Appendix/HPC/pbstorque/submitting-dependency-jobs-using-pbs-torque.md)
    * Miscellaneous
        * [Password-less SSH](Appendix/HPC/password-less-ssh-login.md)
        * [SSH shortcuts](Appendix/HPC/ssh-shortcuts.md)
    * Introduction to UNIX
        * [Unix basics 1](Appendix/unix-basics-1.md)  
        * [Unix basics 2](Appendix/unix-basics-2.md)
        * [Unix basics 3 (grep)](Appendix/unix-basics-3.md)
        * [Unix basics 4 (sed)](Appendix/unix-basics-4.md)
    * XSEDE
        * [What is XSEDE, and how do I get an allocation](Appendix/HPC/xsede/xsede.md)


  * [Other website links of interest](Appendix/OtherLinks.md)
