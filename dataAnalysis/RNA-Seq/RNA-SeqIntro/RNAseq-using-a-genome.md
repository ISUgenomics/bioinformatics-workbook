---
title: RNA Sequence Analysis
layout: single
author: Siva Chudalayandi
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


# RNA-Seq data Analysis

RNA-seq experiments are performed with an aim to comprehend transcriptomic changes in organisms in response to a certain treatment. They are also designed to understand the cause and/or effect of a mutation by measuring the resulting gene expression changes. Thanks to some robust algorithms specifically designed to map short stretches of nucleotide sequences to a genome while being aware of the process of RNA splicing has led to many advances in RNAseq analysis. The overview of RNA-seq analysis is summarized in Fig1.


### Overview
![**Figure 1.**: Overview of the RNAseq workflow](Assets/RNAseq_1.png)

This document will guide you through basic RNAseq analysis, beginning at quality checking of the RNAseq `reads` through to getting the differential gene expression results. We have downloaded an *Arabidopsis* dataset from NCBI for this purpose. Check the <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA348194" target="_blank">BioProject  ⤴</a> page for more information.


# Experimental design

This experiment compares WT and *atrx-1* mutant to analyze how the loss of function  of ATRX chaperone results in changes in gene expression. The ATRX chaperone is a histone chaperone known to be an important player in the regulation of gene expression. RNA was isolated from three WT replicates and three mutant replicates using Trizol. The transcriptome was enriched/isolated using the plant RiboMinus kit for obtaining total RNA. RNA-seq was carried out in Illumina Hiseq 2500. The sequencing reads were generated as paired-end data, hence we have 2 files per replicate.


| Condition | replicate 1 | replicate 2 | replicate 3 |
| --- | --- | --- | --- |
| WT | SRR4420293_1.fastq.gz <br> SRR4420293_2.fastq.gz | SRR4420294_1.fastq.gz <br> SRR4420294_2.fastq.gz | SRR4420295_1.fastq.gz <br> SRR4420295_2.fastq.gz |
| atrx-1 | SRR4420296_1.fastq.gz <br> SRR4420296_2.fastq.gz| SRR4420297_1.fastq.gz <br> SRR4420297_2.fastq.gz| SRR4420298_1.fastq.gz <br> SRR4420298_2.fastq.gz |

# Raw Sequence Data

Generally, if the data is hosted at your local sequencing center, you could download it through a web interface or using `wget` or `curl` commands. However in this example, we will download data hosted on public repositories. There are several public data repositories that host the data.


# 1. Download the data from public repositories

## A) EMBL-EBI's European Nucletide Archive (ENA)
The best option is to download the fastq files directly from the ENA server (e.g., check <a href="https://www.ebi.ac.uk/ena/data/view/PRJNA348194" target="_blank">EBI  ⤴</a>). We can download them directly using `wget` by supplying the download links to each file; for example, in this case:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/005/SRR4420295/SRR4420295_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/005/SRR4420295/SRR4420295_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/006/SRR4420296/SRR4420296_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/006/SRR4420296/SRR4420296_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298_2.fastq.gz

```

## B) National Center for Biotechnology Information (NCBI)
 The data is hosted as SRA (Sequence Read Archives) files on the public archives in NCBI. We can bulk download these using **aspera** high-speed file transfer.

 It is common to install software as environmental modules on your compute cluster. So we begin by loading the relevant modules. The following code expects that you have `sra-toolkit`, GNU `parallel`, and `aspera` installed on your computing cluster.

 <div style="background: mistyrose; padding: 15px; margin-bottom: 20px;">
 <span style="font-weight:800;">WARNING:</span>
 <br><span style="font-style:italic;"> Sometimes, a user might run into <b>perl</b> issues while using edirect. The installed version of perl should support the HTTPS protocol. </span>
 </div><br>

Run the code snippet in the terminal window:

```
module load <path/to/sra-toolkit>
module load <path/to/edirect>
module load <path/to/parallel>

esearch -db sra -query PRJNA348194 | \
  efetch --format runinfo |\
  cut -d "," -f 1 | \
  awk 'NF>0' | \
  grep -v "Run" > srr_numbers.txt

while read line; do echo "prefetch --max-size 100G --transport ascp --ascp-path \"/path/to/aspera/.../ascp|/path.../etc/asperaweb_id_dsa.openssh\" ${line}"; done < srr_numbers.txt > prefetch.cmds
parallel < prefetch.cmds
```
After downloading the SRA files, we convert them to fastq format. We can use the fast-dump command as follows: (this step is slow, so if possible, run these commands using <a href="https://www.gnu.org/software/parallel/" target="_blank">GNU parallel  ⤴</a>). We assume that all SRA files are in a specific folder.
```
module load parallel
INDIR=/path/to/sra/files/        # Folder containing input SRA files for fastq-dump
parallel fastq-dump --split-files --origfmt --gzip" ::: ${INDIR}/*.sra
```

<div style="background: mistyrose; padding: 15px; margin-bottom: 20px;">
<span style="font-weight:800;">WARNING:</span>
<br><span style="font-style:italic;"> Note that <code>fastq-dump</code> runs very slow. </span>
</div><br>


We also need the genome file and associated GTF/GFF file for *Arabidopsis*. These are downloaded directly from <a href="https://www.ncbi.nlm.nih.gov/genome?term=NC_001284&cmd=DetailsSearch" target="_blank">NCBI  ⤴</a>, or <a href="http://plants.ensembl.org/info/website/ftp/index.html" target="_blank">plants Ensembl website  ⤴</a>, or the <a href="https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Gmax" target="_blank">Phytozome website  ⤴</a>. (*^The Phytozome needs logging in and selecting the files.*)

For this tutorial, we downloaded the following files from <a href="https://www.ncbi.nlm.nih.gov/genome/?term=Arabidopsis+thaliana" target="_blank">NCBI  ⤴</a>:
```
Genome Fasta File: GCF_000001735.3_TAIR10_genomic.fna
Annotation file: GCF_000001735.3_TAIR10_genomic.gff

```

We can use `wget` to fetch both the **Genome** (*Fasta*) file and the **Annotation** (*GFF*) file:

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz
```

# 2. Quality Check

We use `fastqc`, a tool that provides a simple way to do quality control checks on raw sequence data coming from high throughput sequencing pipelines (<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">FastQC  ⤴</a>).

*Watch the video below to check the quality of high throughput sequence using FastQC* [source: <a href="https://" target="_blank">BabrahamBioinf  ⤴</a> ]
<iframe width="560" height="315" src="https://www.youtube.com/embed/bz93ReOv87Y" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

It uses various metrics to indicate how your data is. A high-quality Illumina RNAseq file should look something like <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html" target="_blank">this  ⤴</a>. Since there are 6 sets of files (12 files total), we need to run `fastqc` on each one of them. It is convenient to run it in `parallel`.

```
module load fastqc
module load parallel

mkdir fq_out_directory
OUTDIR=fq_out_directory

parallel "fastqc {} -o ${OUTDIR}" ::: *.fastq.gz
```
Because we have in total 6 quality outputs, we will have 6 HTML files and 6 zip files. We can use <a href="http://multiqc.info/" target="_blank">MultiQC  ⤴</a> to aggregate the outputs and get a single HTML file detailing the data quality of all the libraries.

```
cd fq_out_directory
module load python/3           # may need to search for multiqc module
module load py-multiqc         # if not found, try 'module avail multiqc'

multiqc .
[INFO   ]         multiqc : This is MultiQC v0.8
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '.'
[INFO   ]          fastqc : Found 6 reports
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete


```
This will give you a combined HTML file and a folder named **multiqc_data** with containing three files describing the various statistics:
```
ls  
multiqc_data (Folder)
multiqc_report.html

cd multiqc_data
ls

multiqc_fastqc.txt  
multiqc_general_stats.txt  
multiqc_sources.txt
..........
```


You can peruse the complete report or download the plots and view them, for example:

<img src="https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/Assets/fastqc_adapter_content_plot.png" style="width:300px" /> <img src="https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/Assets/fastqc_per_base_n_content_plot.png" style="width:300px" /> <img src="https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/Assets/fastqc_per_base_sequence_quality_plot.png" style="width:300px" />

<!--
![adapter_content](Assets/fastqc_adapter_content_plot.png)
![per_base_n_content](Assets/fastqc_per_base_n_content_plot.png)
![per_base_sequence_quality](Assets/fastqc_per_base_sequence_quality_plot.png)
 -->

If satisfied with the results, proceed with the mapping. If not, then perform quality trimming.  For example, see <a href="http://hannonlab.cshl.edu/fastx_toolkit/" target="_blank">here  ⤴</a>. If the quality is very bad, excluding that sample from the analysis might make more sense.

<!--
USER FEEDBACK:
Please provide some more information about the interpretation of the results or further reading links to learn more about.
What range of the score is acceptable or good? Is it sufficient just to be in the green part of the plot? How to proceed when some positions are much worse than the others?
What is a practical meaning of those 3 parameters (adapter content, per base n content, mean quality scores)? What they affect?
-->

# 3. Mapping reads to the genome #

There are several **mapping** programs available for **aligning RNAseq reads to the genome**. Generic aligners such as `BWA`, `bowtie2`, `BBMap`, etc., are not suitable for mapping RNAseq reads because they are not splice-aware. **RNAseq reads are mRNA reads that only contain exonic regions**, hence mapping them back to the genome requires splitting the individual reads that span an intron. It can be done only by splice-aware mappers. In this tutorial, we will use <a href="https://ccb.jhu.edu/software/hisat2/index.shtml" target="_blank">HISAT2  ⤴</a>, a successor of `Tophat2`.

<!--
USER FEEDBACK:
Why do we want to align reads to the genome? E.g., to identify their positions in the reference genome? or to detect which of reference genes are present in the studied data? or what is the % of coverage between query and reference?
-->


### HiSat2 for mapping

#### Hisat2 indexing

For `HiSat2` mapping, you first need to index the genome and then use the read pairs to map the indexed genome (one set at a time). For indexing the genome, we use the `hisat2-build` command as follows in a SLURM script:

^ *Learn more about <a href="https://bioinformaticsworkbook.org/Appendix/HPC/SLURM/slurm-cheatsheat.html#gsc.tab=0" target="_blank">slurm  ⤴</a> from the hands-on tutorial.*

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=HI_build
#SBATCH --output=HI_build.%j.out
#SBATCH --error=HI_build.%j.err
#SBATCH --mail-user=<user_email_address>
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -o xtrace
cd ${SLURM_SUBMIT_DIR}
scontrol show job ${SLURM_JOB_ID}
ulimit -s unlimited

module load hisat2

GENOME_FNA="/path/to/refrence/GCF_000001735.3_TAIR10_genomic.fna"
GENOME=${GENOME_FNA%.*}                # Drops the .fna extension
hisat2-build ${GENOME_FNA} ${GENOME}

```
Let's go to your working directory, where you downloaded the genome file. Create an empty file using the `touch indexing_hisat2.sh` command, open it using your favorite text editor in the terminal (e.g., `nano` or `vim`) and copy-paste the script from above. Find the `GENOME_FNA=` variable and update the path and filename of the genomic file. Save changes and submit the script into the SLURM queue:

```
sbatch indexing_hisat2.sh
```

<div style="background: mistyrose; padding: 15px; margin-bottom: 20px;">
<span style="font-weight:800;">WARNING:</span>
<br><span style="font-style:italic;"> Note to update the path and filename of the reference genome in the <b>GENOME_FNA</b> variable. If the downloaded file is zipped, decompress it before running the script: <br>
<code>gzip -d {filename}.fna.gz</code>
</span>
</div><br>

<div style="background: #cff4fc; padding: 15px;">
<span style="font-weight:800;">PRO TIP:</span>
<br><span style="font-style:italic;"> Adjust the <b>#SBATCH --time</b> variable, provided in format hours:minutes:seconds. The provided value allocates the time at the computational node. If your setup is too short, the task will be interrupted prefinal. Otherwise, if it is too long, you will have a long wait in the queue. <br>
Setting up <b>#SBATCH --time=1:00:00</b> is sufficient for indexing the genome in this example.</span>
</div><br>

Once complete, you should see several files with the `.ht2` extension. These are the index files.

```
 GCF_000001735.3_TAIR10_genomic.1.ht2
 GCF_000001735.3_TAIR10_genomic.2.ht2
 GCF_000001735.3_TAIR10_genomic.3.ht2
 GCF_000001735.3_TAIR10_genomic.4.ht2
 GCF_000001735.3_TAIR10_genomic.5.ht2
 GCF_000001735.3_TAIR10_genomic.6.ht2
 GCF_000001735.3_TAIR10_genomic.7.ht2
 GCF_000001735.3_TAIR10_genomic.8.ht2
```

<div style="background: #cff4fc; padding: 15px;">
<span style="font-weight:800;">PRO TIP:</span>
<br><span style="font-style:italic;"> If outputs were not generated, check the <b>HI_build.{jobID}.err</b> (and HI_build.{jobID}.out) files to detect the reason for the failure. </span>
</div><br>

 At the mapping step, we simply refer to the index using `GCF_000001735.3_TAIR10_genomic` as described in the next step.


#### Hisat2 Mapping

For mapping, for each set of reads (forward and reverse or R1 and R2), we first set up a `run_hisat2.sh` script as follows:

```
#!/bin/bash

set -o xtrace

# set the reference index:
GENOME="/path/to/refrence/GCF_000001735.3_TAIR10_genomic"

# make an output directory to store the output aligned files
OUTDIR="/path/to/out_dir"
[[ -d ${OUTDIR} ]] || mkdir -p ${OUTDIR}        # If output directory doesn't exist, craete it

p=8 # use 8 threads
R1_FQ="$1" # first argument
R2_FQ="$2" # second argument

# purge and load relevant modules.
module purge

module load hisat2
module load samtools

OUTFILE=$(basename ${R1_FQ} | cut -f 1 -d "_");

hisat2 \
  -p ${p} \
  -x ${GENOME} \
  -1 ${R1_FQ} \
  -2 ${R2_FQ} \
  -S ${OUTDIR}/${OUTFILE}.sam &> ${OUTFILE}.log
samtools view --threads ${p} -bS -o ${OUTDIR}/${OUTFILE}.bam ${OUTDIR}/${OUTFILE}.sam

rm ${OUTDIR}/${OUTFILE}.sam

```

For setup to run with each set of files, we can set a SLURM script (`loop_hisat2.sh`) that loops over each fastq file. Note that this script calls the run_hisat2.sh script for each pair of fastq files supplied as its argument.

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=Hisat2
#SBATCH --output=Hisat2.%j.out
#SBATCH --error=Hisat2.%j.err
#SBATCH --mail-user=<user_email_address>
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -o xtrace
cd ${SLURM_SUBMIT_DIR}
ulimit -s unlimited
scontrol show job ${SLURM_JOB_ID}

for fq1 in *1.*gz; do
  fq2=$(echo ${fq1} | sed 's/1/2/g');
  bash run_hisat2.sh ${fq1} ${fq2};
done &> hisat2_1.log

```

Now submit this job as follows:

```
sbatch loop_hisat2.sh
```

That should create the following files as output:

```
SRR4420298.bam
SRR4420297.bam
SRR4420296.bam
SRR4420295.bam
SRR4420294.bam
SRR4420293.bam
```
### STAR for mapping

#### STAR Index


# 4. Abundance estimation

For quantifying transcript abundance from RNA-seq data, there are many programs available. Two most popular tools include, `featureCounts` and `HTSeq`. We will need a file with aligned sequencing reads (SAM/BAM files generated in previous step) and a list of genomic features (from the GFF file). `featureCounts` is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It also outputs stat info for the overall summarization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons. We can run `featureCounts` on all SAM/BAM files at the same time or individually.

<!--
Why do we want to estimate transcript abundance?
Why and for what is that needed?
-->


#### featureCounts ####

You will need <a href="http://subread.sourceforge.net/" target="_blank">subread  ⤴</a> and `parallel` modules loaded.

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --job-name=Hisat2
#SBATCH --output=Hisat2.%j.out
#SBATCH --error=Hisat2.%j.err
#SBATCH --mail-user=<user_email_address>
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -o xtrace

cd ${SLURM_SUBMIT_DIR}
ulimit -s unlimited
ANNOT_GFF="/PATH/TO/GFF_FILE"          # file should be decompressed, e.g., gzip -d {filename}.gff.gz
INDIR="/Path/TO/BAMFILES/DIR"
OUTDIR=counts
[[ -d ${OUTDIR} ]] || mkdir -p ${OUTDIR}

module purge

module load subread
module load parallel

parallel -j 4 "featureCounts -T 4 -s 2 -p -t gene -g ID -a ${ANNOT_GFF} -o ${OUTDIR}/{/.}.gene.txt {}" ::: ${INDIR}/*.bam
scontrol show job ${SLURM_JOB_ID}

```

Let's go to your working directory, where you downloaded and decompressed the annotation (GFF) file. Create an empty file using the `touch run_feature_counts.sh` command, open it using your favorite text editor in the terminal (e.g., `nano` or `vim`) and copy-paste the script from above. Find the `ANNOT_GFF=` variable and update the path and filename of the annotation file. Also, update the `INDIR=` variable using path to the directory with BAM files from the previous step. Save changes and submit the script into the SLURM queue:

```
sbatch run_feature_counts.sh
```

This creates the following set of files in the specified output folder:

Count Files:

```
SRR4420298.gene.txt
SRR4420293.gene.txt
SRR4420297.gene.txt
SRR4420296.gene.txt
SRR4420295.gene.txt
SRR4420294.gene.txt
```

Each file has a commented line starting with a # which gives the command used to create the count table and the relevant seven columns as follows, for example:

`head SRR4420298.gene.txt`

```
# Program:featureCounts v1.5.2; Command:"featureCounts" "-T" "4" "-s" "2" "-p" "-t" "gene" "-g" "ID" "-a" "/path/to/GCF_000001735.3_TAIR10_genomic.gff" "-o" "/path/to/counts/SRR4420298.gene.txt" "SRR4420298.bam"

Geneid  Chr     Start   End     Strand  Length  SRR4420298.bam
gene0   NC_003070.9     3631    5899    +       2269    13
gene1   NC_003070.9     6788    9130    -       2343    17
gene2   NC_003070.9     11101   11372   +       272     0
gene3   NC_003070.9     11649   13714   -       2066    13
gene4   NC_003070.9     23121   31227   +       8107    32
gene5   NC_003070.9     23312   24099   -       788     0
gene6   NC_003070.9     28500   28706   +       207     0
gene7   NC_003070.9     31170   33171   -       2002    45

```
Additionally, Summary Files are produced. These give the summary of reads that were either ambiguous, multi mapped, mapped to no features, or unmapped among other statistics. We can refer to these to further tweak our analyses etc.

```
SRR4420298.gene.txt.summary
SRR4420293.gene.txt.summary
SRR4420295.gene.txt.summary
SRR4420296.gene.txt.summary
SRR4420294.gene.txt.summary
SRR4420297.gene.txt.summary

```
Using the following command (a combination of paste and awk), we can produce a single count table. This count table could be loaded into R for differential expression analysis.

```
paste <(awk 'BEGIN {OFS="\t"} {print $1,$7}' SRR4420293.gene.txt) \
  <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420294.gene.txt) \
  <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420295.gene.txt) \
  <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420296.gene.txt) \
  <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420297.gene.txt) \
  <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420298.gene.txt) | \
  grep -v '^\#' > At_count.txt
```

You could also edit out the names of the samples to something succinct, for example, S293 instead of SRR4420293.bam.

head At_count.txt

```
Geneid  S293    S294    S295    S296    S297    S298
gene0   11      1       10      28      11      13
gene1   37      3       26      88      22      17
gene2   0       0       0       0       0       0
gene3   6       2       12      40      13      13
gene4   35      6       22      170     53      32
gene5   0       0       0       0       0       0
gene6   0       0       0       0       0       0
gene7   49      15      67      258     83      45
gene8   0       0       0       0       0       0
```

Now we are ready for performing <a href="https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/Differential-Expression-Analysis.html#gsc.tab=0" target="_blank">DGE analysis ⤴</a>!

---
[Table of contents](RNAseq-intro.md)
