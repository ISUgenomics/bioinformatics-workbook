# RNA-Seq data Analysis

RNA-seq experiments are performed with an aim to comprehend transcriptomic changes in organisms in response to a certain treatment. They are also designed to understand the cause and/or effect of a mutation by measuring the resulting gene expression changes. Thanks to some robust algorithms specifically designed to map short stretches of nucleotide sequences to a genome while being aware of the process of RNA splicing has led to many advances in RNAseq analysis. The overview of RNA-seq analysis is summarized in Fig1.


### Overview ###
![**Figure 1.**: Overview of the RNAseq workflow](/assets/RNAseq_1.png)

This document will guide you through basic RNAseq analysis, beginning at quality checking of the RNAseq `reads` through to getting the differential gene expression results. We have downloaded an *Arabidopsis* dataset from NCBI for this purpose. Check the [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA348194) page for more information.


# Experimental design #

This experiment compares WT and *atrx-1* mutant to analyze how loss of function  of ATRX chaperone results in changes in gene expression. The ATRX chaperone is a histone chaperone known to be an important player in regulation of gene expression. RNA was isolated from three WT replicates and three mutant replicates using Trizol. Transcriptome was enriched/isolated using the plant RiboMinus kit for obtaining total RNA. RNA-seq was carried out in Illumina Hiseq 2500. The sequencing reads were generated as paired-end data, hence we have 2 files per replicate.


| Condition | replicate 1 | replicate 2 | replicate 3 |
| --- | --- | --- | --- |
| WT | SRR4420293_1.fastq.gz <br> SRR4420293_2.fastq.gz | SRR4420294_1.fastq.gz <br> SRR4420294_2.fastq.gz | SRR4420295_1.fastq.gz <br> SRR4420295_2.fastq.gz |
| atrx-1 | SRR4420296_1.fastq.gz <br> SRR4420296_2.fastq.gz| SRR4420297_1.fastq.gz <br> SRR4420297_2.fastq.gz| SRR4420298_1.fastq.gz <br> SRR4420298_2.fastq.gz |

# 1. Download the data from public databases #

## NCBI
Generally if the data is hosted at your local sequencing center you could download through a web interface or using `wget` or `curl` commands. In this case, however, we first download the SRA files from the public archives in NCBI in bulk using aspera high speed file transfer. The following code expects that you have sra-toolkit, GNU parallel and aspera installed on your computing cluster. On Ceres, in order to use an installed software, we load the relevant module.

```
module load <path/to/sra-toolkit>
module load <path/to/edirect>
module load <path/to/parallel>
esearch -db sra -query PRJNA276699 | efetch --format runinfo |cut -d "," -f 1 | awk 'NF>0' | grep -v "Run" > srr_numbers.txt
while read line; do echo "prefetch --max-size 100G --transport ascp --ascp-path \"/path/to/aspera/<version>/etc/asperaweb_id_dsa.openssh\" $line"; done<srr_numbers.txt > prefetch.cmds
parallel <prefetch.cmds
```
After downloading the SRA files, we convert it to fastq format. We can use the fast-dump command as follows: (this step is slow and if possible run these commands using [gnu parallel](https://www.gnu.org/software/parallel/)). We assume that all SRA files are in a specific folder.
```
module load parallel
parallel "fastq-dump --split-files --origfmt --gzip" ::: /path/to/SRA/*.sra
```
## EBI
EBI directly hosts the fastq files are available on their server (e.g. check [EBI](https://www.ebi.ac.uk/ena/data/view/PRJNA348194)) we can download them directly using `wget` by supplying the links to each file.

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_1.fastq.gz
```


We also need the genome file and associated GTF/GFF file for for *Arabidopsis*. These are downloaded directly from [NCBI](https://www.ncbi.nlm.nih.gov/genome?term=NC_001284&cmd=DetailsSearch), or [plants Ensembl website](http://plants.ensembl.org/info/website/ftp/index.html) or the [Phytozome website](https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Gmax "Glycine max") (phytozome needs logging in and selecting the files) .

For this tutorial, we downloaded the following files from [NCBI](https://www.ncbi.nlm.nih.gov/genome?term=NC_001284&cmd=DetailsSearch).
```
Genome Fasta File: GCF_000001735.3_TAIR10_genomic.fna
Annotation file: GCF_000001735.3_TAIR10_genomic.gff

```

# 2. Quality Check #

For this we will use `fastqc`, which is a tool that provides a simple way to do quality control checks on raw sequence data coming from high throughput sequencing pipelines ([link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). It provides various metrics to give a indication of how your data is. A high quality illumina RNAseq file should look something like [this](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html). Since there are 6 set of files (12 files total), and we need to run `fastqc` on each one of them. It is convenient to run it in `parallel`.

```
module load fastqc
module load parallel
parallel "fastqc {} -o <fq_out_directory>" ::: *.fastq.gz
```
Because we have a total of 6 quality outputs, we will have 6 html files and 6 zip files. We can use [`multiqc`](http://multiqc.info/) to aggregate the outputs and get a single html file to scan the quality of all the libraries.

```
cd fq_out_directory
module load python_3
multiqc .
[INFO   ]         multiqc : This is MultiQC v0.8
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching '.'
[INFO   ]          fastqc : Found 6 reports
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete


```
This will give you a combined html file and a folder named multiqc_data with containing three files describing the various statistics:
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


You can peruse the complete report or download the plots and view them for example: ![adapter_content](Assets/fastqc_adapter_content_plot.png)

![per_base_n_content](Assets/fastqc_per_base_n_content_plot.png)
![per_base_sequence_quality](Assets/fastqc_per_base_sequence_quality_plot.png)

If satistied with the results, proceed with the mapping. If not, then perform quality trimming. E.g. see [here](http://hannonlab.cshl.edu/fastx_toolkit/). If the quality is very bad it might make more sense to exclude that sample from the analysis.

# 3. Mapping reads to the genome #

There are several mapping programs available for aligning RNAseq reads to the genome. Generic aligners such as BWA, BowTie2, BBMap etc., are not suitable for mapping RNAseq reads because they are not splice aware. RNAseq reads are mRNA reads that only contain exonic regions, hence mapping them back to the genome requires splitting the individual reads that span an intron. This is only done by splice aware mappers. Here for this tutorial, we will use [`HISAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml). HISAT2 is a successor of Tophat2.


### HiSat2 for mapping ###

#### Hisat2 Index ####

For HiSat2 mapping, you need to first index the genome and then use the read pairs to map the indexed genome (one set at a time). For indexing the genome, we use the `hisat2-build` command as follows in a slurm script:

```
#!/bin/bash
#set -o xtrace
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 24:00:00
#SBATCH -J HI_build
#SBATCH -o HI_build.o%j
#SBATCH -e HI_build.e%j
#SBATCH --mail-user=<user_email_address>
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
scontrol show job $SLURM_JOB_ID
ulimit -s unlimited
module use /software/modulefiles/
module load hisat2
GENOME="/path/to/refrence/GCF_000001735.3_TAIR10_genomic.fna"
hisat2-build $GENOME ${GENOME%.*}

```

Once complete, you should see a number of files with `.ht2` extension.  These are the index files.
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
 At the mapping step we simply refer to the index using `GCF_000001735.3_TAIR10_genomic` as described in the next step.


#### Hisat2 Mapping ####

For mapping, each set of reads (forward and reverse or R1 and R2), we first set up a run_hisat2.sh script.
```
#!/bin/bash
set -o xtrace
# set the rerefernce index:
GENOME="/path/to/refrence/GCF_000001735.3_TAIR10_genomic"
# make an output directory to store the output aligned files
mkdir -p /path/to/out_dir
# set that as the output directory
ODIR="/path/to/out_dir"


p=8 # use 8 threads
R1_FQ="$1" # first argument
R2_FQ="$2" # second argument

# purge and load relevant modules.
module purge
module use /software/modulefiles/
module load hisat2


OUTPUT=$(basename ${R1_FQ} |cut -f 1 -d "_");

hisat2 \
  -p ${p} \
  -x ${GENOME} \
  -1 ${R1_FQ} \
  -2 ${R2_FQ} \
  -S $ODIR\/${OUTPUT}.sam &> ${OUTPUT}.log
samtools view --threads 8 -bS -o $ODIR\/${OUTPUT}.bam $ODIR\/${OUTPUT}.sam

rm $ODIR\/${OUTPUT}.sam

```

For setting it up to run with each set of file, we can set a SLURM script that loops over each fastq file. Note that this script calls the run_hisat2.sh script for each pair of fastq file supplied as its argument.
```
#!/bin/bash
set -o xtrace
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 24:00:00
#SBATCH -J Hisat2
#SBATCH -o Hisat2.o%j
#SBATCH -e Hisat2.e%j
#SBATCH --mail-user=csiva@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID

for fq1 in *1.*gz;
do
fq2=$(echo $fq1 | sed 's/1/2/g');
/path/to/run_hisat.sh ${fq1} ${fq2};
done >& hisat2_1.log

```


This should create, following files as output:
```
SRR4420298.bam
SRR4420297.bam
SRR4420296.bam
SRR4420295.bam
SRR4420294.bam
SRR4420293.bam
```
# 2. Abundance estimation #

For quantifying transcript abundance from RNA-seq data, there are many programs available. Two most popular tools include, `featureCounts` and `HTSeq`. We will need a file with aligned sequencing reads (SAM/BAM files generated in previous step) and a list of genomic features (from the GFF file). `featureCounts` is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It also outputs stat info for the overall summarization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons. We can run featureCounts on all SAM/BAM files at the same time or individually.

#### featureCounts ####
You will need [`subread`](http://subread.sourceforge.net/) and `parallel` modules loaded.
```
ANNOT="/path/to/GCF_000001735.3_TAIR10_genomic.gff"
mkdir -p path/to/counts
ODIR="path/to/counts"


module purge
module use /software/modulefiles/
module load subread
module load parallel

parallel -j 4 "featureCounts -T 4 -s 2 -p -t gene -g ID -a $ANNOT -o $ODIR/{.}.gene.txt {}" ::: *.bam

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
Each file has a commented line staring with a # which gives the command used to create the count table and the relevant seven columns as follows, for example:

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
Additionally Summary Files are produced. These give the summary of reads that were either ambiguous, multimapped, mapped to no features or unmapped among other statistics. We can refer to these to further tweak our analyses etc.

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
paste <(awk 'BEGIN {OFS="\t"} {print $1,$7}' SRR4420293.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420294.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420295.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420296.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420297.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420298.gene.txt) | grep -v '^\#' > At_count.txt
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

Now we are ready for performing DGE analysis!

# 3. Differential Gene Expression analysis #

There are many programs that you can use to perform differential expression Some of the popular ones for RNA-seq are [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html),[`edgeR`](http://bioconductor.org/packages/release/bioc/html/edgeR.html), or [`QuasiSeq`](https://cran.r-project.org/web/packages/QuasiSeq/index.html). Here we will demonstrate differential expression using DESeq2.


### Differential Expression with DESeq2 ###
These steps should be done either on RStudio or in R terminal:

```
## RNA-seq analysis with DESeq2
## Largely based on Stephen Turner, @genetics_blog
## https://gist.github.com/stephenturner/f60c1934405c127f09a6

source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")


setwd("~/Atrx/")
dat<-read.table("At_count.txt",header = T,quote = "",row.names = 1)

# Convert to matrix
dat <- as.matrix(dat)
head(dat)

# Assign condition (first three are WT, next three are mutants)

condition <- factor(c(rep("WT",3),rep("Mut",3)))
condition=relevel(condition,ref = "WT")


# Create a coldata frame: its rows correspond to columns of dat (i.e., matrix representing the countData)
coldata <- data.frame(row.names=colnames(dat), condition)

head(coldata)

#            condition
# S293        WT
# S294        WT
# S295        WT
# S296       Mut
# S297       Mut
# S298       Mut


##### DESEq pipeline, first the design and the next step, normalizing to model fitting
dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata,design=~ condition)


dds <- DESeq(dds)

# Plot Dispersions:
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()
```
![qc-dispersions.png](Assets/qc-dispersions.png)
```
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Principal Components Analysis
plotpca(rld)
```
![PCA.png](Assets/PCA.png)
```

# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
```
![Heatmap-Samples](Assets/qc-heatmap-samples.png)
```
# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

```
We observe 204 differentially expressed genes with adjusted p value <= 0.05
|FALSE|TRUE|
|---|---|
|6712|204|

```
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv",quote = FALSE,row.names = F)

```
To get the first few rows do  
`head -4 diffexpr-results.csv`

```
Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,S293,S294,S295,S296,S297,S298
gene32459,4529.272154082,-2.06484914507321,0.288731491703117,-7.15145110390785,8.58653066559819e-13,5.93844460832771e-09,11073.6705141816,5726.45128024,5128.92211910988,1482.46794888894,1473.13299530101,2290.98806677055
gene38073,36575.3870282423,-1.38271999532283,0.199996648356586,-6.91371583816493,4.7212044436934e-12,1.63259249662918e-08,60993.4651420719,43969.9190693109,53655.982748808,21424.5919655829,22323.7590351868,17084.604208493
gene1446,62.2619789946306,-2.86730516347734,0.43036495871904,-6.66249680738816,2.69214169069523e-11,6.20628397761607e-08,186.819432107996,65.9054495168287,73.2220377206551,12.4735356867937,10.8379903446028,24.3134285909079

```



```

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")



## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## This is Stephen Turner's code:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), points(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()
```
![MA plot](Assets/diffexpr-maplot.png)
```
## Plots to Examine Results:

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
```
![Volcano Plot](Assets/diffexpr-volcanoplot.png)
