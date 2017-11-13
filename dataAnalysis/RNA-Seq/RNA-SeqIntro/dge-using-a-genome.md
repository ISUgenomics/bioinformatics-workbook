# RNA-Seq data Analysis

RNA-seq experiments are performed with an aim to comprehend transcriptomic changes in organisms in response to a certain treatment. They are also designed to understand the cause and/or effect of a mutation by measuring the resulting gene expression changes. Thanks to some robust algorithms specifically designed to map short stretches of nucleotide sequences to a genome while being aware of the process of RNA splicing has led to many advances in RNAseq analysis. The overview of RNA-seq analysis is summarized in Fig1.


### Overview ###
![**Figure 1.**: Overview of the RNAseq workflow](/assets/RNAseq_1.png)

This document will guide you through basic RNAseq analysis, beginning at quality checking of the RNAseq `reads` through to getting the differential gene expression results. We have downloaded an *Arabidopsis* dataset from NCBI for this purpose. Check the [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA348194) page for more information.


#Experimental design #

This experiment compares WT and atrx-1 mutant to analyze how loss of function  of ATRX chaperone results in changes in gene expression. RNA was isolated from three WT replicates and three mutant replicates using Trizol. Transcriptome was enriched/isolated using the plant RiboMinus kit for obtaining total RNA. RNA-seq was carried out in Illumina Hiseq 2500. The sequencing reads were generated as paired-end data, hence we have 2 files per replicate.


| Condition | replicate 1 | replicate 2 | replicate 3 |
| --- | --- | --- | --- |
| WT | SRR4420293_1.fastq.gz <br> SRR4420293_2.fastq.gz | SRR4420294_1.fastq.gz <br> SRR4420294_2.fastq.gz | SRR4420295_1.fastq.gz <br> SRR4420295_2.fastq.gz |
| atrx-1 | SRR4420296_1.fastq.gz <br> SRR4420296_2.fastq.gz| SRR4420297_1.fastq.gz <br> SRR4420297_2.fastq.gz| SRR4420298_1.fastq.gz <br> SRR4420298_2.fastq.gz |

# 1. Download the data from NCBI #

Generally if the data is hosted at your local sequencing center you could download through a web interface or using `wget` or `curl` commands. In this case, however, we first download the SRA files from the public archives in NCBI in bulk using aspera high speed file transfer. The following code expects that you have sra-toolkit, GNU parallel and aspera installed on your computing cluster. On Ceres, in order to use an installed software, we load the relevant module.

```
module load <path/to/sra-toolkit>
module load <path/to/edirect>
module load <path/to/parallel>
esearch -db sra -query PRJNA276699 | efetch --format runinfo |cut -d "," -f 1 | awk 'NF>0' | grep -v "Run" > srr_numbers.txt
while read line; do echo "prefetch --max-size 100G --transport ascp --ascp-path \"/path/to/aspera/<version>/etc/asperaweb_id_dsa.openssh\" $line"; done<srr_numbers.txt > prefetch.cmds
parallel <prefetch.cmds
```
After downloading the SRA files, we convert it to fastq format. We can use the fast-dump command as follows: (this step is slow and if possible run these commands using gnu parallel). We assume that all SRA files are in a specific folder.
```
module load parallel
parallel "fastq-dump --split-files --origfmt --gzip" ::: /path/to/SRA/*.sra
```
On the other hand if fastq files are available on a public repository (e.g. [EBI](https://www.ebi.ac.uk/ena/data/view/PRJNA348194)) we can download them directly using wget after copying the links to those files.

```
wget <link/to/fastq.gz>
```


We also need the genome file and associated GTF/GFF file for for Arabidopsis. We download these files directly from the [NCBI](https://www.ncbi.nlm.nih.gov/genome?term=NC_001284&cmd=DetailsSearch), or [plants Ensembl website](http://plants.ensembl.org/info/website/ftp/index.html) or the [Phytozome website](https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Gmax "Glycine max") (phytozome needs logging in and selecting the files) .

We downloaded the following files from [NCBI](https://www.ncbi.nlm.nih.gov/genome?term=NC_001284&cmd=DetailsSearch).
```
Genome Fasta File: GCF_000001735.3_TAIR10_genomic.fna
Annotation file: GCF_000001735.3_TAIR10_genomic.gff

```

# 2. Quality Check #

For this we will use `fastqc`, which is a tool that provides a simple way to do quality control checks on raw sequence data coming from high throughput sequencing pipelines ([link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). It provides various metrics to give a indication of how your data is. A high quality illumina RNAseq file should look something like [this](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html). Since there are 6 set of files (12 files total), and we need to run `fastqc` on each one of them, you can either do a `for` loop or use `parallel` command. We need to submit it as a job to the cluster, but the command should have:

```
module load fastqc
module load parallel
parallel "fastqc {} - <fq_out_directory>" ::: *.fastq.gz
```
Because we have a total of 6 quality outputs, we will have 6 html files and 6 zip files. We can use `multiqc` to aggregate the outputs and get a single html file to scan the quality of all the libraries.

```
cd fq_out_directory
module load python_3
multiqc .
....
....


```
This will give you a combined html file folder with containing three files descbing the various statistics:
```
ls  
........... multiqc_data (Directory)
............. multiqc_report.html
.............
..........
```
![file](/assets/Atrx_multiqc_report.html)

If you change to multiqc_data directory you see these files.
```
cd multiqc_data
ls
multiqc_fastqc.txt  
multiqc_general_stats.txt  
multiqc_sources.txt

```

 The main metrics to check are:
 * Per base sequence quality
 * Adapter content
 * Per base N content

Once you are happy with the results, proceed with the mapping part. If not, then perform quality trimming (see [here](/fastq-quality-trimming.md))

# 3. Mapping reads to the genome #

There are several mapping programs available for aligning RNAseq reads to the genome. Generic aligners such as BWA, BowTie2, BBMap etc., are not suitable for mapping RNAseq reads because they are not splice aware. RNAseq reads are mRNA reads that only contain exonic regions, hence mapping them back to the genome requires splitting the individual read, only done by splice aware mappers. Here for this tutorial, we will use `HiSat2` (derivative of BowTie2 and a successor of Tophat2).


### HiSat2 for mapping ###

#### Hisat2 Index ####

For HiSat2 mapping, you need to first index the genome and then use the read pairs to map the indexed genome (one set at a time). For indexing the genome, `HiSat2` we use the `hisat2-build` command as follows in a slurm script:

```
#!/bin/bash
#set -o xtrace
#SBATCH -p debug74 # optional: to specify the queue
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
GENOME="/project/isu_gif_vrsc/Siva/reference_genomes/GCF_000001735.3_TAIR10_genomic.fna"
hisat2-build $GENOME ${GENOME%.*}

```

Once complete, you should see a number of files with `.ht2` extension. These are the index files. At the mapping step we simply refer to the index.


#### Hisat2 Mapping ####

For mapping, each set of reads (forward and reverse or R1 and R2), we first set up a run_hisast.sh script.
```
#!/bin/bash
set -o xtrace
# set the rerefernce index:
GENOME="/project/isu_gif_vrsc/Siva/reference_genomes/GCF_000001735.3_TAIR10_genomic"
# make an output directory to store the output aligned files
mkdir -p /project/isu_gif_vrsc/Siva/HS_out
# set that as the output directory
ODIR="/project/isu_gif_vrsc/Siva/HS_out"


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

For setting it up to run with each set of file, we can set a SLURM script that loops over each fastq file. Note that this script calls the run_hista2.sh script for each pair of fastq file supplied as its argument.
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
/project/isu_gif_vrsc/Siva/run_hisat.sh ${fq1} ${fq2};
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

For quantifying transcript abundance from RNA-seq data, there are many programs available. Two most popular tools include, `featureCounts` and `HTSeq` and many other tools. We will need a file with aligned sequencing reads (SAM/BAM files generated in previous step) and a list of genomic features (from the GFF file). `featureCounts` is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It also outputs stat info for the overall summarization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons. We can run featureCounts on all SAM/BAM files at the same time or individually.

#### featureCounts ####
You will need `subread` and `parallel` modules loaded.
```
ANNOT="/project/isu_gif_vrsc/Siva/reference_genomes/GCF_000001735.3_TAIR10_genomic.gff"
mkdir -p /project/isu_gif_vrsc/Siva/HS_out/counts
ODIR="/project/isu_gif_vrsc/Siva/HS_out/counts"


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
# Program:featureCounts v1.5.2; Command:"featureCounts" "-T" "4" "-s" "2" "-p" "-t" "gene" "-g" "ID" "-a" "/project/isu_gif_vrsc/Siva/reference_genomes/GCF_000001735.3_TAIR10_genomic.gff" "-o" "/project/isu_gif_vrsc/Siva/HS_out/counts/SRR4420298.gene.txt" "SRR4420298.bam"

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

Again, there are few options here. You can use `DESeq2`,`edgeR`, or `QuasiSeq` (and many more!). Here we will describe how to do this with reference using DESeq2.


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
![qc-dispersions.png](/assets/qc-dispersions.png)
```
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

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


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Plots to Examine Results:

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
```
