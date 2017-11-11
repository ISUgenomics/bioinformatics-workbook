# RNA-Seq data Analysis


This document will guide you through the RNAseq analysis, starting from the quality checking through  getting the differential gene expression results. The next part of the wiki series will guide you through some of the down stream analysis that you can do to the results obtained here. Here is the overview of the RNAseq analysis covered in this tutorial. We have downloaded an Arabidopsis dataset from NCBI for this purpose. Check the [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA348194) page for more information:


### Overview ###
![**Figure 1.**: Overview of the RNAseq workflow](/assets/RNAseq_1.png)


### Experimental design ###

This experiment compares WT and atrx-1 mutant to analyze how ATRX chaperone loss of function results in changes in gene expression. RNA was isolated from three WT replicates and three mutant replicates using Trizol. Transcriptome was enriched/isolated using the plant RiboMinus kit for obtaining total RNA. RNA-seq was carried out in Illumina Hiseq 2500. The sequencing reads were generated as paired-end data, hence we have 2 files per replicate.


| Condition | replicate 1 | replicate 2 | replicate 3 |
| --- | --- | --- | --- |
| WT | SRR4420293_1.fastq.gz <br> SRR4420293_2.fastq.gz | SRR4420294_1.fastq.gz <br> SRR4420294_2.fastq.gz | SRR4420295_1.fastq.gz <br> SRR4420295_2.fastq.gz |
| atrx-1 | SRR4420296_1.fastq.gz <br> SRR4420296_2.fastq.gz| SRR4420297_1.fastq.gz <br> SRR4420297_2.fastq.gz| SRR4420298_1.fastq.gz <br> SRR4420298_2.fastq.gz |

### 1. Download the data from NCBI ###

Generally if the data is hosted at your local sequencing center you could download through a web interface or using `wget` or `curl` commands. In this case, however, we first download the SRA files from the public archives in NCBI in bulk using aspera high speed file transfer.

```
module load <path/to/sra-toolkit>
module load <path/to/edirect>
module load <path/to/parallel>
esearch -db sra -query PRJNA276699 | efetch --format runinfo |cut -d "," -f 1 | awk 'NF>0' | grep -v "Run" > srr_numbers.txt
while read line; do echo "prefetch --max-size 100G --transport ascp --ascp-path \"/path/to/aspera/<version>/etc/asperaweb_id_dsa.openssh\" $line"; done<srr_numbers.txt > prefetch.cmds
parallel <prefetch.cmds
```
After downloading the SRA files, we have to convert it to fastq format. We can use the fast-dump command as follows: (this step is slow and if possible run these commands using gnu parallel). He re we assume that all SRA files are in a specific folder.
```
module load parallel
parallel "fastq-dump --split-files --origfmt --gzip" ::: /path/to/SRA/*.sra
```
On the other hand if fastq files are available on any public repository we can download them directly using wget

```
wget <link/to/fastq.gz>
```


We also need the genome file and associated GTF/GFF file for for Arabidopsis. We will download them directly from the [Phytozome website](https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Gmax "Glycine max") (needs logging in and selecting the files) or [plants Ensembl website](http://plants.ensembl.org/info/website/ftp/index.html). Specifically, you will need;
```
Gmax_275_Wm82.a2.v1.gene.gff3.gz
Gmax_275_v2.0.fa.gz
```
jvggggg

![](/assets/RNAseq_2.png)
**Figure 2:** Files needed for the RNAseq tutoial, genome assembly (unmasked) from the "assembly" directory and gff3 from the "annotation" directory

```
# decompress files
gunzip Gmax_275_Wm82.a2.v1.gene.gff3.gz
gunzip Gmax_275_v2.0.fa.gz
```

### 2. Quality Check ###

For this we will use `fastqc`, which is a tool that provides a simple way to do quality control checks on raw sequence data coming from high throughput sequencing pipelines ([link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). It provides various metrics to give a indication of how your data is. A high qulaity illumina RNAseq file should look something like [this](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html). Since there are 9 set of files (18 files total), and we need to run `fastqc` on each one of them, you can either do a `for` loop or use `parallel` command. We need to submit it as a job to the cluster, but the command should have:

```
module load fastqc
module load parallel
parallel "fastqc {}" ::: *.fastq.gz
```

You will find `.html` files once the job is complete. You can open them using the firefox browser on the HPC (see guide [here](https://gif.biotech.iastate.edu/how-view-files-remote-machine-without-downloading-locally) or download it locally to view them in your local browser. The main metrics to check are:
 * Per base sequence quality
 * Adapter content
 * Per base N content

Once you are happy with the results, proceed with the mapping part. If not, then perform quality trimming (see [here](/fastq-quality-trimming.md))

### 3. Mapping reads to the genome ###

There are several mapping programs available for aligning RNAseq reads back to the genome. Generic aligners such as BWA, BowTie2, BBMap etc., are not suitable for mapping RNAseq reads because they are not splice aware. RNAseq reads are mRNA reads that only contain exoninc regions, hence mapping them back to the genome requires splitting the individual read, only done by splice aware mappers. Here for this tutorial, we will use `HiSat2` (derivative of BowTie2 and a successor of Tophat2).

**Note: you don't have to run all three mapping programs, use any one of the below methods**

#### Using HiSat2 for mapping ####

For HiSat2 mapping, you need to first index the genome and then use the read pairs to map the indexed genome (one set at a time). For indexing the genome, `HiSat2` as is packaged with the `hisat2-build` script. Building index is as follows:

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
Once complete, you should see number of files with `.ht2` extension. These are
the index files.

For mapping, each set of reads (forward and reverse or R1 and R2), we set up a
run script. This script can also be found on our [GitHub page](https://github.com/ISUgenomics/common_analyses/blob/master/runHISAT2.sh).

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

For setting it up to run with each set of file, we can set a SLURM script that
loops over each fastq file:
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
### 2. Abundance estimation ###

For quantifying transcript abundance from RNA-seq data, there are many programs we can use. Two most popular tools inlcude, `featureCounts` and `HTSeq` and . We will need a file with aligned sequencing reads (SAM/BAM files generated in previous step) and a list of genomic features (from the GFF file). `featureCounts` is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. It also outputs stat info for the overall summarization results, including number of successfully assigned reads and number of reads that failed to be assigne due to various reasons. We can run featureCounts on all SAM/BAM files at the same time or individually.


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


Summary Files: These give the summary of reads that were either ambiguous, multimapped, mapped to no features or unammped among other things

```
SRR4420298.gene.txt.summary
SRR4420293.gene.txt.summary
SRR4420295.gene.txt.summary
SRR4420296.gene.txt.summary
SRR4420294.gene.txt.summary
SRR4420297.gene.txt.summary

```
Using the following command, a combination of paste and awk, we can produce a single count table taht can be used to load this data into R for differential expresion

```
paste <(awk 'BEGIN {OFS="\t"} {print $1,$7}' SRR4420293.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420294.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420295.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420296.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420297.gene.txt) <(awk 'BEGIN {OFS="\t"} {print $7}' SRR4420298.gene.txt) | grep -v '^\#' > At_count.txt
```

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

### 5. Differential Gene Expression analysis ###

Again, there are few options here. You can use `edgeR`, `DESeq2`, or `QuasiSeq` (and many more!). Here we will discribe how to do this with reference to the data we have. You can easily modify it to suit your needs (eg., different number of samples/repliates/conditions)

**Note: you don't have to run all three methods, use any one**

#### Option A: edgeR ####
Run the following code on RStudio or R terminal
```
# import data
datain <- read.delim("counts.txt",row.names="Geneid")

# experimental design
DataGroups <- c("CTL", "CTL","CTL","CTL", "TRT", "TRT", "TRT", "TRT")

# load edgeR
library(edgeR)

# create DGE object of edgeR
dgList <- DGEList(counts=datain,group=factor(DataGroups))

# filter data to retain genes that are represented at least 1 counts per million (cpm) in at least 2 samples
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
dgList$samples$lib.size <- colSums(d$counts)

# normalization using TMM method
dgList <- calcNormFactors(dgList, method="TMM")

## data exploration
# MDS plot
png("plotmds.png")
plotMDS(dgList, method="bcv", col=as.numeric(d$samples$group))
dev.off()

# Dispersion estimates
design.mat <- model.matrix(~ 0 + dgList$samples$group)
colnames(design.mat) <- levels(dgList$samples$group)
dgList <- estimateGLMCommonDisp(dgList,design.mat)
dgList <- estimateGLMTrendedDisp(dgList,design.mat, method="power")
dgList <- estimateGLMTagwiseDisp(dgList,design.mat)
png("plotbcv.png")
plotBCV(dgList)
dev.off()

# Differentail expression analysis
fit <- glmFit(dgList, design.mat)
lrt <- glmLRT(fit, contrast=c(1,-1))
edgeR_results <- topTags(lrt, n=Inf)

# plot log2FC of genes and highlight the DE genes
deGenes <- decideTestsDGE(lrt, p=0.05)
deGenes <- rownames(lrt)[as.logical(deGenes)]
png("plotsmear.png")
plotSmear(lrt, de.tags=deGenes)
abline(h=c(-1, 1), col=2)
dev.off()

# save the results as a table
write.table(edgeR_results, file="Results_edgeR.txt")
```

#### Option B: DESeq2 ####
Again, this should be done either on RStudio or in R terminal. Following are the steps

```
## RNA-seq analysis with DESeq2
## Largely based on Stephen Turner, @genetics_blog
## https://gist.github.com/stephenturner/f60c1934405c127f09a6

# Import the data
countdata <- read.table("counts.txt", header=TRUE, row.names=1)

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition (first four are controls, second four and third four contain two different experiments)
(condition <- factor(c(rep("ctl", 4), rep("inf1", 4), rep("inf2", 4))))

# Analysis with DESeq2

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

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

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I (Stephen Turner) like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
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
#### Option C: QuasiSeq ####
Also in RStudio or R terminal. As a first step, save the following code as a file named `QLresultsPvaluePlot.R` using any text editor and place it in the same directory as the count data is in.
```
QLresultsPvaluePlot<-function(QLfit,Strname){
filename=Strname
results<-QL.results(QLfit,Plot=F)
designNum<-dim(results$P.values$QLSpline)[2]
designNames<-colnames(results$P.values$QLSpline)
for (i in 1:designNum){
   print(i)
   if (min(results$P.values$QLSpline[,i])<1 && min(results$P.values$QLSpline[,i])!="NaN"){
       Rnames<-rownames(dataIn)
       if (min(results$Q.values$QLSpline[,i])<10 && min(results$Q.values$QLSpline[,i])!="NaN"){
          if (length(which(results$Q.values$QLSpline[,i]<10))>1){
            meanTrt<-apply(dataIn.norm[which(results$Q.values$QLSpline[,i]<10),which(trt==2)],1,mean)
            meanWt<-apply(dataIn.norm[which(results$Q.values$QLSpline[,i]<10),which(trt==1)],1,mean)
            FoldTrtoverWt <- meanTrt/meanWt
            logTwoFC <-log2(FoldTrtoverWt)
            outData<-cbind(as.matrix(dataIn.norm[which(results$Q.values$QLSpline[,i]<10),]),meanTrt,meanWt,as.matrix(results$P.values$QLSpline[which(results$Q.values$QLSpline[,i]<10),i]),as.matrix(results$Q.values$QLSpline[which(results$Q.values$QLSpline[,i]<10),i]),FoldTrtoverWt,logTwoFC)
            colnames(outData)<-c(colnames(dataIn),"mean_treat","mean_control","Pvalues","Qvalues","fold_change","Log2FC")
            write.table(outData,file=paste(filename,".FulldesignVS.",i,".txt",sep=""))
            }
          if (length(which(results$Q.values$QLSpline[,i]<10))==1){
            outData<-cbind(matrix(dataIn[which(results$Q.values$QLSpline[,i]<10),],1,dim(dataIn)[2]),as.matrix(results$P.values$QLSpline[which(results$Q.values$QLSpline[,i]<10),i]),as.matrix(results$Q.values$QLSpline[which(results$Q.values$QLSpline[,i]<10),i]),(sign(mean(dataIn.norm[which(results$Q.values$QLSpline[,i]<10),which(trt==2)])-mean(dataIn.norm[which(results$Q.values$QLSpline[,i]<10),which(trt==1)])))*mean(dataIn.norm[which(results$Q.values$QLSpline[,i]<10),which(trt==2)])/mean(dataIn.norm[which(results$Q.values$QLSpline[,i]<10),which(trt==1)]))
            colnames(outData)<-c(colnames(dataIn),"Pvalues","Qvalues","fold_change")
            write.table(outData,file=paste(filename,".FullvsDesignVS.",i,".txt",sep=""))
            }
        }
        pdf(file=paste(filename,".",i,".pdf",sep=""),width=5,height=5)
           a<-hist(results$P.values$QLSpline[,i],breaks=seq(0,1,.01),main=paste(Strname,designNames[i]),cex.main=.5)
           b<-a$counts[1]*.75
           bb<-a$counts[1]*.65
           bbb<-a$counts[1]*.55
           text(.5,b,paste("Number of genes qvalue below 0.5 = ",as.character( dim(as.matrix(dataIn[which(results$Q.values$QLSpline[,i]<0.5),i]))[1])),cex=.8)
           text(.5,bb,paste("Number of genes qvalue below 0.3 = ",as.character( dim(as.matrix(dataIn[which(results$Q.values$QLSpline[,i]<0.3),i]))[1])),cex=.8)
           text(.5,bbb,paste("Number of genes qvalue below 0.1 = ",as.character( dim(as.matrix(dataIn[which(results$Q.values$QLSpline[,i]<.1),i]))[1])),cex=.8)
        dev.off()
    }
}
}
```

Next, run these steps on RStudio by setting the work directory to the counts data directory.
```
# set the work directory
setwd("C:/Users/arunk/Google Drive/PostDoc/projects/20170707_RSmith_MosquitoRNAseq/QuassiSeq")
# source the code you just created
source("QLresultsPvaluePlot.R")
# Import the data
uniq<-as.matrix(read.table("counts.txt", header=TRUE,  row.names = 1))

# Check dimensions
cols<-dim(uniq)[2]
# remove the rows with all zero counts
colsm1<-cols
dataIn2<-(uniq[which(rowSums(uniq[,1:cols])>colsm1),])
dataIn3<-dataIn2[which(rowSums(sign(dataIn2[,1:cols]))>1),]
dataIn<-as.matrix(dataIn3)
# normalize data using upperquartile of 0.75
log.offset<-log(apply(dataIn, 2, quantile,.75))
upper.quartiles<-apply(dataIn,2,function(x) quantile(x,0.75))
# calculate scaling factors
scalingFactors<-mean(upper.quartiles)/upper.quartiles
dataIn.norm<-round((sweep(dataIn,2,scalingFactors,FUN="*")))
# standard deviation
sd(dataIn[1,])
sd(dataIn.norm[1,])
# experimental design
trt<-as.factor(c(1,1,1,1,2,2,2,2))
mn<-as.numeric(rep(1,cols))
# QuasiSeq analysis
library(QuasiSeq)
design.list<-vector('list',2)
design.list[[1]]<-model.matrix(~trt)
design.list[[2]]<-mn
log.offset<-log(apply(dataIn, 2, quantile,.75))
fit2<-QL.fit(round(dataIn), design.list,log.offset=log.offset, Model='NegBin')
QLresultsPvaluePlot(fit2,paste("results_",1,2,sep=""))
QLresultsPvaluePlot(fit2,paste("results_",1,2,sep=""))
```
