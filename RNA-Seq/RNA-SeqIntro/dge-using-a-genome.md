# RNA-Seq data Analysis


This wiki will guide you through the RNAseq analysis, starting from the quiality checking till getting the differntial gene expression results. The next part of the wiki series will guide you through some of the down stream analysis that you can do to the results obatined here. Here is the overview of the RNAseq analysis covered in this tutorial.

### Overview ###
![**Figure 1.**: Overview of the RNAseq workflow
](/assets/RNAseq_1.png)

 
### Experimental design ###

This experiment consists of 2 conditions. The first condition is "control" which is mock-infected soybean plants. The second condition is the infected plants with the actual pathogen. Both conditions have four replicates each (total we have 8 pairs of fastq files, 4 pairs belonging to control, and 4 pairs bleonging to pathogen infected soybean). The reads were generated as paired-end data, hence we have 2 files per replicate.


| Condition | Replicate 1 | Replicate 2 | Replicate 3 |  Replicate 4 |
| --- | --- | --- | --- | --- |
| Control | Control.A_R1.fastq.gz <br> Control.A_R2.fastq.gz | Control.B_R1.fastq.gz <br> Control.B_R2.fastq.gz | Control.C_R1.fastq.gz <br> Control.C_R2.fastq.gz |  Control.D_R1.fastq.gz <br> Control.D_R2.fastq.gz  |
| Infected | Infected.A_R1.fastq.gz <br> Infected.A_R2.fastq.gz| Infected.B_R1.fastq.gz <br> Infected.B_R2.fastq.gz| Infected.C_R1.fastq.gz <br> Infected.C_R2.fastq.gz | Infected.D_R1.fastq.gz <br> Infected.D_R2.fastq.gz |

### 1. Download the data ###

For downloading the data, you can use `wget` or `curl` commands, if the data is hosted somewhere. If not, you might have to upload the data to the HPC either using `scp` command or using `rsync` (if data is located locally on your computer), or use `globusURL` to get the data from other computer. Here we will assume that you have the data in our DNA facility (at Iowa State University) and you have access to those files. We will use `wget` command to download them.

```
wget http://upart.biotech.iastate.edu/pub/5_NNNN_01_1_HCC22_956.tar
wget http://upart.biotech.iastate.edu/pub/5_NNNN_01_2_HCC22_956.tar
wget http://upart.biotech.iastate.edu/pub/5_NNNN_01_3_HCC22_956.tar
wget http://upart.biotech.iastate.edu/pub/5_NNNN_01_4_HCC22_956.tar
wget http://upart.biotech.iastate.edu/pub/5_NNNN_01_5_HCC22_956.tar
wget http://upart.biotech.iastate.edu/pub/5_NNNN_01_6_HCC22_956.tar
wget http://upart.biotech.iastate.edu/pub/5_NNNN_01_7_HCC22_956.tar
wget http://upart.biotech.iastate.edu/pub/5_NNNN_01_8_HCC22_956.tar
```

Note that these weblinks are inactive, so if you actually run these commands it will fail as they don't point to any files. Once downloaded, you can untar the archive and you will find the fastq files. To untar:

```
module load parallel
parallel "tar xf {}" ::: *.tar
```

Here we load the `parallel` and then run it effeciently and in parallel on all the tar files. The other option is to run it in a `for` loop, which will take considerable amount of time as it untars one file at a time. After this step, you will have gzipped fastq files. 

We will also need the genome file and associated GTF/GFF file for this wiki. Since the data is for Soybean, we will donwload them directly from the [Phytozome website](https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Gmax "Glycine max") (needs logging in and selecting the files). Specifically, you will need;
```
Gmax_275_Wm82.a2.v1.gene.gff3.gz 
Gmax_275_v2.0.fa.gz
```

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
 
Once you are happy with the results, proceed with the mapping part.

### 3. Mapping reads to the genome ###

There are several mapping programs available for aligning RNAseq reads back to the genome. Generic aligners such as BWA, BowTie2, BBMap etc., are not suitable for mapping RNAseq reads because they are not splice aware. RNAseq reads are mRNA reads that only contain exoninc regions, hence mapping them back to the genome requires splitting the individual read, only done by splice aware mappers. Here for this tutorial, we will use `HiSat2` (derivative of BowTie2), `STAR` aligner and `GSNAP`.

**Note: you don't have to run all three mapping programs, use any one of the below methods**

#### Option A: Use HiSat2 for mapping #### 

For HiSat2 mapping, you need to first index the genome and then use the read pairs to map the indexed genome (one set at a time). For indexing the genome, `HiSat2` as is packaged with the `hisat2-build` script. Building index is as follows:

```
hisat2-build -p 16 Gmax_275_v2.0.fa Gmax_275_v2.0_hisat2
# -p for number of processors
# first argument is the fasta file
# second argument is the base name for indexed files
```
Once complete, you should see number of files with `.ht2l` extension. These are our index files.

For mapping, each set of reads (forward and reverse or R1 and R2), we will set up a run script. This script can also be found on our [GitHub page](https://github.com/ISUgenomics/common_analyses/blob/master/runHISAT2.sh).

```
#!/bin/bash
module load hisat2
module load samtools
DBDIR="/path/containing/HiSat2_index_files"
GENOME="Gmax_275_v2.0_hisat2"
# number of processors to use
p=16
R1_FQ="$1"
R2_FQ="$2"
# output file prefix (uses a portion of the fastq filename, but can be changed if these are not unique)
OUTPUT=$(basename ${R1_FQ} |cut -f 1 -d "_");
# run the HiSat2 aligner
hisat2 \
  -p ${p} \
  -x ${DBDIR}/${GENOME} \
  -1 ${R1_FQ} \
  -2 ${R2_FQ} | \
  -S  ${OUTPUT}.sam &> ${OUTPUT}.log
# convert_gsnap.same file to bam file format
samtools view \
  --threads 16 \
  -b \
  -o ${OUTPUT}.bam \
     ${OUTPUT}.sam
# sort the bam file based on co-ordinates
samtools sort \
  -m 7G \
  -o ${OUTPUT}_sorted.bam \
  -T ${OUTPUT}_temp \
  --threads 16 \
  ${OUTPUT}.bam
```

For setting it up to run with each set of file, we will use this loop:
```
for fastq in *R1.fastq.gz; do
fastq2=$(echo $fastq | sed 's/R1.fastq.gz/R2.fastq.gz/g'); 
echo "./runHISAT2.sh ${fastq} ${fastq2}"; 
done > hisat2.cmds
```

For creating PBS/Slurm submission scripts, use either `makePBSs.py` or `makeSLURMs.py` from the `common_scripts` directory on [GitHub](https://github.com/ISUgenomics/common_scripts). Submit jobs using the for loop.

```
makeSLURMs.py 1 hisat2.cmds
for sub in hisat2_*.sub; do
sbathc $sub; 
done
```

This should create, following files as output:
```
Control.A_sorted.bam
Control.B_sorted.bam
Control.C_sorted.bam
Control.D_sorted.bam
Infected.A_sorted.bam
Infected.B_sorted.bam
Infected.C_sorted.bam
Infected.D_sorted.bam
```

#### Option B: Use STAR for mapping #### 

Again, we need to index the genome first:

```
#!/bin/bash
module load star
FASTA="Gmax_275_v2.0.fa"
GFF="Gmax_275_Wm82.a2.v1.gene.gff3"
DB=$(pwd)
mkdir -p ${DB}/${FASTA%.*}_star
STAR \
  --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir ${DB}/${FASTA%.*}_star
  --genomeFastaFiles ${FASTA}
  --sjdbGTFfile ${GFF}
  --sjdbOverhang 99
```
Once this completes, you should see a folder named `Gmax_275_v2.0_star` with lots of files in it. This is our indexed genome.

For running mapping, we will set up a run script, like we did for HiSat2. You can also find the run script on our [GitHub](https://github.com/ISUgenomics/common_analyses/blob/master/runSTAR.sh)

```
#!/bin/bash
R1="$1"
R2="$2"
OUT=$(basename ${R1} |cut -f 1 -d "_");
GFF="/home/arnstrm/arnstrm/GMAPDB/Gmax_275_Wm82.a2.v1.gene.gff3"
DB="/home/arnstrm/arnstrm/GMAPDB/Gmax_275_v2.0_star"
STAR \
 --runMode alignReads \
 --runThreadN 16 \
 --genomeDir ${DB} \
 --readFilesCommand zcat \
 --outFileNamePrefix ${OUT}_star \
--readFilesIn ${R1} ${R2}
```
For setting it up to run with each set of file, we will use this loop:
```
for fastq in *R1.fastq.gz; do
fastq2=$(echo $fastq | sed 's/R1.fastq.gz/R2.fastq.gz/g'); 
echo "./runSTAR.sh ${fastq} ${fastq2}"; 
done > star.cmds
```

For creating PBS/Slurm submission scripts, use either `makePBSs.py` or `makeSLURMs.py` from the `common_scripts` directory on [GitHub](https://github.com/ISUgenomics/common_scripts). Submit jobs using the for loop.

```
makeSLURMs.py 1 star.cmds
for sub in star_*.sub; do
qsub $sub; 
done
```

This should create, following files as output:
```
Control.A_star.sam
Control.B_star.sam
Control.C_star.sam
Infected1_star.A.sam
Infected1_star.B.sam
Infected1_star.C.sam
Infected2_star.A.sam
Infected2_star.B.sam
Infected2_star.C.sam
```

#### Option C: Use GSNAP for mapping #### 

Final alternative is to use GSNAP aligner. For indexing, use the script below:
```
#!/bin/bash
FASTA="Gmax_275_v2.0.fa"
DB="/work/GIF/arnstrm/GENOMEDB"
module load gmap-gsnap
gmap_build -d ${FASTA%.*}_gsnap -D ${DB} ${FASTA}
```
A directory, named `Gmax_275_v2.0_gsnap` will be created with lots of files in it once the indexing is complete. Next, we will setup a run script for `GSNAP` (as found [here](00))
```
#!/bin/bash
GMAPDB="/work/GIF/arnstrm/GENOMEDB"
DB_NAME="Gmax_275_v2.0_gsnap"
R1="$1"
R2="$2"
OUTFILE=$(basename ${R1} |cut -f 1 -d "_")
gsnap \
   -d ${DB_NAME} \
   -D ${GMAPDB}
   -t 16 \
   -B 5 \
   -N 1
   -m 5 \
   --gunzip \
   --fails-as-input \
   --input-buffer-size=10000000 \
   --output-buffer-size=10000000 \
   -A_gsnap.sam  ${R1} ${R2} > ${OUTFILE}_gsnap_gsnap.sam
```
For setting it up to run with each set of file, we will use this loop:
```
for fastq in *R1.fastq.gz; do
fastq2=$(echo $fastq | sed 's/R1.fastq.gz/R2.fastq.gz/g'); 
echo "./runGSNAP.sh ${fastq} ${fastq2}"; 
done > gsnap.cmds
```

For creating PBS/Slurm submission scripts, use either `makePBSs.py` or `makeSLURMs.py` from the `common_scripts` directory on [GitHub](https://github.com/ISUgenomics/common_scripts). Submit jobs using the for loop.

```
makeSLURMs.py 1 gsnap.cmds
for sub in gsnap_*.sub; do
qsub $sub; 
done
```

This should create, following files as output:
```
Control.A_gsnap.sam
Control.B_gsnap.sam
Control.C_gsnap.sam
Infected1.A_gsnap.sam
Infected1.B_gsnap.sam
Infected1.C_gsnap.sam
Infected2.A_gsnap.sam
Infected2.B_gsnap.sam
Infected2.C_gsnap.sam
```

Once we have the SAM or BAM files generated (from any of the 3 methods above), we can proceed with the abundence estimation.

### 4. Abundence estimation ###

For quantifying transcript abundance from RNA-seq data, there are many programs we can use. Two most popular tools inlcude, `HTSeq` and `featureCounts`. We will need a file with aligned sequencing reads (SAM/BAM files generated in previous step) and a list of genomic features (donwloaded GFF file).

**Note: you don't have to run both these abundence estimation methods, use any one**

#### Option A: HTSeq #### 

As we need to process one SAM/BAM file at a time, we will set up a run script as follows:
```
#!/bin/bash
GFF="Gmax_275_Wm82.a2.v1.gene.gff3"
INFILE="$1"
OUTFILE=$(echo $INFILE | cut -f 1 -d "_")_counts.txt
htseq-count \
   -m intersection-nonempty \
   -t gene \
   -i ID \
      $INFILE $GFF > $OUTFILE
```
```
for sam in *.sam; do
echo "./runHTSEQ.sh ${sam}"; 
done > htseq.cmds
```

For creating PBS/Slurm submission scripts, use either `makePBSs.py` or `makeSLURMs.py` from the `common_scripts` directory on [GitHub](https://github.com/ISUgenomics/common_scripts). Submit jobs using the for loop.

```
makeSLURMs.py 9 htseq.cmds
# as these commands run quickly, we will put them all in one `sub` file
qsub htseq_0.sub; 
done
```

You will have a text file with counts for every SAM/BAM file you provide. Next, we need to merge individual files to make a single file. We can do this using `awk` as follows:
```
awk '{arr[$1]=arr[$1]"\t"$2}END{for(i in arr)print i,arr[i]}' *_counts.txt >> merged_htseq_counts.tsv
mv merged_htseq_counts.tsv counts.txt
```




#### Option B: featureCounts #### 

`featureCounts` is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations. `featureCounts` takes as input SAM/BAM files and an annotation file including chromosomal coordinates of features. It outputs numbers of reads assigned to features (or meta-features). It also outputs stat info for the overall summrization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons (these reasons are included in the stat info). 
We can run this on all SAM/BAM files at the same time.

```
#!/bin/bash
GFF="Gmax_275_Wm82.a2.v1.gene.gff3"
OUTFILE="featureCounts"
module load subread
# package containing featureCounts script
featureCounts \
   -T 16 \
   -p \
   -t gene \
   -a  \
   -o ${OUTFILE}_counts.txt *.bam
```

The ouput will look something like this (__headers might be different, only first 10 lines displayed here__)

```
#Geneid	Chr	Start	End	Strand	Length	Control.A	Control.B	Control.C	Infected1.A	Infected1.B	Infected1.C	Infected2.A	Infected2.B	Infected2.C
Glyma.01G000100.Wm82.a2.v1	Chr01	27355	28320	-	966	24	0	51	93	126	91	121	32	52
Glyma.01G000200.Wm82.a2.v1	Chr01	58975	67527	-	8553	91	1	122	193	214	239	102	111	148
Glyma.01G000300.Wm82.a2.v1	Chr01	67770	69968	+	2199	9	0	12	22	18	21	3	11	18
Glyma.01G000400.Wm82.a2.v1	Chr01	90152	95947	-	5796	169	0	407	480	402	518	502	443	379
Glyma.01G000500.Wm82.a2.v1	Chr01	90289	91197	+	909	0	0	0	0	0	0	0	0	0
Glyma.01G000600.Wm82.a2.v1	Chr01	116094	127845	+	11752	149	4	310	374	402	529	304	352	300
Glyma.01G000700.Wm82.a2.v1	Chr01	143467	155573	+	12107	39	2	78	119	113	129	34	100	107
Glyma.01G000800.Wm82.a2.v1	Chr01	157030	157772	+	743	0	0	0	1	1	0	0	0	0
Glyma.01G000900.Wm82.a2.v1	Chr01	170534	193342	+	22809	240	1	517	759	760	859	462	658	494
```
Since we only need the counts (without the feature information), we will trim this table using `cut` command.
```
cut -f 1,7-15 featureCounts_counts.txt > featureCounts_counts_clean.txt
mv featureCounts_counts_clean.txt counts.txt
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



