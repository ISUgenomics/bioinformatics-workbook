---
title: Kraken2 tutorial
layout: single
author: Rick Masonbrink
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# RSeQC tutorial to identify stranding information in RNAseq


RNA-seq technology has evolved and many datasets are now stranded, meaning strand specificity of origin for each transcript is retained in RNAseq.   
Now comes the conundrum, what if this RNA-seq is not one you generated?  There is always a need to confirm stranding information when not using self-generated data.

## Install RSeQC
```

#Attempt with miniconda did not succeeed.
ml miniconda3
conda create -n rseqc
source activate rseqc
#this did not finish, using pip
conda install -c bioconda rseq

#this works
conda install -c bioconda bedops
#this worked.
pip install --user RSeQC

```
## Generate bam files to assess for stranding information.


#### get the fastq files
```
#get the fastq files
module load sra-toolkit

#unstranded samples
fastq-dump --outdir 01_Align/  --split-files SRR1573504
fastq-dump --outdir 01_Align/  --split-files SRR1573505

#stranded samples
fastq-dump --outdir 01_Align/  --split-files SRR13332812	&
fastq-dump --outdir 01_Align/  --split-files SRR13332813 &
```


#### Get the genome and create the index
```
#download the gene prediction .gff
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff.gz
#download the genome and change the name to something simple.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/167/145/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
gunzip GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna.gz
ln -s GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna Maize.fa

#I am using hisat2 here, but you can use whatever aligner you like.
ml hisat2; hisat2-build Maize.fa Maize
```



### Alignment script
```
At this point we just assume that the reads are unstranded and align.  
If you have a huge number of reads that will take a significant amount of time to align, then take a subset of your fastq files
i.e. head -n 200000 R1.FASTQ >truncated_R1.fastq and head -n 200000 R2.fastq >truncated_R2.fastq

#runHiSat2bam.sh
################################################################################################################################
#note premake the database, hisat2-build genome.fa genome
#Script expects the database to be named the same as the genome file without the last extension (Genome), not (Genome.fa)
#Example run
#sh runHiSat2bam.sh {read1.fastq} {read2.fastq} {Directory of Database, not actual database without / at the end} {genomeName}

PROC=36
R1_FQ="$1"
R2_FQ="$2"
DBDIR="$3"
GENOME="$4"

module load hisat2
hisat2 -p ${PROC}  -x ${GENOME%.*} -1 ${R1_FQ} -2 ${R2_FQ}  -S ${R1_FQ%.*}.sam

module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
mkdir ${R1_FQ%.*}_temp
samtools sort  -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam
################################################################################################################################
```

### Align
```
This places the files in the order that the script requires: (sh runHiSat2bam.sh {read1.fastq} {read2.fastq} {Directory of Database, not actual database without / at the end} {genomeName})
paste <(ls -1 *_1.fastq) <(ls -1 *_2.fastq) |sed 's/\t/ /g' |while read line ; do echo "sh runHiSat2bam.sh "$line" /work/gif/TranscriptomicsWorkshop/remkv6/01_Maize/01_Align Maize.fa";done >Maizealign.sh
```

### Create bed12 file from genic gff
```
#this converts your gff to bed12 format
ml bedops
gff2bed <GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff >maizeGFF.bed12
```


### run RSeQC
```
#assesses the strand information of your alignment
for f in *sorted.bam; do infer_experiment.py --i $f -r maizeGFF.bed12 ;done

```

### Results
```

```
