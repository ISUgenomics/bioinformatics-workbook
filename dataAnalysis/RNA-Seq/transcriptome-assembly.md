# Transcriptome assembly using various methods

## Dataset:
This exercise uses a sample RNAseq data set from GEO database ([GSM1053997](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1053997)).

| Run       | # of Spots | # of Bases | Size    | Published  |
|-----------|------------|------------|---------|------------|
| SRR633726 | 6,571,208  | 1.2G       | 714.8Mb | 2015-07-22 |
| SRR633727 | 6,415,120  | 1.2G       | 641Mb   | 2015-07-22 |

## Download

Data can be downloaded using SRA-toolkit, directly from the NCBI and saved as fastq files. Start an interactive session, load the module and download files

```
salloc -N 1 -n 2 -t 4:00:00
prefetch SRR633726
prefetch SRR633727
fastq-dump --split-files --origfmt --gzip SRR633726
fastq-dump --split-files --origfmt --gzip SRR633727
```
There will be 4 files (paired) generated upon completion of these commands.



## 1. Oases/Velvet
Oases is a de novo transcriptome assembler designed to produce transcripts from short read sequencing technologies, such as Illumina. This has 3 steps, first you run `velveth` to create hash, then you run `velvetg` to make assemblies and finally run the `oases` to generate transcriptome assembly.

Run velveth
```
module load velvet
velveth velvet-oases 21,71,10 \
  -shortPaired -fastq -separate ../R1_pairedout.fastq ../R2_pairedout.fastq
  -shortPaired2 -fastq -separate ../R1_pairedout.fastq ../R2_pairedout.fastq
```
Run velvetg
```
for dir in velvet-oases*; do
  velvetg $dir -read_trkg yes -ins_length 200;
done
```
Run Oases
```
for dir in velvet-oases*; do
  oases $f -min_trans_lgth 100 -ins_length 200;
done
```


## 2. SOAPdenovo-Trans
SOAPdenovo-Trans is a de novo transcriptome assembler basing on the SOAPdenovo framework, adapt to alternative splicing and different expression level among transcripts.The assembler provides a more accurate, complete and faster way to construct the full-length transcript sets.


## 3. Trinity (_de novo_)

For this, you can either merge to create single R1 file  and single R2 file or provide them as separate files. If providing them as separate files:


```
Trinity --seqType fq \
        --max_memory 110G  \
        --left SRR633726_1.fq.gz,SRR633727_1.fq.gz \
        --right SRR633726_2.fq.gz,SRR633727_2.fq.gz \
        --CPU 16
```  

## 4. Trinity (genome guided)

Here, you need to align the reads to the genome and then provide the `bam` file for Trinity. For this exercise, we will use STAR aligner.

Download the reference genome, index and align the short reads (splice aware)

```
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
module load star
STAR \
  --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir TAIR10 \
  --genomeFastaFiles TAIR10_chr_all.fas \

STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir TAIR10 \
  --readFilesCommand zcat \
  --outFileNamePrefix SRR633726 \
  --readFilesIn SRR633726_1.fq.gz SRR633726_2.fq.gz

STAR \
 --runMode alignReads \
 --runThreadN 16 \
 --genomeDir TAIR10 \
 --readFilesCommand zcat \
 --outFileNamePrefix SRR633727 \
 --readFilesIn SRR633727_1.fq.gz SRR633727_2.fq.gz
```

Now, convert the sam files to bam, sort and merge

```
module load samtools
samtools view --threads 16 -b -o SRR633726.bam SRR633726.sam
samtools view --threads 16 -b -o SRR633727.bam SRR633727.sam
samtools sort -m 3G -o SRR633726_sorted.bam -T SRR633726_temp --threads 16 SRR633726.bam
samtools sort -m 3G -o SRR633727_sorted.bam -T SRR633727_temp --threads 16 SRR633727.bam
ls -1 *sorted.bam > bamfiles.txt
samtools merge -b bamfiles.txt --threads 16 merged_SRRs.bam
```

Run trinity.

```
Trinity \
   --max_memory 110G \
   --genome_guided_max_intron 10000 \
   --full_cleanup \
   --output merged_SRRs_trinity \
   --CPU 28 \
   --genome_guided_bam merged_SRRs.bam
```

## 5. IDBA-tran
IDBA-Tran is an iterative De Bruijn Graph De Novo short read assembler for transcriptome. It is purely de novo assembler based on only RNA sequencing reads.


## 6. Cufflinks


Use the merged bam file from Step 4 (Triniyt Genome Guided), for this step.
