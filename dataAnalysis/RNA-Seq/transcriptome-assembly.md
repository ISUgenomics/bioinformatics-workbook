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

```
SRR633726_1.fastq.gz
SRR633726_2.fastq.gz
SRR633727_1.fastq.gz
SRR633727_2.fastq.gz
```

Next, the reference genome. Since we are using _Arabidopsis thaliana_, we will need the genome for this exercise. The Araport genomes are available on [CyVerse data commons](http://datacommons.cyverse.org/browse/iplant/home/araport/public_data).

```
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa TAIR10_genome.fa
```
## Pre-processing

Genome guided transcript assembly program require you to map the raw reads on the genome and provide the bam files. You can use any splice aware mapping programs for this purpose, but in this example, we will just use HiSat2 program.
First, we will index the genome:
```
salloc -N 1 -n 1 -t 4:00:00
module load histat2
hisat2-build TAIR10_genome.fa TAIR10_INDEX
```
Next, we will make a run script `runHiSat2.sh` as follows (**IMPORTANT: please modify the resource requirements like CPU/threads for your machine**). Pay special  attention to the parameters used for mapping: `--dta-cufflinks`, we used this option here to create `XS` tags in the alignment as this will not be included by default (data used in this example is not strand specific). This tag provides the most crucial information for the assembly programs and if it is not included,  they will fail/take long time to finish.

 ```
#!/bin/bash
module load samtools
module load hisat2
CPU=28
R1="$1"
R2=$(echo $R1 |sed 's/_1.fastq.gz/_2.fastq.gz/g')
out=$(echo $R1 |cut -f 1 -d "_")
hisat2 -p $CPU -x TAIR10_INDEX --dta-cufflinks -1 $R1 -2 $R2 > ${out}.sam
samtools view --threads $CPU -b -o ${out}.bam ${out}.sam
samtools sort -o ${out}_sorted.bam -T ${out}_temp --threads $CPU ${out}.bam
```

Now, make the script executable:

```
chmod +x runHiSat2.sh
```

and run mapping
```
./runHiSat2.sh SRR633726_1.fastq.gz
./runHiSat2.sh SRR633727_1.fastq.gz
```

This will result in 2 sam files and 4 bam files (for 2 sets of fastq files). Out of these, we will only be needing 2 bam files ending with `_sorted.bam`

```
SRR633726_sorted.bam
SRR633727_sorted.bam
```


## Genome Guided Transcript assembly:

### 1. Trinity (genome guided)

With the above bam files, an example Trinity run is as follows:

```
Trinity \
   --max_memory 110G \
   --genome_guided_max_intron 10000 \
   --full_cleanup \
   --output merged_SRRs_trinity \
   --CPU 28 \
   --genome_guided_bam SRR633726_sorted.bam,SRR633727_sorted.bam
```
Please modify the options such as `max_memory` or `CPU` suitable for your resources. Unlike other transcript assembly programs (genome guided), Trinity doesn't output a GFF3 formatted file. If you need the GFF3 file, you will need to use gmap-gsnap program to align the reads back to genome and obtain GFF3 file.

```
gmap_build  -D $(pwd -P) -d TAIR10_gmap TAIR10_genome.fa
gmap -D $(pwd -P) -d TAIR10_gmap -B 4 -t 16 -f gff3_match_cdna Trinity_GG.fasta > Trinity_GG.gff3
```
Here, `Trinity_GG.fasta` being the output of Trinity program.


### 2. Class2

Using the bam files, you can merge them to generate single bam file:

```
module load samtools
ls *_sorted.bam > bamfiles.txt
samtools merge -b bamfiles.txt --threads 16 merged_SRR.bam
```
then, run Class2

```
module load class2
bam=merged_SRR.bam
out=merged_SRR
perl $CLASS2_HOME/run_class.pl \
   -a $bam \
   -o ${out}_class.gtf \
   -p 28 \
   --verbose \
   --clean
```
Change `-p` (number of threads) suitable for your cluster.

### 3. Cufflinks

The merged BAM file from the above step can also be used for this program. For running Cullflinks:

```
module load boost/1.50.0
module load gcc/4.8.4
module load eigen/3.2.8
module load samtools/0.1.19
module load python/2.7.11_gcc
module load cufflinks/2.2.1
bam="merged_SRR.bam"
out=$(basename $bam |cut -f 1 -d "_")
cufflinks \
   --output-dir $out \
   --num-threads 28 \
   --verbose \
   --no-update-check \
   $bam
```

### 4. StringTie

Again, merged bam file is necessary.

```
module load stringtie/1.3.3
bam="merged_SRR.bam"
out=$(basename $bam |cut -f 1 -d "_")
stringtie \
   ${bam} \
   -p 28 \
   -v \
   -o ${out}_stringtie.gtf
```


## De novo Transcript assembly:
### 1. Oases/Velvet
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


### 2. SOAPdenovo-Trans
SOAPdenovo-Trans is a de novo transcriptome assembler basing on the SOAPdenovo framework, adapt to alternative splicing and different expression level among transcripts. The assembler provides a more accurate, complete and faster way to construct the full-length transcript sets.

make a run script `runSOAPtrans.sh` as follows:
```
#!/bin/bash
module load soapdenovo-trans
srr=$1
kmer=61
SOAPdenovo-Trans-127mer all -s config.txt -o ${srr}_${kmer} -K ${kmer} -p 16
done
```
make it executable:
```
chmod +x runSOAPtrans.sh
```

and run it as:

```
./runSOAPtrans.sh SRR633726
./runSOAPtrans.sh SRR633727
```

### 3. Trinity (_de novo_)

For this, you can either merge to create single R1 file  and single R2 file or provide them as separate files. If providing them as separate files:


```
Trinity --seqType fq \
        --max_memory 110G  \
        --left SRR633726_1.fq.gz,SRR633727_1.fq.gz \
        --right SRR633726_2.fq.gz,SRR633727_2.fq.gz \
        --CPU 16
```  
