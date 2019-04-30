---
title: Introduction to Braker2
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Tutorial of how to run Braker2 gene prediction pipeline

I have outlined how to perform a gene annotation of the soybean cyst nematode (*Heterodera glycines*) using Braker2.1.

Here is the original publication of Braker1 which provides some of the background necessary to understand how Braker2 works.
[Braker1 pub](https://academic.oup.com/bioinformatics/article/32/5/767/1744611)


Here is a link to a powerpoint pdf that provides some explanation of how Braker2 works.
[Braker2.1](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=2ahUKEwjs97T3kqLeAhUJ7YMKHecCDaYQFjAAegQICRAC&url=https%3A%2F%2Fpag.confex.com%2Fpag%2Fxxvi%2Frecordingredirect.cgi%2Foid%2FRecording2918%2Fpaper30184_1.pdf&usg=AOvVaw26ajX_Y-Ml3G7M-k_ozx-R)

Here is a the user guide for setting different parameters in Braker2.1.
[Braker2.1 User Guide](https://academic.oup.com/bioinformatics/article/32/5/767/1744611)

##  General steps before running Braker2.1
There are more the few steps that will need to be performed before you can run the formal braker.pl script.

### Prerequisite data for prediction
1. Gather all transcriptional data for your organism, (ESTs, RNA-seq, Transcripts, Isoseq) (required)
2. Gather all proteins for your organism
3. Gather repeat libraries for genome masking
4. Get your genome (required)

### Prepare repeat masking and genome alignments
1. Create repeats library and softmask your genome with RepeatMasker
2. Trim and align your RNA-seq data to your genome
3. Align transcripts, ESTs, or other long transcriptional reads to genome with splice aware aligner
4. Convert Sam to Bam and sort
5. Run Braker2


## Create repeat database and mask genome
Repeats can magnificently affect gene predictions, so it is best to mask repeats in the genome prior to prediction. I will leave simple repeats unmasked, as they are frequently part of genes.
For further information on how to run repeatmodeler and repeatmasker a tutorial has already been created in the Bioinformatics workbook [Repeats Tutorial](https://github.com/ISUgenomics/bioinformatics-workbook/blob/master/dataAnalysis/ComparativeGenomics/RepeatModeler_RepeatMasker.md)
```
#run Repeatmodeler
module load GIF2/repeatmodeler/1.0.8
module load GIF2/perl/5.22.1
BuildDatabase -name SCNDB -engine ncbi -pa 16 genome738sl.polished.mitoFixed.fa
RepeatModeler -database SCNDB -engine ncbi -pa 16

#Run Repeatmasker
module load GIF2/repeatmasker/4.0.6
ln -s RM_3245.WedMay21605262018/consensi.fa.classified
RepeatMasker -pa 16 -gff -nolow -lib consensi.fa.classified genome738sl.polished.mitoFixed.fa
```
## Download RNA-seq
```
#/work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm
vi Fastq2download.list
############################
SRR6230579
SRR6230580
SRR6230581
SRR6230582
SRR6230583
SRR6230584
SRR6230585
SRR6230586
SRR6230587
############################
less Fastq2download.list |while read line ; do echo "fastq-dump --split-files --origfmt --gzip "$line" &"; done >sraDownloader.sh
sh sraDownloader.sh


ln -s ../01_Alignment2Genome/738GenomeDB/
ln -s ../../CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/genome738sl.polished.mitoFixed.fa


```

## Align RNA-seq
At the GIF, a very common procedure is to align RNA-seq data for subsequent use, so having a script read to go for this makes things a bit easier each time.  

This script "runAlignConvertSortStats.sh" will trim a set of paired reads, align them to a reference, convert the sam file to bam, sort the bam file, and then collect alignment statistics.
```
runAlignConvertSortStats.sh #note: I premade the database so it would not be made again and again
###############################################################################
#!/bin/bash

#You must provide the following. Note variable DBDIR does not need a "/" at the end.
# sh runAlignConvertSortStats.sh 28 sequence_1.fastq sequence_2.fastq /work/GIF/remkv6/files genome.fa


PROC=16
R1_FQ="$1"
R2_FQ="$2"
DBDIR="$3"
GENOME="$4"

#Trim with trimgalore!
module load trimgalore
module load py-setuptools/35.0.2-py2-hqrsjff
trim_galore --paired ${R1_FQ} ${R2_FQ}


#Align with Hisat2!
#If aligning more than one set of RNAseq, build the database outside the script
module load hisat2
#hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p ${PROC} -x ${GENOME%.*} -1 ${R1_FQ%%.*}_val_1.fq.gz -2 ${R2_FQ%%.*}_val_2.fq.gz -S ${R1_FQ%.*}.sam


#Convert sam to bam, and sort!
module load samtools
samtools view --threads ${PROC} -b -o ${R1_FQ%.*}.bam ${R1_FQ%.*}.sam
samtools sort -m 7G -o ${R1_FQ%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads ${PROC} ${R1_FQ%.*}.bam


#Check alignment quality with Picard!
module load GIF2/java/1.8.0_25-b17
module load picard/2.17.0-ft5qztz
java -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/picard-2.17.0-ft5qztzntoymuxiqt3b6yi6uqcmgzmds/bin/picard.jar CollectAlignmentSummaryMetrics  REFERENCE_SEQUENCE=${DBDIR}/${GENOME} INPUT=${R1_FQ%.*}_sorted.bam OUTPUT=${R1_FQ%.*}.bam_alignment.stats
```

Set up the files for running with the runAlignConvertSortStats.sh script.
```
paste <(ls -1 *1.fastq.gz ) <(ls -1 *2.fastq.gz) |awk '{print "sh runAlignConvertSortStats.sh",$1,$2,"/work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa"}' >align.sh

#align.sh
################################################################################
sh runAlignConvertSortStats.sh SRR6230579_1.fastq.gz SRR6230579_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
sh runAlignConvertSortStats.sh SRR6230580_1.fastq.gz SRR6230580_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
sh runAlignConvertSortStats.sh SRR6230581_1.fastq.gz SRR6230581_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
sh runAlignConvertSortStats.sh SRR6230582_1.fastq.gz SRR6230582_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
sh runAlignConvertSortStats.sh SRR6230583_1.fastq.gz SRR6230583_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
sh runAlignConvertSortStats.sh SRR6230584_1.fastq.gz SRR6230584_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
sh runAlignConvertSortStats.sh SRR6230585_1.fastq.gz SRR6230585_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
sh runAlignConvertSortStats.sh SRR6230586_1.fastq.gz SRR6230586_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
sh runAlignConvertSortStats.sh SRR6230587_1.fastq.gz SRR6230587_2.fastq.gz /work/GIF/remkv6/Baum/03_GlandRepeat738/02_AlignmentWholeWorm genome738sl.polished.mitoFixed.fa
################################################################################
```


## Align Transcripts/EST/Isoseq to genome
Alignment of transcripts, ESTs, etc. is another common procedure, so creating a reusable script makes this process much quicker.

This script "runGmap.sh" will create a gmap database, align your sequences to a reference with gmap, create a sam file, convert the sam file to a bam file, and then sort the bam.


```
#!/bin/bash

#Makes a database and  aligns your sequences into reference.
#sh runGmap.sh <database name> <folder of database file ending with a "/"> <Fasta file> <query file>

#examples
#sh run_gmap.sh red_abalone_02Jun2017_5fUJu /work/GIF/remkv6/Serb/03_DavideGMAP/ red_abalone_02Jun2017_5fUJu.fasta DavideQuerydna.fasta
#sh run_gmap.sh  m.yessoensisGenome /work/GIF/remkv6/Serb/03_DavideGMAP/ DavideQuerydna.fasta
#sh run_gmap.sh Crassostreagigasgenome /work/GIF/remkv6/Serb/03_DavideGMAP/ Crassostreagigasgenome.fa DavideQuerydna.fasta


module load gmap-gsnap/2018-03-25-qa3kh3t
module load samtools
dbname=$1
dbloc=$2
dbfasta=$3
query=$4

#Build the database!
#build the gmap database out of the script if multiple alignments are to be done
#gmap_build -d $dbname  -D $dbloc $dbfasta

#Align the transcripts!
gmap -D $dbloc -d $dbname -B 5 -t 16  --input-buffer-size=1000000 --output-buffer-size=1000000 -f samse  $query >${dbname%.*}.${query%.*}.sam

#Convert the sam to bam
samtools view --threads 16 -b -o ${dbname%.*}.${query%.*}.bam ${dbname%.*}.${query%.*}.sam

#sort the bam file
samtools sort -m 7G -o ${query%.*}_sorted.bam -T ${R1_FQ%.*}_temp --threads 16 ${query%.*}.bam
```


Set up the files for running with runGmap.sh.  
```
sh runGmap.sh genome /work/GIF/remkv6/USDA/19_braker/ genome738sl.polished.mitoFixed.fa consensus_isoforms.fasta
sh runGmap.sh genome /work/GIF/remkv6/USDA/19_braker/ genome738sl.polished.mitoFixed.fa H.glycinesEST.fasta
```

Create the gmap database, as no need to create it twice.
```
gmap_build  -d genome genome738sl.polished.mitoFixed.fa
```

## Run Braker

Just two prerequisites before running braker

1. Get the geneMark key to your home directory.  For all not-for profit institutions, a geneMark key can be downloaded [here](http://exon.gatech.edu/GeneMark/license_download.cgi).  YES, even if your system already has GeneMark installed, you will need the key in your home directory.  
```
mv gm_key_64 ~/.gm_key
```
2. Make sure your augustus/config/species folder is writeable.  
```
#Here is how to check it for my augustus installation
ls -ld /work/GIF/software/programs/augustus/3.3.2/config/species/
drwxrwsrwx. 111 arnstrm its-hpc-condo-severin 111 Oct 17 17:02 /work/GIF/software/programs/augustus/3.3.2/config/species/
```
If you cannot get write permissions to this folder, then you will need to copy the entire augustus folder installation to your current directory.  There you will have to specify in braker where the augusutus config file is in your braker command (--AUGUSTUS_CONFIG_PATH=/path/to/config/folder)

Braker2 can use both transcriptional data and protein alignments to predict genes.  However, only intron/exon borders are improved using the proteins.
```
module load GIF/braker/2.1.0
braker.pl --species=Hglycines --genome=genome738sl.polished.mitoFixed.fa.masked --bam=genome.consensus_isoformsSorted.bam,genome.H.glycinesEST_sorted.bam,SRR6230579_1.fastq_sorted.bam,SRR6230580_1.fastq_sorted.bam,SRR6230581_1.fastq_sorted.bam,SRR6230582_1.fastq_sorted.bam,SRR6230583_1.fastq_sorted.bam,SRR6230584_1.fastq_sorted.bam,SRR6230585_1.fastq_sorted.bam,SRR6230586_1.fastq_sorted.bam,SRR6230587_1.fastq_sorted.bam --prot_seq=H.glycinesProtsNCBI.fasta --prg=exonerate

```


## Troubleshooting Braker errors
If braker fails for some reason, there are multiple places to look for errors

```
#The directory I ran braker
/work/GIF/remkv6/USDA/19_braker

#stdout
braker_0.e384359
#stderr
braker_0.o384359
#braker.log file
/work/GIF/remkv6/USDA/19_braker/braker/Hglycines/braker.log
#genemark stdout
/work/GIF/remkv6/USDA/19_braker/braker/Hglycines/GeneMark-ET.stdout
#genemark stderr
/work/GIF/remkv6/USDA/19_braker/braker/Hglycines/errors/filterGenemark.stderr
#exonerate log
/work/GIF/remkv6/USDA/19_braker/braker/Hglycines/startAlign_exonerate.log
```

You will inevitably see this error if you ever have to rerun braker in the same directory.  However, the error script provides the solution. Add "--useexisting" to your braker script
```
ERROR: in file /work/GIF/software/programs/braker/2.1.0/braker.pl at line 2858
/work/GIF/software/programs/augustus/3.3.2/config/species/Hglycines already exists. Choose another species name, delete this directory or use the existing species with the option --useexisting. Be aware that existing parameters will then be overwritten during training.
```

If braker ever fails due to time constraints or in a later braker analysis, you can just resubmit your same script (perhaps with --useexisting) and braker will start where it left off.

---
[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)
