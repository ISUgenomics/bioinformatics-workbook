---
layout: redirected
sitemap: false
permalink: /dataAnalysis/GenomeAssembly/Hybrid/Scaffolding_with_HiC_Juicer
redirect_to: /dataAnalysis/GenomeAssembly/Hybrid/Juicer_Juicebox_3dDNA_pipeline
---

---
title: How to analyze Hi-C data with Juicer and scaffold your genome assembly using 3D-DNA
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

#  Assess your needs

This tutorial was conducted with a 160Mb genome and 1 billion Hi-C reads, using juicer's CPU scripts on an HPC with 16 procs and 128Gb ram.

Things you will need to get started.  
1.  Hi-C reads in fastq format
2.  Genome
3.  The restriction enzyme used for your HiC data
4.  If you want to run hiccups(optional), you'll need a GPU node
5.  Depending on your genome size and amount of repetitive content, you may want to create a black list to prevent juicer from running forever on the dedup step. The blacklist will remove reads in these highly repetitive areas from the merged_sort.txt output from juicer.

# Software Dependencies
Most of these are pretty common among HPC for bioinformatics.  I was lucky and didnt have to install anything.

```
blast -- I used 2.7.1
bedtools -- I used 2.27.1
samtools -- I used 1.8-6
bwa -- I used 0.7.17
gnutls -- I used 3.5.13
jdk (java development kit)1.7-1.8 -- I used 1.8
bioawk -- I used 1.0-3
lastz -- I used 1.03.73
python -- I used 2.7.15
parallel -- I used 20170322
```

## Decide if you want a black list to get rid of reiterated simple repeats that kill juicer at the dedup step
If not, move to the next step (Initial setup of juicer scripts).
```
makeblastdb -in MysteryGenome.fasta -dbtype nucl -out Genome.DB

blastn -db Genome.DB -dust no -num_threads 16 -outfmt 6 -query CentromereAndTelomereRepeats.fasta -evalue 100 -num_alignments 100000  -out Repeats2Genome.blastout

#this gets the coordinates, respective of strand, removes duplicates, and merges overlapping coordinates
less Centromere2Genome.blastout |awk '{if($10>$9){print $2,$9,$10} else {print $2,$10,$9}}' |tr " " "\t" |sort -k1,1V -k2,3n |uniq >CentromereUnmerge.bed ;bedtools merge -d 1000 -i CentromereUnmerge.bed >Centromeres.bed

Have this ready before you go to the deduplication stage

```

## Initial setup of juicer scripts
```
#my starting directory
#/09_JuicerScaff/04_HicReads

git clone https://github.com/theaidenlab/juicer.git
cd juicer/
ln -s CPU/scripts scripts
cd scripts/
wget http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.6_jcuda.0.8.jar
ln -s juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
cd ..
```
## softlink and index your reference
```
#/09_JuicerScaff/04_HicReads/juicer

mkdir references
cd references
ln -s ../../../MysteryGenome.fasta
module load bwa
bwa index MysteryGenome.fasta
cd ..
```

### Predict the fragment sizes from a restriction enzyme digest for your genome.
```
#09_JuicerScaff/04_HicReads/juicer

mkdir restriction_sites
cd restriction_sites/
module load python/2.7.15-ief5zfp
python ../misc/generate_site_positions.py MboI MaskedMysteryGenome.fasta /09_JuicerScaff/04_HicReads/juicer/references/MysteryGenome.fasta
#created "MysteryGenome.fasta_MboI.txt"
cd ..
```

### Softlink all fastq files to fastq folder
```
#/09_JuicerScaff/04_HicReads/juicer

mkdir fastq
cd fastq
for f in /work/GIF/remkv6/Baum/01_DovetailScaffolding/02_DovetailFastq/CP4477_hic_hiseq/*gz ; do ln -s $f;done
cd ..
```
### Create a chromosome sizes file
```
module load bioawk
bioawk -c fastx '{print $name"\t"length($seq)}' MysteryGenome.fasta >chrom.sizes
```

### Run juicer
```
#/09_JuicerScaff/04_HicReads/juicer

module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
bash scripts/juicer.sh  -y restriction_sites/MaskedMysteryGenome.fasta_MboI.txt -z references/MaskedMysteryGenome.fasta -p chrom.sizes
```
For me, this ran for about 2 days before juicer started generating a merged_nodups.txt in the aligned folder. This is when I KILL THE PROCESS, as it takes more than 4 days to deduplicate my billion reads.  
Adding a -S dedup allows the juicer script to restart the deduplication process, but all previous progress is erased.

Here is my workaround to get the dedup to run in parallel using the cpu scripts, which enable submission to multiple nodes.

### Dedup workaround
```
#/09_JuicerScaff/04_HicReads/02_juicerUnmasked/aligned
#remove the unfinished files from aligned/ folder
rm dups.txt; rm merged_nodups.txt;rm opt_dups.txt

# This step is only relevant if you are creating a blacklist of simple repeats
less merged_sort.txt |awk '{if($3<$7) {print $2"\t"$3"\t"$7+150"\t"$0} else {print $2"\t"$7"\t"$3+150"\t"$0}}'|bedtools intersect -v -a - -b ../../02_InvestigateDupRegions/FilterRepeatscaffold4_5570000-5610000.bed |cut -f 4- >BlacklistedMerged_sort.txt
mv merged_sort.txt Merged_SortOld.txt
mv BlacklistedMerged_sort.txt merged_sort.txt



#split the merged_sort.txt file (I split by 64 to make it faster and run on 4 nodes x 16 processes
split -d  --additional-suffix=merged_sort.txt -n l/64  merged_sort.txt &

#This sets up the native directory structure of juicer, without copying unecessary files
for f in x*txt; do mkdir $f.dir; cd $f.dir; cp ../../* .; ln -s ../../fastq/; ln -s ../../misc/ ;ln -s ../../debug/; ln -s ../../references/;ln -s ../../restriction_sites/;ln -s ../../scripts/;ln -s ../../splits/; mkdir aligned; mv ../$f aligned/merged_sort.txt; cd ../   ;done
#creates 64 folders like this
###############################################################################################
x01merged_sort.txt.dir
x02merged_sort.txt.dir
x03merged_sort.txt.dir
x04merged_sort.txt.dir
x05merged_sort.txt.dir
x06merged_sort.txt.dir
x07merged_sort.txt.dir
x08merged_sort.txt.dir
x09merged_sort.txt.dir
x10merged_sort.txt.dir
... etc
##############################################################################################

# sample: file structure of x01merged_sort.txt.dir
##############################################################################################
LICENSE  README.md  aligned  chrom.sizes  debug  fastq  juiceIt_0.sub  juicer.sh  misc  references  restriction_sites  scripts  splits
##############################################################################################
# sample: contents of x01merged_sort.txt.dir/aligned/
##############################################################################################
 merged_sort.txt
##############################################################################################
# This takes every directory that was made from the above script, moves to that directory, and creates execution script for juicer to run in parallel.  Dont forget -S dedup.
for f in x*dir; do echo "cd "$f"; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes";done >>dedup.sh
#sample output from the above script
###############################################################
cd x00merged_sort.txt.dir; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes
cd x01merged_sort.txt.dir; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes
cd x02merged_sort.txt.dir; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes
cd x03merged_sort.txt.dir; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes
cd x04merged_sort.txt.dir; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes
cd x05merged_sort.txt.dir; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes
cd x06merged_sort.txt.dir; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes
cd x07merged_sort.txt.dir; bash scripts/juicer.sh -S dedup -y restriction_sites/MysteryGenome.fasta_MboI.txt -z references/MysteryGenome.fasta -p chrom.sizes
...etc
###############################################################
# I submitted these in parallel to four nodes
```
# Two choices moving forward
#### Choice 1: Create your merged_nodups.txt file and generate your .hic files
This choice is faster and will give you the assembly files you want more quickly.  However this path will not give you all of the standard juicer output.  
```
################################################################################
#/09_JuicerScaff/04_HicReads/juicer/aligned
#Rename your old merged_sort.txt
mv merged_sort.txt Round1MergedSort.txt
cat x*dir/aligned/merged_nodups.txt |sort --parallel=16 -k2,2d -k6,6d > merged_nodups.txt

#Beware, if you kill juicer early before the .hic file is completed, you may find lines of binary in your merged_nodups.txt files.


# generate a hic file manually with the deduplicated reads
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
java -Xmx2g -jar ../scripts/juicer_tools.jar pre merged_nodups.txt merged_nodups.hic ../chrom.sizes
```
#### Choice 2: Rerun juicer on the merged_nodups.txt reads to completely remove duplicates. (about 24hrs for me)
The program will then generate all the associated file output with this method, but does not dramatically affect the number of duplicates found.

To be more careful with removing all duplicates, you should run juicer -S dedup on your newly merged_nodups.txt file, and rename it to  merged_sort.txt. Then rerun through juicer.
The file is 82GB smaller this time around and only takes a day to run without getting stuck on duplicates
This removes only 115KB of duplicates the second run through.
```
#/09_JuicerScaff/04_HicReads/juicer/aligned
#Rename your old merged_sort.txt
mv merged_sort.txt Round1MergedSort.txt

#Concatenate all split merged_nodups.txt files
cat x*dir/aligned/merged_nodups.txt |sort --parallel=16 -k2,2d -k6,6d > merged_sort.txt

#To be careful to remove all duplicates, you should run juicer -S dedup on your newly generated merged_sort.txt.
This removes only 15KB of duplicates the second run through.
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
bash scripts/juicer.sh  -S dedup -y restriction_sites/MaskedMysteryGenome.fasta_MboI.txt -z references/MaskedMysteryGenome.fasta -p chrom.sizes

```

# Run the 3D DNA pipeline to generate a scaffolded genome assembly that can be manipulated in juicebox
```
/09_JuicerScaff/04_HicReads/01_JuiceBox
git clone https://github.com/theaidenlab/3d-dna.git
cd 3d-dna/

module load GIF2/lastz/1.03.73
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
module load python
module load parallel/20170322-36gxsog
bash run-asm-pipeline.sh -m haploid /09_JuicerScaff/04_HicReads/juicer/references/MysteryGenome.fasta /09_JuicerScaff/04_HicReads/juicer/aligned/merged_nodups.txt

runs in abot 8-9hrs 16cpu.  140gb bam merged_nodups.txt file.
```

## How to use JuiceBox

Download and install Juicebox

[Juicebox Download](https://github.com/aidenlab/Juicebox/wiki/Download)

There are some sources of information via tutorial videos of genome assembly corrections using Juicebox.

1. [Introduction to Juicebox for Hi-C exploration](https://www.youtube.com/watch?v=xjNXyeUSfZM)
2. [JuiceBox tutorial for genome assembly correction](
https://www.youtube.com/watch?v=Nj7RhQZHM18)
3. [A helpful lightning speed correction of the barrel medic genome](https://www.youtube.com/watch?v=IMmVp8FodmY)

### Uploading files to Juicebox
Juicer uses the .assembly and .hic files generated from 3D-DNA and juicer/3D-DNA.
![MysteryGenome loaded into JuiceBox](https://isugenomics.github.io/bioinformatics-workbook/assets/JuiceBox.png)

This is what a genome assembly looks like from a population of highly heterozygous individuals.  It is not nearly as neat as the barrel medic genome, but lots of scaffolding opportunities are possible.

## Sources

* [Host website](https://github.com/aidenlab/juicer) -- Has all of the software and tutorials hosted here.
* [Aiden lab forum](http://aidenlab.org/forum.html) -- Very helpful if you have a problem with their software.
* [Original publication of mosquito genome assembly using Hi-C](http://science.sciencemag.org/content/356/6333/92)

---

* [Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
