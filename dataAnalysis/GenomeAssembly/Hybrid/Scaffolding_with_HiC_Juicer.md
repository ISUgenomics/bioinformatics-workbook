---
title: How to analyze Hi-C data with Juicer, scaffold your assembly with HiC using 3D-DNA
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

#  Assess your needs
```
This tutorial was conducted with a 160Mb genome and 1 billion Hi-C reads, using juicer's CPU scripts on an HPC with 16 procs and 128Gb ram. There are more than a few things you will need to get started.  
1.  Hi-C reads in fastq format
2.  Genome
3.  The restriction enzyme used for your HiC data
4.  If you want to run hiccups(optional), you'll need a GPU node
5.  Depending on your genome size and amount of repetitive content, you may want to create a black list to prevent juicer from running forever on the dedup step.  The blacklist will remove reads in these highly repetitive areas from the merged_sort.txt output from juicer. Essentially simple repeats are just evil for this step
```
# Software Dependencies  of this tutorial
```
Most of these are pretty common among HPC for bioinformatics.  I was lucky and didnt have to install anything.
blast
bedtools
samtools
bwa
gnutls/3.5.13
jdk (java development kit) 1.8
bioawk
lastz
python
parallel
```

## Decide if you want a black list to get rid of reiterated simple repeats that kill juicer at the dedup step
If not, move to the next step (Initial setup of juicer scripts).
```
makeblastdb -in MisAssFixed.Pilon.fasta -dbtype nucl -out Genome.DB

blastn -db Genome.DB -dust no -num_threads 16 -outfmt 6 -query CentromereAndTelomereRepeats.fasta -evalue 100 -num_alignments 100000  -out Repeats2Genome.blastout

#this gets the coordinates, respective of strand, removes duplicates, and merges overlapping coordinates
less Centromere2Genome.blastout |awk '{if($10>$9){print $2,$9,$10} else {print $2,$10,$9}}' |tr " " "\t" |sort -k1,1V -k2,3n |uniq >CentromereUnmerge.bed ;bedtools merge -d 1000 -i CentromereUnmerge.bed >Centromeres.bed

Have this ready before you go to the deduplication stage

```

### Initial setup of juicer scripts
```
#my starting directory
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads

git clone https://github.com/theaidenlab/juicer.git
cd juicer/
ln -s CPU/scripts scripts
cd scripts/
wget http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.6_jcuda.0.8.jar
ln -s juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
cd ..
```
### softlink and index your reference
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer

mkdir references
cd references
ln -s ../../../MisAssFixed.Pilon.fasta
module load bwa
bwa index MisAssFixed.Pilon.fasta
cd ..
```

### Predict the fragment sizes from a restriction enzyme digest for your genome.
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer

mkdir restriction_sites
cd restriction_sites/
module load python/2.7.15-ief5zfp
python ../misc/generate_site_positions.py MboI MaskedMisAssFixed.Pilon.fasta /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer/references/MisAssFixed.Pilon.fasta
#created "MisAssFixed.Pilon.fasta_MboI.txt"
cd ..
```

### Softlink all fastq files to fastq folder
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer

mkdir fastq
cd fastq
for f in /work/GIF/remkv6/Baum/01_SCNDovetailScaffolding/02_DovetailFastq/CP4477_hic_hiseq/*gz ; do ln -s $f;done
cd ..
```
### Create a chromosome sizes file
```
module load bioawk
bioawk -c fastx '{print $name"\t"length($seq)}' MisAssFixed.Pilon.fasta >chrom.sizes
```

### Run juicer
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer

module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
bash scripts/juicer.sh  -y restriction_sites/MaskedMisAssFixed.Pilon.fasta_MboI.txt -z references/MaskedMisAssFixed.Pilon.fasta -p chrom.sizes
```
For me, this ran for about 2 days before juicer started generating a merged_nodups.txt in the aligned folder. This is when I KILL THE PROCESS, as it takes more than 4 days to deduplicate my billion reads.  
Adding a -S dedup allows the juicer script to restart the deduplication process, but all previous progress is erased.

I honestly could never get the SLURM scripts to submit to my system, and thus was forced to work with the cpu scripts.  Here is my workaround to get the dedup to run in parallel, and enable submission to multiple nodes.

### Dedup workaround
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/02_juicerUnmasked/aligned
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
# sample: contents of x01merged_sort.txt.dir
##############################################################################################
 merged_sort.txt
##############################################################################################
# This takes every directory that was made from the above script, moves to that directory, and creates execution script for juicer to run in parallel.  Dont forget -S dedup.
for f in x*dir; do echo "cd "$f"; bash scripts/juicer.sh -S dedup -y restriction_sites/MisAssFixed.Pilon.fasta_MboI.txt -z references/MisAssFixed.Pilon.fasta -p chrom.sizes";done >>dedup.sh

# I submitted these in parallel to four nodes
```
### Create your merged_nodups.txt file and generate your .hic and .assembly files
```
#/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer/aligned
#This took hours to concatenate.
cat x*dir/aligned/merged_nodups.txt > merged_nodups.txt

# generate a hic file manually with the deduplicated reads
module load bwa
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
java -Xmx2g -jar ../scripts/juicer_tools.jar pre merged_nodups.txt merged_nodups.hic ../chrom.sizes
```


# Run the 3D DNA pipeline to generate a scaffolded genome assembly that can be manipulated in juicebox
```
/work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/01_JuiceBox
git clone https://github.com/theaidenlab/3d-dna.git
cd 3d-dna/

module load GIF2/lastz/1.03.73
module load gnutls/3.5.13-7a3mvfy
#must be 1.8 jdk
module load jdk/8u172-b11-rnauqmr
module load python
module load parallel/20170322-36gxsog
bash run-asm-pipeline.sh -m haploid /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer/references/MisAssFixed.Pilon.fasta /work/GIF/remkv6/Baum/04_Dovetail2Restart/04_GapFilling/09_JuicerScaff/04_scnHicReads/juicer/aligned/merged_nodups.txt

runs in abot 8-9hrs 16cpu.  140gb bam merged_nodups.txt file.
```
