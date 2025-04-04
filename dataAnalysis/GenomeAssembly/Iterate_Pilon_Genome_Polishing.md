---
title: Multiple rounds of Pilon polishing/gap-filling on a genome assembly
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# How to iterate multiple rounds of Pilon polishing/gap-filling on a genome assembly

After assembling a genome, polishing is usually done to correct small assembly errors. Typically, multiple rounds of polishing is not recommended, due to declining BUSCO scores. However, in certain situations, multiple rounds of polishing may be warranted.  In this case, my reads were obtained from thousands of individuals (population assembly). Exacerbating the problem, soybean cyst nematodes do not like homozygosity, so these individuals (and reads) were fairly diverse. The only time I consider multiple rounds of polishing acceptable is when this considerable diversity in reads exists.   

### Prerequisite software

* Hisat2 --  a quick aligner
* Samtools -- who doesnt need to convert sam to bam
* Pilon -- my polisher of choice due to its overwhelmingly informative output.  


### Setting up your directory

My working directory
```
/work/gif/remkv6/USDA/03_LoopingPolisher
```
Softlink your genome
```
ln -s /work/gif/remkv6/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta
```

Softlink the reads
```
ln -s /work/gif/archiveNova/lane2/CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz
ln -s /work/gif/archiveNova/lane2/CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz
```

### Create a run script for pilon
```
#runPilon.sh
###############################################################################
#!/bin/bash
```
Create unix variables
```
DIR="$1"
GENOME="$2"
R1_FQ="$3"
R2_FQ="$4"
```

Index the genome and perform short read mapping using Hisat2
```
module load hisat2
hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p 36 -x ${GENOME%.*} -1 $R1_FQ -2 $R2_FQ -S ${GENOME%.*}.${R1_FQ%.*}.sam
```
Convert your sam to bam, sort, and index.
Note that I use only a part of my available threads to ensure I do not run into RAM usage issues.  This is the safest approach when iterating.
```
module load samtools
samtools view --threads 12 -b -o ${GENOME%.*}.${R1_FQ%.*}.bam ${GENOME%.*}.${R1_FQ%.*}.sam
mkdir Samtemp
samtools sort  -o ${GENOME%.*}.${R1_FQ%.*}_sorted.bam -T Samtemp --threads 16 ${GENOME%.*}.${R1_FQ%.*}.bam
samtools index ${GENOME%.*}.${R1_FQ%.*}_sorted.bam
```
I am running pilon here using the max memory available to me, using a temp folder I created above. Pilon can also have issues with RAM usage, so in this case I cut my thread usage to half of those available. My chunk size is also smaller than default, again catering to potential ram issues that would affect iteration
```
module load pilon/1.22-s7zrot6
java -Xmx200g -Djava.io.tmpdir=Samtemp -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/pilon-1.22-s7zrot6o5yqjh6oxpdxsxcdiswpjioyy/bin/pilon-1.22.jar  --genome ${GENOME} --frags  ${GENOME%.*}.${R1_FQ%.*}_sorted.bam --output ${GENOME%.*}.pilon --outdir ${DIR} --changes --fix all --threads 18  --chunksize 30000
```

The mapping files can get out of hand and use all of your storage, so when I have the loop working, I delete the sam and bam files after each round.
```
rm *sam
rm *bam
```
##### runPilon without comments
<details>
  <summary>Click to see content</summary>
  <pre>
#!/bin/bash
DIR="$1"
GENOME="$2"
R1_FQ="$3"
R2_FQ="$4"
module load hisat2
hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p 36 -x ${GENOME%.*} -1 $R1_FQ -2 $R2_FQ -S ${GENOME%.*}.${R1_FQ%.*}.sam

module load samtools
samtools view --threads 12 -b -o ${GENOME%.*}.${R1_FQ%.*}.bam ${GENOME%.*}.${R1_FQ%.*}.sam
mkdir Samtemp
samtools sort  -o ${GENOME%.*}.${R1_FQ%.*}_sorted.bam -T Samtemp --threads 16 ${GENOME%.*}.${R1_FQ%.*}.bam
samtools index ${GENOME%.*}.${R1_FQ%.*}_sorted.bam

module load pilon/1.22-s7zrot6
java -Xmx200g -Djava.io.tmpdir=Samtemp -jar /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/pilon-1.22-s7zrot6o5yqjh6oxpdxsxcdiswpjioyy/bin/pilon-1.22.jar  --genome ${GENOME} --frags  ${GENOME%.*}.${R1_FQ%.*}_sorted.bam --output ${GENOME%.*}.pilon --outdir ${DIR} --changes --fix all --threads 18  --chunksize 30000

  rm *sam
  rm *bam
  </pre>
  </details>

### How to iterate the runPilon.sh script
Typically I would run one round of polishing like this
```
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz
```

However, to iterate we need to put this into a loop. Lets run 20 rounds of polishing
```
for f in {01..20}; do echo "sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome"${f}".fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome"${f}".pilon.fasta" ;done |paste - <(for f in {02..21}; do echo $line" SCNgenome"${f}".fasta";done) |sed 's/\t/ /g' >PolishLoop.sh


```

###### PolishLoop.sh
<details>
  <summary>Click to see content</summary>
  <pre>

sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome01.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome01.pilon.fasta  SCNgenome02.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome02.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome02.pilon.fasta  SCNgenome03.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome03.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome03.pilon.fasta  SCNgenome04.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome04.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome04.pilon.fasta  SCNgenome05.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome05.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome05.pilon.fasta  SCNgenome06.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome06.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome06.pilon.fasta  SCNgenome07.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome07.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome07.pilon.fasta  SCNgenome08.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome08.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome08.pilon.fasta  SCNgenome09.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome09.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome09.pilon.fasta  SCNgenome10.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome10.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome10.pilon.fasta  SCNgenome11.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome11.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome11.pilon.fasta  SCNgenome12.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome12.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome12.pilon.fasta  SCNgenome13.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome13.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome13.pilon.fasta  SCNgenome14.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome14.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome14.pilon.fasta  SCNgenome15.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome15.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome15.pilon.fasta  SCNgenome16.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome16.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome16.pilon.fasta  SCNgenome17.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome17.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome17.pilon.fasta  SCNgenome18.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome18.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome18.pilon.fasta  SCNgenome19.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome19.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome19.pilon.fasta  SCNgenome20.fasta
sh runPilon.sh /work/gif/remkv6/USDA/03_LoopingPolisher SCNgenome20.fasta CrazyWackyScienceDNA-1_S1_L002_R1_001.fastq.gz CrazyWackyScienceDNA-1_S1_L002_R2_001.fastq.gz; mv SCNgenome20.pilon.fasta  SCNgenome21.fasta

</pre>
</details>

### Evaluate the results of polishing the genome in a loop
Since we set the --changes flag in our Pilon runs, we know the number of changes made in each round
```
43577 SCNgenome01.pilon.changes
 4396 SCNgenome02.pilon.changes
 1241 SCNgenome03.pilon.changes
  496 SCNgenome04.pilon.changes
  308 SCNgenome05.pilon.changes
  186 SCNgenome06.pilon.changes
  105 SCNgenome07.pilon.changes
  117 SCNgenome08.pilon.changes
   74 SCNgenome09.pilon.changes
   41 SCNgenome10.pilon.changes
   36 SCNgenome11.pilon.changes
   35 SCNgenome12.pilon.changes
   30 SCNgenome13.pilon.changes
   31 SCNgenome14.pilon.changes
   25 SCNgenome15.pilon.changes
   21 SCNgenome16.pilon.changes
   23 SCNgenome17.pilon.changes
   21 SCNgenome18.pilon.changes
   22 SCNgenome19.pilon.changes
   36 SCNgenome20.pilon.changes

```
How much has the size of the genome changed?
About 600kb of sequence was removed from the genome after 20 rounds.
```
#using new_Assemblathon.pl script to get info on assemblies.
for f in *fasta; do ~/common_scripts/new_Assemblathon.pl $f;done >Assemblathons.txt
paste <(for f in *fasta; do ls $f;done ) <(grep "Total size of scaffolds" Assemblathons.txt) |sed 's/   //g' |less

SCNgenome01.fasta        Total size of scaffolds  157982452
SCNgenome02.fasta        Total size of scaffolds  157667082
SCNgenome03.fasta        Total size of scaffolds  157577436
SCNgenome04.fasta        Total size of scaffolds  157526821
SCNgenome05.fasta        Total size of scaffolds  157495502
SCNgenome06.fasta        Total size of scaffolds  157456603
SCNgenome07.fasta        Total size of scaffolds  157444744
SCNgenome08.fasta        Total size of scaffolds  157434711
SCNgenome09.fasta        Total size of scaffolds  157427429
SCNgenome10.fasta        Total size of scaffolds  157427423
SCNgenome11.fasta        Total size of scaffolds  157426682
SCNgenome12.fasta        Total size of scaffolds  157426677
SCNgenome13.fasta        Total size of scaffolds  157425936
SCNgenome14.fasta        Total size of scaffolds  157422293
SCNgenome15.fasta        Total size of scaffolds  157422298
SCNgenome16.fasta        Total size of scaffolds  157422294
SCNgenome17.fasta        Total size of scaffolds  157422298
SCNgenome18.fasta        Total size of scaffolds  157421197
SCNgenome19.fasta        Total size of scaffolds  157420463
SCNgenome20.fasta        Total size of scaffolds  157407943
SCNgenome21.fasta        Total size of scaffolds  157407206
```
Has the total N content diminished with each round?
```
paste <(ls -l *fasta) <(for f in *fasta; do sed 's/N/\nN/g'  <$f |grep -c "N" ;done)
SCNgenome01.fasta  1681753
SCNgenome02.fasta  1674378
SCNgenome03.fasta  1670502
SCNgenome04.fasta  1669372
SCNgenome05.fasta  1668709
SCNgenome06.fasta  1668603
SCNgenome07.fasta  1668416
SCNgenome08.fasta  1668416
SCNgenome09.fasta  1668367
SCNgenome10.fasta  1668055
SCNgenome11.fasta  1668008
SCNgenome12.fasta  1667816
SCNgenome13.fasta  1667674
SCNgenome14.fasta  1667674
SCNgenome15.fasta  1667620
SCNgenome16.fasta  1667620
SCNgenome17.fasta  1667530
SCNgenome18.fasta  1666495
SCNgenome19.fasta  1666495
SCNgenome20.fasta  1666303
SCNgenome21.fasta  1666303

```
### Options to convert to gap filling or other types of Pilon fixes
By making a small change in the runPilon.sh script, you can also create a looping gap filler that will continually build sequence at gap edges until they are complete (in theory). Other fun options are available.

```
--fix fixlist
   A comma-separated list of categories of issues to try to fix:
     "snps": try to fix individual base errors;
     "indels": try to fix small indels;
     "gaps": try to fill gaps;
     "local": try to detect and fix local misassemblies;
     "all": all of the above (default);
     "bases": shorthand for "snps" and "indels" (for back compatibility);
     "none": none of the above; new fasta file will not be written.
   The following are experimental fix types:
     "amb": fix ambiguous bases in fasta output (to most likely alternative);
     "breaks": allow local reassembly to open new gaps (with "local");
     "circles": try to close circlar elements when used with long corrected reads;
     "novel": assemble novel sequence from unaligned non-jump reads.
```

### References

* https://github.com/broadinstitute/pilon/wiki


[Back to the Assembly and Annotation Index page](../GenomeAnnotation/annotation_and_assembly_index.md)
