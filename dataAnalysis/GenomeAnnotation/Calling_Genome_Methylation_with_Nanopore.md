---
title: Nanopore Methylation with Nanopolish and Promoter Assessments
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


# A pipeline for assessing methylation in a genome using nanopore sequencing

DNA methylation plays a role in chromatin compaction and gene expression. Finding which genes occupy methylated regions of the genome, may provide some insight into gene regulation. Here we will use nanopolish to identify methylation sites in a genome with nanopore reads.   

### Prerequisite software and data
* samtools -- standard install on most HPCs
* minimap2 -- a quick aligner for long reads
* nanopolish -- I use a conda environment below to install.
* nanopore fast5 (raw) and fastq data

### Install nanopolish
```
/work/gif/remkv6/USDA/01_NanoporeMethylation

git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish/
make
#python dependencies to use programs from the scripts folder
pip install -r scripts/requirements.txt --user
export PATH="/work/gif/remkv6/USDA/01_NanoporeMethylation/nanopolish:$PATH"
```

### Create file structure requested by Nanopolish

The structure is to have two separate folders.  One for your fastq reads and one for the corresponding fast5 data.  
```
/work/gif/remkv6/USDA/01_NanoporeMethylation

mkdir ooze_1_1fast5   ;cd ooze_1_1fast5  ; for f in /work/gif/archiveNova/03_RevoltingBlobMethylationNanopore/1577/1577/1577/20201119_1838_X1_FAO51552_10ff60f0/fast5_pass/barcode03/* ; do ln -s $f ; done ;cd ../
mkdir ooze_1_1fastq   ;cd ooze_1_1fastq  ; for f in /work/gif/archiveNova/03_RevoltingBlobMethylationNanopore/1577/1577/1577/20201119_1838_X1_FAO51552_10ff60f0/fastq_pass/barcode03/* ; do ln -s $f ; done ;cd ../
mkdir ooze_1_2fast5   ;cd ooze_1_2fast5  ; for f in /work/gif/archiveNova/03_RevoltingBlobMethylationNanopore/1577/1577/1577_2/20201123_1600_X1_FAO68129_332eba1a/fast5_pass/barcode03/*  ; do ln -s $f ; done ;cd ../
mkdir ooze_1_2fastq   ;cd ooze_1_2fastq  ; for f in /work/gif/archiveNova/03_RevoltingBlobMethylationNanopore/1577/1577/1577_2/20201123_1600_X1_FAO68129_332eba1a/fastq_pass/barcode03/*  ; do ln -s $f ; done ;cd ../
mkdir ooze_2_1fast5   ;cd ooze_2_1fast5  ; for f in /work/gif/archiveNova/03_RevoltingBlobMethylationNanopore/1577/1577/1577/20201119_1838_X1_FAO51552_10ff60f0/fast5_pass/barcode04/* ; do ln -s $f ; done ;cd ../
mkdir ooze_2_1fastq   ;cd ooze_2_1fastq  ; for f in /work/gif/archiveNova/03_RevoltingBlobMethylationNanopore/1577/1577/1577/20201119_1838_X1_FAO51552_10ff60f0/fastq_pass/barcode04/* ; do ln -s $f ; done ;cd ../
mkdir ooze_2_2fast5   ;cd ooze_2_2fast5  ; for f in /work/gif/archiveNova/03_RevoltingBlobMethylationNanopore/1577/1577/1577_2/20201123_1600_X1_FAO68129_332eba1a/fast5_pass/barcode04/*  ; do ln -s $f ; done ;cd ../
mkdir ooze_2_2fastq   ;cd ooze_2_2fastq  ; for f in /work/gif/archiveNova/03_RevoltingBlobMethylationNanopore/1577/1577/1577_2/20201123_1600_X1_FAO68129_332eba1a/fastq_pass/barcode04/*  ; do ln -s $f ; done ;cd ../
```
Once the fast5 and fastq data is softlinked, you need to concatenate all of your fastq data for each sample
```
cd ooze_1_1fastq; cat *fastq >ooze1_1.fastq; cd ../ooze_1_2fastq; cat *fastq >ooze1_2.fastq ; cd ooze_2_1fastq; cat *fastq >ooze2_1.fastq; cd ../ooze_2_2fastq; cat *fastq >ooze2_2.fastq
```

### Index with nanopolish
```
/work/gif/remkv6/USDA/01_NanoporeMethylation
ml miniconda3; source activate nanopolish
nanopolish index -d ooze_1_1fast5/ ooze_1_1fastq/ooze1_1.fastq &
nanopolish index -d ooze_1_2fast5/ ooze_1_2fastq/ooze1_2.fastq &
nanopolish index -d ooze_2_1fast5/ ooze_2_1fastq/ooze2_1.fastq &
nanopolish index -d ooze_2_2fast5/ ooze_2_2fastq/ooze2_2.fastq &
```

### Align nanopore reads from each sample to your genome
My working directory
```
/work/gif/remkv6/USDA/01_NanoporeMethylation
```
Softlink your genome
```
ln -s /work/gif/remkv6/04_DovetailRevoltingBlobGenome/49_RenameChromosomes/01_Transfer2Box/RevoltingBlobgenome.fasta
```
Create the temp files needed to sort after alignment
```
for f in *fastq/*fastq; do mkdir ${f%.*}"tmp";done
```
Execute alignment of nanopore reads, sort and index
```
for f in *fastq; do echo "ml minimap2;ml samtools;minimap2 -a -x map-ont RevoltingBlobgenome.fasta "$f"| samtools sort -T "${f%.*}"tmp -o "${f%.*}".sorted.bam ;samtools index "${f%.*}".sorted.bam";done >align.sh
```

### Alignment stats

If alignment is poor, you will not get a good assessment of the methylation changes in your sample. Mine was not great, I should tweak the parameters to get a better alignment.
```
/work/gif/remkv6/USDA/01_NanoporeMethylation

for f in *fastq/*bam; do samtools flagstat  $f >>alignmentStats.txt;done
paste <(grep "total" alignmentStats.txt) <(grep "%" alignmentStats.txt) |less
```
Yes this is a poor alignment, but for demonstration purposes, I will continue.  I'd prefer this % to be near 70 at least.
```
450024 + 0 in total (QC-passed reads + QC-failed reads) 214902 + 0 mapped (47.75% : N/A)
276447 + 0 in total (QC-passed reads + QC-failed reads) 126140 + 0 mapped (45.63% : N/A)
488875 + 0 in total (QC-passed reads + QC-failed reads) 267888 + 0 mapped (54.80% : N/A)
623642 + 0 in total (QC-passed reads + QC-failed reads) 326584 + 0 mapped (52.37% : N/A)
```

### Assess different types of methylation for each sample
I create a script below that will assess all detectable types of methylation possible with nanopolish.     

A list the script below uses.
```
vi types
###########
cpg
gpc
###########
```
Create scripts to create the run scripts for all different types of methylation.  
```
less types |while read line; do echo "for f in *bam; do echo \"ml miniconda3; source activate nanopolish;  export PATH=\\\"work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:\\\$PATH\\\" ;nanopolish call-methylation -q "$line" -r \"\${f%%.*}\".fastq -b \"\$f\" -t 36 --progress -g RevoltingBlobgenome.fasta >\"\${f%..*}\"_methylation_"${line%.*}".tsv \";done >>AllMethylationScripts.sh" ;done |less
```
AllMethylationScripts.sh
<details>
  <summary>Click to see content</summary>
  <pre>

for f in *fastq/*bam; do echo "ml miniconda3; source activate nanopolish;  export PATH=\"work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:\$PATH\" ;nanopolish call-methylation -q cpg -r "${f%%.*}".fastq -b "$f" -t 36 --progress -g RevoltingBlobgenome.fasta >"${f%..*}"_methylation_cpg.tsv ";done >>AllMethylationScripts.sh
for f in *fastq/*bam; do echo "ml miniconda3; source activate nanopolish;  export PATH=\"work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:\$PATH\" ;nanopolish call-methylation -q gpc -r "${f%%.*}".fastq -b "$f" -t 36 --progress -g RevoltingBlobgenome.fasta >"${f%..*}"_methylation_gpc.tsv ";done >>AllMethylationScripts.sh

#######################################################################################################################################################

The scripts above generate these.
###################################################################################################################################################################
ml miniconda3; source activate nanopolish;  export PATH="work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:$PATH" ;nanopolish call-methylation -q cpg -r ooze_1_1fastq/ooze1_1.fastq -b ooze_1_1fastq/ooze1_1.sorted.bam -t 36 --progress -g RevoltingBlobgenome.fasta >ooze_1_1fastq/ooze1_1.sorted.bam_methylation_cpg.tsv
ml miniconda3; source activate nanopolish;  export PATH="work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:$PATH" ;nanopolish call-methylation -q cpg -r ooze_1_2fastq/ooze1_2.fastq -b ooze_1_2fastq/ooze1_2.sorted.bam -t 36 --progress -g RevoltingBlobgenome.fasta >ooze_1_2fastq/ooze1_2.sorted.bam_methylation_cpg.tsv
ml miniconda3; source activate nanopolish;  export PATH="work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:$PATH" ;nanopolish call-methylation -q cpg -r ooze_2_1fastq/ooze2_1.fastq -b ooze_2_1fastq/ooze2_1.sorted.bam -t 36 --progress -g RevoltingBlobgenome.fasta >ooze_2_1fastq/ooze2_1.sorted.bam_methylation_cpg.tsv
ml miniconda3; source activate nanopolish;  export PATH="work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:$PATH" ;nanopolish call-methylation -q cpg -r ooze_2_2fastq/ooze2_2.fastq -b ooze_2_2fastq/ooze2_2.sorted.bam -t 36 --progress -g RevoltingBlobgenome.fasta >ooze_2_2fastq/ooze2_2.sorted.bam_methylation_cpg.tsv
ml miniconda3; source activate nanopolish;  export PATH="work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:$PATH" ;nanopolish call-methylation -q gpc -r ooze_1_1fastq/ooze1_1.fastq -b ooze_1_1fastq/ooze1_1.sorted.bam -t 36 --progress -g RevoltingBlobgenome.fasta >ooze_1_1fastq/ooze1_1.sorted.bam_methylation_gpc.tsv
ml miniconda3; source activate nanopolish;  export PATH="work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:$PATH" ;nanopolish call-methylation -q gpc -r ooze_1_2fastq/ooze1_2.fastq -b ooze_1_2fastq/ooze1_2.sorted.bam -t 36 --progress -g RevoltingBlobgenome.fasta >ooze_1_2fastq/ooze1_2.sorted.bam_methylation_gpc.tsv
ml miniconda3; source activate nanopolish;  export PATH="work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:$PATH" ;nanopolish call-methylation -q gpc -r ooze_2_1fastq/ooze2_1.fastq -b ooze_2_1fastq/ooze2_1.sorted.bam -t 36 --progress -g RevoltingBlobgenome.fasta >ooze_2_1fastq/ooze2_1.sorted.bam_methylation_gpc.tsv
ml miniconda3; source activate nanopolish;  export PATH="work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/:$PATH" ;nanopolish call-methylation -q gpc -r ooze_2_2fastq/ooze2_2.fastq -b ooze_2_2fastq/ooze2_2.sorted.bam -t 36 --progress -g RevoltingBlobgenome.fasta >ooze_2_2fastq/ooze2_2.sorted.bam_methylation_gpc.tsv
</pre>
</details>


### Get the methylation frequencies
```
/work/gif/remkv6/USDA/01_NanoporeMethylation
for f in *fastq/*tsv;do ln -s $f;done
for f in *tsv; do /work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/01_Alignment/nanopolish/scripts/calculate_methylation_frequency.py $f >${f%.*}_frequency.tsv; done

#some error in the script that this corrects
for f in *_frequency.tsv ; do sed -i 's/num_CG_in_group/num_motifs_in_group/g' $f;done

#makes part of the script below
paste <(ls -1 Ooze1*gpc_frequency.tsv |tr "\n" " " |sed 's/^/<(cat /g' |sed 's/$/)/g' ) <(ls -1 Ooze2*gpc_frequency.tsv |tr "\n" " " |sed 's/^/<(cat /g' |sed 's/$/)/g')    |while read line; do echo "python3 /work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/methylation_example/compare_methylation.py "$line;done |less



#had to manually add the awk command, as it was difficult to code in.
python3 /work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/methylation_example/compare_methylation.py <(cat Ooze1_1.sorted.bam_methylation_gpc_frequency.tsv Ooze1_2.sorted.bam_methylation_gpc_frequency.tsv |awk '{if (NR==1) {print} else if(substr($1,1,1)=="C") {print}else {next}}') <(cat Ooze2_1.sorted.bam_methylation_gpc_frequency.tsv Ooze2_2.sorted.bam_methylation_gpc_frequency.tsv |awk '{if (NR==1) {print} else if(substr($1,1,1)=="C") {print}else {next}}') >gpcCompare.tsv

python3 /work/gif/remkv6/04_DovetailRevoltingBlobGenome/57_NanoporeMethylation/methylation_example/compare_methylation.py <(cat Ooze1_1.sorted.bam_methylation_cpg_frequency.tsv Ooze1_2.sorted.bam_methylation_cpg_frequency.tsv |awk '{if (NR==1) {print} else if(substr($1,1,1)=="C") {print}else {next}}') <(cat Ooze2_1.sorted.bam_methylation_cpg_frequency.tsv Ooze2_2.sorted.bam_methylation_cpg_frequency.tsv |awk '{if (NR==1) {print} else if(substr($1,1,1)=="C") {print}else {next}}') >cpgCompare.tsv


```

### Identify promoter sequences upstream of genes
DNA methylation typically affects the expression of genes by methylating promoter sequences.  The easiest way to get this is with your genome, is to use bedtools flank with your gene prediction gff to assess methylation 2kb upstream of gene (a proxied promoter).  

Softlink my gene prediction/annotation
```
ln -s /work/gif/remkv6/04_DovetailRevoltingBlobGenome/49_RenameChromosomes/01_Transfer2Box/OrderedRevoltingBlobGenePredictions.gff3
```

Get the "contig_name length"
```
ln -s /work/gif/remkv6/04_DovetailRevoltingBlobGenome/49_RenameChromosomes/01_Transfer2Box/RevoltingBlobgenome.fasta
ml bioawk; bioawk -c fastx '{print $name,length($seq)}' RevoltingBlobgenome.fasta >chrom.sizes
```

Make promoter coordinates in a bed format that bedtools will accept.
```
ml bedtools2; bedtools flank -i <(awk '$3=="gene"' OrderedRevoltingBlobGenePredictions.gff3)  -g chrom.sizes -l 2000 -r 0 -s > genes.2kb.promoters.bed
```

### Compare methylation among samples with promoter sequences
This is some testing based on at least 5 aligned reads to a region, and distinguishes differences between samples based on columns 3 and 5, and assesses how many times a promoter is methylated

##### CPG methylation

```
less cpgCompare.tsv |awk 'NR>1 &&$2>5 && $4>5 && $3>$5' |sort -k3,3nr -k5,5n |sed 's/:/\t/1' |sed 's/-/\t/1' |tr "\t" " " |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |bedtools intersect -wo -a - -b genes.2kb.promoters.bed |cut -f 13 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |sort -k1,1nr  >Sample_1CPG.tsv
less cpgCompare.tsv |awk 'NR>1 &&$2>5 && $4>5 && $3<$5' |sort -k3,3nr -k5,5n |sed 's/:/\t/1' |sed 's/-/\t/1' |tr "\t" " " |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |bedtools intersect -wo -a - -b genes.2kb.promoters.bed |cut -f 13 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |sort -k1,1nr  >Sample_2CPG.tsv
```

 <details>
  <summary>Sample_1CPG.tsv Click to see content</summary>
  <pre>
  14 BLOB19558
13 BLOB19559
13 BLOB19560
12 BLOB14134
12 BLOB19557
 8 BLOB19569
 6 BLOB06502
 6 BLOB07680
 6 BLOB19562
 4 BLOB04612
 4 BLOB14593
 3 BLOB20119
 3 BLOB21251
 2 BLOB01050
 2 BLOB01732
 2 BLOB16387
 2 BLOB16388
 2 BLOB19309
 2 BLOB19568
 1 BLOB00857
 1 BLOB02492
 1 BLOB12086
 1 BLOB14594
 1 BLOB14601
 1 BLOB16386
 1 BLOB19561
 1 BLOB20672
 1 BLOB21356

</pre>
</details>

  <details>
    <summary>Sample_2CPG.tsv Click to see content</summary>
    <pre>
    26 BLOB19557
23 BLOB14134
21 BLOB19559
20 BLOB19560
19 BLOB19558
18 BLOB19562
14 BLOB06502
 7 BLOB19569
 6 BLOB20119
 5 BLOB14594
 5 BLOB21251
 4 BLOB01050
 4 BLOB14593
 4 BLOB19568
 3 BLOB00857
 3 BLOB01732
 3 BLOB16388
 3 BLOB20672
 2 BLOB07680
 2 BLOB16386
 2 BLOB16387
 2 BLOB19561
 1 BLOB04612
 1 BLOB09973
 1 BLOB09974
 1 BLOB13689
 1 BLOB14595
 1 BLOB19309
</pre>
</details>

##### GPC methylation
```
less gpcCompare.tsv |awk 'NR>1 &&$2>5 && $4>5 && $3<$5' |sort -k3,3nr -k5,5n |sed 's/:/\t/1' |sed 's/-/\t/1' |tr "\t" " " |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |bedtools intersect -wo -a - -b genes.2kb.promoters.bed |cut -f 13 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |sort -k1,1nr  >Sample_1GPC.tsv

less gpcCompare.tsv |awk 'NR>1 &&$2>5 && $4>5 && $3>$5' |sort -k3,3nr -k5,5n |sed 's/:/\t/1' |sed 's/-/\t/1' |tr "\t" " " |sed 's/ /\t/1' |sed 's/ /\t/1' |sed 's/ /\t/1' |bedtools intersect -wo -a - -b genes.2kb.promoters.bed |cut -f 13 |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq -c |sort -k1,1nr  >Sample_2GPC.tsv
```
<details>
  <summary>Sample_1GPC.tsv Click to see content</summary>
  <pre>
  10 BLOB19569
   7 BLOB19557
   7 BLOB19560
   7 BLOB19562
   6 BLOB04612
   6 BLOB19558
   5 BLOB02447
   5 BLOB14593
   5 BLOB19559
   5 BLOB20119
   3 BLOB01050
   3 BLOB01732
   3 BLOB16386
   3 BLOB19568
   2 BLOB01731
   2 BLOB09973
   2 BLOB09974
   2 BLOB14594
   2 BLOB14595
   2 BLOB16388
   2 BLOB19309
   2 BLOB21251
   2 BLOB22228
   1 BLOB00081
   1 BLOB00082
   1 BLOB06502
   1 BLOB12086
   1 BLOB16387
   1 BLOB21356
</pre>
</details>
  <details>
    <summary>Sample_2GPC.tsv Click to see content</summary>
    <pre>
    27 BLOB19559
27 BLOB19560
26 BLOB19558
20 BLOB19557
18 BLOB19562
11 BLOB04612
10 BLOB20119
 7 BLOB14593
 7 BLOB16386
 7 BLOB21251
 6 BLOB16387
 6 BLOB16388
 5 BLOB01050
 5 BLOB02447
 5 BLOB06502
 5 BLOB13689
 5 BLOB14134
 5 BLOB19569
 3 BLOB01732
 3 BLOB19309
 3 BLOB19561
 3 BLOB19568
 2 BLOB09973
 2 BLOB09974
 2 BLOB14594
 2 BLOB14595
 2 BLOB20672
 1 BLOB02492
 1 BLOB02922
 1 BLOB18593
 1 BLOB18594
 1 BLOB20967
 1 BLOB22377
 1 BLOB22378
</pre>
</details>




### References
* https://github.com/jts/nanopolish
* https://nanopolish.readthedocs.io/en/latest/

[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)
