---
title: RSeQC tutorial
layout: single
author: Rick Masonbrink
author_profile: true
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
fastq-dump --outdir 01_Align/  --split-files SRR1573504 &
fastq-dump --outdir 01_Align/  --split-files SRR1573505 &

#stranded samples
fastq-dump --outdir 01_Align/  --split-files SRR13332812&
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
<details>
  <summary>Maizealign.sh content</summary>
  <pre>

```
sh runHiSat2bam.sh SRR13332812_1.fastq SRR13332812_2.fastq /work/gif/TranscriptomicsWorkshop/remkv6/01_Maize/01_Align Maize.fa
sh runHiSat2bam.sh SRR13332813_1.fastq SRR13332813_2.fastq /work/gif/TranscriptomicsWorkshop/remkv6/01_Maize/01_Align Maize.fa
sh runHiSat2bam.sh SRR1573504_1.fastq SRR1573504_2.fastq /work/gif/TranscriptomicsWorkshop/remkv6/01_Maize/01_Align Maize.fa
sh runHiSat2bam.sh SRR1573505_1.fastq SRR1573505_2.fastq /work/gif/TranscriptomicsWorkshop/remkv6/01_Maize/01_Align Maize.fa
```
</pre>
</details>

### Create bed12 file from genic gff
```
#this converts your gff to bed12 format
ml bedops
gff2bed <GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.gff >maizeGFF.bed12
```
<details>
  <summary>maizeGFF.bed12 content</summary>
  <pre>

```
NC_001666.2     0       140384  NC_001666.2:1..140384   .       +       RefSeq  region  .       ID=NC_001666.2:1..140384;Dbxref=taxon:4577;Is_circular=true;Name=Pltd;gbkey=Src;genome=chloroplast;mol_type=genomic DNA
NC_001666.2     0       140384  id-NC_001666.2:1..140384        .       +       RefSeq  sequence_feature        .       ID=id-NC_001666.2:1..140384;Note=[JLA]%3B junction IRA-LSC (circular molecule);gbkey=misc_feature
NC_001666.2     88      1150    cds-NP_043004.1 .       -       RefSeq  CDS     0       ID=cds-NP_043004.1;Parent=gene-ZemaCp002;Dbxref=UniProtKB/Swiss-Prot:P48183,Genbank:NP_043004.1,GeneID:845199;Name=NP_043004.1;gbkey=CDS;gene=psbA;locus_tag=ZemaCp002;product=photosystem II protein D1;protein_id=NP_043004.1;transl_table=11
NC_001666.2     88      1150    gene-ZemaCp002  .       -       RefSeq  gene    .       ID=gene-ZemaCp002;Dbxref=GeneID:845199;Name=psbA;gbkey=Gene;gene=psbA;gene_biotype=protein_coding;locus_tag=ZemaCp002
NC_001666.2     1385    1420    exon-ZemaCt121-2        .       -       RefSeq  exon    .       ID=exon-ZemaCt121-2;Parent=rna-ZemaCt121;Dbxref=GeneID:845256;gbkey=tRNA;gene=trnK;locus_tag=ZemaCt121;product=tRNA-Lys
NC_001666.2     1385    1420    id-ZemaCt121-2  .       -       RefSeq  exon    .       ID=id-ZemaCt121-2;Parent=gene-ZemaCt121;Dbxref=GeneID:845256;exon_number=2;gbkey=exon;gene=trnK;locus_tag=ZemaCt121;number=2
NC_001666.2     1385    3946    gene-ZemaCt121  .       -       RefSeq  gene    .       ID=gene-ZemaCt121;Dbxref=GeneID:845256;Name=trnK;gbkey=Gene;gene=trnK;gene_biotype=tRNA;locus_tag=ZemaCt121
NC_001666.2     1385    3946    rna-ZemaCt121   .       -       RefSeq  tRNA    .       ID=rna-ZemaCt121;Parent=gene-ZemaCt121;Dbxref=GeneID:845256;gbkey=tRNA;gene=trnK;locus_tag=ZemaCt121;product=tRNA-Lys
NC_001666.2     1420    3909    id-ZemaCt121    .       -       RefSeq  intron  .       ID=id-ZemaCt121;Parent=gene-ZemaCt121;Dbxref=GeneID:845256;exon_number=1;gbkey=intron;gene=trnK;locus_tag=ZemaCt121;number=1
NC_001666.2     1673    3215    cds-NP_043005.2 .       -       RefSeq  CDS     0       ID=cds-NP_043005.2;Parent=gene-ZemaCp003;Dbxref=UniProtKB/Swiss-Prot:P48190,Genbank:NP_043005.2,GeneID:845178;Name=NP_043005.2;gbkey=CDS;gene=matK;locus_tag=ZemaCp003;product=maturase K;protein_id=NP_043005.2;transl_table=11
NC_001666.2     1673    3215    gene-ZemaCp003  .       -       RefSeq  gene    .       ID=gene-ZemaCp003;Dbxref=GeneID:845178;Name=matK;gbkey=Gene;gene=matK;gene_biotype=protein_coding;locus_tag=ZemaCp003
NC_001666.2     3909    3946    exon-ZemaCt121-1        .       -       RefSeq  exon    .       ID=exon-ZemaCt121-1;Parent=rna-ZemaCt121;Dbxref=GeneID:845256;gbkey=tRNA;gene=trnK;locus_tag=ZemaCt121;product=tRNA-Lys
NC_001666.2     3909    3946    id-ZemaCt121-1  .       -       RefSeq  exon    .       ID=id-ZemaCt121-1;Parent=gene-ZemaCt121;Dbxref=GeneID:845256;exon_number=1;gbkey=exon;gene=trnK;locus_tag=ZemaCt121;number=1
NC_001666.2     4490    4708    cds-NP_043006.1 .       -       RefSeq  CDS     2       ID=cds-NP_043006.1;Parent=gene-ZemaCp004;Dbxref=UniProtKB/Swiss-Prot:P27723,Genbank:NP_043006.1,GeneID:845232;Name=NP_043006.1;gbkey=CDS;gene=rps16;locus_tag=ZemaCp004;product=ribosomal protein S16;protein_id=NP_043006.1;transl_table=11
NC_001666.2     4490    4708    id-ZemaCp004-2  .       -       RefSeq  exon    .       ID=id-ZemaCp004-2;Parent=gene-ZemaCp004;Dbxref=GeneID:845232;exon_number=2;gbkey=exon;gene=rps16;locus_tag=ZemaCp004;number=2
NC_001666.2     4490    5604    gene-ZemaCp004  .       -       RefSeq  gene    .       ID=gene-ZemaCp004;Dbxref=GeneID:845232;Name=rps16;gbkey=Gene;gene=rps16;gene_biotype=protein_coding;locus_tag=ZemaCp004
NC_001666.2     4708    5564    id-ZemaCp004    .       -       RefSeq  intron  .       ID=id-ZemaCp004;Parent=gene-ZemaCp004;Dbxref=GeneID:845232;exon_number=1;gbkey=intron;gene=rps16;locus_tag=ZemaCp004;number=1
NC_001666.2     5564    5604    cds-NP_043006.1 .       -       RefSeq  CDS     0       ID=cds-NP_043006.1;Parent=gene-ZemaCp004;Dbxref=UniProtKB/Swiss-Prot:P27723,Genbank:NP_043006.1,GeneID:845232;Name=NP_043006.1;gbkey=CDS;gene=rps16;locus_tag=ZemaCp004;product=ribosomal protein S16;protein_id=NP_043006.1;transl_table=11
NC_001666.2     5564    5604    id-ZemaCp004-1  .       -       RefSeq  exon    .       ID=id-ZemaCp004-1;Parent=gene-ZemaCp004;Dbxref=GeneID:845232;exon_number=1;gbkey=exon;gene=rps16;locus_tag=ZemaCp004;number=1
NC_001666.2     6772    6844    exon-ZemaCt122-1        .       -       RefSeq  exon    .       ID=exon-ZemaCt122-1;Parent=rna-ZemaCt122;Dbxref=GeneID:845265;gbkey=tRNA;gene=trnQ;locus_tag=ZemaCt122;product=tRNA-Gln
NC_001666.2     6772    6844    gene-ZemaCt122  .       -       RefSeq  gene    .       ID=gene-ZemaCt122;Dbxref=GeneID:845265;Name=trnQ;gbkey=Gene;gene=trnQ;gene_biotype=tRNA;locus_tag=ZemaCt122
NC_001666.2     6772    6844    rna-ZemaCt122   .       -       RefSeq  tRNA    .       ID=rna-ZemaCt122;Parent=gene-ZemaCt122;Dbxref=GeneID:845265;gbkey=tRNA;gene=trnQ;locus_tag=ZemaCt122;product=tRNA-Gln
NC_001666.2     7198    7384    cds-NP_043007.1 .       +       RefSeq  CDS     0       ID=cds-NP_043007.1;Parent=gene-ZemaCp005;Dbxref=UniProtKB/Swiss-Prot:P48188,Genbank:NP_043007.1,GeneID:845208;Name=NP_043007.1;gbkey=CDS;gene=psbK;locus_tag=ZemaCp005;product=photosystem II protein K;protein_id=NP_043007.1;transl_table=11
NC_001666.2     7198    7384    gene-ZemaCp005  .       +       RefSeq  gene    .       ID=gene-ZemaCp005;Dbxref=GeneID:845208;Name=psbK;gbkey=Gene;gene=psbK;gene_biotype=protein_coding;locus_tag=ZemaCp005
NC_001666.2     7774    7885    cds-NP_043008.1 .       +       RefSeq  CDS     0       ID=cds-NP_043008.1;Parent=gene-ZemaCp006;Dbxref=UniProtKB/Swiss-Prot:P09970,Genbank:NP_043008.1,GeneID:845206;Name=NP_043008.1;gbkey=CDS;gene=psbI;locus_tag=ZemaCp006;product=photosystem II protein I;protein_id=NP_043008.1;transl_table=11
NC_001666.2     7774    7885    gene-ZemaCp006  .       +       RefSeq  gene    .       ID=gene-ZemaCp006;Dbxref=GeneID:845206;Name=psbI;gbkey=Gene;gene=psbI;gene_biotype=protein_coding;locus_tag=ZemaCp006
NC_001666.2     8011    8099    exon-ZemaCt123-1        .       -       RefSeq  exon    .       ID=exon-ZemaCt123-1;Parent=rna-ZemaCt123;Dbxref=GeneID:845270;gbkey=tRNA;gene=trnS;locus_tag=ZemaCt123;product=tRNA-Ser
NC_001666.2     8011    8099    gene-ZemaCt123  .       -       RefSeq  gene    .       ID=gene-ZemaCt123;Dbxref=GeneID:845270;Name=trnS;gbkey=Gene;gene=trnS;gene_biotype=tRNA;locus_tag=ZemaCt123
NC_001666.2     8011    8099    rna-ZemaCt123   .       -       RefSeq  tRNA    .       ID=rna-ZemaCt123;Parent=gene-ZemaCt123;Dbxref=GeneID:845270;gbkey=tRNA;gene=trnS;locus_tag=ZemaCt123;product=tRNA-Ser
NC_001666.2     9082    10144   cds-NP_043009.1 .       +       RefSeq  CDS     0       ID=cds-NP_043009.1;Parent=gene-ZemaCp007;Dbxref=UniProtKB/Swiss-Prot:P48184,Genbank:NP_043009.1,GeneID:845202;Name=NP_043009.1;gbkey=CDS;gene=psbD;locus_tag=ZemaCp007;product=photosystem II protein D2;protein_id=NP_043009.1;transl_table=11
NC_001666.2     9082    10144   gene-ZemaCp007  .       +       RefSeq  gene    .       ID=gene-ZemaCp007;Dbxref=GeneID:845202;Name=psbD;gbkey=Gene;gene=psbD;gene_biotype=protein_coding;locus_tag=ZemaCp007
NC_001666.2     10091   11513   cds-NP_043010.1 .       +       RefSeq  CDS     0       ID=cds-NP_043010.1;Parent=gene-ZemaCp008;Dbxref=UniProtKB/Swiss-Prot:P48187,Genbank:NP_043010.1,GeneID:845201;Name=NP_043010.1;Note=CP43;gbkey=CDS;gene=psbC;locus_tag=ZemaCp008;product=photosystem II 44 kDa protein;protein_id=NP_043010.1;transl_table=11
NC_001666.2     10091   11513   gene-ZemaCp008  .       +       RefSeq  gene    .       ID=gene-ZemaCp008;Dbxref=GeneID:845201;Name=psbC;gbkey=Gene;gene=psbC;gene_biotype=protein_coding;locus_tag=ZemaCp008
NC_001666.2     11668   11756   exon-ZemaCt124-1        .       -       RefSeq  exon    .       ID=exon-ZemaCt124-1;Parent=rna-ZemaCt124;Dbxref=GeneID:845271;gbkey=tRNA;gene=trnS;locus_tag=ZemaCt124;product=tRNA-Ser
NC_001666.2     11668   11756   gene-ZemaCt124  .       -       RefSeq  gene    .       ID=gene-ZemaCt124;Dbxref=GeneID:845271;Name=trnS;gbkey=Gene;gene=trnS;gene_biotype=tRNA;locus_tag=ZemaCt124
NC_001666.2     11668   11756   rna-ZemaCt124   .       -       RefSeq  tRNA    .       ID=rna-ZemaCt124;Parent=gene-ZemaCt124;Dbxref=GeneID:845271;gbkey=tRNA;gene=trnS;locus_tag=ZemaCt124;product=tRNA-Ser
NC_001666.2     12016   12205   cds-NP_043011.1 .       +       RefSeq  CDS     0       ID=cds-NP_043011.1;Parent=gene-ZemaCp009;Dbxref=UniProtKB/Swiss-Prot:Q33300,Genbank:NP_043011.1,GeneID:1466358;Name=NP_043011.1;Note=YCF9;gbkey=CDS;gene=psbZ;locus_tag=ZemaCp009;product=photosystem II protein Z;protein_id=NP_043011.1;transl_table=11
NC_001666.2     12016   12205   gene-ZemaCp009  .       +       RefSeq  gene    .       ID=gene-ZemaCp009;Dbxref=GeneID:1466358;Name=psbZ;gbkey=Gene;gene=psbZ;gene_biotype=protein_coding;locus_tag=ZemaCp009
NC_001666.2     12508   12579   exon-ZemaCt125-1        .       +       RefSeq  exon    .       ID=exon-ZemaCt125-1;Parent=rna-ZemaCt125;Dbxref=GeneID:845248;gbkey=tRNA;gene=trnG;locus_tag=ZemaCt125;product=tRNA-Gly
NC_001666.2     12508   12579   gene-ZemaCt125  .       +       RefSeq  gene    .       ID=gene-ZemaCt125;Dbxref=GeneID:845248;Name=trnG;gbkey=Gene;gene=trnG;gene_biotype=tRNA;locus_tag=ZemaCt125
NC_001666.2     12508   12579   rna-ZemaCt125   .       +       RefSeq  tRNA    .       ID=rna-ZemaCt125;Parent=gene-ZemaCt125;Dbxref=GeneID:845248;gbkey=tRNA;gene=trnG;locus_tag=ZemaCt125;product=tRNA-Gly
NC_001666.2     12985   12986   id-NC_001666.2:12986..12986     .       +       RefSeq  sequence_conflict       .       ID=id-NC_001666.2:12986..12986;gbkey=conflict;replace=;zero_length_insertion=True
NC_001666.2     12998   12999   id-NC_001666.2:12999..12999     .       +       RefSeq  sequence_conflict       .       ID=id-NC_001666.2:12999..12999;gbkey=conflict;replace=a;zero_length_insertion=True
NC_001666.2     13072   13146   exon-ZemaCt126-1        .       -       RefSeq  exon    .       ID=exon-ZemaCt126-1;Parent=rna-ZemaCt126;Dbxref=GeneID:845280;gbkey=tRNA;gene=trnfM;locus_tag=ZemaCt126;product=tRNA-Met
NC_001666.2     13072   13146   gene-ZemaCt126  .       -       RefSeq  gene    .       ID=gene-ZemaCt126;Dbxref=GeneID:845280;Name=trnfM;gbkey=Gene;gene=trnfM;gene_biotype=tRNA;locus_tag=ZemaCt126
NC_001666.2     13072   13146   rna-ZemaCt126   .       -       RefSeq  tRNA    .       ID=rna-ZemaCt126;Parent=gene-ZemaCt126;Dbxref=GeneID:845280;gbkey=tRNA;gene=trnfM;locus_tag=ZemaCt126;product=tRNA-Met
NC_001666.2     13244   13292   exon-ZemaCt127-2        .       -       RefSeq  exon    .       ID=exon-ZemaCt127-2;Parent=rna-ZemaCt127;Dbxref=GeneID:845249;gbkey=tRNA;gene=trnG;locus_tag=ZemaCt127;product=tRNA-Gly
NC_001666.2     13244   13292   id-ZemaCt127-2  .       -       RefSeq  exon    .       ID=id-ZemaCt127-2;Parent=gene-ZemaCt127;Dbxref=GeneID:845249;exon_number=2;gbkey=exon;gene=trnG;locus_tag=ZemaCt127;number=2
NC_001666.2     13244   14013   gene-ZemaCt127  .       -       RefSeq  gene    .       ID=gene-ZemaCt127;Dbxref=GeneID:845249;Name=trnG;gbkey=Gene;gene=trnG;gene_biotype=tRNA;locus_tag=ZemaCt127
NC_001666.2     13244   14013   rna-ZemaCt127   .       -       RefSeq  tRNA    .       ID=rna-ZemaCt127;Parent=gene-ZemaCt127;Dbxref=GeneID:845249;gbkey=tRNA;gene=trnG;locus_tag=ZemaCt127;product=tRNA-Gly
NC_001666.2     13292   13990   id-ZemaCt127    .       -       RefSeq  intron  .       ID=id-ZemaCt127;Parent=gene-ZemaCt127;Dbxref=GeneID:845249;exon_number=1;gbkey=intron;gene=trnG;locus_tag=ZemaCt127;number=1
NC_001666.2     13990   14013   exon-ZemaCt127-1        .       -       RefSeq  exon    .       ID=exon-ZemaCt127-1;Parent=rna-ZemaCt127;Dbxref=GeneID:845249;gbkey=tRNA;gene=trnG;locus_tag=ZemaCt127;product=tRNA-Gly
NC_001666.2     13990   14013   id-ZemaCt127-1  .       -       RefSeq  exon    .       ID=id-ZemaCt127-1;Parent=gene-ZemaCt127;Dbxref=GeneID:845249;exon_number=1;gbkey=exon;gene=trnG;locus_tag=ZemaCt127;number=1
NC_001666.2     14497   14707   cds-NP_043012.1 .       -       RefSeq  CDS     0       ID=cds-NP_043012.1;Parent=gene-ZemaCp010;Dbxref=UniProtKB/TrEMBL:Q33301,Genbank:NP_043012.1,GeneID:1466359;Name=NP_043012.1;Note=ORF69;gbkey=CDS;locus_tag=ZemaCp010;product=hyp
etc...
```
</pre>
</details>



### run RSeQC
```
#assesses the strand information of your alignment
for f in *sorted.bam; do infer_experiment.py --i $f -r maizeGFF.bed12 ;done
```

<details>
  <summary>Printed content of the above for loop</summary>
  <pre>

```
infer_experiment.py --i SRR13332812_1_sorted.bam -r maizeGFF.bed12
infer_experiment.py --i SRR13332813_1_sorted.bam -r maizeGFF.bed12
infer_experiment.py --i SRR1573504_1Gene_sorted.bam -r maizeGFF.bed12
infer_experiment.py --i SRR1573505_1Gene_sorted.bam -r maizeGFF.bed12
```
</pre>
</details>

### Results for unstranded reads
```

[E::idx_find_and_load] Could not retrieve index file for 'SRR1573504_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.5082
Fraction of reads explained by "1++,1--,2+-,2-+": 0.2406
Fraction of reads explained by "1+-,1-+,2++,2--": 0.2512
[E::idx_find_and_load] Could not retrieve index file for 'SRR1573505_1Gene_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.3724
Fraction of reads explained by "1++,1--,2+-,2-+": 0.3053
Fraction of reads explained by "1+-,1-+,2++,2--": 0.3222
```

### Results for stranded reads
```
[E::idx_find_and_load] Could not retrieve index file for 'SRR13332812_1_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.4494
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0349
Fraction of reads explained by "1+-,1-+,2++,2--": 0.5157
[E::idx_find_and_load] Could not retrieve index file for 'SRR13332813_1_sorted.bam'
Reading reference gene model maizeGFF.bed12 ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.4496
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0331
Fraction of reads explained by "1+-,1-+,2++,2--": 0.5174
```

We see that the stranded reads have an RF orientation and can make adjustments in the aligner to align these reads using strand information.  

### References

[RSeQC Sourceforge](http://rseqc.sourceforge.net/#use-pip3-to-install-rseqc-v3-0-0-or-newer)
