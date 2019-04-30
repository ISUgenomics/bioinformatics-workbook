---
title: Genome Annotation
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# MEME and FIMO tutorial
This tutorial will go through how to identify a set of dna motifs that among a list of sequences. Once the motifs have been identified, they can be used to search a larger database of sequences.  This is really useful when trying to find patterns of conserved sequences in large databases of sequences.

Here I will try to identify dna sequences that are conserved in the upstream sequences of genes involved in seed dormancy.

## Collect sequences for motif finding
```
#/work/GIF/remkv6/USDA/18_MEMESeedDormancyFIMO
wget https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
#gene upstream sequences for arabidopsis TAIR10
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_blastsets/upstream_sequences/TAIR10_upstream_1000_translation_start_20101028

#how many genes involved in ending seed dormancy? (substr just removes isoform data)
less ATH_GO_GOSLIM.txt |grep "seed dormancy" |awk '{print $3}' |sort|uniq|awk '{print substr($1,0,9)}' |sort|uniq|wc
    569     569    5287

#index the fasta
module load cdbfasta
cdbfasta TAIR10_upstream_1000_translation_start_20101028

#grab only those genes that are involved in seed dormancy
less ATH_GO_GOSLIM.txt |grep "seed dormancy" |awk '{print $3}' |sort|uniq|awk '{print substr($1,0,9)}' |sort|uniq|cdbyank TAIR10_upstream_1000_translation_start_20101028.cidx >SeedDormancy.fasta

#Lets establish a background sequence file of utrs that establish the background frequency of dinucleotides
fasta-get-markov -m 1 TAIR10_upstream_1000_translation_start_20101028 >1orderMarkovPromoters
```

## Run MEME to identify DNA motifs

Here is the meme manual, which discusses the many different options that can be used to identify motifs from a set of fasta sequences.

http://meme-suite.org/doc/meme.html?man_type=web

This meme run of 569 x 1000bp upstream sequences involved in seed dormancy took 15.5 hrs with 16 processors.

```
#this one compensates for the sequences being DNA, allows zero or one motif (zoops), a max of 5 motifs, both strands (revcomp), a larger maximum memory available, a maximum motif size of 100bp, a max number of sites equal to 500, and a markov background file of dinucleotide frequencies for the promoter (-bfile 1orderMarkovPromoters)

module unuse /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module use /opt/rit/modules
module load meme
meme SeedDormancy.fasta -oc SeedDormancyMemeLg -dna -mod zoops -nmotifs 5 -revcomp -maxsize 1000000 -maxsites 500  -maxw 100 -bfile 1orderMarkovPromoters
```
MEME provides three output files that are useful, and txt file, xml file, and an html file.  The txt is easily parsable, the xml file is useful for input, and the html file provides a nice output to peruse by eye.
[HTML Output](https://isugenomics.github.io/bioinformatics-workbook/assets/memeTut.html)
## Scan for motifs in the upstream sequences for each gene, FIMO

Here is a link to the readme for each option, including those that I did not use below.

http://meme-suite.org/doc/fimo.html

This essentially finds the motifs in the promoter fasta file, using the same background file for dinucleotide frequency that was used for MEME.

This run took a few minutes to complete with a single thread.
```
#/work/GIF/remkv6/USDA/18_MEMESeedDormancyFIMO

fimo -oc RevisedSeedDormancyFIMO --bgfile 1orderMarkovPromoters SeedDormancyMemeLg/meme.xml  TAIR10_upstream_1000_translation_start_20101028
```
[HTML Output](https://isugenomics.github.io/bioinformatics-workbook/assets/fimoTut.html)

## Look to see which genes have these motifs
```
#How many genes are there in the arabidpsis genome TAIR10
grep -c ">" ../TAIR10_upstream_1000_translation_start_20101028
   27416

#How many of these genes have each motif?
less fimo.txt |awk '$2=="MEME-1"{print $3}' |sort|uniq|wc
     16719   16719  167190
less fimo.txt |awk '$2=="MEME-2"{print $3}' |sort|uniq|wc
     16017   16017  160170
less fimo.txt |awk '$2=="MEME-3"{print $3}' |sort|uniq|wc
     17044   17044  170440
less fimo.txt |awk '$2=="MEME-4"{print $3}' |sort|uniq|wc
     11257   11257  112570
less fimo.txt |awk '$2=="MEME-5"{print $3}' |sort|uniq|wc
      9345    9345   93450

#How many genes have all 5 of these motifs
less fimo.txt |awk '{print $3}' |sort|uniq -c|awk '$1==5' |wc
      2011    4022   36198
#How many genes have all 4/5 of these motifs
less fimo.txt |awk '{print $3}' |sort|uniq -c|awk '$1==4' |wc
      2217    4434   39906
#How many genes have all 3/5 of these motifs
less fimo.txt |awk '{print $3}' |sort|uniq -c|awk '$1==3' |wc
      2441    4882   43938
#How many genes have all 2/5 of these motifs
less fimo.txt |awk '{print $3}' |sort|uniq -c|awk '$1==2' |wc
      2447    4894   44046
#How many genes have all 1/5 of these motifs      
less fimo.txt |awk '{print $3}' |sort|uniq -c|awk '$1==1' |wc
      2222    4444   39999   

#How many genes in the arabidopsis goslim functional annotation have at all five of these motifs?
less fimo.txt |awk '{print $3}' |sort|uniq -c|awk '$1==5' |awk '{print $2}' |grep -f - ../ATH_GO_GOSLIM.txt |awk '{print $1}' |sort|uniq|wc
   2234    2234   22323

#How many of the above are seed dormancy genes?
less fimo.txt |awk '{print $3}' |sort|uniq -c|awk '$1==5' |awk '{print $2}' |grep  -f - ../ATH_GO_GOSLIM.txt |grep "seed dormancy" |awk '{print $1}' |sort|uniq|wc
       33      33     330

```
Well, this isn't so promising for identifying a commonality among promoters for genes involved in seed dormancy, but it does demonstrate how to find motifs and search for motifs using statistically robust methodology with MEME and FIMO.

---
[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)
