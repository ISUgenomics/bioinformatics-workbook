# MEME and FIMO tutorial
This tutorial will go through how to identify a set of dna motifs that among a list of sequences. Once the motifs have been identified, they can be used to search a larger database of sequences.  This is really useful when trying to find patterns of conserved sequences in large databases of sequences.

Here I will try to identify dna sequences that are conserved in the upstream sequences of genes involved in seed dormancy.

## Collect sequences for motif finding
```
#/work/GIF/remkv6/USDA/18_MEMESeedDormancyFIMO
wget https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas


#how many genes involved in ending seed dormancy? (substr just removes isoform data)
less ATH_GO_GOSLIM.txt |grep "seed dormancy" |awk '{print $3}' |sort|uniq|awk '{print substr($1,0,9)}' |sort|uniq|wc
    569     569    5287


#This essentially searches for genes with the annotation "seed dormancy", extracts the name from the gff, then uses the name to extract the position of each gene, such that 1000bp can be extracted 5' from the end of genes.
less ATH_GO_GOSLIM.txt |grep "seed dormancy" |awk '{print $3}' |sort|uniq|awk '{print substr($1,0,9)}' |sort|uniq|grep -f - TAIR10_GFF3_genes.gff |awk '$3=="gene"' |awk '{if($6=="-"){print substr($1,4,4),$5,$5+1000} else {print substr($1,4,4),$4-1000,$4}}' |awk '{print "samtools faidx  TAIR10_chr_all.fas "$1":"$2"-"$3}' >extractSeedDormancyPromoters.sh



#I used positions to find the 1000 bp upstream sequences, but found the curated sequences online later.
module load samtools
samtools faidx TAIR10_chr_all.fas
sh extractSeedDormancyPromoters.sh >SeedDormancy.fasta


### Revised approach using curated upstream regions
module load cdbfasta
cdbfasta TAIR10_upstream_1000_translation_start_20101028
less ATH_GO_GOSLIM.txt |grep "seed dormancy" |awk '{print $3}' |sort|uniq|awk '{print substr($1,0,9)}' |sort|uniq|cdbyank TAIR10_upstream_1000_translation_start_20101028.cidx >SeedDormancy.fasta

#Lets establish a background sequence file of utrS that are not seed dormancy related
less ATH_GO_GOSLIM.txt |grep "seed dormancy" |awk '{print $3}' |sort|uniq|awk '{print substr($1,0,9)}' |sort|uniq|grep -v -f - <(grep ">" TAIR10_upstream_1000_translation_start_20101028) |awk '{print $1}' |sed 's/>//g' |cdbyank TAIR10_upstream_1000_translation_start_20101028.cidx >AllOther5UTR.fasta

```

## Run MEME


![Meme Manual](http://meme-suite.org/doc/meme.html?man_type=web)
```
#this one compensates for the sequences being DNA, allows zero or one motif, a max of 5 motifs, both strands, a larger maximum memory available, and a max number of sites equal to 453.
echo "module load meme" >meme.sh
echo "meme SeedDormancy.fasta -oc SeedDormancyMeme -dna -mod zoops -nmotifs 5 -revcomp -maxsize 100000000 -maxsites 453 " >>meme.sh

#This one allows the identification of larger motifs, 3 motifs max, and a max site of 500, minimum width of 50bp, and max width of 200bp
meme SeedDormancy.fasta -oc SeedDormancyMemeLg -dna -mod zoops -nmotifs 3 -revcomp -maxsize 1000000 -maxsites 500 -minw 50 -maxw 200

```

## Scan for motifs in the upstream sequences for each gene

Here is a link to the readme for each option.

![FIMO Manual](http://meme-suite.org/doc/fimo.html)

```
#/work/GIF/remkv6/USDA/18_MEMESeedDormancyFIMO

#found the curated gene upstream sequences for arabidopsis TAIR10
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_blastsets/upstream_sequences/TAIR10_upstream_1000_translation_start_20101028


#FIMO on the small motifs
/opt/rit/app/meme/4.12.0/bin/fimo -oc smallMotifs SeedDormancyMeme/meme.xml  TAIR10_upstream_1000_translation_start_20101028


#FIMO on the larger motifs
/opt/rit/app/meme/4.12.0/bin/fimo SeedDormancyMemeLg/meme.xml  TAIR10_upstream_1000_translation_start_20101028

```

## Look to see which genes have these motifs

```
### The Small Motifs
less fimo.txt |awk '$2=="MEME-1"{print $3}' |sort|uniq|wc
  14990   14990  149900
less fimo.txt |awk '$2=="MEME-2"{print $3}' |sort|uniq|wc
  16939   16939  169390
less fimo.txt |awk '$2=="MEME-3"{print $3}' |sort|uniq|wc
  19426   19426  194260
less fimo.txt |awk '$2=="MEME-4"{print $3}' |sort|uniq|wc
  10764   10764  107640
less fimo.txt |awk '$2=="MEME-5"{print $3}' |sort|uniq|wc
  19826   19826  198260

###The Large Motifs
#How many genes have motif 1
less fimo.txt |awk '$2=="MEME-1"{print $3}' |sort|uniq|wc
  14262   14262  142620
#How many genes have motif 2  
less fimo.txt |awk '$2=="MEME-2"{print $3}' |sort|uniq|wc
  14480   14480  144800
#How many genes have motif 3
less fimo.txt |awk '$2=="MEME-3"{print $3}' |sort|uniq|wc
   4333    4333   43330
```
