---
title: Kraken metagenomics on nematodes
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

#  Here I will try to see what kind of bacteria and viruses lie within the Tylenchida nematode RNAseq.


## Build appropriate kraken2 database
```
#/work/GIF/remkv6/USDA/21_kraken

#softlink the fastq
mkdir fastq
for f in ../../../Baum/05_738NewAnalyses/04_NewRNAseqOnlineAlignment4Jbrowse/01_trim/*fastq.gz; do ln -s $f;done



#download taxonomy data from ncbi
module load GIF/kraken2
#kraken2-build --download-taxonomy --db PlantViral
#This should have worked, but NCBI removed a couple files from their ftp site (est and gss)
Had to modify the download_taxonomy.sh script to skip downloading these files.


vi download_taxonomy.sh
########
Changed from
TAXONOMY_DIR="$KRAKEN2_DB_NAME/taxonomy"
if [ -z "$KRAKEN2_PROTEIN_DB" ]

To
for subsection in gb wgs
if [ -z "Viral" ]
########
sh download_taxonomy.sh
```

#### To add your genomes to the kraken database, you will have to look up the taxonomy ID and add this to each fasta header. i.e. (>sequence"|kraken:taxid|390850)
```
I found the taxids here: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
bioawk -c fastx '{print ">"$name"|kraken:taxid|6326\n"$seq}' B.xylophilus.fa >B.xylophilusTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|166010\n"$seq}' D.destructor.fa >D.destructorTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|166011\n"$seq}' D.dipsaci.fa >D.dipsaciTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|1517492\n"$seq}' G.ellingtonae.fa >G.ellingtonaeTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|36090\n"$seq}' G.pallida.fa >G.pallidaTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|31243\n"$seq}' G.rostochiensis.fa >G.rostochiensisTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|6304\n"$seq}' M.arenaria.fa >M.arenariaTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|390850\n"$seq}' M.enterolobii.fa >M.enterolobiiTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|298850\n"$seq}' M.floridensis.fa >M.floridensisTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|189291\n"$seq}' M.graminicola.fa >M.graminicolaTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|6305\n"$seq}' M.hapla.fa >M.haplaTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|6306\n"$seq}' M.incognita.fa >M.incognitaTax.fa
bioawk -c fastx '{print ">"$name"|kraken:taxid|6303\n"$seq}' M.javanica.fa >M.javanicaTax.fa
bioawk -c fastx '{print ">"$name"kraken:taxid|51029\n"$seq}' H.glycines.fa >H.glycinesTax.fa

```

### We are now ready to set the database to download and build
```
module use /work/GIF/software/modules
module load GIF/kraken2
module load blast-plus
kraken2-build --download-library viral --db NematodeViral
kraken2-build --download-library bacteria --db NematodeViral
kraken2-build --add-to-library B.xylophilusTax.fa -db NematodeViral
kraken2-build --add-to-library D.destructorTax.fa -db NematodeViral
kraken2-build --add-to-library D.dipsaciTax.fa -db NematodeViral
kraken2-build --add-to-library G.ellingtonaeTax.fa -db NematodeViral
kraken2-build --add-to-library G.pallidaTax.fa -db NematodeViral
kraken2-build --add-to-library G.rostochiensisTax.fa -db NematodeViral
kraken2-build --add-to-library H.glycinesTax.fa -db NematodeViral
kraken2-build --add-to-library M.arenariaTax.fa -db NematodeViral
kraken2-build --add-to-library M.enterolobiiTax.fa -db NematodeViral
kraken2-build --add-to-library M.floridensisTax.fa -db NematodeViral
kraken2-build --add-to-library M.graminicolaTax.fa -db NematodeViral
kraken2-build --add-to-library M.haplaTax.fa -db NematodeViral
kraken2-build --add-to-library M.incognitaTax.fa -db NematodeViral
kraken2-build --add-to-library M.javanicaTax.fa -db NematodeViral
kraken2-build --build --db NematodeViral --threads 16

#this took 5:26:32 on an hpc with 16 threads and 128GB ram.  
```

## Download RNASEQ samples
```
module load sra-toolkit
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files SRR3162514 #Globodera ellingtonae
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files SRR7775195 #Globodera pallida
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files ERR202492 #Globodera pallida
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files ERR202487 #Globodera rostochiensis
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files SRR7943144 #Ditylenchus destructor
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files DRR141214 #Bursaphelenchus xylophilus
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files ERR790020 # Meloidogyne javanica
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files SRR8691582 #Meloidogyne incognita
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files SRR2389452 #Meloidogyne graminicola
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files ERR790021 #Meloidogyne arenaria
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files SRR6269845 #Heterodera glycines
fastq-dump --outdir 04_DownloadedRNAseq/ --gzip --split-files SRR6269844 & #Heterodera glycines
```

### Create Kraken scripts
```
paste <(ls -1 */*R1*gz) <(ls -1 */*R1*gz) |while read line; do echo "kraken2 -db Plant --threads 16 --report "$line".report --gzip-compressed  --unclassified-out "${line%.*}"unclassified#.fq --classified-out "${line%.*}"classified#.fq --paired "$line" > "${line%.*} ;done |awk '{print $1,$2,$3,$4,$5,$6,$8,$9,$10,$12,$13,$15,$16,$17,$18,$19,$21".Kraken.out"}' >kraken.sh


#kraken.sh
##############################################################################################################################################################################
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/DRR141214_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/DRR141214_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/DRR141214_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/DRR141214_1.fastq.gz 04_DownloadedRNAseq/DRR141214_2.fastq.gz > 04_DownloadedRNAseq/DRR141214_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/ERR202487_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/ERR202487_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/ERR202487_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/ERR202487_1.fastq.gz 04_DownloadedRNAseq/ERR202487_2.fastq.gz > 04_DownloadedRNAseq/ERR202487_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/ERR202492_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/ERR202492_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/ERR202492_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/ERR202492_1.fastq.gz 04_DownloadedRNAseq/ERR202492_2.fastq.gz > 04_DownloadedRNAseq/ERR202492_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/ERR790020_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/ERR790020_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/ERR790020_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/ERR790020_1.fastq.gz 04_DownloadedRNAseq/ERR790020_2.fastq.gz > 04_DownloadedRNAseq/ERR790020_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/SRR3162514_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/SRR3162514_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/SRR3162514_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/SRR3162514_1.fastq.gz 04_DownloadedRNAseq/SRR3162514_2.fastq.gz > 04_DownloadedRNAseq/SRR3162514_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/SRR7775195_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/SRR7775195_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/SRR7775195_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/SRR7775195_1.fastq.gz 04_DownloadedRNAseq/SRR7775195_2.fastq.gz > 04_DownloadedRNAseq/SRR7775195_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/SRR7943144_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/SRR7943144_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/SRR7943144_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/SRR7943144_1.fastq.gz 04_DownloadedRNAseq/SRR7943144_2.fastq.gz > 04_DownloadedRNAseq/SRR7943144_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/SRR7943144_2_val_2.fq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/SRR7943144_2_val_2.fqunclassified#.fq --classified-out 04_DownloadedRNAseq/SRR7943144_2_val_2.fqclassified#.fq --paired 04_DownloadedRNAseq/SRR7943144_1_val_1.fq.gz 04_DownloadedRNAseq/SRR7943144_2_val_2.fq.gz > 04_DownloadedRNAseq/SRR7943144_2_val_2.fq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/DRR141214_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/DRR141214_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/DRR141214_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/DRR141214_1.fastq.gz 04_DownloadedRNAseq/DRR141214_2.fastq.gz > 04_DownloadedRNAseq/DRR141214_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/ERR202487_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/ERR202487_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/ERR202487_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/ERR202487_1.fastq.gz 04_DownloadedRNAseq/ERR202487_2.fastq.gz > 04_DownloadedRNAseq/ERR202487_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/ERR202492_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/ERR202492_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/ERR202492_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/ERR202492_1.fastq.gz 04_DownloadedRNAseq/ERR202492_2.fastq.gz > 04_DownloadedRNAseq/ERR202492_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/ERR790020_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/ERR790020_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/ERR790020_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/ERR790020_1.fastq.gz 04_DownloadedRNAseq/ERR790020_2.fastq.gz > 04_DownloadedRNAseq/ERR790020_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/SRR3162514_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/SRR3162514_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/SRR3162514_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/SRR3162514_1.fastq.gz 04_DownloadedRNAseq/SRR3162514_2.fastq.gz > 04_DownloadedRNAseq/SRR3162514_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/SRR7775195_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/SRR7775195_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/SRR7775195_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/SRR7775195_1.fastq.gz 04_DownloadedRNAseq/SRR7775195_2.fastq.gz > 04_DownloadedRNAseq/SRR7775195_2.fastq.Kraken.out
kraken2 -db NematodeViral --threads 16 --report 04_DownloadedRNAseq/SRR7943144_2.fastq.gz.report --gzip-compressed --unclassified-out 04_DownloadedRNAseq/SRR7943144_2.fastqunclassified#.fq --classified-out 04_DownloadedRNAseq/SRR7943144_2.fastqclassified#.fq --paired 04_DownloadedRNAseq/SRR7943144_1.fastq.gz 04_DownloadedRNAseq/SRR7943144_2.fastq.gz > 04_DownloadedRNAseq/SRR7943144_2.fastq.Kraken.out
etc...
################################################################################################################################################################################add the correct modules and submit to hpc node

module use /work/GIF/software/modules
module load GIF/kraken2
module load perl
```

### Summarize kraken data output
```
#generates list of only those with .01% of reads and removed nematodes from the report

for f in *report; do echo "awk '\$1>0 && \$3>10' "$f" |uniq|sort -k1,1nr |grep -v \"Meloidogyne\" |grep -v \"Heterodera\" |grep -v \"Globodera\" |grep -v \"Bursaphelenchus\" |grep -v \"Ditylenchus\" >"$f".summary" ;done >summarizer.sh
sh summarizer.sh


I took these files, added the species name to the fifth column, removed those entries that had fewer than 100 reads allocated, kept only genera, species, and subspecies, and then concatenated all files for a network in cytoscape.

```
![Kraken](../../assets/KrakenNetwork.png)
Large green hexagons are the source species RNASEQ, red diamonds are viruses, and triangles are bacteria present in two or more species.  

---
[Table of contents](compGenomics_index.md)
