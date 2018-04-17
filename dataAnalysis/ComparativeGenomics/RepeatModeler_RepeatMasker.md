# Repeatmodeler is a repeat-identifying software that can provide a list of repeat family sequences to mask repeats in a genome with RepeatMasker.
Things to consider with this software is that it can take a long time with large genomes (>1Gb==>96hrs on a 16 cpu node).  You also need to set the correct parameters in repeatmodeler so that you get repeats that are not only grouped by family, but are also annotated.

Repeatmodeler http://www.repeatmasker.org/RepeatModeler/
RepeatMasker http://www.repeatmasker.org/RMDownload.html

### Get your genome and unzip
```
wget ftp://ftp.bioinfo.wsu.edu/species/Gossypium_raimondii/JGI_221_G.raimondii_Dgenome/assembly/G.raimondii_JGI_221_v2.0.assembly.fasta.gz
gunzip G.raimondii_JGI_221_v2.0.assembly.fasta.gz
```

### Run repeatmodeler
This can take longer than 96 hours on one node with 16threads if the genome is larger than 1Gb
```
module load GIF2/repeatmodeler/1.0.8
# make up a name for your database, choose your search engine, the number of threads, and the genome file
BuildDatabase -name  G.raimondii_JGI_221_v2.0.assembly.DB -engine ncbi -pa 16 G.raimondii_JGI_221_v2.0.assembly.fasta

#Now run repeatmodeler with the database made above
RepeatModeler -database  G.raimondii_JGI_221_v2.0.assembly.DB -engine ncbi -pa 16
```



### Run RepeatMasker
This takes about 24-48 hours to finish on a genome less than 1Gb.
```
module load repeatmasker/4.0.6
#I moved to a different directory, so I softlinked my classified file.  Make sure you use the consensi.fa.classified file, or your repeats will just be masked by repeatmasker, but unannotated.

ln -s RM_10210.MonJan151105342018/consensi.fa.classified

RepeatMasker -pa 16 -gff -lib consensi.fa.classified G.raimondii_JGI_221_v2.0.assembly.fasta
```
