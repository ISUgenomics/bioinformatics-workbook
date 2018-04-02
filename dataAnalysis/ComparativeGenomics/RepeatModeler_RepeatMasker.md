# Repeatmodeler is a repeat-identifying software that can provide a list of repeat family sequences that can then be used to mask repeats in a genome with RepeatMasker.
Things to consider with this software is that it can take a long time with large genomes (>1Gb==>96hrs on a 16 cpu node).  You also need to set the correct parameters in repeatmodeler so that you get repeats that are not only grouped by family, but are also annotated.

Repeatmodeler http://www.repeatmasker.org/RepeatModeler/
RepeatMasker http://www.repeatmasker.org/RMDownload.html

### Get your genome
```
wget ftp://ftp.bioinfo.wsu.edu/species/Gossypium_raimondii/JGI_221_G.raimondii_Dgenome/assembly/G.raimondii_JGI_221_v2.0.assembly.fasta.gz

gunzip G.raimondii_JGI_221_v2.0.assembly.fasta.gz
```

### Run repeatmodeler
```
module load GIF2/repeatmodeler/1.0.8
BuildDatabase -name  G.raimondii_JGI_221_v2.0.assembly.DB -engine ncbi -pa 16 G.raimondii_JGI_221_v2.0.assembly.fasta
RepeatModeler -database  G.raimondii_JGI_221_v2.0.assembly.DB -engine ncbi -pa 16
```



### Run RepeatMasker
```
module load repeatmasker/4.0.6
ln -s RM_10210.MonJan151105342018/consensi.fa.classified
RepeatMasker -pa 16 -gff -lib consensi.fa.classified G.raimondii_JGI_221_v2.0.assembly.fasta
```
