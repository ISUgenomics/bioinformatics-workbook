---
title: Genome Repeats Identification
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Repeatmodeler is a repeat-identifying software that can provide a list of repeat family sequences to mask repeats in a genome with RepeatMasker.
Things to consider with this software is that it can take a long time with large genomes (>1Gb==>96hrs on a 16 cpu node).  You also need to set the correct parameters in repeatmodeler so that you get repeats that are not only grouped by family, but are also annotated.

Repeatmodeler http://www.repeatmasker.org/RepeatModeler/
RepeatMasker http://www.repeatmasker.org/RMDownload.html

### Get your genome and unzip
```
#I will be using the Araport11 Arabidopsis genome
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
TAIR10_chr_all.fas
```

### Run repeatmodeler
This can take longer than 96 hours on one node with 16threads if the genome is larger than 1Gb.  For arabidopsis it took 9.5 hours
```
# make up a name for your database, choose your search engine, the number of threads, and the genome file
module load GIF2/repeatmodeler/1.0.8
BuildDatabase -name TAIR10_chr_all.DB -engine rmblast TAIR10_chr_all.fas
RepeatModeler -database TAIR10_chr_all.DB -engine ncbi -pa 16
```



### Run RepeatMasker
This can take takes about 24-48 hours to finish on a genome over 1Gb.  However the arabidopsis run below took 14 minutes with 16 threads.
```
#I moved to a different directory, so I softlinked my classified file.  Make sure you use the consensi.fa.classified file, or your repeats will just be masked by repeatmasker, but unannotated.

#make sure you softlink the classified file, otherwise you will not get a table of classified elements after the run.
ln -s RM_3245.WedMay21605262018/consensi.fa.classified


#This will produce a gff for the repeat mapping, a masked fasta file, and a table summarizing the repeats found in the genome.
module load GIF2/repeatmasker/4.0.6
RepeatMasker -pa 16 -gff -lib consensi.fa.classified TAIR10_chr_all.fas
```

This is what I find in Arabidopsis.
```
==================================================
file name: TAIR10_chr_all.fas
sequences:             7
total length:  119667750 bp  (119482146 bp excl N/X-runs)
GC level:         36.06 %
bases masked:   17878420 bp ( 14.94 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
SINEs:              381        91996 bp    0.08 %
      ALUs            0            0 bp    0.00 %
      MIRs            0            0 bp    0.00 %

LINEs:             2505      1306247 bp    1.09 %
      LINE1        2407      1277415 bp    1.07 %
      LINE2           0            0 bp    0.00 %
      L3/CR1          0            0 bp    0.00 %

LTR elements:      6875      7287909 bp    6.09 %
      ERVL            0            0 bp    0.00 %
      ERVL-MaLRs      0            0 bp    0.00 %
      ERV_classI      0            0 bp    0.00 %
      ERV_classII     0            0 bp    0.00 %

DNA elements:      5916      2718611 bp    2.27 %
     hAT-Charlie    435       163286 bp    0.14 %
     TcMar-Tigger     0            0 bp    0.00 %

Unclassified:      7261      3594123 bp    3.00 %

Total interspersed repeats: 14998886 bp   12.53 %


Small RNA:          442       115046 bp    0.10 %

Satellites:        1064       981082 bp    0.82 %
Simple repeats:   35831      1435203 bp    1.20 %
Low complexity:    9032       443907 bp    0.37 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element


The query species was assumed to be homo
RepeatMasker version open-4.0.6 , default mode

run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensi.fa.classified"
RepBase Update 20160829, RM database version 20160829

```
Now there is also a GFF that can be used for many other genomic comparisons.
```
##gff-version 2
##date 2018-05-03
##sequence-region TAIR10_chr_all.fas
1       RepeatMasker    similarity      1       115     13.1    +       .       Target "Motif:A-rich" 1 107
1       RepeatMasker    similarity      1066    1097    10.0    +       .       Target "Motif:(C)n" 1 32
1       RepeatMasker    similarity      1155    1187    17.1    +       .       Target "Motif:(TTTCTT)n" 1 33
1       RepeatMasker    similarity      4291    4328     8.4    +       .       Target "Motif:(AT)n" 1 38
1       RepeatMasker    similarity      5680    5702     9.3    +       .       Target "Motif:(T)n" 1 23
1       RepeatMasker    similarity      8669    8699     0.0    +       .       Target "Motif:(CT)n" 1 31
1       RepeatMasker    similarity      9961    10030   20.7    +       .       Target "Motif:(AT)n" 1 70
1       RepeatMasker    similarity      10814   10885   28.7    +       .       Target "Motif:(AT)n" 1 71
1       RepeatMasker    similarity      11915   11960   12.0    +       .       Target "Motif:(ATC)n" 1 46
1       RepeatMasker    similarity      11985   12001    0.0    +       .       Target "Motif:(GAA)n" 1 17
1       RepeatMasker    similarity      12875   12917   16.5    +       .       Target "Motif:(TTCTTG)n" 1 44
1       RepeatMasker    similarity      12985   12997   21.8    +       .       Target "Motif:(TGGTTTT)n" 1 36
1       RepeatMasker    similarity      12998   13021    8.9    +       .       Target "Motif:(T)n" 1 24
1       RepeatMasker    similarity      13346   13368    4.5    +       .       Target "Motif:(AG)n" 1 23
1       RepeatMasker    similarity      15436   15469   10.8    +       .       Target "Motif:(ATT)n" 1 32
1       RepeatMasker    similarity      15872   15895    4.3    +       .       Target "Motif:(TA)n" 1 25
1       RepeatMasker    similarity      16804   16838    0.0    +       .       Target "Motif:(TTG)n" 1 31
1       RepeatMasker    similarity      17009   17256    7.3    +       .       Target "Motif:rnd-5_family-1313" 1203 1449
1       RepeatMasker    similarity      17256   17735    6.9    +       .       Target "Motif:rnd-5_family-1313" 1230 1808
1       RepeatMasker    similarity      17817   18078   17.5    +       .       Target "Motif:rnd-4_family-1372" 2109 2389
1       RepeatMasker    similarity      18100   18642    3.0    +       .       Target "Motif:rnd-5_family-843" 1 549
1       RepeatMasker    similarity      18661   18731    2.8    +       .       Target "Motif:rnd-5_family-843" 1156 1226
1       RepeatMasker    similarity      20510   20557   19.0    +       .       Target "Motif:(AT)n" 1 50
1       RepeatMasker    similarity      23109   23167   17.3    +       .       Target "Motif:GA-rich" 1 52
1       RepeatMasker    similarity      34438   34456    5.6    +       .       Target "Motif:(TTG)n" 1 19
1       RepeatMasker    similarity      37736   37777    0.0    +       .       Target "Motif:(GAA)n" 1 41
1       RepeatMasker    similarity      37790   37825    6.9    +       .       Target "Motif:GA-rich" 1 34
1       RepeatMasker    similarity      41361   41385    0.0    +       .       Target "Motif:(TA)n" 1 25
1       RepeatMasker    similarity      41549   41567    5.5    +       .       Target "Motif:(TA)n" 1 19
1       RepeatMasker    similarity      41716   41760   30.0    +       .       Target "Motif:A-rich" 1 45
1       RepeatMasker    similarity      42444   42488   18.2    +       .       Target "Motif:(T)n" 1 45
1       RepeatMasker    similarity      43659   43713   30.3    +       .       Target "Motif:(AATTTT)n" 1 54
1       RepeatMasker    similarity      46523   46564   22.0    +       .       Target "Motif:(TTC)n" 1 42
1       RepeatMasker    similarity      46832   46864   18.1    +       .       Target "Motif:A-rich" 1 32
1       RepeatMasker    similarity      47067   47290   17.8    -       .       Target "Motif:rnd-5_family-3415" 1427 1673
1       RepeatMasker    similarity      49387   49415   15.3    +       .       Target "Motif:(TTCTT)n" 1 30
1       RepeatMasker    similarity      50433   50515   29.6    +       .       Target "Motif:(CAG)n" 1 83
1       RepeatMasker    similarity      53368   53434   20.4    +       .       Target "Motif:(TAATTTG)n" 1 69
1       RepeatMasker    similarity      55677   55919    7.9    +       .       Target "Motif:rnd-5_family-7947" 1 242
1       RepeatMasker    similarity      55880   56020   20.0    +       .       Target "Motif:rnd-5_family-2172" 1 117
1       RepeatMasker    similarity      56021   56236   14.9    +       .       Target "Motif:rnd-5_family-1066" 747 946
1       RepeatMasker    similarity      56237   56293   20.0    +       .       Target "Motif:rnd-5_family-2172" 118 165
1       RepeatMasker    similarity      56311   56576    6.8    +       .       Target "Motif:rnd-5_family-4425" 1129 1398
1       RepeatMasker    similarity      59061   59095   19.6    +       .       Target "Motif:A-rich" 1 35
1       RepeatMasker    similarity      61867   61893   16.6    +       .       Target "Motif:(T)n" 1 27
1       RepeatMasker    similarity      62347   62372    0.0    +       .       Target "Motif:(AAG)n" 1 26
1       RepeatMasker    similarity      62392   62803   16.2    +       .       Target "Motif:rnd-5_family-2965" 1 531
1       RepeatMasker    similarity      63507   63531    0.0    +       .       Target "Motif:(TCTTTC)n" 1 25
1       RepeatMasker    similarity      64859   64872   27.7    +       .       Target "Motif:A-rich" 1 37
#truncated file for visualization of gff
```

---
[Table of contents](Repeats_index.md)
