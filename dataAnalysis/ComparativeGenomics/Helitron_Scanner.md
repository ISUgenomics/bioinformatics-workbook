---
title: Genome Repeats Identification
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Finding Helitrons in your genome using HelitronScanner

Helitrons are rolling-circle transposons that are hard to identify, due to their lack of terminal repeats and target site duplications. Because Helitrons have poor sequence conservation, Helitronscanner was developed to identify DNA motifs conserved in helitrons, and then uses the motifs to identify similar regions in your genome.

Running Helitronscanner is essentially a four step process: identify 5' end, identify 3' end, pair ends, and create helitron fasta.
```
#get the software
wget http://bo.csam.montclair.edu/du/assets/filesdb/HelitronScanner.zip
#get your genome
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
```

#### Identify 5' and 3' heltron sequences based on homology to known.
```
#scan for the 3' end using the default training set (1MB slices and 16 threads)
java -jar /work/GIF/remkv6/USDA/12_HelitronScanner/HelitronScanner/HelitronScanner.jar scanTail -lf TrainingSet/tail.lcvs  -bs 1000000 -g TAIR10_chr_all.fas  -th 16 -o tail.helitronscanner.out
#scan for the 5' end using the default training set(1MB slices and 16 threads)
java -jar /work/GIF/remkv6/USDA/12_HelitronScanner/HelitronScanner/HelitronScanner.jar scanHead -lf TrainingSet/head.lcvs  -bs 1000000 -g TAIR10_chr_all.fas  -th 16 -o head.helitronscanner.out
```

#### Pair the helitron ends
```
#you can mess with the -hlr option to increase or decrease the min and max length of predictions
java -jar /work/GIF/remkv6/USDA/12_HelitronScanner/HelitronScanner/HelitronScanner.jar pairends -hs head.helitronscanner.out -ts tail.helitronscanner.out -hlr 200:20000 -o paired.helitrons
```

Output from pairing
```
>1
20615158:20626081 [17:21]
>2
10010305:10020455 [16:16]
>3
13257658:13262521 [16:20]
>4

>5
13501701:13511415 [24:18] 14061194:14080301 [17:16] 25256642:25265336 [23:18]
```
#### Create the fasta sequenes for each helitrons
```
java -jar /work/GIF/remkv6/USDA/12_HelitronScanner/HelitronScanner/HelitronScanner.jar draw -p paired.helitrons -g TAIR10_chr_all.fas -o draw_helitrons_pure --pure
```

How many and how big are the helitrons identified?
```
less draw_helitrons_pure.hel.fa|bioawk -c fastx '{print length($seq)}' |summary.sh
Total:  63,457
Count:  6
Mean:   10,576
Median: 9,933
Min:    4,864
Max:    19,108
#So 6 with a mean length of 10.5kb.
```
