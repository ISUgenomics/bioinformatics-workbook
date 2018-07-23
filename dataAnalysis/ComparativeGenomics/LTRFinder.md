# LTR finder, annotate LTR elements in a genome, and obtain GFF for comparisons

LTR-finder is a software focused on extensively charcaterizing LTR retrotransposons in a genome.  The algorithm was developed around identifying the common structural features of LTR retrotransposons.  Specifically, LTRfinder is quite good at identifying pairs of LTR's in the genome, and further annotates the intervening sequences for protein coding domains using Prosite. This software will by no means give a complete assessement of repeats in a genome (i.e. repeatmodeler,REPET,etc), but a refined look at LTR retroelements.<br/>
[LTR-finder publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1933203/) <br/>
[Software Homepage](http://tlife.fudan.edu.cn/ltr_finder/) <br/>
[Software Download](https://code.google.com/archive/p/ltr-finder/source) <br/>
Note, PROSITE is needed for LTR retroelement annotation, but not necessary <br/>

### Get your genome
```
#I will be using the Araport11 Arabidopsis genome
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
TAIR10_chr_all.fas
 ```

### LTR finder parameters
```
Usage  : [options] <INPUT_FASTA_FILE>
         -o NUM     gap open penalty, default is 3
         -t NUM     gap extension penalty, default is 1
         -e NUM     gap end penalty, default is 1
         -m NUM     match score, default is 2
         -u NUM     unmatch score, default is -2
         -D NUM     Max distance between 5'&3'LTR, default is 20000
         -d NUM     Min distance between 5'&3'LTR, default is 1000
         -L NUM     Max length of 5'&3'LTR, default is 3500
         -l NUM     Min length of 5'&3'LTR, default is 100
         -p NUM     min length of exact match pair, default is 20
         -g NUM     Max gap between joined pairs, default is 50
         -G NUM     Max gap between RT sub-domains, default is 2
         -T NUM     Min sub-domains found in a RT domain, default is 4
         -j NUM     Threshold for join new sequence in existed alignment
                    new alignment similarity higher than this will be joined,
                    default is 0.70
         -J NUM     Threshold for split existed alignment to two part
                    new alignment similarity lower than this will be split,
                    set this threshold lower than -j, means turn it off,
                    default is 0.90
         -S NUM     output Score limit, default is 6.00, [0,10]
         -M NUM     min LTR similarity threshold, default is 0.00, [0,1]
         -B NUM     Boundary alignment sharpness threshold, higher one.
                     one of the two edge's sharpness must higher than
                     this threshold, default is 0.400, [0,1]
         -b NUM     Boundary alignment sharpness threshold, lower one.
                     both of the two edge's sharpness must higher than
                     this threshold, default is 0.400, [0,1]
         -r NUM     PBS detecting threshold, min tRNA match length: 14, [1,18]
         -w NUM     output format: [0]-full, 1-summary, 2-table.
         -O NUM     output alignment length(only affect -w0), default is 40
         -P STR     SeqIDs, will only calculate matched SeqID
                      POSIX style regular express is supported.
         -s filename      tRNA sequence file(FASTA format)
         -f filename      data file used to draw figure
         -a ps_scan_dir   Use ps_scan to predict protein domain
         -x         Output in html format
         -E         LTR must have edge signal
                    (at least two of PBS,PPT,TSR)
         -C         detect Centriole, delete highly repeat regions
         -F 01string      Filter to choose desired result,default is 0
                     10000000000 5'-LTR must have TG
                     01000000000 5'-LTR must have CA
                     00100000000 3'-LTR must have TG
                     00010000000 3'-LTR must have CA
                     00001000000 TSR must be found
                     00000100000 PBS must be found
                     00000010000 PPT must be found
                     00000001000 RT domain muse be found
                     00000000100 Integrase core must be found
                     00000000010 Integrase c-term must be found
                     00000000001 RNase H must be found
                              -h         help
```

### Simple run of LTR finder on default settings with output as summary and table
```
/work/GIF/software/programs/ltrfinder/1.0.5/ltr_finder TAIR10_chr_all.fas -w 1 2 >LTRfinder.out
```
Process the ltrfinder output and convert to bed format for genomic comparison.
```
#remove the line with "Sequence"; get the chromosome and location; get rid of grep artifact "--" and empty lines; if character 1 from column 1 is "[" then print the second column (chromosome number) else print only the relevant locational information($3,$5,$8); get rid of "Strand:"; put chromosome and coordinates on same line; add the missing tab in the first line; get rid of all first tabs; convert spaces to tab for standard bed format.
grep -v "Sequence" LTRfinder.out |grep "Location" -B 1 |sed 's/--//g' |sed '/^$/d' |awk '{if(substr($1,1,1)=="[") {print $2} else {print $3,$5,$8}}' |sed 's/Strand://g' |tr "\n" "\t" |sed 's/+/+\n/g' |sed 's/-/-\n/g' |awk '{if(NR==1) {print "\t"$0} else {print $0}}' |sed 's/\t//1' |tr " " "\t" >LTRfinder.bed
```
