# How to Search for Positive/Negative/Neutral Selection with Codeml

Have you ever wanted to know which genes across a set of genomes are undergoing positive (diversifying), neutral, or negative (purifying) selection? Here is a tutorial that utilizes multiple genomes and gene predictions to assess this phenomenon.  

I modeled this analysis after a previously published work that demonstrated positive selection among E. coli genes --> https://genome.cshlp.org/content/17/9/1336/T3.expansion.html


### Prerequisite software
All of these were already installed on my HPC.
```
cufflinks/2.2.1
cdbfasta/2017-03-16
orthofinder/2.5.2
diamond/2.0.4
clustalo/1.2.4
pal2nal.v14
```


### Download your datasets
```
20_TestDataset/

GCF_000007405.1 -- shigella flexneri
GCF_001721125.1 -- E. coli -- Food-borne Pathogen Omics Research Center, FORC
GCF_003697165.2 -- E. coli -- University of Arkansas for Medical Sciences
GCF_013357365.1 -- E. coli -- USDA-ARS, Western Regional Research Center
GCF_000005845.2 -- E. coli -- Univ of Wisconsin

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/697/165/GCF_003697165.2_ASM369716v2/GCF_003697165.2_ASM369716v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/697/165/GCF_003697165.2_ASM369716v2/GCF_003697165.2_ASM369716v2_genomic.gff.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/405/GCF_000007405.1_ASM740v1/GCF_000007405.1_ASM740v1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/405/GCF_000007405.1_ASM740v1/GCF_000007405.1_ASM740v1_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/721/125/GCF_001721125.1_ASM172112v1/GCF_001721125.1_ASM172112v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/721/125/GCF_001721125.1_ASM172112v1/GCF_001721125.1_ASM172112v1_genomic.gff.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/357/365/GCF_013357365.1_ASM1335736v1/GCF_013357365.1_ASM1335736v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/357/365/GCF_013357365.1_ASM1335736v1/GCF_013357365.1_ASM1335736v1_genomic.gff.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_translated_cds.faa.gz
```

### Extract proteins and transcripts from genome and gff
I always create my own transcripts and proteins from the gff, as in rare occurrences proteins and transcripts do not concur.
```
20_TestDataset/

module load cufflinks

gffread  -g GCF_000007405.1_ASM740v1_genomic.fna GCF_000007405.1_ASM740v1_genomic.gff  -x GCF_000007405.1_Transcripts.fasta -y GCF_000007405.1_Proteins.fasta
gffread -g GCF_001721125.1_ASM172112v1_genomic.fna GCF_001721125.1_ASM172112v1_genomic.gff -x GCF_001721125.1_Transcripts.fasta -y GCF_001721125.1_Proteins.fasta
gffread -g GCF_003697165.2_ASM369716v2_genomic.fna GCF_003697165.2_ASM369716v2_genomic.gff -y GCF_003697165.2_Proteins.fasta -x GCF_003697165.2_Transcripts.fasta
gffread -g GCF_013357365.1_ASM1335736v1_genomic.fna GCF_013357365.1_ASM1335736v1_genomic.gff -y GCF_013357365.1_Proteins.fasta -x GCF_013357365.1_Transcripts.fasta
gffread -g GCF_000005845.2_ASM584v2_genomic.fna GCF_000005845.2_ASM584v2_genomic.gff -x GCF_000005845.2_Transcripts.fasta -y  GCF_000005845.2_Proteins.fasta
```

### Format the fasta files so they are easier to use
You need to make sure that the proteins from each species have a unique name.  i.e. species one and species two cannot have a protein named g123456.1.   
```
20_TestDataset/

mkdir 01_CleanProteins
#This removes periods and everything outside of the first column of the fasta header.
for f in *Proteins.fasta; do sed  '/^[^>]/s/\.//g' $f |awk '{print $1}' >01_CleanProteins/Clean${f};done
```

### Run orthofinder to classify orthologous groups
If you have a really large dataset, you can set orthofinder to use pre-computed BLAST/Diamond results
```
20_TestDataset/
ml diamond;ml orthofinder; orthofinder -a 36 -t 36 -f 01_CleanProteins
```

### Extract gene family proteins and transcripts
```
20_TestDataset/

# make a giant fasta file containing all proteins from every species
cat 01_CleanProteins/*fasta >AllCleanProteins.fasta
ml cdbfasta; cdbfasta AllCleanProteins.fasta

#make a giant fasta file containing all transcripts from every species
for f in *_Transcripts.fasta; do awk '{print $1}' $f >Clean${f};done
ml cdbfasta; cdbfasta AllCleanTranscripts.fasta

# extract proteins and transcripts by looping through the tree files and formatting to get gene names.  This can be painful at times, as orthofinder does add the fasta file name to the gene names.
for f in 01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/*tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanProteins.fasta.cidx >${f}_Proteins.fasta;done

# Do the same for the transcripts, thus it is much easier if your protein and transcript names match
for f in 01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/*; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanTranscripts.fasta.cidx >${f}_Transcripts.fasta;done
```

### Generate protein alignments with clustalo
```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

ml clustalo/1.2.4-iy2mf6y
for f in *_Proteins.fasta; do echo "ml clustalo;clustalo -i "$f" -o "${f%.*}"alignment";done >generateAlignment.sh
```

### run pal2nal

The suggestion is to use pal2nal with "-nogap", though codeml has a way to deal with these. However, if you do not have sequence in your pal2nal output, likely you had gaps along the entire sequence across all of the species analyzed.
```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -zxvf pal2nal.v14.tar.gz

for f in *alignment; do echo "pal2nal.v14/pal2nal.pl  "$f" "${f%_*}"_Transcripts.fasta -output paml -nogap >"${f%_*}"_pal2nal" ; done >pal2nal.sh
```

### Understanding Codeml parameters
```
There are lots of options that can be tweaked to test for positive selection. There are many more that I did not mention, and I only mention the ones for which I have a rudimentary understanding.

runmode = 0 means I am specifying user trees
seqtype = 1 means I am using codon sequences
CodonFreq = 0 means I am assuming equal codon frequencies (codon substitution model)
model = 0 whether ω should be allowed to vary among branches in the tree
NSsites = 0 1 2 7 8 whether ω should be allowed to vary among sites , several of these can be specified at once.
icode = 0  is specifies that I want to use the universal genetic code
fix_kappa = 0 -- initial transition/transversion ratio value set to 0
kappa = 0 -- assumption starting value for expectation of transitions/transversions
fix_omega = 1 -- fixes omega to test for significance between current and foreground branch
omega = 0.001 -- not 100% sure here, but likely the initial value set for omega in current branch
RateAncestor = 0 -- ancestral sequence not constructed
clock = 0  -- assumption of a molecular clock is not allowed, rates are free to vary from branch to branch
Small_Diff = 5e-08 -- supposed to use multiple values to make sure your results are not sensitive to this number uses 1e-8 to 1e-5.
Malpha = 0 -- when substitution rates are assumed to vary from site to site, Malpha specifies whether one gamma distribution will be applied across all sites
method = 0 -- nucleotide substitution model
ndata = 1 -- number of separate datasets
```

### Prepare Codeml control files and run Codeml
```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

#generate the control files appropriately named, but on one line
for f in OG*pal2nal; do echo "seqfile = ../"$f" \t treefile = ../"${f%_*}"\t outfile = "${f%.*}".out  \t verbose = 1 \t runmode = 0 \t seqtype = 1 \t CodonFreq = 0 \t model = 0 \t NSsites = 0 1 2 7 8 \t icode = 0 \t RateAncestor = 0 \t clock = 0 \t  \t method = 0 \t ndata = 1" > ${f%.*}"codeml.ctl";done

#move orthgroup named files to new folders, rename the control files to codeml.ctl, as codeml only works on codeml.ctl
for f in OG*ctl; do mkdir ${f}dir; mv $f ${f}dir/codeml.ctl; done

#transform the single line of tab separated values to a newline separated document
sed -i 's/\\t/\n/g' *dir/codeml.ctl

#run codeml -- Creates an sh script that can be separated to submit to multiple nodes.  
for f in *ctldir; do echo "ml paml/4.9h-m5kkb2e; cd "$f"; codeml; cd ../ ";done >IterateCodeml.sh
```

### Which orthogroups are demonstrating positive selection
To avoid digging into each folder, I extract the log-likelihood for each model tested for each orthgroup.
* ntime is the number of branch lengths
* np is the number of parameters
* The third value is your log likelihood value

```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

#this prints the file name and grabs the values we want (np and Log Liklihood)
grep --with-filename "lnL"  OG*_treecodeml.ctldir/*out |awk '{print $1,$4,$5}' |sed 's/)://g' |tr " " "\t" |awk '{ if (NR % 5 == 0) {print $0"#"} else {print}}'  |tr "\n" "\t" |sed 's/#/\n/g' |awk '{print $1,$2,$3,$5,$6,$8,$9,$11,$12,$14,$15}' |less

#example output
OG0000004_treecodeml.ctldir/OG0000004_tree.out:lnL(ntime: 61 -235.745781 62 -235.745804 64 -235.745781 62 -235.753920 64 -235.753943
OG0000008_treecodeml.ctldir/OG0000008_tree.out:lnL(ntime: 55 -918.513428 56 -917.211691 58 -917.211691 56 -916.214085 58 -916.214672
OG0000009_treecodeml.ctldir/OG0000009_tree.out:lnL(ntime: 55 -740.121154 56 -719.787244 58 -719.530605 56 -720.923654 58 -719.128922
```

### Create separate files for each test, as the number of parameters varys among and between models.
```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/
#M0 vs M1 -- Support for neutral selection can be identified if M1 provides a better fit than M0
grep --with-filename "lnL"  OG*_treecodeml.ctldir/*out |awk '{print $1,$4,$5}' |sed 's/)://g' |tr " " "\t" |awk '{ if (NR % 5 == 0) {print $0"#"} else {print}}'  |tr "\n" "\t" |sed 's/#/\n/g' |awk '{print $1,$2,$3,$5,$6,$8,$9,$11,$12,$14,$15}' |awk '{print  $1,($^C$2),(2*($5+(-1*$3))) }' |  tr " " "\t" >M1M2ChiSquares.tab

#M1 vs M2 -- Support for positive selection can be identified if M2 provides a better fit than M1
grep --with-filename "lnL"  OG*_treecodeml.ctldir/*out |awk '{print $1,$4,$5}' |sed 's/)://g' |tr " " "\t" |awk '{ if (NR % 5 == 0) {print $0"#"} else {print}}'  |tr "\n" "\t" |sed 's/#/\n/g' |awk '{print $1,$2,$3,$5,$6,$8,$9,$11,$12,$14,$15}' |awk '{if ($6>$4) {print  $1,($6-$4),(2*($7+(-1*$5)))}else {print "FLIP"$1,($4-$6),(2*($5+(-1*$7)))} }'  |  tr " " "\t" >M1M2ChiSquares.tab

#M7 vs M8 --The M8–M7 comparison offers a more stringent test for positive selection
grep --with-filename "lnL"  OG*_treecodeml.ctldir/*out |awk '{print $1,$4,$5}' |sed 's/)://g' |tr " " "\t" |awk '{ if (NR % 5 == 0) {print $0"#"} else {print}}'  |tr "\n" "\t" |sed 's/#/\n/g' |awk '{print $1,$2,$3,$5,$6,$8,$9,$11,$12,$14,$15}' |awk '{if($10>$8){print  $1, ($10-$8), (2*($11+(-1*$9)))} else {print $1,($8-$10), (2*($9+(-1*$11)))} }' |  tr " " "\t" >M7M8ChiSquares.tab
```
### Create a chi-squared test database of critical values for significance
We need the significance values for the chi square tests, so with the paml module loaded type: chi2
Then we select the degres of freedom column with the 0.0500 header name and copy it to a chi2lists file (tab separated here).  
chi2lists
```
1       3.8415
2       5.9915
3       7.8147
4       9.4877
5       11.0705
6       12.5916
7       14.0671
8       15.5073
9       16.9190
10      18.3070
11      19.6751
12      21.0261
13      22.3620
14      23.6848
15      24.9958
```

### Tests for neutral selection comparing models M0 and M1a
```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

Is model M1 better than null(M0)? In how many orthogroups?  --neutral selection
less M0M1ChiSquares.tab |awk '{print $2}' |while read line; do grep -w "$line" chi2lists ;done |paste M0M1ChiSquares.tab - |awk '$3>$5' |wc
     91     455    6998

less M0M1ChiSquares.tab |awk '{print $2}' |while read line; do grep -w "$line" chi2lists ;done |paste M0M1ChiSquares.tab - |awk '$3>$5' |sed 's/:/\t/g' |cut -f 1 |sed 's/FLIP//g' |while read line; do grep --with-filename -A 200 "Model 1" $line |grep "\*" - |grep -v "branch" - |grep -v "Positively" -;done > NeutralSelectedSitesM1.txt

#number of orthogroups with significant log ratio threshold and Naive Empirical Bayes (NEB) analysis
less NeutralSelectedSitesM1.txt |awk '{print $1}' |sort|uniq|wc
    89      89    4272

#number of neutrally selected sites among the 89 orthogroups
less NeutralSelectedSitesM1.txt |wc
       1319    7241  110317
```


### Positive selection tests comparing models M1a and M2a
```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

#How many orthogroups had significant chi squared tests for log likelihood?
less M1M2ChiSquares.tab |awk '{print $2}' |while read line; do grep -w "$line" chi2lists ;done |paste M1M2ChiSquares.tab - |awk '$3>$5' |wc
105     525    8500

less M1M2ChiSquares.tab |awk '{print $2}' |while read line; do grep -w "$line" chi2lists ;done |paste M1M2ChiSquares.tab - |awk '$3>$5' |sed 's/:/\t/g' |cut -f 1 |sed 's/FLIP//g' |while read line; do grep --with-filename -A 200 "Model 2" $line |grep "\*" - |grep -v "branch" - |grep -v "Positively" -;done > PositivelySelectedSitesM2.txt

#number of orthogroups with significant log ratio threshold and Naive Empirical Bayes (NEB) analysis
less PositivelySelectedSitesM2.txt |awk '{print $1}' |sort |uniq|wc
      101     101    4848

#number of positively selected sites among the 101 orthogroups
wc PositivelySelectedSitesM2.txt
  1572   8636 131512 PositivelySelectedSitesM2.txt
```



### Positive selection tests comparing models M7 and M8 -- a more stringent site model
```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

How many orthogroups had significant chi squared tests for log likelihood?
less M7M8ChiSquares.tab |awk '{print $2}' |while read line; do grep -w "$line" chi2lists ;done |paste M0M1ChiSquares.tab - |awk '$3>$5' |sed 's/:/\t/g' |cut -f 1 |wc
   109     109    5123

#positively selected sites -- only 106 comparisons had sites with significance
less M7M8ChiSquares.tab |awk '{print $2}' |while read line; do grep -w "$line" chi2lists ;done |paste M0M1ChiSquares.tab - |awk '$3>$5' |sed 's/:/\t/g' |cut -f 1 |while read line; do grep --with-filename -A 500 "M7" $line |grep "\*" - |grep -v "branch" - |grep -v "Positively" - ;done >PositivelySelectedSitesM8.txt

#number of orthgroups with a significant log ratio threshold and Naive Emperical Bayes (NEB) analysis
less PositivelySelectedSitesM8.txt |awk '{print $1}' |sort|uniq|wc
   106     106    5088

#number of positively selected sites among the 106 orthogroups
wc PositivelySelectedSitesM8.txt
  1560   8888 131827 PositivelySelectedSitesM8.txt
```

### How many of the M0 models show positive selection across the entire orthogroup

Some of these values are ~999.000, which usually means your Ds was 0, so they can be tossed.
```
20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

grep "omega" *dir/*out |awk '$4>1 && $4<998' >PositivelySelectedOrthogroupsM0.txt

wc PositivelySelectedOrthogroupsM0.txt
  81  324 5837 PositivelySelectedOrthogroupsM0.txt

```
### Results compared with previous a previous study
Note that this is not an exact comparison, as I used different genome's from the same species. Here I compare my results to the published results demonstrating genes under positive selection in Escherichia coli.
https://genome.cshlp.org/content/17/9/1336/T2.expansion.html

Create a text list gene names from the linked table above -- previouslypublished.txt
```
fluA
eaeH
nmpC
ubiF
ompF
ompA
ycgV
yciD
rzpR
tra8_2
yddk
pqqL
nohA
purR
yeeU
yeeV
mdtC
ompC
insA_6
rfaC
wecD
lamB
tra8_3
```
Now lets extract the gene names from the positively selected orthogroups -- including M0
```
#grabs all the orthogroup names, and uses them to concatenate all the orthogroup trees
cat PositivelySelected* |sed 's/_/\t/g' |awk '{print $1}' |sort |uniq |sed 's/$/_tree.txt/g' |tr "\n" " " |sed 's/^/cat /g' |sed 's/$/ >AllPositiveTress.txt/g' >concatenatePositiveSelectionTrees.sh
sh concatenatePositiveSelectionTrees.sh


#need to concatenate all of the gff's, as some genes in the orthogroups may not have the annotations
cat ../../../../*gff >../../../../Allgffs.gff

#this grabs the functional gene names from the concatenated gffs.
less AllPositiveTress.txt |sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2 |while read line; do grep $line ../../../../Allgffs.gff ;done |awk '$3=="gene"' |sed 's/Name=/\t/g' |sed 's/;/\t/2' |cut -f 10 |sort |uniq >E.coliGenesUnderPositiveSelection.txt

#using grep to adjust for differences in capitalization
grep -f  E.coliGenesUnderPositiveSelection.txt previouslypublished.txt |wc
    23      23     121

All 23 positively selected genes in the previously published table were confirmed here. An addition 180 orthogroups also found to demonstrate positive selection in this analysis.

```

Now lets extract the gene names from the positively selected orthogroups -- Excluding M0
```
#grabs all the orthogroup names, and uses them to concatenate all the orthogroup trees
cat PositivelySelectedSitesM2.txt PositivelySelectedSitesM8.txt |sed 's/_/\t/g' |awk '{print $1}' |sort |uniq |sed 's/$/_tree.txt/g' |tr "\n" " " |sed 's/^/cat /g' |sed 's/$/ >AllSignificantPositiveTrees.txt/g' >concatenateSignificantPositiveSelectionTrees.sh
sh concatenateSignificantPositiveSelectionTrees.sh


#this grabs the functional gene names from the concatenated gffs.
less AllSignificantPositiveTrees.txt |sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2 |while read line; do grep $line ../../../../Allgffs.gff ;done |awk '$3=="gene"' |sed 's/Name=/\t/g' |sed 's/;/\t/2' |cut -f 10 |sort |uniq >E.coliGenesUnderSignificantPositiveSelection.txt

#using grep to adjust for differences in capitalization
grep -f  E.coliGenesUnderSignificantPositiveSelection.txt previouslypublished.txt |wc
    23      23     121

All 23 positively selected genes in the previously published table were confirmed here. An addition 103 orthogroups also found to demonstrate positive selection with significance in this analysis.
```

### Further interpretation
```
Of the 3,996 orthogroups surveyed, 89 orthogroups were confirmed to have neutral selection with model M1, 101 orthogroups had positive selection with model M2, 106 orthogroups had positive selection with model M8, and 81 were found to have positive selection across the entire orthogroup with model M0.

Next steps would be to survey your significant orthogroups to identify which areas of positive/neutral selection are significant to your project. Do these positively selected proteins play a role in some adaptive function?  Do they correspond to protein-protein binding domains, protein-DNA binding domains, protein-RNA binding domains? Are these proteins anomalous annotation errors or do they have a functional definition?  There are lots of things to explore, your mind is the only limit.
```




### References

* The PAML manual was tremendously helpful in increasing my understanding -- http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf
* An earlier version of the paml manual that has some useful descriptions of input files -- https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.332.818&rep=rep1&type=pdf
* I endlessly queried biostars to help improve my understanding --  codeml site:www.biostars.org
* This tutorial for pairwise comparisons is what got me running initially -- https://www.protocols.io/view/introduction-to-calculating-dn-ds-ratios-with-code-qhwdt7e
* Here is another set of tutorials that was helpful in understanding the value to the statistical method --  https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjY2ear9IH2AhUnDzQIHc4VBaQQFnoECAIQAQ&url=https%3A%2F%2Fs3-eu-west-1.amazonaws.com%2Fpfigshare-u-files%2F5999646%2Fworked_example.pdf&usg=AOvVaw2fyo5IUrZVUOS5fk0FGr_l

* This paper was crucial in my understanding of the log ratio tests https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5015

* A codeml tutorial that helps in the understanding the output of codeml and a detailed explanation of each model -- https://lorenzogatti.me/2016_FiPS_Tutorials/solutions_tutorial02.html


### Disclaimer
I am a biologist with bioinformatics expertise, not a statistician. I do not consider myself an expert in running or interpreting codeml. That being said, I spent lots of time reading to understand this material and feel it is accurate.    
