# How to Search for Positive/Negative/Neutral Selection with Codeml

Have you ever wanted to know which genes across a set of genomes are undergoing positive (diversifying) or negative (purifying) selection?  Here is a tutorial that utilizes multiple genomes and gene predictions to assess this phenomenon.  

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
mkdir 01_CleanProteins
#This removes periods and everything outside of the first column of the fasta header.
for f in *Proteins.fasta; do sed  '/^[^>]/s/\.//g' $f |awk '{print $1}' >01_CleanProteins/Clean${f};done
```

### Run orthofinder to classify orthologous groups
If you have a really large dataset, you can set orthofinder to use pre-computed BLAST/Diamond results
```
/work/gif/remkv6/Macintosh/01_AphidEffector/20_TestDataset
ml diamond;ml orthofinder; orthofinder -a 36 -t 36 -f 01_CleanProteins
```

### Extract gene family proteins and transcripts
```
/work/gif/remkv6/Macintosh/01_AphidEffector/20_TestDataset

# make a giant fasta file containing all proteins from every species
cat 01_CleanProteins/*fasta >AllCleanProteins.fasta
ml cdbfasta; cdbfasta AllCleanProteins.fasta

#make a giant fasta file containing all transcripts from every species
for f in *_Transcripts.fasta; do awk '{print $1}' $f >Clean${f};done
ml cdbfasta; cdbfasta AllCleanTranscripts.fasta

# extract proteins and transcripts by looping through the tree files and formatting to get gene names.  This can be painful at times, as orthofinder does add the fasta file name to the gene names.
for f in 01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/*tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanProteins.fasta.cidx >${f}_Proteins.fasta;done

for f in 01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/*; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanTranscripts.fasta.cidx >${f}_Transcripts.fasta;done
```

### Generate protein alignments with clustalo
```
/work/gif/remkv6/Macintosh/01_AphidEffector/20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees

ml clustalo/1.2.4-iy2mf6y
for f in *_Proteins.fasta; do echo "ml clustalo;clustalo -i "$f" -o "${f%.*}"alignment";done >generateAlignment.sh
```

### run pal2nal

The suggestion is to use pal2nal with "-nogap", though codeml has a way to deal with these. However, if you do not have sequence in your pal2nal output, likely you had gaps along the entire sequence across all of the species analyzed.
```
/work/gif/remkv6/Macintosh/01_AphidEffector/20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees

wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -zxvf pal2nal.v14.tar.gz

for f in *alignment; do echo "pal2nal.v14/pal2nal.pl  "$f" "${f%_*}"_Transcripts.fasta -output paml -nogap >"${f%_*}"_pal2nal" ; done >pal2nal.sh
```
### Understanding Codeml parameters
```
There are lots of options that can be tweaked to test for positive selection. There are many more that I did not mention.

runmode = 0 means I am specifying user trees
seqtype = 1 means I am using codon sequences
CodonFreq = 0 means I am assuming equal codon frequencies (codon substitution model)
model = 0 a single dn/ds ratio for each branch in my tree
NSsites = 2 3 test for overall selection and discrete values across branches, several of these can be specified at once.
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
### run Codeml
```
/work/gif/remkv6/Macintosh/01_AphidEffector/20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees

#generate the control files appropriately named, but on one line
for f in OG*pal2nal; do echo "seqfile = ../"$f" \t treefile = ../"${f%_*}"\t outfile = "${f%.*}".out \t noisy = 9 \t verbose = 1 \t runmode = 0 \t seqtype = 1 \t CodonFreq = 0 \t model = 0 \t NSsites = 0 1 2 3  \t icode = 0 \t fix_kappa = 0 \t kappa = 0  \t fix_omega = 0 \t omega = 0.5 \t  RateAncestor = 0 \t clock = 0 \t  Small_Diff = 5e-06 \t Malpha = 0 \t method = 0 \t ndata = 1" > ${f%.*}"codeml.ctl";done

#move orthgroup named files to new folders, codeml only works on codeml.ctl
for f in OG*ctl; do mkdir ${f}dir; mv $f ${f}dir/codeml.ctl; done

#get it out of tab format
sed -i 's/\\t/\n/g' *dir/codeml.ctl

#run codeml -- separate this script and submit to multiple nodes.  
for f in *ctldir; do echo "ml paml/4.9h-m5kkb2e; cd "$f"; codeml; cd ../ ";done >IterateCodeml.sh
```

### Which orthogroups are demonstrating positive selection
To avoid digging into each folder, I extract the dn/ds(omega) value for each orthgroup.
```
/work/gif/remkv6/Macintosh/01_AphidEffector/20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees

grep "omega" */*out |awk '$4>1' |wc
     27     108    1959
Of the 3,996 orthogroups surveyed, 27 were found to have positive selection.
```


### References
```
The PAML manual was tremendously helpful in increasing my understanding -- http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf
I endlessly queried biostars to help improve my understanding -- getSE do in codeml site:www.biostars.org
This tutorial for pairwise comparisons is what got me running initially -- https://www.protocols.io/view/introduction-to-calculating-dn-ds-ratios-with-code-qhwdt7e
```
