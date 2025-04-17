---
title: Extracting and Assemble the Mitochondrial Genome from Nuclear-targeted Nanopore Reads 
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

In this tutorial, we'll walk through extracting and assembling the mitochondrial genome of Dodonaea viscosa using Oxford Nanopore reads.
We'll perform quality filtering, subset the reads, assemble multiple draft assemblies, then circularize and polish the best assembly.


# Setup: Organize your workspace

Where are your nuclear-targeted or other long reads?
```
#both nanopore runs filtered for quality, adapters trimmed, and 5k length cutoff
ln -s /work/gif3/sharu/viscosa/03_Porechop/trimmed_5kplus.fq
```
Where is your nuclear genome?
```
MitoRemovedpolished_3.fasta
```

In what location is this analyis being performed?
```
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly
```

# Find related species with an assembled mitochondrial genome

To aid in assembly and annotation, find mitochondrial genomes from related species:

  1. Search the NCBI Nucleotide database for "Sapindaceae mitochondrion genome."

  2. Filter for plants and mitochondria.

  3. Download sequences for related species. Some examples in my case:

* Xanthoceras sorbifolium
* Nephelium lappaceum
* Litchi chinensis
* Sapindus mukorossi 
  
Once you find a few assemblies of related species, you can visualize how they are related by plotting 'genus species' names into the NCBI Taxonomy Browser tool: Common Tree

https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi 


# Identifying organellar genomic reads from nuclear-targeted reads

##### Install Minimap2

Creates environment named minimap2, and calls bioconda to install minimap2. Activate environment to access minimap.
```
micromamba create -n minimap2 -c bioconda minimap2
micromamba activate minimap2
```

### Map reads to related species' mitochondrial genome
```
sh runMinimap.sh trimmed_5kplus.fq RelatedNCBIMitochondrialGenomeSequences.fasta
```
runMinimap.sh
```
##############################################################################
#!/bin/bash
query=$1
target=$2
outname="${query%.*}_${target%.*}_minimap2.paf"
module load minimap2
minimap2 -x map-ont -t 36 $target $query > ${outname}
##############################################################################
```

##### Install Seqtk

We will be using Seqtk to extract fastq sequences. This code creates an environment named seqtk, and calls bioconda to install seqtk. Activate environment to access seqtk.
```
micromamba create -n seqtk -c bioconda seqtk
micromamba activate seqtk
```

### Extract reads from mapping with at least 2kb alignment length using Seqtk

The length of the alignment filter is to get rid of any nuclear integrated mitochondrial genomic fragments. If your genome assembly has very few or small NUMTs or NUMPs then, this is likley the easiest way to obtain an organelle-targeted set of reads from a nuclear-targeted set of reads.
```
less trimmed_5kplus.fq_RelatedNCBIMitochondrialGenomeSequences_minimap2.paf |awk -F"\t" '($4-$3)>2000 {print $1}' |uniq |seqtk subseq trimmed_5kplus.fq.gz - >MitoNanopore2k.fastq
```

### As an alternative to the step above, we can map the reads to the nuclear genome and extract only those that did not map. The thought here is that if there is not too contamination in the reads, we will be able to eliminate all reads that are nuclear and keep only low-quality reads, organelle reads, and contamination reads to use for organelle genome assembly.

```
#softlink genome, reads and mapping script
ln -s ../../03_DononeaViscosaDeposition/MitoRemovedpolished_3.fasta
ln -s ../trimmed_5kplus.fq
ln -s ../runMinimap.sh

#map reads to nuclear genome
echo "sh runMinimap.sh trimmed_5kplus.fq MitoRemovedpolished_3.fasta" >Map2Nucleus.sh

#get names of all reads in the dataset
awk 'NR % 4 == 1' trimmed_5kplus.fq |sed 's/^@//g'  |awk '{print $1}' >AllTrimmedReadNamesFormatted.list

cat AllTrimmedReadNamesFormatted.list <(awk -F"\t" '$12==60 && $10>500' {print $1}' trimmed_5kplus_MitoRemovedpolished_3_minimap2.paf|uniq ) |sort |uniq -c |awk '$1==1 {print $2}' >OrganelleReads.list
seqtk subseq trimmed_5kplus.fq OrganelleReads.list >OrganelleReads.fastq
```

### As an alternative to the previous two steps, we can identify reads that have mitochondrial genes using miniprot and assemble those reads.

```
ml micromamba; micromamba activate miniprot

#convert fastq to fasta
micromamba activate seqtk; seqtk seq -a trimmed_5kplus.fq >trimmed_5kplus.fasta

# map proteins to the reads
echo "miniprot -t 64 -S trimmed_5kplus.fasta NamedSmukorossiProteins.fasta >SmukrossiProteins2OrganelleReads.paf " >miniprotToTrimmed5kplus.sh

#Grab all read names that have a mitochondrial protein alignment
less SmukrossiProteins2OrganelleReads.paf |cut -f 6 |sort |uniq >MiniprotOrganelleReads.list

#extract the reads with seqtk
seqtk subseq trimmed_5kplus.fq MiniprotOrganelleReads.list >MiniprotOrganelleReads.fastq

#note that this generated a very small number of reads that I could not split to sufficient depth. I ran this example without splitting before assembly.
```
##### Install Trycycler

Trycycler provides positional splitting of a fastq file, leading to a higher-quality assembly consensus.
```
micromamba create -c bioconda -c conda-forge -n trycycler trycycler
micromamba activate trycycler
```

### Create subsets of reads from different positions in the fastq file using Trycycler

Here trycycler will generate 12 fastq files in the read_subsetNuclearClean and read_subsets folders, that are close to equal in size and without any read overlap between samples
```
trycycler subsample --reads OrganelleReads.fastq --out_dir read_subsetNuclearClean  -t 36
trycycler subsample --reads MitoNanopore2k.fastq --out_dir read_subsets -t 36
```
Alternatively generate random subsets with Seqtk by changing the seed number.
```
# By modifying '-s' and your output filename, you can generate different subsets of 10% of the reads.   
seqtk sample -s100 MitoNanopore2k.fastq 0.1 > subset.fastq
```

# Installing software for creating assemblies from read subsets

After choosing your method to split the fastq files, we will need to use these to generate multiple independent assemblies.  Here I will be using Unicycler, Miniasm, Raven, CANU, and FLYE to generate different assemblies.

##### Install Unicycler assembler

```
git clone https://github.com/rrwick/Unicycler.git

#create environment
python -m venv Unicycler
# activate the environment
source Unicycler/bin/activate
# install Unicycler
python3 setup.py install --prefix=/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly
```

##### Install FLYE assembler

```
micromamba -n flye -c bioconda::flye
micromamba activate flye
```

##### Install Miniasm assembler 

```
 micromamba -n miniasm -c bioconda::miniasm
 micrombamba activate miniasm
```

##### Install Raven assembler
```
micromamba create -n raven -c bioconda -c conda-forge raven-assembler
micromamba activate raven
```

##### Install CANU
```
micromamba create -n canu-env -c bioconda -c conda-forge canu
micromamba activate canu-env
```

##### Install Circlator to circularize assemblies at the origin for better comparison
```
micromamba create -n circlator-env python=3.7 -y
micromamba activate circlator-env
micromamba install -c bioconda -c conda-forge circlator
```

# Generate assemblies of read subsets

**Generate FLYE Assemblies**
```
for f in *fastq; do echo "ml micromamba; micromamba activate flye; flye --nano-raw $f --threads 36 --out-dir ${f%.*}FlyeOut --genome-size 800000 --meta --scaffold" -m 2000 ; done >flye.sh
```

**Generate Miniasm Assemblies**
```
for f in *fastq; do echo "sh AssembleMitoMiniasm.sh "$f" "${f%.*}"_MiniasmOut";done >miniasmAssemblies.sh
```

AssembleMitoMiniasm.sh
```bash
#############################################################################################################################
#!/bin/bash

# Ensure script stops on errors
set -e

# Check if correct number of arguments are provided
if [[ $# -ne 2 ]]; then
    echo "Usage: sh assemble_mito.sh <reads.fastq> <output_prefix>"
    exit 1
fi

# Assign command-line arguments to variables
READS="$1"
PREFIX="$2"

# Create an output directory
OUTDIR="${PREFIX}_output"
mkdir -p "$OUTDIR"

# Check if input file exists
if [[ ! -f "$READS" ]]; then
    echo "Error: File '$READS' not found!"
    exit 1
fi

# Convert FASTQ to FASTA and save in the output directory
echo "Converting FASTQ to FASTA..."
awk 'NR%4==1 {print ">" substr($0,2)} NR%4==2 {print}' "$READS" > "$OUTDIR/${PREFIX}_reads.fasta"

# Run minimap2 for read overlap
echo "Running minimap2..."
ml micromamba;micromamba activate minimap2;minimap2 -x ava-ont -t 36 -k 19 -w5  "$OUTDIR/${PREFIX}_reads.fasta" "$OUTDIR/${PREFIX}_reads.fasta" > "$OUTDIR/${PREFIX}.paf"

# Run miniasm
echo "Running miniasm..."
ml micromamba; micromamba activate miniasm;miniasm;miniasm -f "$OUTDIR/${PREFIX}_reads.fasta" "$OUTDIR/${PREFIX}.paf" -s 500 > "$OUTDIR/${PREFIX}.gfa"

# Convert GFA to FASTA
echo "Converting GFA to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' "$OUTDIR/${PREFIX}.gfa" > "$OUTDIR/${PREFIX}.final.fasta"

echo "Assembly complete! Output files are in $OUTDIR"
#################################################################################################################################
```

**Generate CANU 2.2 assemblies**
```
for f in sample*fastq; do echo "ml micromamba; micromamba activate canu-env; canu -p mito -d "${f%.*}"CanuOut genomeSize=550k -nanopore-raw "$f ;done >canu.sh
```

**Generate Raven assemblies**
```
 for f in sample*fastq; do echo "ml micromamba; micromamba activate raven;mkdir "${f%.*}"Raven ; cd "${f%.*}"Raven ; ln -s ../"$f"; raven "$f" >"${f%.*}"Assembly.fasta" ;done  >raven.sh   
```

**Generate Unicycler Assemblies**
```
# current location
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/
source bin/activate

# create assembly commands for each split
for f in *fastq; do echo "unicycler -l "$f" -t 36 --mode bold -o "${f%.*}"Unicycler --threads 36 --racon_path /work/gif3/masonbrink/Heuther/09_Unicycler/Unicycler/racon/build/bin/racon" ;done >UnicyclerRuns.sh
```

# Evaluation of Assemblies -- find largest scaffold in each assembly with close to the expected number of genes

### Acquire sapindaceae mitochondrial genes
```
#Downloaded protein sequences from annotated Sapindus mukorossi mitochondrion
vi SmukorossiProteins.fasta
grep -c ">" SmukorossiProteins.fasta
41

#Rename proteins by protein name, not accession
awk '{if(NF>1) {print ">"$2} else {print $0}}' SmukorossiProteins.fasta >NamedSmukorossiProteins.fasta
```
**Get all of the genomes into one folder with unique file names -- Reads mapped to related mitochondrial assembly**
```
#grabs all the unicycler assemblies and renames them by folder name
 for f in ../sample*Unicycler/assembly.fasta; do echo "ln -s "$f ; echo $f |sed 's|/|\t|g' |cut -f 2| awk '{print $1".fasta"}' ;done |tr "\n" " " |sed 's/ln -s/\nln -s/g' >SoftlinkUnicycler.sh
sh SoftlinkUnicycler.sh

#grabs all flye assemblies and renames by folder name
for f in ../read_subsets/sample*FlyeOut/assembly.fasta; do echo "ln -s "$f ; echo $f |sed 's|/|\t|g' |cut -f 3| awk -F"\t" '{print $0".fasta"}' ;done |tr "\n" " " |sed 's/ln -s/\nln -s/g' >SoftlinkFlye.sh

#softlinks all the miniasm assemblies which already have unique names
for f in ../read_subsets/*MiniasmOut_output/*final.fasta ; do ln -s $f;done

#softlinks and renames all canu assemblies 
for f in ../read_subsets/*CanuOut/mito.contigs.fasta; do echo "ln -s "$f ; echo $f |sed 's|/|\t|g' |cut -f 3| awk -F"\t" '{print $0".fasta"}' ;done |tr "\n" " " |sed 's/ln -s/\nln -s/g' >CanuSoftlink.sh

#softlinks and renames all Raven assemblies, which already have unique names
for f in ../read_subsets/*Raven/*Assembly.fasta ;do ln -s $f;done
```
**Get all of the genomes into one folder with unique file names -- Reads that could not map to an organelle-less nuclear genome**
```
#grabs all the unicycler assemblies and renames them by folder name
 for f in ../read_subsetNuclearClean/*Unicycler/assembly.fasta; do echo "ln -s "$f ; echo $f |sed 's|/|\t|g' |cut -f 3| awk '{print $1"NuclearClean.fasta"}' ;done |tr "\n" " " |sed 's/ln -s/\nln -s/g' >SoftlinkUnicycler2.sh
sh SoftlinkUnicycler.sh

#grabs all flye assemblies and renames by folder name
for f in ../read_subsetNuclearClean/*FlyeOut/assembly.fasta; do echo "ln -s "$f ; echo $f |sed 's|/|\t|g' |cut -f 3| awk -F"\t" '{print $0"NuclearClean.fasta"}' ;done |tr "\n" " " |sed 's/ln -s/\nln -s/g' >SoftlinkFlye2.sh

#softlinks all the miniasm assemblies which do not have unique names, and needed renamed this time
 for f in ../read_subsetNuclearClean/*MiniasmOut_output/*final.fasta; do echo "ln -s "$f ; echo $f |sed 's|/|\t|g' |cut -f 3| awk -F"\t" '{print $0"NuclearClean.fasta"}' ;done |tr "\n" " " |sed 's/ln -s/\nln -s/g' >SoftlinkMiniasm2.sh

#softlinks and renames all canu assemblies # all of these failed due to memory exhaustion
for f in ../read_subsetsNuclearClean/*CanuOut/mito.contigs.fasta; do echo "ln -s "$f ; echo $f |sed 's|/|\t|g' |cut -f 3| awk -F"\t" '{print $0"NuclearClean.fasta"}' ;done |tr "\n" " " |sed 's/ln -s/\nln -s/g' |less

#softlinks and renames all Raven assemblies so they do not overwrite the ones above
for f in ../read_subsets/*Raven/*Assembly.fasta ;do ln -s $f ${f%.*}NuclearClean.fasta;done
```

### Align mitochondrial proteins to assemblies

```
#Align proteins to the genomes.
for f in *.fasta; do miniprot $f NamedSmukorossiProteins.fasta >${f%.*}.genes; done

#how many of these 109 assemblies (not all shown here) had most of the genes on how many contigs and how large was the longest contig with mapping genes? Note I forced only assemblies with at least 35 unique mitochondrial gene alignments to print below
for f in *genes; do paste <(ls -1 $f) <(cut -f 1 $f |sort|uniq|wc -l) <(cut -f 6 $f|sort|uniq|wc -l)  <(cut -f 7 $f |sort -k1,1nr |awk '{print $1}' ) ;done |awk '$2>35' |less                                     

1kClean60qualAssembly.genes     38      3       339819
2kRawAssembly.genes     38      1       536045
3kRawAssembly.genes     38      1       536052
5kCleanAssembly.genes   38      2       463238
assembly.genes  38      1       535684
CleanReads3kAssembly.genes      38      1       535684
CleanReads5kAssembly.genes      38      2       462995
Flye5kRawAssembly.genes 38      4       251593
FlyeMiniprotAssembly.genes      38      4       275765
RawReads5kAssembly.genes        37      2       372648
sample_01FlyeOut.genes  36      2       347764
sample_02Assembly.genes 39      4       470606
sample_02FlyeOut.genes  38      4       180048
sample_03Assembly.genes 39      3       463133
sample_04Assembly.genes 39      4       449804
sample_04FlyeOut.genes  38      5       257488
sample_05Assembly.genes 39      5       412423
sample_05FlyeOut.genes  37      4       205704
sample_06Assembly.genes 39      4       412399
sample_06FlyeOut.genes  37      3       225437
sample_07Assembly.genes 39      3       528267
sample_07FlyeOut.genes  38      2       381058
sample_08Assembly.genes 40      6       316564
sample_08FlyeOut.genes  38      5       123612
sample_09Assembly.genes 39      4       522057
sample_09FlyeOut.genes  38      5       229821
sample_10Assembly.genes 39      5       463042
sample_10FlyeOut.genes  36      4       128482
sample_11Assembly.genes 39      6       277342
sample_11FlyeOut.genes  38      6       163011
sample_12Assembly.genes 40      6       412414
sample_12FlyeOut.genes  36      5       252031
UCManualInterventionRemoveNUMTReads.genes       38      2       462795
```
So it looks like I commonly find 38-39 genes in these assemblies, but only two had the expected 40 genes. I ran a few assemblies with varying levels of filtering using all of the reads too, so those are included and some happen to be really good assemblies. The fewest contigs with the most genes seems to occur with assembly.fasta, which is an assembly using all trimmed reads that were able to get a 3kb alignment to the related Smukorossi mitochondrial genome. However, thee are a couple assemblies that have all 40 genes, which we may use to incorporate into the singular contig at a later step. 

### Best Mitochondrial assembly refinement

```
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/assemblies

#rename the assembly so I do not lose it
cp ../MitoNanopore3kRawUnicycler/assembly.fasta MitoNanopore3kRawUnicyclerassembly.fasta

#just a standard blast to self
makeblastdb -in MitoNanopore3kRawUnicyclerassembly.fasta -dbtype nucl
blastn -query MitoNanopore3kRawUnicyclerassembly.fasta -db MitoNanopore3kRawUnicyclerassembly.fasta -outfmt 6 -num_threads 36 -out MitoNanopore3kRawUnicyclerassembly.blastout
```

**Blast Results: truncated to show the regions that are large and identical**
```
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
1       1       100.000 535684  0       0       1       535684  1       535684  0.0     9.892e+05
1       1       99.851  43743   25      21      298590  342311  46555   2832    0.0     80382
1       1       99.851  43743   25      21      2832    46555   342311  298590  0.0     80382
1       1       99.969  31916   0       5       503771  535683  427604  459512  0.0     58874
1       1       99.969  31916   0       5       427604  459512  503771  535683  0.0     58874
1       1       99.984  6318    0       1       450654  456971  170106  163790  0.0     11660
1       1       99.984  6318    0       1       163790  170106  456971  450654  0.0     11660
1       1       99.968  6318    0       2       526823  533139  170106  163790  0.0     11655
1       1       99.968  6318    0       2       163790  170106  533139  526823  0.0     11655
1       1       99.640  2502    3       5       343015  345511  2501    1       0.0     4566
1       1       99.640  2502    3       5       1       2501    345511  343015  0.0     4566

So there are three duplicate regions that can be condensed to two regions. Since two are very close, I will just merge them
 1-2501, 2832-46555, 503771-535683
```

**Trim regions from the end of the assembly that are duplicate**
```
ml samtools
#indexes the fasta file
samtools faidx MitoNanopore3kRawUnicyclerassembly.fasta
#These coordinates were used to extract the fasta from a blast to self identifying duplicated regions of the genome
echo "1:46555-503771" |samtools faidx MitoNanopore3kRawUnicyclerassembly.fasta -r - >TrimmedCleanMito.fasta
```

**Check for the presence of all genes in the trimmed assembly**
```
#map the proteins without 
miniprot -S -j 0 TrimmedCleanMito.fasta NamedSmukorossiProteins.fasta >GenesInTrimmedCleanMito.paf
```

**Can we add the missing two genes that were found in other assemblies?**
```
This concatenates the names of the genes in each assembly, finds the ones that are only in one assembly, and extracts those scaffolds
cat <(less ../assemblies/GenesInTrimmedCleanMito.paf |cut -f 1|sort|uniq ) <(cut -f 1 sample_08Assembly.genes |cut -f 1|sort|uniq) |sort|uniq -c |awk '$1==1{print $2}' |grep -w -f - sample_08Assembly.genes |cut -f 6 |cdbyank sample_08Assembly.fasta.cidx >MissingMitoScaffolds.fasta

#create database and run a standard self BLASTn
makeblastdb -in TrimmedCleanMito.fasta -dbtype nucl
blastn -db TrimmedCleanMito.fasta -query ../Assemblies/MissingMitoScaffolds.fasta -outfmt 6 -num_threads 8 -out MissingMitoScaffolds2TrimmedCleanMito.blastout                                                    

So there are two genes that were recognized as missing from the assembly
orf119
rps1

However, these genes map to two different 90kb scaffolds in both of the assemblies that were able to align these genes.  Aligning these 90kb scaffolds to the TrimmedCleanMito.fasta assembly and the Sapindus mukrossi mitochondrial genome did not provide any hits, suggesting that these two genes may be NUMT genes integrated into the nuclear genome and can safely be ignored.
```

**Circularize with circlator**
```
#This reorients the fasta to start at the origin
micromamba activate circlator
circlator fixstart TrimmedCleanMito.fasta ReOriented
```

**Align reads to the circularized mitochondrial assembly** 
```
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler

ln -s ../../runMinimapNbamSort.sh
#use minimap to align reads to the genome to produce a sorted bam file.
sh runMinimapNbamSort.sh MitoNanopore3kRaw.fastq ReOriented.fasta

#runMinimap.sh
##############################################################################
#!/bin/bash
query=$1
target=$2
outname="${query%.*}_${target%.*}_minimap2.sam"
module load minimap2
minimap2 -x asm5 -a -t 36 $target $query > ${outname}

ml samtools;samtools view --threads 36 -b -o ${outname%.*}.bam ${outname}
samtools sort  -o ${outname%.*}_sorted.bam -T TEMP --threads 36 ${outname%.*}.bam
samtools index ${outname%.*}_sorted.bam
##############################################################################
```

**Run Pilon on the genome with the aligned reads**
```
#run pilon with the sorted bam file and the genome to ensure the assembly is of high quality
echo " ml pilon ; java -Xmx190g -Djava.io.tmpdir=TEMP -jar  /opt/rit/el9/20230413/app/linux-rhel9-x86_64_v3/gcc-11.2.1/pilon-1.22-iojt4j5x62smcab6j4mjn77ejfqimebo/bin/pilon-1.22.jar --genome ReOriented.fasta --unpaired MitoNanopore3kRaw_ReOriented_minimap2_sorted.bam --outdir /work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/assemblies/  --changes --fix all --threads 36   ">pilon.sh 
```
**Results of assembly and polishing**
```
                                         Number of scaffolds          1
                                     Total size of scaffolds     457217
                                            Longest scaffold     457217
                                           Shortest scaffold     457217
                                 Number of scaffolds > 1K nt          1 100.0%
                                Number of scaffolds > 10K nt          1 100.0%
                               Number of scaffolds > 100K nt          1 100.0%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size     457217
                                        Median scaffold size     457217
                                         N50 scaffold length     457217
                                          L50 scaffold count          1
                                         n90 scaffold length     457217
                                          L90 scaffold count          1
                                                 scaffold %A      27.60
                                                 scaffold %C      22.62
                                                 scaffold %G      22.30
                                                 scaffold %T      27.48
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0
```

# Find the chloroplast assemblies in our mitochondrial assemblies

### Best Chloroplast assembly

With all of the assemblies you we've made trying to assemble the mitochondria, we most likely assembled a chloroplast genome too, as they are much simpler than plant mitochondria. Lets align the previously published chloroplast genome to the assemblies. 

### Find the chloroplast assembly

```
#download the chloroplast genome 
NC_036099.1 Dodonaea viscosa chloroplast, complete genome -- 159,375bp

#map the published chloroplast assembly to the assemblies we've made
for f in *fasta; do sh runMinimap.sh DviscosaChloroplast.fasta $f;done

Here we can see which contigs were of the appropriate size and complete via how much genome they aligned
cat *minimap2.paf |awk '$7>150000 && $7<200000'|less

This will grab only those minimap alignments that are in the same orientation as the published assembly. There were way too many, so filtered for even better assemblies
cat *minimap2.paf |awk '$4>157000 && $3<20'|less 

In my case, I had 46 likely complete chloroplast assemblies, but only four were in the same orientation and started near the same base
NC_036099.1     159375  5       159366  -       tig00000001     194930  8636    167738  139848  159617  60      tp:A:P  cm:i:25896      s1:i:139710     s2:i:80342      dv:f:0.0081     rl:i:0
NC_036099.1     159375  5       158489  +       Utg177960       268827  56996   215246  138870  158751  60      tp:A:P  cm:i:25672      s1:i:138733     s2:i:108725     dv:f:0.0090     rl:i:962
NC_036099.1     159375  5       159366  +       tig00000001     190469  30999   190097  139817  159615  60      tp:A:P  cm:i:25895      s1:i:139681     s2:i:72177      dv:f:0.0081     rl:i:0
NC_036099.1     159375  5       159366  -       Utg177424       252033  5114    147499  140664  159651  60      tp:A:P  cm:i:25994      s1:i:138502     s2:i:71811      dv:f:0.0086     rl:i:960

Column 7 tells us the length of the contig the chloroplast genome hit to, and it appears that either this genome is much larger than the published version or it contains a false amount of duplication.  The easiest way to identify if this is the case is to do a self-blastn to the contig
```
**BLASTn of Utg177960 to Utg177960**
```
#index and extract the chloroplast contig
cdbfasta sample_06Assembly.fasta
echo "Utg177960" |cdbyank sample_06Assembly.fasta.cidx >MyChloroplast.fasta

#create the blast database
makeblastdb -in MyChloroplast.fasta -dbtype nucl
#run the self blast
blastn -db MyChloroplast.fasta -query MyChloroplast.fasta -outfmt 6 -out Mychloroplast2mychloroplast.blastout

#How much of the assembly is identical? This shows the top best hit which is everything and then the next best hits which are reciprocal best hits
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
Utg177960       Utg177960       100.000 268827  0       0       1       268827  1       268827  0.0     4.964e+05
Utg177960       Utg177960       99.980  56113   0       7       159146  215251  1       56109   0.0     1.036e+05
Utg177960       Utg177960       99.980  56113   0       7       1       56109   159146  215251  0.0     1.036e+05
Utg177960       Utg177960       99.983  52625   1       7       215252  267872  144613  197233  0.0     97123
Utg177960       Utg177960       99.983  52625   1       7       144613  197233  215252  267872  0.0     97123

So there are two duplicate regions from 1-56109 and 215252-267872, lets extract the remaining sequence
echo "Utg177960 56109 215252" |cdbyank sample_06Assembly.fasta.cidx -R |bioawk -c fastx '{print $name,length($seq)}' >UnpolishedMyChloroplastDeduplicated.fasta
```
**Circularize with circlator**
```
ml micromamba; micromamba activate circlator
circlator fixstart UnpolishedMyChloroplastDeduplicated.fasta MyReOrientedChloroplastGenome
```
**Polish chloroplast assembly**
```
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/ChloroHuge

ln -s ../runMinimapNbamSort.sh
sh runMinimapNbamSort.sh MitoNanopore3kRaw.fastq ReOrientedChloroplastGenome.fasta
echo " ml pilon ; java -Xmx190g -Djava.io.tmpdir=TEMP -jar  /opt/rit/el9/20230413/app/linux-rhel9-x86_64_v3/gcc-11.2.1/pilon-1.22-iojt4j5x62smcab6j4mjn77ejfqimebo/bin/pilon-1.22.jar --genome MyReOrientedChloroplastGenome.fasta --unpaired MitoNanopore3kRaw_ReOrientedChloro_minimap2_sorted.bam --outdir /work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/  --changes --fix all --threads 36  ">pilonAllChloro.sh
```

**Results of assembly and polishing**
```
new_Assemblathon.pl MyReOrientedChloroplastGenomePilon.fasta

---------------- Information for assembly 'MyReOrientedChloroplastGenome.fasta' ----------------


                                         Number of scaffolds          1
                                     Total size of scaffolds     159144
                                            Longest scaffold     159144
                                           Shortest scaffold     159144
                                 Number of scaffolds > 1K nt          1 100.0%
                                Number of scaffolds > 10K nt          1 100.0%
                               Number of scaffolds > 100K nt          1 100.0%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size     159144
                                        Median scaffold size     159144
                                         N50 scaffold length     159144
                                          L50 scaffold count          1
                                         n90 scaffold length     159144
                                          L90 scaffold count          1
                                                 scaffold %A      30.65
                                                 scaffold %C      19.17
                                                 scaffold %G      18.73
                                                 scaffold %T      31.45
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0
```

