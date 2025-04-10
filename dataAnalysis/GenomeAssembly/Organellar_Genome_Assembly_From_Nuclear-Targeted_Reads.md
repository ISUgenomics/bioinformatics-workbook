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
  
Once you find a few assemblies of related species, you can visualize how they are related by plotting 'genus species' names into the NCBI Tasonomy Browser tool: Common Tree

https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi 


# Identifying organellar genomic reads from nuclear-targeted reads

##### Install Minimap2
Creates environment named minimap2, and calls bioconda to install minimap2. Activate environment to access minimap.
```
micromamba create -n minimap2 -c bioconda minimap2
micromamba activate minimap2
```



### Map reads to mitochondrial genomes of related species
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

### Install Trycycler

Trycycler provides positional splitting of a fastq file, leading to a higher-quality assembly consensus.
```
micromamba create -c bioconda -c conda-forge -n trycycler trycycler
micromamba activate trycycler
```

### Create subsets of reads from different positions in the fastq file using Trycycler
```
trycycler subsample --reads MitoNanopore2k.fastq --out_dir read_subsets  -t 36
```
Alternatively generate random subsets with Seqtk.
```
# By modifying '-s' and your output filename, you can generate different subsets of 10% of the reads.   
seqtk sample -s100 MitoNanopore2k.fastq 0.1 > subset.fastq
```

# Installing software for creating assemblies from read subsets

After choosing your method to split the fastq files, we will need to use these to generate multiple independent assemblies.  Here I will be using Unicycler, Miniasm, and FLYE to generate different assemblies.

### Install Unicycler assembler

```
git clone https://github.com/rrwick/Unicycler.git

#create environment
python -m venv Unicycler
# activate the environment
source Unicycler/bin/activate
# install Unicycler
python3 setup.py install --prefix=/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly
```

### Install FLYE assembler

```
micromamba -n flye -c bioconda::flye
micromamba activate flye
```

### Install Miniasm assembler 

```
 micromamba -n miniasm -c bioconda::miniasm
 micrombamba activate miniasm
```

### Install Circlator to circularize assemblies at the origin for better comparison
```
micromamba create -n circlator-env python=3.7 -y
micromamba activate circlator-env
micromamba install -c bioconda -c conda-forge circlator
```
# Generate assemblies of read subsets

**Generate FLYE Assemblies**

```
#current location
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/
#softlinks fastq files in read_subsets/ folder to the current folder
for f in read_subsets/*fastq; do ln -s $f;done


for f in *fastq; do echo "ml micromamba; micromamba activate flye; flye --nano-raw $f --threads 36 --out-dir ${f%.*}FlyeOut --genome-size 800000 --meta --scaffold" -m 2000 ; done >flye.sh
```

**Generate Miniasm Assemblies**

```
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/read_subsets
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
###################################################################################################################################
```

**Generate Unicycler Assemblies**
```
# current location
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/
source bin/activate

# create assembly commands for each split
for f in *fastq; do echo "unicycler -l "$f" -t 36 --mode bold -o "${f%.*}"Unicycler --threads 36 --racon_path /work/gif3/masonbrink/Heuther/09_Unicycler/Unicycler/racon/build/bin/racon" ;done >UnicyclerRunAll.sh
```

# Evaluation of Assemblies -- find largest scaffold in each assembly with close to the expected number of genes

```


#Extracted proteins from Mitos2 annotation of a previous mitochondrial genome
less sequence.txt |sed 's/>lcl|.*\[gene=/>/g' |sed 's/].*//g' >NamedMitoProteins.fasta

#Align proteins to the genomes -- should be 38 total. Chloroplast has a previously published genome for comparison, so this was not needed for the chloroplast

for f in *assembly.fasta; do miniprot $f NamedMitoProteins.fasta >${f%.*}.genes; done

for f in *genes; do paste <(ls -1 $f) <(cut -f 1 $f |sort|uniq|wc) ;done

#Grab size of the largest scaffold in the assembly
for f in *fasta; do  paste <(ls  -1 *fasta) <(new_Assemblathon.pl $f |grep "Longest scaffold");done
```

### Best Mitochondrial assembly 
```
cp ../MitoNanopore3kRawUnicycler/assembly.fasta MitoNanopore3kRawUnicyclerassembly.fasta

makeblastdb -in MitoNanopore3kRawUnicyclerassembly.fasta -dbtype nucl
blastn -query MitoNanopore3kRawUnicyclerassembly.fasta -db MitoNanopore3kRawUnicyclerassembly.fasta -outfmt 6 -num_threads 36 -out MitoNanopore3kRawUnicyclerassembly.blastout

samtools faidx MitoNanopore3kRawUnicyclerassembly.fasta
#These coordinates were used to extract the fasta from a blast to self identifying duplicated regions of the genome
echo "1:41000-503000" |samtools faidx MitoNanopore3kRawUnicyclerassembly.fasta -r - >TrimmedCleanMito.fasta
```
**Circularize with circlator**
```
#This reorients the fasta to start at the origin
micromamba activate circlator
circlator fixstart TrimmedCleanMito.fasta ReOriented
```

**Polish the mitochondrial assembly** 
```
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler

ln -s ../../runMinimapNbamSort.sh
sh runMinimapNbamSort.sh MitoNanopore3kRaw.fastq ReOriented.fasta
echo " ml pilon ; java -Xmx190g -Djava.io.tmpdir=TEMP -jar  /opt/rit/el9/20230413/app/linux-rhel9-x86_64_v3/gcc-11.2.1/pilon-1.22-iojt4j5x62smcab6j4mjn77ejfqimebo/bin/pilon-1.22.jar --genome ReOriented.fasta --unpaired    MitoNanopore3kRaw_ReOriented_minimap2_sorted.bam --outdir /work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/assemblies/  --changes --fix all --threads 36   ">pilon.sh 
```
**Results of assembly and polishing**
```
                                         Number of scaffolds          1
                                     Total size of scaffolds     462647
                                            Longest scaffold     462647
                                           Shortest scaffold     462647
                                 Number of scaffolds > 1K nt          1 100.0%
                                Number of scaffolds > 10K nt          1 100.0%
                               Number of scaffolds > 100K nt          1 100.0%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size     462647
                                        Median scaffold size     462647
                                         N50 scaffold length     462647
                                          L50 scaffold count          1
                                         n90 scaffold length     462647
                                          L90 scaffold count          1
                                                 scaffold %A      27.62
                                                 scaffold %C      22.63
                                                 scaffold %G      22.27
                                                 scaffold %T      27.48
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0
```

### Best Chloroplast assembly
```
cd sample_11Unicycler/

#grabbed the longest contig and did an online blast to see duplication positions to remove.

#extracts just the chloroplast contig sections that are not duplicated
vi List
####################################################################################################################################
1:1-67868
1:94909-158707
####################################################################################################################################
samtools faidx assembly.fasta -r List >Chloroplast.fasta
#Connect the two contigs
vi Chloroplast.fasta
#rewrap the fasta lines
tr "\n" "\t" <Chloroplast.fasta |sed 's/\t/#/1' |sed 's/\t//g' |sed 's/#/\n/g' >MyChloroplast.fasta
```
**Circularize with circlator**
```
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler
ln -s sample_11Unicycler/MyChloroplast.fasta
#reorient the fasta to start at the origin
ml micromamba
micromamba activate circlator
circlator fixstart MyChloroplast.fasta ReOrientedChloro
```
**Polish chloroplast assembly**
```
/work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/ChloroHuge]

ln -s ../runMinimapNbamSort.sh
sh runMinimapNbamSort.sh MitoNanopore3kRaw.fastq ReOrientedChloro.fasta
echo " ml pilon ; java -Xmx190g -Djava.io.tmpdir=TEMP -jar  /opt/rit/el9/20230413/app/linux-rhel9-x86_64_v3/gcc-11.2.1/pilon-1.22-iojt4j5x62smcab6j4mjn77ejfqimebo/bin/pilon-1.22.jar --genome ReOrientedChloro.fasta --unpaired MitoNanopore3kRaw_ReOrientedChloro_minimap2_sorted.bam --outdir /work/gif3/masonbrink/USDA/04_MitochondrialIsolationAndAssembly/Unicycler/  --changes --fix all --threads 36  ">pilonAllChloro.sh
```

**Results of assembly and polishing**
```
                                         Number of scaffolds          1
                                     Total size of scaffolds     131997
                                            Longest scaffold     131997
                                           Shortest scaffold     131997
                                 Number of scaffolds > 1K nt          1 100.0%
                                Number of scaffolds > 10K nt          1 100.0%
                               Number of scaffolds > 100K nt          1 100.0%
                                 Number of scaffolds > 1M nt          0   0.0%
                                Number of scaffolds > 10M nt          0   0.0%
                                          Mean scaffold size     131997
                                        Median scaffold size     131997
                                         N50 scaffold length     131997
                                          L50 scaffold count          1
                                         n90 scaffold length     131997
                                          L90 scaffold count          1
                                                 scaffold %A      32.01
                                                 scaffold %C      18.20
                                                 scaffold %G      18.74
                                                 scaffold %T      31.06
                                                 scaffold %N       0.00
                                         scaffold %non-ACGTN       0.00
                             Number of scaffold non-ACGTN nt          0
```

