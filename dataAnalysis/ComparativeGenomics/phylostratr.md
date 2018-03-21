#  Determining evolutionary origin of all genes in a genome

`phylostratr` is an R package for generating phylostratigraphy of all genes in a genome. You will need complete predicted proteins translated from the transcripts for a given genome.


## Installation

Following instructions will work on any `apt` based package manager systems such as Ubuntu or Mint. If you have any other system, you can use the `dockerfile` and run them as container (in the future).

```
# install R
sudo apt install r-base-core
# this may or may not be the latest version
# check R version
R --version
# if less than 3.4.2, continue with updating R
```

### Update R

First, you need to know what Ubuntu you are running. For that:
```
lsb_release -a
# note down the codename from the output, eg: artful
sudo su
echo "deb http://www.stats.bris.ac.uk/R/bin/linux/ubuntu artful/" >> /etc/apt/sources.list
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
apt-get update
apt-get upgrade
R --version
# check R version again to make sure you have the latest version
```
You can proceed, if you have the R version > 3.4.3

```
R version 3.4.4 (2018-03-15) -- "Someone to Lean On"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under the terms of the
GNU General Public License versions 2 or 3.
For more information about these matters see
http://www.gnu.org/licenses/.
```

### Install dependencies

Again, since this is a apt based system, you can easily install them as follows:
```
sudo apt-get install virtualbox-guest-dkms
sudo apt-get install -y libcurl-dev
sudo apt-get install libssl-dev
sudo apt-get install libxml2-dev
```
For R packages, you need to first get to the R terminal

```
R
```
From there, install packages:
```
install.packages("curl")
install.packages("devtools")
install.packages("taxize")
install.packages("reshape2)
install.packages("taxizedb")
install.packages("dplyr")
install.packages("readr")
install.packages("phylostratr")
install.packages("magrittr")
install.packages("devtools)
# phylostratr
install_github("arendsee/phylostratr")
```

With this, you will be set to run `phylostratr` analysis.

### Running `phylostratr`

We will run `phylostratr` in 3 parts. In the first part, we will find the taxonomy id for the species of interest from NCBI and choose all the species for determining the gene age. Second, we will run the BLAST (all possible 1 vs. 1). We will do this in a HPC cluster to speed up the process. Last, we will finish running the `phylostratr`, writing results in a csv file and save plots as pdf.

#### 1. Setting up focal species and downloading protein sequences

We will use _Arabidopsis thaliana_ as an example. To find the Taxonomy id, you can look up the scientific name on [NCBI taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy). The ID is `3702`

```
library(reshape2)
library(taxizedb)
library(dplyr)
library(readr)
library(phylostratr)
library(magrittr)
focal_taxid <- '3702'
strata <-
# Get stratified relatives represented in UniProt
uniprot_strata(focal_taxid, from=2) %>%
# Select a diverse subset of 5 or fewer representatives from each stratum.
strata_apply(f=diverse_subtree, n=5, weights=uniprot_weight_by_ref()) %>%
# Use a prebuilt set of prokaryotic species
use_recommended_prokaryotes %>%
# Add yeast and human proteomes
add_taxa(c('4932', '9606')) %>%
# Download proteomes, storing the filenames
uniprot_fill_strata

# Replace the UniProt Arabidopsis thaliana proteome with the Araport11 proteome.
# You will need to replace the filename with the path to your own file. You can
# use the UniProt genes, however UniProt includes multiple versions of each gene,
# raising the total to 89135 genes. Since A. thaliana is the focal species, it
# is important to be very explicit about version of the proteome.
strata@data$faa[['3702']] <- 'Araport11_genes.201606.pep.fasta'
```

You will find the `uniprot-seqs` directory in your workdir, with many files identified as `<taxaid>.faa`. We will take this for performing BLAST on HPC.

#### 2. BLAST on HPC


We can run the BLAST directly inside the  `faa` directory. First you need to generate BLASTdb files for each `faa` file. This can be done by running `makeblastdb` on all of them:

```
# request a compute node
salloc -N 1 -n 16 -t 1:00:00
# load the blast module
module load ncbi-blast
# run makeblastdb
for faa in *.faa; do
makeblastdb -in $faa -dbtype prot -parse_seqids -out ${faa%.*};
done
```
This will take few minutes to complete. Obtain the path for this folder
```
pwd -P
```
save the path for this folder for the next step

We will then set up a `runBLASTp.sh` script as follows (with the path to blastdb from last step):

```
#!/bin/bash
daname="$2"
infile="$1"
outfile="$(basename "${infile%.*}")_${daname}.out"
database="/work/LAS/path/todb/${daname}"
module load ncbi-blast
blastp \
 -query "${infile}" \
 -db "${database}" \
 -out "${outfile}" \
 -num_threads 2 \
 -outfmt "6 qseqid sseqid qstart qend sstart send evalue score"
```
Generate commands for each blast run. Assuming `3702.faa` as your focal taxa-id, you can do this:
```
for faa in *.faa; do
echo "./runBLASTp.sh 3702.faa ${faa%.*}";
done > blastp.cmds
```

Now use the `makeSLURMp.py` script to generate parallel jobs. Count the number of lines in `blastp.cmds` and divide it by 10 (since 10 is the max number of jobs you can run on Condo), use that as second argument:

```
makeSLURMp.py NN blastp.cmds
```
the script can be found [here](https://github.com/ISUgenomics/common_scripts/blob/master/makeSLURMp.py)

To submit all the jobs, you can use the loop again.

```
for sub in blast_*sub; do
sbatch $sub;
done
```
The results will be written to a tab delimited file, named as `3702_NNNNN.out`, where NNNN is the taxa-id. This requires further modification to make it work with `phylostratr`. You will need to add the taxa-id as the last column and header for each field.

```
for out in *.out; do
taxaid=$(echo ${out%.*} |cut -f 2 -d "_");
echo -e "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tscore\tstaxid" > ${taxaid}.tab
awk '{print $0"\t"$'{taxaid}'}' ${out} >> ${taxaid}.tab;
done
```
Finally, move the results to the workdir so the `phylostratr` can see the results

```
mv *.tab ../
```

#### 3. Generate phylostratigraphy results

Finally, load the blast results, startify  and write results/plots.

```
results <- strata_blast(strata, blast_args=list(nthreads=8)) %>%
  strata_besthits %>%
  merge_besthits
phylostrata <- stratify(results)
write.csv(phylostrata, "phylostrata_table.csv")
tabled <- table(stratify(results)$mrca_name)
write.csv(tabled, "phylostrata_stats.csv")
# Get metadata. Note, this will take awhile
strata <- strata %>%
# for each species, add map of UniProt IDs to PFAM domains
strata_uniprot_pfam_map %>%
# for each species, add proteome summary data
add_proteome_stats %>%
# for each species, add list of organelle encoded UniProt proteins
add_organelle_proteins
prot <- proteome_stats_table(strata)
strata2 <- strata_convert(strata, target='all', to='name')
# You can explore the proteome stats interactively with plotly:
ggplotly(plot_proteome_stats(strata2))
ggplotly(plot_proteome_lengths(strata2))
proteome_report_table(strata2) %>%
merge(get_phylostrata_map(strata2)) %>%
arrange(-ps) %>%
select(species, ps, N, min, q25, median, q75, max) %>%
kable()
# generate heatmaps
plot_heatmaps(results, heatmaps.pdf, focal_id = 381124)
plot(strata)
```
