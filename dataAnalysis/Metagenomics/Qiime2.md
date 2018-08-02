---
title: "Amplicon analysis with QIIME2"
excerpt: "An example workflow using QIIME2 version 2017.7"
layout: single
author: "Adam Rivers"

---
By Adam Rivers, Designed from the official [QIIME2 tutorials](https://docs.qiime2.org/2017.7/tutorials/)
{% include toc %}


**Why QIIME 2?** There are a number of great software packages for general amplicon analysis.
Some examples include [Mothur](https://www.mothur.org/), [Phyloseq](https://joey711.github.io/phyloseq/), [Dada2](https://benjjneb.github.io/dada2/), [UPARSE](http://www.drive5.com/uparse/) and [QIIME 1](http://qiime.org/).
The most widely used software may be QIIME 1. QIIME 1 is a collection of
custom tools and wrappers around other software that makes it easy to customize amplicon
analysis, but that flexibility sometimes makes it hard to track the provenance
of data or be sure you are doing the right thing. QIIME 2 has a very
different model for data analysis that wraps data and information about that data
into one object, which addresses some of the prior shortcomings. QIIME 2 also
incorporates a major advance that has happened in the last year: the use of
exact "Sequence Variants" (SV) rather than "Operational Taxonomic Units" (OTU).
Finally QIIME 2 still has a great development team behind it and is poised to
become one of the primary amplicon analysis methods.  For that reason we are
teaching QIIME 2 while it is still in its pre-release stage (Don't you feel hi-tech?).
{: .notice--info}


# Data sets
For this tutorial We will be looking at data from this paper:
> Castrillo, G., Teixeira, P.J.P.L., Paredes, S.H., Law, T.F., de Lorenzo, L., Feltcher, M.E., Finkel, O.M., Breakfield, N.W., Mieczkowski, P., Jones, C.D., Paz-Ares, J., Dangl, J.L., 2017. Root microbiota drive direct integration of phosphate stress and immunity. Nature 543, 513–518. [doi:10.1038/nature21417](https://dx.doi.org/10.1038/nature21417)

The authors knocked out the function of different genes involved in the phosphate stress response of *Arabidopsis thaliana*. In two experiments they grew the mutants in natural soil, sequenced the 16S V4 region and found that the disruption of phosphate regulation genes altered the microbial community of the plants.  We will be analyzing this data in several different ways and comparing our results.

If you are running this tutorial on the Ceres computer cluster the the data is available at:
 ```bash
 # Fastq files
 /project/microbiome_workshop/amplicon/data/dangldatasubsampled

 # Metadata mapping file
  /project/microbiome_workshop/amplicon/data/mapping.txt

 # A Manifest file containing the names, locations and
 # orientation of the read files for import to QIIME2
   /project/microbiome_workshop/amplicon/data/manifest.csv
 ```

 The data has been modified from the archived data files by combining technical replicates, removing samples not used in the first experiment in the paper and subsampling each sample down to 10,000 reads to speed up analysis during this tutorial. A total of 146 samples will be analyzed.

 The the whole unprocessed dataset can be downloaded from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/data/view/PRJEB15671). Be sure to download the "submitted files" not the "processed files" or the filename will not match with the metadata file.


# Connecting to Ceres

Ceres is the computer cluster for the USDA Agricultural Research Service's SCInet computing environment. From Terminal or Putty (for Windows users) create a secure shell connection to Ceres
```bash
ssh -o TCPkeepAlive=yes -o ServeraliveInterval=20 -o ServerAliveCountMAx=100 <user.name>@login.scinet.science
```
Once you are logged into Ceres you can request access to an interactive node.
In a real analysis you would create a script that runs all the commands in sequence and submit the script through a program called Sbatch, part of the computer's job scheduling software named Slurm.

To request access to an interactive node:
```bash
# Request access to one node of the cluster
# using the queue "short"
# to see available queues use the command "sinfo"
salloc -p short -N 1 -n 40 -t 06:00:00

# Set up Xvfb graphics software so that QIIME2 can generate figures.
# Xvfb creates a virtual X session. This needs to be run in the background.
# If you encounter an error saying the session already exists, then use a
# different session number, e.g. "Xvfb :2"
Xvfb :1 -screen 0 800x600x16 &

# Add a Display variable to the local environment variables
export DISPLAY=:1.0

# Load the QIIME2 module
module load qiime2/gcc/64/2.2017.7
```
When you are done at the end of the tutorial end your session like this.
```bash
# to log off shut down the graphics window
killall Xvfb
# sign out of the allocated node
exit
# sign out of Ceres head node
exit
```

# Set up your working directory

```bash
# In your homespace or other desired location, make a
# directory and move into it
mkdir qiime2-phosphate-tutorial
cd qiime2-phosphate-tutorial
```

# Understanding QIIME2 files
QIIME2 uses two different file types that contain the data and metadata from an analysis: ```.qza``` files are data files while ```.qzv``` files are visualizations. The authors of QIIME2 call these data files "data artifacts" to indicate that they are objects containing data and metadata about an experiment. It is really just a zip file containing a specially formatted directory with data and metadata. You can see what type of data is contained in a data file with the command ```qiime tools peek filename.qza```. All of QIIME2 files can be viewed using an online browser that is available at [https://view.qiime2.org](https://view.qiime2.org). ```.qza``` files will contain basic info (name, universally unique identifier, data type and data format) as well ad a graph of data provenance.```.qzv``` files will contain all of that and graphic visualizations.  The raw data in these files can be accessed using the command ```qiime tools export```.

# Import your paired-end sequences
For this project the reads were sequences using Illumina paired-end, 250 base pair reads with forward and reverse reads in separate files. The fastq is imported in to a QIIME2 data artifact ending in ```.qza```

```bash
time qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /project/microbiome_workshop/amplicon/data/manifest.csv \
  --output-path demux.qza \
  --source-format PairedEndFastqManifestPhred33
```
Time to run: 2 minutes

What's this `time` thing? You can add the `time` command to any command line task
 to see how long it took to run.
 {: .notice--info}

Output:
* ```demux.qza```

# Examine the quality of the data
We can view the characteristics of the dataset and the quality scores of the data by creating a QIIME2 visualization artifact.

```bash
time qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
 ```
 Time to run: 1 minute

 Output:
 * ```demux.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdemux.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/demux.qzv)

This will create a visualization file. You can download the file to your local computer. From a new terminal window on your local computer copy the file:

```bash
scp <user.name>@login.scinet.science:/path/to/data .
```

Now you can view the file on your local computer using the [QIIME2 visualization server](https://view.qiime2.org).  Alternatively you can view the precomputed file on that server using the button above.

When viewing the data look for the point in the forward and reverse reads where quality scores decline below 25-30. We will need to trim reads to this point to create high quality sequence variants.

# Selecting Sequence Variants

The process of selecting sequence variants is the core processing step in amplicon analysis. This takes the place of "OTU picking" a method of clustering similar data together that was the common method for dealing with sequencing errors until last year.  Three different methods have been published to select sequence variants, [Dada2](https://dx.doi.org/10.1038/nmeth.3869) uses and statistical error correction model, [Deblur](https://dx.doi.org/10.1128/mSystems.00191-16) takes an information theoretic approach and [UNOISE2](https://doi.org/10.1101/081257) applies a heuristic. Each of these methods attempt to remove or correct reads with sequencing errors and then remove chimeric sequences originating from different DNA templates.
For the next step you can select either the Dada2 method or the Deblur method.

**A note on parallel processing**. Both Dada2 and Deblur can independently process each sample. This becomes very important as experiments grow to thousands of samples. In this tutorial we are taking advantage of " course grained parallelism" provided by the multiple cores in our processor. However, Deblur and Dada2 (after the error model learning step) can select sequence variants in a sample totally independently, making the problem "embarrassingly parallel." This means that samples can be split into many different files and processed on arbitrarily many machines simultaneously.
{: .notice--info}

Sequence variant selection is the slowest step in the tutorial. For that reason it is best to submit this step using the SLURM Sbatch scheduler.

## Option 1: Dada2 (Slower)

```bash
time qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --o-table table-dada2 \
  --o-representative-sequences rep-seqs-dada2 \
  --p-trim-left-f 9 \
  --p-trim-left-r 9 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 200 \
  --p-n-threads 40 \
  --p-n-reads-learn 200000
```
To submit this command using Sbatch:
```bash
sbatch /project/microbiome_workshop/amplicon/example/qiime2-phosphate-tutorial/dada2.sh
```
Time to run: 35 minutes

Output:
* ```rep-seqs-dada2.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Frep-seqs-dada2.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/rep-seqs-dada2.qza)
* ```table-dada2.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftable-dada2.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/table-dada2.qzv)

## Option 2: Deblur (Faster)
Deblur only uses forward reads at this time. You could get around this by merging your data with an outside tool like [BBmerge](http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) then importing your data as single ended. For simplicity, in this tutorial we will just use the forward reads.

```bash
 time qiime deblur denoise-16S \
   --i-demultiplexed-seqs demux.qza \
   --p-trim-length 220 \
   --output-dir deblurresults \
   --p-jobs-to-start 36
```

To submit this command using Sbatch:
```bash
sbatch /project/microbiome_workshop/amplicon/example/qiime2-phosphate-tutorial/deblur.sh
```
Time to run: 4 minutes

Output:
* ```deblurresults/representative_sequences.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Frepresentative_sequences.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/representative_sequences.qza)
* ```deblurresults/stats.qza```
[View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Fstats.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/stats.qza)
* ```deblurresults/table.qza```
[View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fdeblurresults%2Ftable.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/deblurresults/table.qza)

Okay, we have just done the hard part of amplicon sequence analysis.  At this point we have our BIOM count table, the representative sequence variants and a stats file for Deblur.

We have just called sequence variants two different ways. In a real workflow you would only use one method.  From here on out we will use the output of dada2 only: ```table-dada2.qza```.  

# Adding metadata and examining count tables

```bash
time qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv \
  --m-sample-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt
  ```
  Time to run: 30 seconds

 Output:
 * ```table-dada2.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftable-dada2.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/table-dada2.qzv)

# Phylogenetics
There are a number of diversity metrics like unifrac distance that require the construction of a phylogenetic tree.

## Multiple sequence alignment
First Mafft is used to align the sequences

```bash
time qiime alignment mafft \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza
```
Time to run: 1 minute

Output:
* ```aligned-rep-seqs.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Faligned-rep-seqs.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/aligned-rep-seqs.qza)

## Masking sites
Some sites in the alignment are not phylogenetically informative. These sites are masked.

```bash
time qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza
```
Time to run: 1 minute

Output:
* ```masked-aligned-rep-seqs.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fmasked-aligned-rep-seqs.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/masked-aligned-rep-seqs.qza)

## Creating a tree
Fastree is used to generate a phylogenetic tree from the masked alignment.
```bash
time qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza
```
Time to run: 1 minute

Output:
* ```unrooted-tree.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Funrooted-tree.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/unrooted-tree.qza)

## Midpoint rooting
Fastree creates an unrooted tree. We can root the tree at it's midpoint with this command:
```bash
time qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
 ```
Time to run: 5 seconds

Output:
* ```rooted-tree.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Frooted-tree.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/rooted-tree.qza)

# Taxonomic analysis
Sequence variants are of limited usefulness by themselves. Often we are interested in what kinds of organisms are present in our sample, not just the diversity of the sample.  To identify these sequence variants two things are needed:  a reference database and an algorithm for identifying the sequence using the database.

The primary databases are:

Database | Description | License
---------|-------------|--------
[Greengenes](http://greengenes.secondgenome.com/) | A curated database of archaea and bacteria - static since 2013 | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US)
[Silva](https://www.arb-silva.de/) | The most up-to-date and extensive  database of prokaryotes and eukaryotes, several versions | [Free academic](https://www.arb-silva.de/silva-license-information) / [Paid commercial license](http://www.ribocon.com/silva_licenses)
[The RDP database](https://rdp.cme.msu.edu/) | A large collection of archaeal bacterial and fungal sequences | [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/deed.en_US)
[UNITE](https://unite.ut.ee/) | The primary database for fungal ITS and 28S data | Not stated

There are several methods of taxonomic classification available. The most commonly used classifier is the [RDP classifier](https://rdp.cme.msu.edu/classifier/classifier.jsp). Other software includes [SINTAX](http://www.drive5.com/usearch/manual/cmd_sintax.html) and [16S classifier](http://metabiosys.iiserb.ac.in/16Sclassifier/). We will be using the QIIME2's built-in naive Bayesian classifier (which is built on Scikit-learn but similar to RDP), noting that the method, while fast and powerful, has a tendency  [over-classify](http://www.drive5.com/usearch/manual/tax_err.html) reads.

There are two steps to taxonomic classification: [training the classifier](https://docs.qiime2.org/2017.7/tutorials/feature-classifier/) (or using a [pre-trained](https://docs.qiime2.org/2017.7/data-resources/) dataset) and classifying the sequence variants.  Generally it is best to train the classifier on the exact region of the 16S, 18S or ITS you sequenced.
For this tutorial we will be using a classifier model trained on the Silva 99% database trimmed to the V4 region.

```bash
time qiime feature-classifier classify-sklearn \
  --i-classifier  /project/microbiome_workshop/amplicon/data/taxonomy/gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.qza
```
Time to run: 4 minutes

Output:
* ```taxonomy.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxonomy.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxonomy.qza)

```bash
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
  ```
Time to run: 1 second

* ```taxonomy.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxonomy.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxonomy.qzv)


Create a bar plot visualization of the taxonomy data:
```bash
qiime taxa barplot \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt \
  --o-visualization taxa-bar-plots.qzv
```
Time to run: 1 minute

Output:
* ```taxa-bar-plots.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxa-bar-plots.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxa-bar-plots.qzv)


Looking at the the ```taxonomy.qzv``` file using https://view/qiime2.org We can see the data presented at different taxonomic levels and grouped by different experimental factors. If we drill down to taxonomic level 5 something looks a bit odd. There's lots of "Rickettsiales;f__mitochondria".  This is really  plant mitochondrial contamination. Some of these samples also have chloroplast contamination.  This kind of Taxonomic filtering isn't available in QIIME2 yet but it can be be done manually.

## Filtering contaminants
By viewing the ```taxonomy.qzv``` file in the browser we can easily search for the sequence ID's that do or do not match "mitochondria or chloroplast". We can also do this via the command line:

First create a list of all sequence variant ids that are not mitochondria or chloroplasts.

```bash
# Export taxonomy data to tabular format
qiime tools export --output-dir taxonomy-export taxonomy.qza

# search for matching lines with grep then select the id column
grep -v -i "mitochondia|chloroplast|Feature" taxonomy-export/taxonomy.tsv | cut  -f 1 > no-chloro-mito-ids.txt
```

Convert our Qiime data artifact to the underlying [biom](http://biom-format.org/) file
```bash
# Export data to biom format
qiime tools export --output-dir dada2-table-export table-dada2.qza
# Move into the directory
cd dada2-table-export

# Convert the HDF5 biom file to a tsv biom file
biom subset-table \
  --input-hdf5-fp feature-table.biom \
  --axis observation \
  --ids ../no-chloro-mito-ids.txt \
  --output-fp feature-table-subset.biom

# Create a new QIIME2 data artifact with the filtered Biom file
qiime tools import \
  --input-path feature-table-filtered.biom \
  --output-path ../table-dada2-filtered.qza \
  --type FeatureTable[Frequency]

cd ..
```
Time to run: 2 minutes

Output:
* ```table-dada2-filtered.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftable-dada2-filtered.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/table-dada2-filtered.qza)


Since we have altered the qza file we can create a new bar plots:

```bash
qiime taxa barplot \
  --i-table table-dada2-filtered.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt \
  --o-visualization taxa-bar-plots-filtered.qzv

```
Time to run: 1 minute

Output:
* ```taxa-bar-plots-filtered.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Ftaxa-bar-plots-filtered.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/taxa-bar-plots-filtered.qzv)


# Alpha and Beta diversity analysis
For mostly historical reasons one of the first questions that amplicon sequencing was used for was to look at within sample and between sample ecological diversity alpha and beta diversity. Diversity metrics and collector's curves were generated for for samples to understand how environmental factors affected diversity.

```bash
time qiime diversity core-metrics \
  --i-table table-dada2-filtered.qza \
  --i-phylogeny rooted-tree.qza \
  --p-sampling-depth 4000 \
  --output-dir core-diversity
```

Output (all files are in the directory ```core-diversity```):
* ```bray_curtis_distance_matrix.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fbray_curtis_distance_matrix.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/bray_curtis_distance_matrix.qza)
* ```bray_curtis_pcoa_results.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fbray_curtis_pcoa_results.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/bray_curtis_pcoa_results.qza)
* ```evenness_vector.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fevenness_vector.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/evenness_vector.qza)
* ```faith_pd_vector.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Ffaith_pd_vector.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/faith_pd_vector.qza)
* ```jaccard_distance_matrix.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fjaccard_distance_matrix.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/jaccard_distance_matrix.qza)
* ```jaccard_pcoa_results.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fjaccard_pcoa_results.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/jaccard_pcoa_results.qza)
* ```observed_otus_vector.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fobserved_otus_vector.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/observed_otus_vector.qza)
* ```shannon_vector.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fshannon_vector.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/shannon_vector.qza)
* ```unweighted_unifrac_distance_matrix.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Funweighted_unifrac_distance_matrix.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/unweighted_unifrac_distance_matrix.qza)
* ```unweighted_unifrac_pcoa_results.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Funweighted_unifrac_pcoa_results.qza)  [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/unweighted_unifrac_pcoa_results.qza)
* ```weighted_unifrac_distance_matrix.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fweighted_unifrac_distance_matrix.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/weighted_unifrac_distance_matrix.qza)
* ```weighted_unifrac_pcoa_results.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcore-diversity%2Fweighted_unifrac_pcoa_results.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/core-diversity/weighted_unifrac_pcoa_results.qza)

 Each of these output files can be examined using other ```qiime diversity``` subcommands.


# Ordination
Ordination is a dimensionality reduction technique that enables the visualization of sample differences. QIIME has a plugin called emperor that calculates a Bray-Curtis dissimilarity matrix and uses principal coordinates analysis (PCoA). you could also export the pcoa data and plot it yourself in the package of your choice.

 ```bash
 qiime emperor plot --i-pcoa core-diversity/bray_curtis_pcoa_results.qza --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping.txt --o-visualization pcoa-visualization.qzv
 ```
 Time to run: 15 seconds

 Output:
 * ```pcoa-visualization.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fpcoa-visualization.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/pcoa-visualization.qzv)


# Differential abundance of sequence variants
If you are doing an experimental manipulation rather than just observing an environment you will likely have an experimental design with treatments and want to know which bacteria respond to these treatments. The best way to go this is an active area of applied statistics research.

The problem is challenging for several reasons:
* The data is compositional so the abundance of each taxa affects every other taxa.
* The data is over-dispersed count data, fitting (arguably) a negative binomial model.
* The data is sparse and some 0's mean a taxa is not present while other zeros mean an organism is present at a level below the limit of detection for the sequences sampled.
* Each library is sampled to a different depth so the issue of how to standardize the data comes up (simply dividing by the sum does not work).

**To rarefy or not to rarefy?** It is a common but controversial practice to downsample count data to the lowest count in your dataset to get around the issue of differential sequencing depth. In their paper titled "Waste not want not, why rarifying microbiome data is inadmissible" [McMurdie et al. (2014)](https://doi.org/10.1371/journal.pcbi.1003531) point out that this is a large waste of data and statistical power, and advocate for using differential expression software like DESeq2 that uses special normalizations and a negative binomial distribution to model data. The software uses a generalized linear model so it has a very flexible experimental design interface.  [Weiss et al. (2017)](https://doi.org/10.1186/s40168-017-0237-y) argue that the assumptions underling both the normalization and the distribution used by DEseq2 and other normalization methods are inappropriate for microbiome data. [Ancom](https://dx.doi.org/10.3402%2Fmehd.v26.27663) uses a zero inflated Gaussian model but only allows for simple two-way comparisons not richer statistical models. Gneiss [(paper)](http://doi.org/10.1128/mSystems.00162-16), [(tutorial)](https://forum.qiime2.org/t/gneiss-tutorial/932) is currently the only compositional method available in QIIME2 (support for ANCOM was dropped).  The only thing everyone agrees on is that agrees you can't just do a T-test.
{: .notice--warning}



## Gneiss analysis
Gneiss applies a method for compositional data borrowed from geology to "sidestep" the question of the absolute changes of sequences and "instead look at the balance between particular subsets of the microbial community," [(Morton et al. 2017)](https://doi.org/10.1128/mSystems.00162-16). For more on the concept of balances see this [post](https://github.com/biocore/gneiss/blob/master/ipynb/balance_trees.ipynb). Sequence variants are hierarchically clustered based on environmental gradients or co-occurrence. Then the isometric log ratios are calculated and compared for subsets of taxa. While individual taxa cannot be compared, groups of taxa responding to environmental effects can be compared. Statistical tests of for differences can also be applied.

Since we don't have a particular environmental gradient that is structuring lets start by doing a hierarchical clustering based on co-occurrence patterns.

First add pseudocounts to each cell in the matrix. This is done so that log transformations can be taken across the table.

```bash
time qiime gneiss add-pseudocount \
    --i-table table-dada2-filtered.qza \
    --p-pseudocount 1 \
    --o-composition-table composition.qza
```
Time to run: 6 seconds

Output:
* ```composition.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fcomposition.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/composition.qza)


Perform [Ward's agglomerative clustering](https://arxiv.org/abs/1111.6285)
```bash    
time qiime gneiss correlation-clustering \
    --i-table composition.qza\
    --o-clustering hierarchy.qza
```
Time to run: 5 minutes

Output:
* ```hierarchy.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fhierarchy.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/hierarchy.qza)

A tree has now been generated that can be used for making comparisons of sample groups.

Calculate the isometric log transforms on each internal node of the tree

```bash
time qiime gneiss ilr-transform \
    --i-table composition.qza \
    --i-tree hierarchy.qza \
    --o-balances balances.qza
```
Time to run: 15 seconds

Output:
* ```balances.qza``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fbalances.qza) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/balances.qza)

The balances are normally distributed an can now be analyzed using mixed linear
models  We can perform a regression on the three categorical data types, Genotype, Fraction (soil or endophytic compartment) or Soil).  Themodel explains about 10% of the total variation at all nodes of the trees. This is typical for these complex experiments.  The amount that can be explained increases as we move up the covariance tree. Overall the most predictive factor is Genotype which is encouraging.

```bash
time qiime gneiss ols-regression \
    --p-formula "Genotype+Soil+Fraction" \
    --i-table balances.qza \
    --i-tree hierarchy.qza \
    --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping3.txt \
    --o-visualization regression_summary.qzv
```

Output:
* ```regression_summary.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fregression_summary.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/regression_summary.qzv)


One of the assumptions if the ordinary least squares model is that the fixed factors are random, in other words the authors randomly arrived at the genotypes they knocked out. Of course that's not true, the genotypes were selected because they had an impact on the phosphorus stress response.  They are Fixed factors. A mixed linear model can account for fixed and random factors and effects. Gneiss offerers a linear mixed model regression too but the interface seems to be in development  so there is not much I can say about it but we can try it now. Statistical modeling is done by the [statsmodels](http://www.statsmodels.org/stable/mixed_linear.html) python package.

```bash
qiime gneiss lme-regression \
  --p-formula "Genotype" \
  --i-table balances.qza \
  --i-tree hierarchy.qza \
  --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping3.txt \
  --p-groups Soil \
  --o-visualization linear_mixed_effects_model.qzv \
```

Output:
* ```linear_mixed_effects_model.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Flinear_mixed_effects_model.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/linear_mixed_effects_model.qzv)


We can look at the most statistically significant balances and examine what taxa make up those partitions.

```bash
qiime gneiss balance-taxonomy \
    --i-balances balances.qza \
    --i-tree hierarchy.qza \
    --i-taxonomy taxonomy.qza \
    --p-taxa-level 2 \
    --p-balance-name 'y0' \
    --m-metadata-file /project/microbiome_workshop/amplicon/data/mapping3.txt  \
    --m-metadata-category Genotype \
    --o-visualization y0_taxa_summary.qzv
```


Output:
* ```y0_taxa_summary.qzv``` [View](https://view.qiime2.org/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2FMicrobiome-workshop%2Fassets%2Fqiime%2Fy0_taxa_summary.qzv) \| [Download](https://usda-ars-gbru.github.io/Microbiome-workshop/assets/qiime/y0_taxa_summary.qzv)


In this case the y0 balance is a split between samples that have plants in them and raw soil. It makes sense that this is the largest effect.  What happens if you run balance y2 or decrease the taxonomic level?

# Other analyses

This tutorial covered a range of analyses that can be done with microbiome data but there are other types on analyses that can be done too.

* Functional analysis - Several packages attempt to impute function from taxonomy including [PiCrust](https://picrust.github.io/picrust/), [Tax4fun](http://tax4fun.gobics.de/), [Piphillin](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0166104)
* Inferring ecological interaction networks -[SPIEC-EASI](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004226). Co-Variance is an issue, but SPIEC-EASI attempts to model conditional independence.
* Data management tools - [Qiita](https://qiita.ucsd.edu/), and from ARS scientist Dan Manter [myPhyloDB](http://www.myphylodb.org/)
* Set analysis - From ARS Scientist Devin Coleman-Derr [MetaComet](https://probes.pw.usda.gov/MetaCoMET/)
* Other general analysis tools - [Mothur](https://www.mothur.org/) and the R-based [Phyloseq](https://joey711.github.io/phyloseq/)
