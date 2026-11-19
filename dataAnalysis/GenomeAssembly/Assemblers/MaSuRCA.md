---
title: "Hybrid Genome Assembly with MaSuRCA: Illumina Paired-End and Long Reads"
layout: single
author: Rick Masonbrink
author_profile: true
last_modified_at: 2026-11-19
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

## Background

Main software website links:

- [MaSuRCA GitHub repository and manual](https://github.com/alekseyzimin/masurca)
- [MaSuRCA 4.1.4 release and installation](https://github.com/alekseyzimin/masurca/releases/tag/v4.1.4)
- [MaSuRCA paper: The MaSuRCA genome assembler](https://doi.org/10.1093/bioinformatics/btt476)
- [Hybrid assembly with the MaSuRCA mega-reads algorithm](https://doi.org/10.1101/gr.213405.116)

MaSuRCA combines accurate Illumina reads with long reads to assemble a genome. It first builds super-reads from the Illumina data, uses these to construct corrected long sequences called mega-reads, and then assembles those sequences. This tutorial uses **Illumina paired-end reads and PacBio long reads**, with an alternative configuration for Oxford Nanopore reads. No mate-pair libraries are needed.

The instructions target **MaSuRCA 4.1.4**, the latest release listed when this page was revised on November 19, 2026. Use the [release archive](https://github.com/alekseyzimin/masurca/releases/tag/v4.1.4) together with the [MaSuRCA manual and installation instructions](https://github.com/alekseyzimin/masurca).

This workflow is particularly useful for projects with existing Illumina data and older, less accurate long reads. The PacBio reads used here are continuous long reads (CLR), not HiFi reads. For a project with HiFi data, also consider the workbook's [hifiasm tutorial](https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/Assemblers/hifiasm_tutorial.html#gsc.tab=0).

## Assembly workflow

A standard hybrid run performs these operations internally:

1. Correct Illumina read errors and construct super-reads.
2. Use the super-reads to build corrected long sequences called mega-reads.
3. Assemble the mega-reads using the selected backend; this example uses Flye.

MaSuRCA first generates an `assemble.sh` script from a configuration file. Running that script performs the assembly and saves intermediate files that can support restarting an interrupted run.

## How to run MaSuRCA

### Parameter considerations

| Parameter or input | Purpose |
| --- | --- |
| `PE = pe 350 150 R1 R2` | Names the library `pe`, gives its fragment mean and standard deviation, and specifies the two FASTQ files. |
| `PACBIO` | Supplies the combined PacBio CLR file. |
| `GRAPH_KMER_SIZE = auto` | Lets MaSuRCA choose the super-read graph k-mer size. |
| `USE_LINKING_MATES = 0` | Uses the long-read hybrid configuration; this setting concerns links derived from paired-end reads. |
| `USE_GRID = 0` | Keeps the pipeline inside the single allocated node. |
| `LHE_COVERAGE = 25` | Limits long-read input to approximately 25× of the longest reads. |
| `MEGA_READS_ONE_PASS = 0` | Uses two passes of mega-read construction. |
| `NUM_THREADS = 64` | Matches the CPU allocation below. |
| `JF_SIZE = 2400000000` | Uses 20 times the 120 Mb planning estimate for the Jellyfish hash size; this is not a RAM limit. |
| `FLYE_ASSEMBLY = 1` | Assembles corrected mega-reads with MaSuRCA's bundled Flye backend. |

Start with the hybrid settings below and adjust them using read coverage, read lengths, genome size, and available resources. MaSuRCA's manual specifies at least 50× Illumina coverage for the hybrid workflow and advises against the Flye backend below 15× long-read coverage.

Fragment mean and standard deviation refer to the DNA fragments, not individual read lengths. The 350 bp mean and 75 bp standard deviation used below are provisional estimates for this exercise.

### Running MaSuRCA

Check the module name and version on your cluster:

```bash
module avail masurca
# If your cluster installation is out of date like the one below, proceed with the local installation below. 
module load masurca/4.1.0
command -v masurca
masurca --help
```

If a suitable module is unavailable or outdated, use the release installation below which does not require root privileges. 

```bash
mkdir 10_masurca
cd 10_masurca

curl -fL --retry 3 \
  -o MaSuRCA-4.1.4.tar.gz \
  https://github.com/alekseyzimin/masurca/releases/download/v4.1.4/MaSuRCA-4.1.4.tar.gz

tar -xzf MaSuRCA-4.1.4.tar.gz
cd MaSuRCA-4.1.4
./install.sh > install.log 2>&1
cd ../

export PATH="MaSuRCA-4.1.4/bin:$PATH"
command -v masurca
masurca --help

#make sure you have Flye in your path also, which requires python 2.7.
ml miniconda
conda create -n masurca_runtime -c conda-forge python=2.7
conda activate masurca_runtime
export PATH=MaSuRCA-4.1.4/Flye/bin/:$PATH

flye --version 
```

The developers advise against Bioconda installations of MaSuRCA because dependency conflicts can break the pipeline. Run the assembly with its bundled tools, outside environments that override them.

### Acquire sample data

We will retain *Arabidopsis thaliana* Landsberg erecta (Ler) as the example. The runs below belong to **BioSample SAMEA5751963** in **BioProject PRJEB31147**, providing a consistent sample source for both technologies.

| Run | Data type | Archived bases | Approximate coverage using 120 Mb |
| --- | --- | ---: | ---: |
| [ERR3624574](https://trace.ncbi.nlm.nih.gov/Traces?run=ERR3624574) | Illumina paired-end | 9.70 Gb | 81× |
| [ERR3415825](https://www.ebi.ac.uk/ena/browser/view/ERR3415825) | PacBio CLR subreads | 3.72 Gb | 31× |
| [ERR3415826](https://trace.ncbi.nlm.nih.gov/Traces?run=ERR3415826) | PacBio CLR subreads | 4.56 Gb | 38× |

Together, the two PacBio runs provide approximately 69× coverage. Coverage is calculated as total sequenced bases divided by the assumed genome size. 

#### Download the reads

Download the runs using SRA Toolkit and extract FASTQ files. Run these commands from the project directory. 

```bash
ml sratoolkit/3.2.1
mkdir fastq
cd fastq
prefetch ERR3624574 ERR3415825 ERR3415826 --max-size u
mkdir tmp
fasterq-dump ERR3624574/ERR3624574.sra  --split-files -t tmp -e 8 
fasterq-dump ERR3415825/ERR3415825.sra  -t tmp -e 8 
fasterq-dump ERR3415826/ERR3415826.sra  -t tmp -e 8 
```

Combine the two long-read runs because MaSuRCA expects one long-read input file.

```bash
cat ERR3415825.fastq ERR3415826.fastq >Ler_pacbio.fastq
rm ERR3415825.fastq ERR3415826.fastq
#pigz does gzip with multiple threads
pigz -p 24 *fastq
cd ../
```

Supply the original paired-end FASTQ files. MaSuRCA performs its own Illumina error correction; routine trimming or external error correction can affect its super-read construction. 

### Genome assembly with *A. thaliana* Illumina and PacBio reads

#### Create the configuration file

Move to the assembly directory and generate `sr_config.txt`.

**sr_config.txt**
```text
# REQUIRED: Begin the input-data section.
DATA

# REQUIRED: Illumina library definition.
# Fields: unique library ID, mean fragment length, standard deviation, R1, R2.
# The 350 bp mean and 150 bp SD are assumptions; replace with measured estimates.
# These numbers are not the sum of your Illumina read lengths, it is the fragment size that was sequenced + the intervening unsequenced part of the fragment.
PE = pe 350 150 fastq/ERR3624574_1.fastq.gz fastq/ERR3624574_2.fastq.gz

# REQUIRED FOR THIS HYBRID WORKFLOW: PacBio reads combined into one file.
# Long reads are optional for MaSuRCA generally.
PACBIO = fastq/Ler_pacbio.fastq.gz

# REQUIRED: End the input-data section.
END

# REQUIRED: Begin the assembly-parameter section.
PARAMETERS

# OPTIONAL: Automatically select the super-read graph k-mer size.
GRAPH_KMER_SIZE = auto

# OPTIONAL: Recommended hybrid setting with sufficient long-read coverage.
# set this to 0 if you have more than 15x coverage by long reads, 1 for illumina-only assemblies.
USE_LINKING_MATES = 0

# OPTIONAL: Run within the current machine or allocated compute node.
# Setting 1 enables distributed mega-read construction to multiple nodes and requires GRID_ENGINE and GRID_QUEUE settings.
# SGE versus SLURM is selected separately with GRID_ENGINE.
USE_GRID = 0

# OPTIONAL: Retain approximately 25x coverage of the longest long reads.
LHE_COVERAGE = 25

# OPTIONAL: Use two mega-read passes; 1 selects the faster one-pass approach.
MEGA_READS_ONE_PASS = 0

# REQUIRED: Number of CPU threads; match your scheduler allocation.
NUM_THREADS = 64

# REQUIRED: Jellyfish hash size, not a memory limit in bytes.
# This value is 20 times the assumed 120 Mb genome size.
JF_SIZE = 2400000000

# Enable Flye assembly of corrected mega-reads since we have more than 15x coverage of long reads.
FLYE_ASSEMBLY = 1

# REQUIRED: End the assembly-parameter section.
END
```


MaSuRCA's manual specifies at least 50× Illumina coverage for this hybrid approach and advises against the Flye backend below 15× long-read coverage. 

### Execute MaSuRCA assembly 

```
#This command is quick and uses your config generated above to create the assemble script below.
masurca sr_config.txt > masurca.configure.log 2>&1

#Run inside an allocated compute node, or include this command in your scheduler job script.
bash assemble.sh > masurca.assembly.log 2>&1
```

## Expected files generated during assembly

For the PacBio hybrid assembly above, look for the files below. Exact intermediate names can vary with the configuration and software version.

| File or directory | What it contains |
| --- | --- |
| `sr_config.txt` | Input paths and assembly parameters. |
| `assemble.sh` | Generated script that runs the assembly pipeline. |
| `masurca.configure.log` | Configuration checks and script-generation messages. |
| `masurca.assembly.log` | Assembly progress, diagnostic messages, and completion information. |
| `pe.cor.fa` | Error-corrected Illumina reads, an intermediate product. |
| `work1/` | Intermediate files used to construct super-reads. |
| `mr.*` | Mega-read construction outputs and intermediate files. |
| `flye.mr.*/assembly.scaffolds.fasta` | Final scaffold FASTA for the selected Flye backend. |

For the Flye backend, MaSuRCA documents the final scaffold FASTA as `assembly.scaffolds.fasta` within a `flye.mr.*` directory. Locate it and check the final messages in the assembly log:

```bash
find . -type f -path './flye.mr.*/assembly.scaffolds.fasta' -print
tail -n 60 masurca.assembly.log
```

## Assembly modes

### Illumina paired-end plus Oxford Nanopore

For a project with matching Nanopore and Illumina reads, replace the `PACBIO` line in sr_config.txt with:

```text
NANOPORE = /absolute/path/to/nanopore.fastq.gz
```

Leave the paired-end input in place. 

### Multiple Illumina paired-end libraries

For multiple paired-end libraries, add a `PE = ` line for each with a unique two-character prefix (i.e. pe, pa). Do not concatenate libraries with different fragment distributions into one library entry.


## Check the assembly statistics

### Install seqkit if not available
```bash
#load your HPC's seqkit module 
ml seqkit/2.4.0 

#if your HPC does not have seqkit installed, use miniconda to install your own seqkit
ml miniconda
conda create -n seqkit -c conda-forge -c bioconda seqkit
conda activate seqkit
```

### seqkit stats -a output definitions

| Column | Meaning |
| --- | --- |
| `file` | Input filename. |
| `format` | Sequence format, such as FASTA or FASTQ. |
| `type` | Sequence type, such as DNA. |
| `num_seqs` | Number of sequences; here, the number of assembly records. |
| `sum_len` | Total sequence length, including gap characters. |
| `min_len` | Shortest sequence length. |
| `avg_len` | Mean sequence length. |
| `max_len` | Longest sequence length. |
| `Q1`, `Q2`, `Q3` | First quartile, median, and third quartile of sequence lengths. |
| `sum_gap` | Gap-character count under SeqKit's gap definition; not the count of ambiguous `N` bases. |
| `N50` | Sequence length at which sequences of that length or longer account for at least half of the total span. |
| `Q20(%)`, `Q30(%)` | Base-quality summaries for FASTQ reads; not assembly accuracy measurements for a FASTA file. |
| `GC(%)` | Percentage of bases that are G or C. |

### Arabidopsis Illumina plus PacBio dataset

Use SeqKit for basic sequence statistics and BUSCO for conserved gene completeness. Install evaluation programs separately from MaSuRCA, and activate that environment only after the assembly finishes:

```bash
 seqkit stats -a flye.mr.67.17.15.0.02/assembly.scaffolds.fasta 
 ```


| file | format | type | num_seqs | sum_len | min_len | avg_len | max_len | Q1 | Q2 | Q3 | sum_gap | N50 | Q20(%) | Q30(%) | GC(%) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| flye.mr.67.17.15.0.02/assembly.scaffolds.fasta | FASTA | DNA | 169 | 119,259,814 | 517 | 705,679.4 | 13,791,725 | 8,095 | 52,500 | 272,857 | 0 | 8,149,967 | 0 | 0 | 36.05 |


## Errors and troubleshooting

| Problem | What to check |
| --- | --- |
| Configuration fails | Verify absolute paths, unique library prefixes, section-ending `END` lines, and installation dependencies. Read `masurca.configure.log`. |
| Job runs out of memory | Inspect scheduler records, request more memory, and review the Jellyfish hash allocation. More threads do not resolve insufficient RAM. |
| Job hits its time limit | Check which stage stopped and retain intermediates so completed work can be reused. |
| A tool or Perl module cannot be found | Check the MaSuRCA installation log and executable paths; avoid conflicting environments. |
| Corrected reads or final assembly are unexpectedly sparse | Check usable coverage, read lengths, sample identity, and the first failed stage in the logs. |
| No final FASTA appears | Inspect the log for a failed upstream step; intermediate files alone do not indicate success. |

After a timeout or interruption, first make sure the earlier job is no longer running. Preserve its logs, address the cause, and resubmit from the same assembly directory. Regenerating `assemble.sh` allows MaSuRCA to account for existing intermediates. Inspect any failed or partial outputs before restarting and use stage-specific guidance from the developers rather than deleting files indiscriminately. If you change the input reads or assembly strategy, start a separate assembly directory.

## Summarized interpretation

The MaSuRCA hybrid assembly contains 119.3 Mb in 169 scaffolds, with an N50 of 8.15 Mb and a longest scaffold of 13.79 Mb. This total is broadly consistent with the approximately 120 Mb represented in published Arabidopsis thaliana Landsberg erecta (Ler) assemblies. Jiao and Schneeberger (2020) reported 118.5 Mb of chromosome-anchored sequence for Ler, although their assemblies also included unanchored, predominantly repetitive sequences. Consequently, the slightly larger total here does not independently indicate assembly inflation, and neither assembly span should be interpreted as a complete measurement of the physical genome.

These statistics indicate substantial assembly continuity but do not establish completeness or accuracy. BUSCO should be used to evaluate conserved gene completeness, while duplicated BUSCOs may warrant investigation for retained haplotypes or other redundant sequence. Merqury can complement BUSCO by comparing k-mers from Illumina reads with the assembly to estimate consensus accuracy (QV) and sequence completeness, while highlighting potential missing or duplicated sequence without requiring a reference genome. Contamination screening and examination of read coverage, GC content, and taxonomic assignments can help distinguish nuclear sequences from organellar sequences and potential contaminants. Small scaffolds should therefore be evaluated before removal. Read mapping and comparison with the published Ler assembly can further identify unsupported joins, missing regions, and possible structural differences. 

- [Published Ler assembly: Jiao and Schneeberger (2020)](https://doi.org/10.1038/s41467-020-14779-y)


## Further reading

- [MaSuRCA 4.1.4 release](https://github.com/alekseyzimin/masurca/releases/tag/v4.1.4)
- [MaSuRCA manual and configuration reference](https://github.com/alekseyzimin/masurca)
- [Zimin et al. (2013): The MaSuRCA genome assembler](https://doi.org/10.1093/bioinformatics/btt476)
- [Zimin et al. (2017): Hybrid assembly with the mega-reads algorithm](https://doi.org/10.1101/gr.213405.116)
- [Jiao and Schneeberger (2020): Chromosome-level assemblies of multiple Arabidopsis genomes](https://doi.org/10.1038/s41467-020-14779-y)
- [SeqKit documentation](https://bioinf.shenwei.me/seqkit/usage/)


[Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
