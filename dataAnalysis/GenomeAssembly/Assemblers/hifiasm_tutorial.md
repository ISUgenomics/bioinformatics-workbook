---
title: "Genome Assembly with hifiasm: HiFi Reads, Hi-C Phasing, and Trio Binning"
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

## Background

Main software website links:

- [hifiasm GitHub repository and quick start](https://github.com/chhylp123/hifiasm)
- [hifiasm manual](https://hifiasm.readthedocs.io/en/latest/)
- [hifiasm parameter reference](https://hifiasm.readthedocs.io/en/latest/parameter-reference.html)
- [hifiasm paper: Haplotype-resolved de novo assembly using phased assembly graphs](https://doi.org/10.1038/s41592-020-01056-5)

hifiasm is a de novo genome assembler originally developed for accurate PacBio HiFi long reads. It can use sequence differences between homologous chromosomes to separate haplotypes during assembly. Additional Hi-C reads or parental sequencing data can improve haplotype resolution.

This tutorial uses PacBio HiFi reads from one diploid individual. Recent versions also support ONT R10 assembly through a separate `--ont` mode; the commands below are for HiFi data. PacBio CLR/subreads should not be substituted for HiFi reads.

## Assembly workflow

A standard run performs these operations internally:

1. Correct residual read errors and identify read overlaps.

2. Construct and simplify assembly graphs.

3. Generate primary contigs and partially phased contigs.

Run these operations with one command. hifiasm also saves binary intermediate files that can accelerate subsequent runs.

## How to run hifiasm

### Parameter considerations

| Parameter or input | Purpose |
| --- | --- |
| `-o sample.asm` | Set the output prefix; this is not a directory-creation option. |
| `-t 32` | Use 32 CPU threads. Match this to the available CPU allocation. |
| `reads.hifi.fastq.gz` | Supply HiFi reads in FASTA or FASTQ format, optionally gzip-compressed. |
| Multiple input files | List files from the same sample after the options. |
| `--hg-size 100m` | Optionally provide an estimated haploid genome size to help coverage inference. It is not required. |
| `-f0` | Disable the initial Bloom filter for small genomes; do not apply this automatically to large genomes. |
| `-l0` | Disable haplotig purging; the hifiasm FAQ recommends this for inbred or highly homozygous genomes. |
| `--primary` | Request primary and alternate assemblies. A standard HiFi-only run already produces a primary assembly. |
| `-i` | Ignore cached binary intermediates and recompute from input reads. |

Start with defaults for a heterozygous diploid sample. Coverage, heterozygosity, repeats, and genome size affect resource use and assembly quality. The [FAQ](https://hifiasm.readthedocs.io/en/latest/faq.html) gives approximately 13× coverage per haplotype as a usual minimum, rather than a guarantee of a complete assembly.

### Running hifiasm

Check the module name and version on your cluster:

```bash
module avail hifiasm
# Load the module used for this tutorial.
module load hifiasm/0.25.1
hifiasm --version
hifiasm -h

# Alternatively install hifiasm, yak, and seqkit
conda create -n hifiasm -c conda-forge -c bioconda hifiasm yak seqkit
conda activate hifiasm
```

### Acquire sample data

| Data type | Run accession | Verified description | Amount |
|---|---|---|---:|
| **Arabidopsis thaliana PacBio reads** | [ERR6210723](https://www.ncbi.nlm.nih.gov/sra/?term=ERR6210723) | *A. thaliana* Col-0 HiFi reads used for the Col-CEN assembly; project **PRJEB46164** | ~14.65 Gb |
| **Arabidopsis Hi-C** | [SRR1504819](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1504819) | *A. thaliana* Col-0; paired-end HindIII Hi-C, HiSeq 2000; project **PRJNA253621** | ~21.94 Gb |
| **Screwworm offspring PacBio HiFi** | [SRR33955207](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33955207) | *C. hominivorax* F1 male from a **J-06 female × Yaviza male** cross; PacBio Revio; project **PRJNA1237801** | ~88.18 Gb |
| **Screwworm maternal Illumina** | [SRR33955209](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33955209) | *C. hominivorax* **J-06 female parent**; paired-end NextSeq 2000; project **PRJNA1237801** | ~23.51 Gb |
| **Screwworm paternal Illumina** | [SRR33955208](https://www.ncbi.nlm.nih.gov/sra/?term=SRR33955208) | *C. hominivorax* **Yaviza male parent**; paired-end NextSeq 2000; project **PRJNA1237801** | ~20.99 Gb |

Amounts in the table are total sequenced bases, not compressed download sizes.

Create the output and temporary directories before extracting reads:

```bash
mkdir -p fastq tmp
```

```bash
ml sratoolkit/3.2.1
# Arabidopsis thaliana data
prefetch --max-size u ERR6210723 SRR1504819
fasterq-dump ERR6210723/ERR6210723.sra -O fastq -t tmp -e 8 --progress
fasterq-dump SRR1504819/SRR1504819.sra --split-files -O fastq -t tmp -e 8 --progress
# Cochliomyia hominivorax (New World screwworm)
prefetch --max-size u SRR33955207 SRR33955209 SRR33955208
fasterq-dump SRR33955207/SRR33955207.sra -O fastq -t tmp -e 8 --progress
fasterq-dump SRR33955209/SRR33955209.sra --split-files -O fastq -t tmp -e 8 --progress
fasterq-dump SRR33955208/SRR33955208.sra --split-files -O fastq -t tmp -e 8 --progress
```

`prefetch` normally applies a download-size limit. If a run exceeds that limit, rerun the same command with an appropriate `--max-size` value (or `--max-size u` for no limit). 

From the project directory, compress the extracted files with:

```bash
#reduce the space that these reads occupy
cd fastq
pigz -p 64 *.fastq
#alternatively if pigz is not available 
gzip *.fastq
cd ..
```

`pigz` writes one gzip file per input and removes each uncompressed original after successful compression. The compression command runs in the foreground; `&` would place a command in the background. Use no more threads than your allocation provides.

### Genome assembly with *A. thaliana* PacBio HiFi

Use a separate output prefix for each independent experiment. Run the commands from the project directory so the relative `fastq/` paths resolve correctly. Match `-t` to the CPUs available to your job. To add additional HiFi files, list their paths separated by spaces after fastq/ERR6210723.fastq.gz and before 2> HIFIonly/sample.log.

In these commands, `2>` redirects diagnostic messages and progress information to a log file. It is Bash syntax, not a hifiasm option. 

#### Assembly with one HiFi dataset

```bash
mkdir -p HIFIonly
hifiasm -o HIFIonly/sample.asm -t 64 \
    fastq/ERR6210723.fastq.gz \
    2> HIFIonly/sample.log
```

## Expected files generated during assembly

For the HiFi-only commands above, look for these files. Names can change with assembly mode and software version.

| File ending | What it contains |
|---|---|
| `.ec.bin` | Cached error-corrected reads, reused when restarting assembly. |
| `.ovlp.source.bin`, `.ovlp.reverse.bin` | Cached read-overlap information, reused to avoid recalculating overlaps. |
| `.bp.r_utg.gfa` | Raw unitig graph—an intermediate assembly retaining haplotype variation. |
| `.bp.p_utg.gfa` | Processed unitig graph, with small bubbles simplified. |
| `.bp.p_ctg.gfa` | Primary contig assembly; can contain a mixture of parental haplotypes. |
| `.bp.hap1.p_ctg.gfa`, `.bp.hap2.p_ctg.gfa` | Two **partially phased** contig assemblies. |
| `.noseq.gfa` | Corresponding graph without DNA sequences, useful for visualization. |
| `.lowQ.bed` | Coordinates of regions flagged as low quality. |
| `sample.log` | Assembly progress, diagnostic messages, and statistics. |

For routine use of one representative genome, convert the primary `.bp.p_ctg.gfa` file to FASTA. The `.bp.hap1` and `.bp.hap2` files are partially phased and do not establish parental origin.

Additional options to consider: `--primary` requests a primary assembly plus alternate contigs; `-l0` disables haplotig purging and also produces primary/alternate outputs. 

## Assembly modes


### HiFi plus Hi-C

Provide paired Hi-C reads using `--h1` for R1 and `--h2` for R2:

```bash
mkdir HIFI_HIC
hifiasm -o HIFI_HIC/sample.hic.asm -t 64 \
    --h1 fastq/SRR1504819_1.fastq.gz \
    --h2 fastq/SRR1504819_2.fastq.gz \
    fastq/ERR6210723.fastq.gz 2> HIFI_HIC/sample.log
```

The haplotype outputs use `.hic.hap1.p_ctg.gfa` and `.hic.hap2.p_ctg.gfa`. 

### Trio assembly with screwworm (*C. hominivorax*)

With parental short reads, first create parental k-mer databases using [yak](https://github.com/lh3/yak), then assemble the offspring's HiFi reads:

```bash
#in case you did not install yak above
conda install -c conda-forge -c bioconda yak
```

For a single end reads per parent, the general syntax is:

```bash
# Generic single-file examples; skip these two commands in this tutorial for the paired-end dataset below.
yak count -k31 -b37 -t64 -o paternal.yak paternal.fastq.gz
yak count -k31 -b37 -t64 -o maternal.yak maternal.fastq.gz
```

For this tutorial assembly, supply the paired-end reads as follows:

```bash
mkdir TrioBinning
cd TrioBinning
# Female parent
yak count -k31 -b37 -t64 -o female.yak \
    <(zcat ../fastq/SRR33955209_1.fastq.gz ../fastq/SRR33955209_2.fastq.gz) \
    <(zcat ../fastq/SRR33955209_1.fastq.gz ../fastq/SRR33955209_2.fastq.gz)
# Male parent
yak count -k31 -b37 -t64 -o male.yak \
    <(zcat ../fastq/SRR33955208_1.fastq.gz ../fastq/SRR33955208_2.fastq.gz) \
    <(zcat ../fastq/SRR33955208_1.fastq.gz ../fastq/SRR33955208_2.fastq.gz)
hifiasm -o MaleFemaleScrewworm.asm -t 64 \
    -1 male.yak -2 female.yak \
    ../fastq/SRR33955207.fastq.gz 2> MaleFemaleScrewworm.asm.log
```

The two identical streams supplied to each `yak count` command are intentional: with `-b37`, yak reads the sequences in two passes. Each stream includes both mates from one parent; the first-pass counts are cleared before the second pass, so the final counts are not doubled. The `<(...)` syntax is Bash process substitution.

Outputs include `.dip.hap1.p_ctg.gfa` and `.dip.hap2.p_ctg.gfa`, corresponding to the parental databases passed with `-1` and `-2`. Here, hap1 is paternal (Yaviza; SRR33955208) and hap2 is maternal (J-06; SRR33955209). Both assemblies represent inherited haplotypes of the same male offspring, rather than separate male and female individuals.

## Check the assembly statistics

Run each conversion in the directory containing its GFA file. Use a sequence-containing GFA, not a `.noseq.gfa` file. 

Before annotation, evaluate total length, contig count, N50, gene completeness, contamination, and read support. Where possible, assess consensus accuracy and retained or missing sequence with independent k-mer data.

### seqkit stats -a output definitions
| Column | Meaning |
|---|---|
| `file` | Input filename. |
| `format` | Sequence file format, such as FASTA or FASTQ. |
| `type` | Sequence type, such as DNA, RNA, or protein. |
| `num_seqs` | Number of sequences; for an assembly FASTA, usually the number of contigs or scaffolds. |
| `sum_len` | Total length of all sequences, including gap characters. |
| `min_len` | Length of the shortest sequence. |
| `avg_len` | Mean sequence length: `sum_len / num_seqs`. |
| `max_len` | Length of the longest sequence. |
| `Q1` | First quartile of sequence lengths: approximately 25% of sequences are this length or shorter. |
| `Q2` | Median sequence length: approximately half the sequences are this length or shorter. |
| `Q3` | Third quartile of sequence lengths: approximately 75% of sequences are this length or shorter. |
| `sum_gap` | Total number of gap characters, such as `-`; this is not a count of ambiguous `N` bases. |
| `N50` | Sequence length at which sequences of that length or longer account for at least 50% of the total sequence length. |
| `Q20(%)` | Percentage of bases with Phred quality scores ≥20; only applicable to FASTQ data. |
| `Q30(%)` | Percentage of bases with Phred quality scores ≥30; only applicable to FASTQ data. |
| `GC(%)` | Percentage of bases that are G or C. |

### Arabidopsis HiFi only dataset

Run in `HIFIonly/`.

```bash
awk -F '\t' '$1 == "S" {print ">" $2; print $3}'  sample.asm.bp.p_ctg.gfa >sample.asm.bp.p_ctg.fasta
seqkit stats -a sample.asm.bp.p_ctg.fasta
```

Example output:

```text

file                       format  type  num_seqs      sum_len  min_len    avg_len     max_len      Q1      Q2      Q3  sum_gap         N50  Q20(%)  Q30(%)  GC(%)
sample.asm.bp.p_ctg.fasta  FASTA   DNA      1,225  184,939,908    2,478  150,971.4  26,129,794  27,398  32,688  47,438        0  12,446,441       0       0  36.74
```

### Arabidopsis HiFi plus Hi-C datasets

Run in `HIFI_HIC/`.

```bash
awk -F '\t' '$1 == "S" {print ">" $2; print $3}' sample.hic.asm.hic.p_ctg.gfa >sample.hic.asm.hic.p_ctg.fasta
seqkit stats -a sample.hic.asm.hic.p_ctg.fasta 
```

Example output:

```text
file                            format  type  num_seqs      sum_len  min_len    avg_len     max_len      Q1      Q2      Q3  sum_gap         N50  Q20(%)  Q30(%)  GC(%)
sample.hic.asm.hic.p_ctg.fasta  FASTA   DNA      1,110  181,041,676    2,478  163,100.6  26,129,794  28,103  33,688  48,726        0  12,446,441       0       0  36.75
```

### Screwworm trio: paternal and maternal haplotype

Run in `fastq/`, where the trio assembly was created.

```bash
awk -F '\t' '$1 == "S" {print ">" $2; print $3}' MaleFemaleScrewworm.asm.dip.hap1.p_ctg.gfa >PaternalScrewworm.asm.dip.hap1.p_ctg.fasta
awk -F '\t' '$1 == "S" {print ">" $2; print $3}' MaleFemaleScrewworm.asm.dip.hap2.p_ctg.gfa >MaternalScrewworm.asm.dip.hap2.p_ctg.fasta

seqkit stats -a *fasta
```

Example output:

```text
file                                          format  type  num_seqs      sum_len  min_len      avg_len     max_len        Q1      Q2       Q3  sum_gap         N50  Q20(%)  Q30(%)  GC(%)
PaternalScrewworm.asm.dip.hap1.p_ctg.fasta  FASTA   DNA        171  452,923,884   23,090  2,648,677.7  45,033,906  47,451.5  98,990  610,554        0  17,497,353       0       0  27.83
MaternalScrewworm.asm.dip.hap2.p_ctg.fasta  FASTA   DNA        181  467,724,117   14,474  2,584,111.1  65,710,413    33,490  72,368  508,348        0  18,032,292       0       0   27.6
```

Evaluate each haplotype separately. Concatenating both haplotypes into one reference can create apparent gene duplication and ambiguous mapping.

## Errors and troubleshooting

| Problem | What to check |
| --- | --- |
| `hifiasm: command not found` | Check the module name and confirm that `hifiasm --version` works. |
| Job killed or out of memory | Check the job log and recorded memory usage, then request more memory if needed. For small genomes, consider `-f0`. |
| Job reaches its time limit | Request more time. Rerunning with the same prefix can reuse completed binary intermediates. |
| Primary assembly is unexpectedly large | Check contamination, the genome-size estimate, and the inferred homozygous coverage before changing purging settings. |
| Assembly is fragmented | Review read quality, coverage, repeats, and sample complexity. |
| Expected `.bp.p_ctg.gfa` is absent | Check the selected mode, version, log, and actual output names before converting files. |
| New reads or correction settings are being tested | Use a fresh prefix or `-i` to avoid reusing inappropriate intermediates. |

Avoid changing `--hom-cov` or purging parameters simply to force the assembly to an expected size. Base changes on read coverage and independent quality evidence.

## Further reading

- [hifiasm quick start and current features](https://github.com/chhylp123/hifiasm)
- [HiFi-only assembly](https://hifiasm.readthedocs.io/en/latest/pa-assembly.html)
- [Output file descriptions](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html)
- [Frequently asked questions](https://hifiasm.readthedocs.io/en/latest/faq.html)
- [Original hifiasm paper](https://doi.org/10.1038/s41592-020-01056-5)

[Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)