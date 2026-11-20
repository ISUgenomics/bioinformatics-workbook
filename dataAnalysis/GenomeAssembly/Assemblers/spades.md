---
title: "Metagenome Assembly with SPAdes 4.3: Illumina and Hybrid Assembly"
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

## Background

Main software website links:

- [SPAdes GitHub repository](https://github.com/ablab/spades)
- [SPAdes manual](https://ablab.github.io/spades/)
- [SPAdes 4.3.0 release](https://github.com/ablab/spades/releases/tag/v4.3.0)
- [metaSPAdes paper: A New Versatile Metagenomic Assembler](https://doi.org/10.1101/gr.213959.116)
- [SPAdes hybrid assembly documentation](https://ablab.github.io/spades/hybrid.html)

metaSPAdes is the metagenomic assembly workflow included in SPAdes. It assembles reads from microbial communities, where different organisms and strains can have very different sequencing depths. Long reads can provide additional information for resolving repeated sequences and connections in the assembly graph.

This tutorial uses **SPAdes 4.3.0** to compare two assemblies of the **NWC2 natural whey culture metagenome from Gruyère cheese**:

| Assembly | Input | Program and mode |
| --- | --- | --- |
| Illumina only | Paired-end Illumina reads | SPAdes 4.3.0 with `--meta` |
| Hybrid | The same Illumina pairs plus Nanopore reads | SPAdes 4.3.0 with `--meta --nanopore` |

The SPAdes manual describes hybrid metagenomic assembly as **experimental**. We will evaluate whether adding long reads improves this dataset rather than assume it produces a better assembly. Closely related strains can remain collapsed or generate ambiguous graph paths.

## How to run metaSPAdes

### Parameter considerations

| Parameter or input | Required? | Purpose |
| --- | --- | --- |
| `--meta` | Required with `spades.py` for this workflow | Select the metagenomic pipeline; `metaspades.py` selects it automatically. |
| `-1 R1.fastq.gz` | Required here | Supply forward paired-end reads. |
| `-2 R2.fastq.gz` | Required here | Supply matching reverse reads in the same order. |
| `-o IlluminaOnly` | Required | Set the output directory for this experiment. |
| `--nanopore reads.fastq.gz` | Required for the hybrid comparison | Supply supplementary Nanopore reads. |
| `-t 64` | Optional; explicitly set here | Match the allocated CPU threads. |
| `-m 128` | Optional; explicitly set here | Set SPAdes' memory limit in GB, below the job's RAM allocation. |
| `-k` | Optional; omitted here | Override automatic k-mer selection; values must be odd, below 128, and compatible with short-read lengths. |
| `--frugal` | Optional; omitted here | Save memory with potentially altered assembly quality or contiguity. |
| `--continue` | Optional restart instruction | Resume with the previous settings; specify only `-o` alongside it. |

### Running SPAdes

Check to see what modules are already available on your HPC. 
```bash
module --default avail spades seqkit fastp fastplong datasets unzip quast sra
```

Install the software that is missing from your HPC, but use modules with current software versions when available. Most HPC's have a conda module, but if not see here to install miniconda locally. [Miniconda Installation](https://www.anaconda.com/docs/getting-started/installation#linux)
```bash
#quast has incompatibilites with the other software, so create another conda environment.
module load miniconda
conda create -n quast -c conda-forge -c bioconda python=3.11 quast
conda activate quast

conda deactivate 
conda create -n metagenome_assembly -c conda-forge -c bioconda seqkit fastp fastplong ncbi-datasets-cli unzip sra-tools spades=4.3.0
conda activate metagenome_assembly
```

### Acquire sample data

We will use **NWC2**, a natural whey starter culture used in Gruyère cheese production, deposited by Agroscope under [BioProject PRJNA477604](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA477604).

| Data type | NCBI SRA run | Instrument | Archived sequenced bases | Approximate compressed FASTQ size |
| --- | --- | --- | ---: | ---: |
| Illumina paired-end | [SRR7585899](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7585899) | MiSeq | 454,935,620 | 0.32 GB across both mates |
| Oxford Nanopore | [SRR7585900](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7585900) | MinION | 1,385,178,778 | 1.38 GB |

The study recovered four dominant bacterial genomes: *Streptococcus thermophilus*, *Lactobacillus delbrueckii*, and two *L. helveticus* strains, along with plasmids and phages. 

Download the Illumina and nanopore reads.
```bash
mkdir fastq
cd fastq
mkdir tmp
prefetch --max-size u SRR7585899 SRR7585900
fasterq-dump SRR7585899/SRR7585899.sra --split-files -t tmp -e 16 
fasterq-dump SRR7585900/SRR7585900.sra -t tmp -e 8 

#compress for efficient storage
pigz -p 24 *fastq

#return to project directory
cd ..
```

**Coverage varies by organism and strain.** Do not divide all sequenced bases by the genome size of a single bacterium and report that as metagenome coverage.

## Filter and trim sequence reads

### Prepare the Illumina reads

Adapter readthrough and other technical sequences can disrupt metagenome assembly. Inspect raw-read quality and use fastp to remove detected adapters and filter poor-quality pairs:

```bash
mkdir logs 

fastp \
  -i fastq/SRR7585899_1.fastq.gz \
  -I fastq/SRR7585899_2.fastq.gz \
  -o fastq/SRR7585899_R1.clean.fastq.gz \
  -O fastq/SRR7585899_R2.clean.fastq.gz \
  --detect_adapter_for_pe \
  --thread 24 \
  --html fastq/illumina.fastp.html \
  --json fastq/illumina.fastp.json \
  > logs/fastp.log 2>&1
```

The 150 bp minimum is a tutorial filtering choice, not a validated optimum. Review the before/after reports and retained bases. Both assemblies use exactly the same retained Illumina pairs; orphaned or discarded reads are not added back.

### Prepare the Nanopore reads 

```bash
fastplong \
  -i fastq/SRR7585900.fastq.gz \
  -o fastq/SRR7585900.trimmed.fastq.gz \
  --qualified_quality_phred 9 \
  --disable_length_filtering \
  --html nanopore.fastplong.html \
  --json nanopore.fastplong.json
```

Fastplong documented defaults classify bases below Q15 as unqualified and reject reads with more than 40% unqualified bases. With older Nanopore data this can remove substantial amounts of data, potentially reducing useful coverage from less abundant organisms. This command attempts to detect and trim Nanopore adapters, removes reads with more than the default 40% of bases below Q9, disables length filtering, and saves the processed reads with HTML and JSON quality reports.

### Evaluate trimming results
```bash
seqkit stats -a fastq/*fastq.gz
```

| file | format | type | num_seqs | sum_len | min_len | avg_len | max_len | Q1 | Q2 | Q3 | sum_gap | N50 | N50_num | Q20(%) | Q30(%) | AvgQual | GC(%) | sum_n |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| fastq/SRR7585899_1.fastq.gz | FASTQ | DNA | 1,410,764 | 226,409,627 | 32 | 160.5 | 301 | 83 | 141 | 238 | 0 | 220 | 403,410 | 97.73 | 94.01 | 25.95 | 45.19 | 152,409 |
| fastq/SRR7585899_2.fastq.gz | FASTQ | DNA | 1,410,764 | 228,525,993 | 32 | 162 | 301 | 83 | 141 | 243 | 0 | 224 | 397,748 | 91.39 | 83.25 | 20.46 | 45.26 | 152,810 |
| fastq/SRR7585899_R1.clean.fastq.gz | FASTQ | DNA | 1,391,400 | 222,535,706 | 31 | 159.9 | 301 | 83 | 141 | 236 | 0 | 218 | 398,643 | 98.32 | 94.78 | 27.87 | 45.21 | 18 |
| fastq/SRR7585899_R2.clean.fastq.gz | FASTQ | DNA | 1,391,400 | 221,618,947 | 31 | 159.3 | 301 | 83 | 140 | 234 | 0 | 216 | 397,893 | 93.3 | 85.43 | 21.86 | 45.2 | 0 |
| fastq/SRR7585900.fastq.gz | FASTQ | DNA | 407,027 | 1,385,178,778 | 138 | 3,403.2 | 118,642 | 586 | 1,177 | 3,179 | 0 | 9,562 | 35,181 | 16.7 | 0 | 8.45 | 39.51 | 0 |
| fastq/SRR7585900.trimmed.fastq.gz | FASTQ | DNA | 321,533 | 1,128,723,866 | 143 | 3,510.4 | 102,363 | 607 | 1,235 | 3,341 | 0 | 9,704 | 28,423 | 17.96 | 0 | 8.75 | 39.35 | 0 |

## Assembly modes


### Metagenome assembly with Illumina reads

Run from the project directory inside a compute allocation with 64 CPUs and 128 GB RAM. These CPU/RAM settings are starting allocations from what is available in my HPC. More threads do not necessarily produce proportional speed improvements.

```bash
spades.py --meta \
  -1 fastq/SRR7585899_R1.clean.fastq.gz \
  -2 fastq/SRR7585899_R2.clean.fastq.gz \
  -t 64 -m 128 \
  -o IlluminaOnly \
  > logs/IlluminaOnly.log 2>&1
```

`-o` creates the output directory. `> ... 2>&1` is shell redirection of progress and diagnostic messages, not a SPAdes option. SPAdes also writes its own `spades.log`.


### Hybrid metagenome assembly: Illumina plus Nanopore

Keep the Illumina inputs and resource settings unchanged. Add the Nanopore file and use a new output directory:

```bash
spades.py --meta \
  -1 fastq/SRR7585899_R1.clean.fastq.gz \
  -2 fastq/SRR7585899_R2.clean.fastq.gz \
  --nanopore fastq/SRR7585900.trimmed.fastq.gz \
  -t 64 -m 128 \
  -o Hybrid \
  > logs/Hybrid.log 2>&1
```

This is a separate hybrid metagenomic assembly, not a polishing step applied to the Illumina assembly. `--meta` retains metagenomic settings; `--nanopore` adds long-read evidence. 

Long reads may bridge repeat-associated breaks, but strain mixtures and uneven coverage still complicate interpretation. In particular, NWC2 contains two related *L. helveticus* strains. A longer assembly does not establish that those strains were correctly separated.


## Expected files generated during assembly

For both assemblies, look inside the directory supplied with `-o`:

| File or directory | What it contains |
| --- | --- |
| `contigs.fasta` | Final contig sequences. |
| `scaffolds.fasta` | Scaffold sequences, potentially including gaps represented by Ns. |
| `assembly_graph_with_scaffolds.gfa` | Assembly graph and scaffold paths, normally in GFA 1.2 format. |
| `assembly_graph.fastg` | Assembly graph in FASTG format. |
| `contigs.paths`, `scaffolds.paths` | Reported sequence paths through the graph. |
| `corrected/` | Read-correction outputs/intermediates when correction is enabled, as in these commands. |
| `spades.log` | Progress, warnings, and final output locations. |
| `params.txt`, `input_dataset.yaml` | Recorded parameters and input-library descriptions. |
| `K*/` | Intermediate files for individual k-mer iterations. |

Use the final contig/scaffold FASTA files for evaluation, not FASTA files from the `K*/` directories. The `_cov_` value in SPAdes headers is k-mer coverage, not a direct per-base coverage estimate.


## Check the assembly statistics

### seqkit stats -a output definitions

| Column | Meaning |
| --- | --- |
| `file` | Input filename. |
| `format`, `type` | File format and sequence type. |
| `num_seqs` | Number of sequences, usually contigs or scaffolds for an assembly. |
| `sum_len` | Total sequence length. |
| `min_len`, `avg_len`, `max_len` | Minimum, mean, and maximum sequence lengths. |
| `Q1`, `Q2`, `Q3` | First quartile, median, and third quartile of sequence lengths. |
| `sum_gap` | Count under SeqKit's gap-character definition; not an ambiguous-N count. |
| `N50` | Length at which sequences of that length or longer account for at least half of the total span. |
| `Q20(%)`, `Q30(%)` | FASTQ base-quality summaries; not consensus accuracy measurements for assembly FASTA files. |
| `GC(%)` | Percentage of bases that are G or C. |

### Compare the two assemblies

Run from the project directory:

```bash
seqkit stats -a IlluminaOnly/contigs.fasta IlluminaOnly/scaffolds.fasta Hybrid/contigs.fasta Hybrid/scaffolds.fasta
```

| file | format | type | num_seqs | sum_len | min_len | avg_len | max_len | Q1 | Q2 | Q3 | sum_gap | N50 | N50_num | Q20(%) | Q30(%) | AvgQual | GC(%) | sum_n |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IlluminaOnly/contigs.fasta | FASTA | DNA | 4,225 | 6,149,075 | 56 | 1,455.4 | 53,777 | 98 | 125 | 860 | 0 | 8,288 | 192 | 0 | 0 | 0 | 41.65 | 0 |
| IlluminaOnly/scaffolds.fasta | FASTA | DNA | 4,177 | 6,153,654 | 56 | 1,473.2 | 53,777 | 98 | 122 | 808 | 0 | 8,650 | 184 | 0 | 0 | 0 | 41.62 | 4,800 |
| Hybrid/contigs.fasta | FASTA | DNA | 3,215 | 6,342,402 | 56 | 1,972.8 | 134,674 | 102 | 120 | 830.5 | 0 | 12,820 | 103 | 0 | 0 | 0 | 41.64 | 0 |
| Hybrid/scaffolds.fasta | FASTA | DNA | 3,077 | 6,606,546 | 56 | 2,147.1 | 275,586 | 102 | 112 | 631 | 0 | 19,683 | 74 | 0 | 0 | 0 | 39.99 | 261,035 |

### Compare with published community genomes using MetaQUAST

Use multiple matching community references rather than a single bacterial chromosome. NCBI lists these NWC2 assemblies:

| Published genome | BioSample | Assembly |
| --- | --- | --- |
| *Streptococcus thermophilus* NWC_2_1 | SAMN09476686 | GCA_049952445.1 |
| *Lactobacillus delbrueckii* subsp. *lactis* NWC_2_2 | SAMN09476687 | GCA_049952475.1 |
| *Lactobacillus helveticus* NWC_2_3 | SAMN09476688 | GCA_003814305.1 |
| *Lactobacillus helveticus* NWC_2_4 | SAMN09476689 | GCA_003814325.1 |

Download the versioned genome FASTAs using NCBI Datasets, keeping one FASTA per genome assembly:

```bash
mkdir reference
datasets download genome accession \
  GCA_049952445.1 GCA_049952475.1 GCA_003814305.1 GCA_003814325.1 \
  --include genome --filename reference/Genomes.zip

unzip -q reference/Genomes.zip -d reference

#rename genomes by accession only 
for accession in GCA_049952445.1 GCA_049952475.1 GCA_003814305.1 GCA_003814325.1; do
  cp reference/ncbi_dataset/data/"$accession"/*_genomic.fna \
     reference/"$accession".fna
done

#get genome stats on downloaded accessions
seqkit stats -a reference/*.fna
```

| file | format | type | num_seqs | sum_len | min_len | avg_len | max_len | Q1 | Q2 | Q3 | sum_gap | N50 | N50_num | Q20(%) | Q30(%) | AvgQual | GC(%) | sum_n |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| reference/GCA_003814305.1.fna | FASTA | DNA | 2 | 2,232,981 | 22,170 | 1,116,490.5 | 2,210,811 | 22,170 | 1,116,490.5 | 2,210,811 | 0 | 2,210,811 | 1 | 0 | 0 | 0 | 37.51 | 0 |
| reference/GCA_003814325.1.fna | FASTA | DNA | 3 | 2,230,130 | 22,170 | 743,376.7 | 2,177,422 | 26,354 | 30,538 | 1,103,980 | 0 | 2,177,422 | 1 | 0 | 0 | 0 | 37.38 | 0 |
| reference/GCA_049952445.1.fna | FASTA | DNA | 1 | 1,971,439 | 1,971,439 | 1,971,439 | 1,971,439 | 1,971,439 | 1,971,439 | 1,971,439 | 0 | 1,971,439 | 1 | 0 | 0 | 0 | 38.96 | 0 |
| reference/GCA_049952475.1.fna | FASTA | DNA | 2 | 2,278,067 | 8,888 | 1,139,033.5 | 2,269,179 | 8,888 | 1,139,033.5 | 2,269,179 | 0 | 2,269,179 | 1 | 0 | 0 | 0 | 48.73 | 0 |

Run a comparison of your assembled contigs with the published genomes using metaquast
```bash
conda deactivate
conda activate quast
mkdir evaluation

metaquast.py IlluminaOnly/contigs.fasta Hybrid/contigs.fasta \
  -r reference \
  --labels IlluminaOnly,Hybrid --min-contig 500 \
  -t 64 -o evaluation/metaquast
```

### Quast output

| Column | Describes | Explanation |
|---|---|---|
| Genome | Reference | Accession of the published reference assembly. |
| Species | Reference | Organism represented by the reference assembly. |
| # fragments | Reference | Number of sequences in the reference assembly, including chromosomes and any included plasmids. |
| Length, bp | Reference | Total length of all sequences in the reference assembly. |
| Mean genome fraction, % | Assembly evaluation | Percentage of the reference covered by assembly alignments, averaged across the evaluated assemblies included in the summary. |
| # misassembled blocks | Assembly evaluation | Number of aligned assembly blocks flagged as misassembled relative to the reference; this does not indicate errors in the reference itself. |

| Genome | Species | Reference sequences | Reference length (bp) | Mean genome fraction (%) | # misassembled blocks |
|---|---|---:|---:|---:|---:|
| GCA_003814305.1 | *Lactobacillus helveticus* | 2 | 2,232,981 | 51.062 | 29 |
| GCA_003814325.1 | *Lactobacillus helveticus* | 3 | 2,230,130 | 80.771 | 41 |
| GCA_049952445.1 | *Streptococcus thermophilus* | 1 | 1,971,439 | 97.325 | 30 |
| GCA_049952475.1 | *Lactobacillus delbrueckii* subsp. *lactis* | 2 | 2,278,067 | 90.496 | 103 |

| Statistic | IlluminaOnly | Hybrid |
|---|---:|---:|
| **Alignment-based statistics** | | |
| Genome fraction (%) | 80.048 | 78.843 |
| Duplication ratio | 1.005 | 1.146 |
| Largest alignment | 53,777 | 93,041 |
| Total aligned length | 5,512,189 | 5,760,928 |
| NA50 | 8,569 | 11,718 |
| NA90 | 1,356 | 1,616 |
| auNA | 12,591 | 22,865 |
| LA50 | 175 | 109 |
| LA90 | 821 | 630 |
| **Misassemblies** | | |
| # misassemblies | 19 | 82 |
| # relocations | 14 | 70 |
| # translocations | 0 | 0 |
| # inversions | 0 | 0 |
| # interspecies translocations | 5 | 12 |
| # misassembled contigs | 19 | 49 |
| Misassembled contigs length | 225,927 | 795,141 |
| # possibly misassembled contigs | 14 | 15 |
| # possible misassemblies | 18 | 23 |
| # local misassemblies | 2 | 9 |
| # scaffold gap ext. mis. | 0 | 0 |
| # scaffold gap loc. mis. | 0 | 0 |
| # unaligned mis. contigs | 3 | 5 |
| **Unaligned** | | |
| # fully unaligned contigs | 115 | 66 |
| Fully unaligned length | 197,975 | 152,614 |
| # partially unaligned contigs | 15 | 18 |
| Partially unaligned length | 40,266 | 119,331 |
| **Per-base quality** | | |
| # mismatches per 100 kbp | 131.98 | 383.31 |
| # mismatches | 7,275 | 22,082 |
| # indels per 100 kbp | 12.26 | 131.89 |
| # indels | 676 | 7,598 |
| # indels (≤ 5 bp) | 594 | 7,417 |
| # indels (> 5 bp) | 82 | 181 |
| Indels length | 2,368 | 16,756 |
| # N's per 100 kbp | 0 | 0 |
| # N's | 0 | 0 |
| **Statistics without reference** | | |
| # contigs | 1,306 | 951 |
| # contigs (≥ 0 bp) | 4,225 | 3,215 |
| # contigs (≥ 1,000 bp) | 984 | 751 |
| # contigs (≥ 5,000 bp) | 332 | 303 |
| # contigs (≥ 10,000 bp) | 138 | 144 |
| # contigs (≥ 25,000 bp) | 25 | 43 |
| # contigs (≥ 50,000 bp) | 3 | 19 |
| Largest contig | 53,777 | 134,674 |
| Total length | 5,761,340 | 6,044,727 |
| Total length (≥ 0 bp) | 6,149,075 | 6,342,402 |
| Total length (≥ 1,000 bp) | 5,528,341 | 5,897,503 |
| Total length (≥ 5,000 bp) | 3,966,574 | 4,792,872 |
| Total length (≥ 10,000 bp) | 2,591,671 | 3,635,553 |
| Total length (≥ 25,000 bp) | 877,432 | 2,070,864 |
| Total length (≥ 50,000 bp) | 157,715 | 1,283,636 |
| N50 | 8,911 | 14,368 |
| N90 | 1,764 | 2,572 |
| auN | 13,101 | 26,421 |
| L50 | 169 | 92 |
| L90 | 732 | 481 |
| **Similarity statistics** | | |
| # similar correct contigs | 55 | 55 |
| # similar misassembled blocks | 3 | 3 |


Inspect the combined and per-reference reports for genome fraction, aligned block continuity, mismatches/indels, and misassembly calls. These published genomes were themselves recovered and polished from sequencing data, so they are comparison references rather than independent ground truth. Closely related strains can share alignments, and structural differences require verification with alignments. Map Illumina and Nanopore reads back to each assembly to investigate depth variation and questionable joins.

The four assembly downloads do not necessarily represent every phage or rare member of the original sample. Do not classify every unaligned contig as contamination. Check which plasmid records are included in each downloaded assembly before interpreting plasmid recovery.

## Errors and troubleshooting

| Problem | What to check |
| --- | --- |
| `spades.py: command not found` | Check that the module is loaded or export spades to path. |
| Mode conflict | Use `--meta`; remove `--isolate` and `--careful`. |
| Coverage cutoff rejected | Remove `--cov-cutoff`; it is not supported in metaSPAdes. |
| Job killed or memory limit exceeded | Inspect scheduler records; increase RAM and keep `-m` below the allocation. |
| Read pairs do not match | Check extraction filenames, pair order, and preprocessing reports. |
| Hybrid assembly improves little | Inspect long-read lengths, usable depth, quality, and strain complexity. |
| Larger N50 but lower reference recovery | Check for missing low-abundance sequences or collapsed/chimeric paths. |

To resume an interrupted run with the same settings, specify only `--continue` and `-o`:

```bash
spades.py --continue -o Hybrid > logs/Hybrid.resume.log 2>&1
```

To change supported resource settings while restarting, use `--restart-from last`; do not repeat input files or mode flags:

```bash
spades.py --restart-from last -t 64 -m 128 -o Hybrid \
  > logs/Hybrid.restart.log 2>&1
```

Confirm that the earlier job has stopped before restarting. Use a new output directory if you change the input reads or the intended experiment. 

## Summarized interpretation

This comparison tests whether Nanopore links improve an Illumina-based assembly of the same natural whey community. The short-read run provides a baseline, while the hybrid run adds information that may help resolve repeat-associated graph paths. Evaluate contig count and N50 together with recovery of each published genome, strain-specific sequence, read support, and possible structural errors. Assembly span should reflect a mixture of genomes and mobile elements, not the expected size of one bacterium. The four dominant bacterial genomes provide useful comparisons, but their successful recovery does not indicate recovery of the complete bacterial community. 

Adding Nanopore reads improved assembly contiguity. Among contigs at least 500 bp long, N50 increased from 8,911 to 14,368 bp, maximum contig length increased from 53,777 to 134,674 bp, and contig count decreased from 1,306 to 951. Including all contigs, total assembly span increased from 6.15 to 6.34 Mb. However, reference genome fraction decreased from 80.05% to 78.84%, while duplication ratio increased from 1.005 to 1.146 and reported misassemblies increased from 19 to 82. Mismatches increased from 131.98 to 383.31 per 100 kbp, and indels increased from 12.26 to 131.89 per 100 kbp. These results indicate improved contiguity accompanied by reduced reference recovery and greater disagreement with the references, warranting investigation of sequence accuracy, redundant sequence, and strain resolution. Hybrid scaffolds contained 261,035 Ns compared with 4,800 in the Illumina-only scaffolds, indicating that some additional scaffold continuity reflects gaps rather than resolved sequence. Although strain differences can contribute to reference disagreements, these results do not support calling the hybrid assembly better overall.
## Further reading

- [SPAdes 4.3.0 release](https://github.com/ablab/spades/releases/tag/v4.3.0)
- [metaSPAdes paper](https://doi.org/10.1101/gr.213959.116)
- [SPAdes options and metagenomic mode](https://ablab.github.io/spades/running.html)
- [SPAdes hybrid assembly](https://ablab.github.io/spades/hybrid.html)
- [SPAdes output descriptions](https://ablab.github.io/spades/output.html)
- [Somerville et al. (2019): Long-read based assembly of low-complexity metagenomes](https://doi.org/10.1186/s12866-019-1500-0)
- [NWC2 BioSample SAMN09580370](https://www.ncbi.nlm.nih.gov/biosample/SAMN09580370)
- [NWC2 BioProject PRJNA477604](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA477604)
- [Illumina run SRR7585899](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7585899)
- [Nanopore run SRR7585900](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7585900)
- [NCBI prefetch and fasterq-dump instructions](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
- [fastp documentation](https://github.com/OpenGene/fastp)
- [fastplong documentation](https://github.com/OpenGene/fastplong?utm_source=chatgpt.com)
- [QUAST and MetaQUAST documentation](https://github.com/ablab/quast)

[Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
