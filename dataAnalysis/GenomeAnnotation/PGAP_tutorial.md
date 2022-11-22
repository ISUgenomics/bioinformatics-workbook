---
title: NCBI Prokaryotic Genome Annotation Pipeline (PGAP)
layout: single
author: Sharu Paul
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

## Tutorial for NCBI PGAP
The <b>NCBI Prokaryotic Genome Annotation Pipeline (PGAP)</b> is designed to annotate bacterial and archaeal genomes.


```
https://github.com/ncbi/pgap
```

## Installation

PGAP is available as a docker image and NCBI [Quickstart](https://github.com/ncbi/pgap/wiki/Quick-Start)
page has instructions for using PGAP with docker.
I used PGAP on Nova, which has singularity as a module. Singularity can be used
with docker containers.

```
module load singularity
```

Pull the installation file:

```
singularity pull https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py
```

Install PGAP

```
export PGAP_INPUT_DIR=/PATH/TO/YOUR/work/directory
./pgap.py --update
```

And test

```
./pgap.py -r -o mg37_results YOUR/PATH/TO/test_genomes/MG37/input.yaml
```

This may take a while to finish depending on node/memory used.

### Errors
I ran into two main errors in setup and running PGAP.

#### Out of memory
This happens when the node used is low on memory. According to NCBI, PGAP needs 32 GB of memory, better to have more.

#### Out of disk space
Change the installation directory to your work directory otherwise PGAP will install
in home directory (which usually has low memory if using HPC) by default and run out of space.
Alternatively, create a `.pgap` directory in your work directory and softlink it in your home directory.
You might still run into this error on HPC if you have not [set up your home directory](http://datascience.101workbook.org/06-IntroToHPC/00-HOME-DIRECTORY/00-setting-up-home-directory).

## Run PGAP
Three files are needed to run PGAP; Assembly fasta file, metadata YAML file, and a input YAML file. You can make the YAML files using a text editor.
For details on input files check NCBI's [input files](https://github.com/ncbi/pgap/wiki/Input-Files) page.

I used PGAP for annotation of a genome of an unknown species in genus Spirochaeta. Here is an example run:
1. Spirochaete.fasta
2. Metadata YAML file (#data_submol.yaml):

```
organism:
    genus_species: 'Spirochaeta'
```

Additional information can be added to this metadata file. Check NCBI's [input files](https://github.com/ncbi/pgap/wiki/Input-Files) page for an example.
3. Input YAML file (#Input.yaml):

```
fasta:
    class: File
    location: /PATH/TO/GENOME/ASSEMBLY/FILE/Spirochaete.fasta
submol:
    class: File
    location: /PATH/TO/data_submol.yaml
```

Once you have the input files ready, run PGAP using following command. This is an example of a run without additional options, check more options using `./pgap.py -h` command.

```
./pgap.py -r -o output_directory_name Input.yaml
```

This can take hours depending on size of the genome and the memory allocated.
Best to run the job using slurm on HPC, i.e., create a job script and submit.


[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)
