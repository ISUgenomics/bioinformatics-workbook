---
title: "Introduction to NextFlow"
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

## Learning Objectives

1. What is nextflow?
2. Installing nextflow
3. Basic nextflow commands to run a workflow

## Introduction

[Nextflow]() is a workflow framework that can be used by a bioinformatician to integrate all of her/his/their bash/python/perl/other scripts into a one cohesive pipeline that are portable,reproducible, scalable and checkpointed. Nextflow is its own DSL (Domain Specific Language) that extends a language called [groovy](https://groovy-lang.org) which is an extension of Java that has language feature similarity to Python, Ruby and Smalltalk.

## Nextflow features

* **Fast protyping** -- let's you write a computational pipeline from smaller tasks
* **Reproducibility** -- supports Docker and Singularity containers
* **Portable** -- can run locally, Slurm, SGE, PBS, and cloud (Google, Kubernetes and AWS)
* **Unified parallelism** -- can process chunks through the entire pipeline (QC -> align -> call snps)
* **Continuous checkpoints** -- each chunk and process it goes through is checkpointed
* **Stream oriented** -- promotes programming approach extending Unix pipes model.

## Prerequisites


* [How to install Java 8 on Mac](https://stackoverflow.com/questions/24342886/how-to-install-java-8-on-mac)

```
Note: Oracle Java 8/9/10 is no longer available for public download (license change).
First install and update brew from Terminal:
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew tap homebrew/cask-versions
brew update
NEW as of June 2019
To install the JDKs from AdoptOpenJDK:
brew tap adoptopenjdk/openjdk
... then one of the following:
brew cask install adoptopenjdk8
brew cask install adoptopenjdk9
brew cask install adoptopenjdk10
brew cask install adoptopenjdk11
```

* NextFlow


## Installation

* Linux

  * ```
  curl -s https://get.nextflow.io | bash
  ```

* Mac
  * [Install Xcode from Mac App Store](https://apps.apple.com/us/app/xcode/id497799835?mt=12)
  * [Install Java](https://www.java.com/en/download/)

  * ```
  brew install nextflow
  ```

* PC
  * Install the Unix subsystem
  * Install Java 8 or above
    * ```sudo apt get install openjdk-8-jdk```
  * Install Nextflow
    * ```
  curl -s https://get.nextflow.io | bash
  ```



## Running a nextflow workflow


Nextflow has many commands but we are going to focus on the run, pull and drop commands

To execute any of these commands, we type `nextflow ` then the command

```bash
Commands:
  clean         Clean up project cache and work directories
  clone         Clone a project into a folder
  config        Print a project configuration
  console       Launch Nextflow interactive console
  drop          Delete the local copy of a project
  help          Print the usage help for a command
  info          Print project and system runtime information
  kuberun       Execute a workflow in a Kubernetes cluster (experimental)
  list          List all downloaded projects
  log           Print executions log and runtime info
  pull          Download or update a project
  run           Execute a pipeline project
  self-update   Update nextflow runtime to the latest available version
  view          View project script file(s)
```

## nextflow run

`nextflow run` can be used on a nextflow `main.nf` script or can be run directly from github repository.  Let's try this using a nextflow blast script developed by Iowa State's Genome Informatic's Facility (ISUGIF) by calling it directly from the github repo.

**Command: Showing the usage statement**

```
nextflow run isugifNF/blast --help
```

**Output**

The output is a general usage statement that can be used in place of `--help` to run blast.

```
N E X T F L O W  ~  version 20.07.1
Launching `isugifNF/blast` [backstabbing_franklin] - revision: 11f393fd09 [master]
NOTE: Your local project version looks outdated - a different revision is available in the remote repository [89887cd5b9]
Usage:
      The typical command for running the pipeline is as follows:
      nextflow run parallelBLAST.nf --query QUERY.fasta --genome GENOME.fasta -profile local
      nextflow run parallelBLAST.nf --query QUERY.fasta --dbDir "blastDatabaseDirectory" --dbName "blastPrefixName" -profile local

      Mandatory arguments:
       --query                        Query fasta file of sequences you wish to BLAST
       --genome                       Genome from which BLAST databases will be generated
       or
       --query                        Query fasta file of sequences you wish to BLAST
       --dbDir                        BLAST database directory (full path required)
       --dbName                       Prefix name of the BLAST database
       -profile                       Configuration profile to use. Can use multiple (comma separated)
                                      Available: test, condo, ceres, local, nova

       Optional arguments:
       --outdir                       Output directory to place final BLAST output
       --outfmt                       Output format ['6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen frames salltitles qcovs']
       --options                      Additional options for BLAST command [-evalue 1e-3]
       --outfileName                  Prefix name for BLAST output [blastout]
       --threads                      Number of CPUs to use during blast job [16]
       --chunkSize                    Number of fasta records to use when splitting the query fasta file
       --app                          BLAST program to use [blastn;blastp,tblastn,blastx]
       --help                         This usage statement.

```



**Command: running the test example**

This nextflow script happens to have a test dataset to show it is functioning properly and is perfect to show how nextflow operates.

```
nextflow run isugifNF/blast -profile test,local
```

**OUTPUT**

The command will produce the following output which tells us that it was run locally and executed two processes `software_check` and `runBlast`.  The processes run were determined by the nextflow script from [isugifNF/blast](www.github.com/isugifNF/blast).

```
N E X T F L O W  ~  version 20.07.1
Launching `isugifNF/blast` [spontaneous_hawking] - revision: 11f393fd09 [master]
NOTE: Your local project version looks outdated - a different revision is available in the remote repository [89887cd5b9]
executor >  local (2)
executor >  local (2)
[f3/9756f3] process > software_check [100%] 1 of 1 ✔
[ea/7c08f2] process > runBlast (1)   [100%] 1 of 1 ✔
WARN: Task runtime metrics are not reported when using macOS without a container engine
```

This script produces a work directory which is default for all nextflow scripts and an `out_dir` directory which can be specified with `--outdir` (see help above)

Here is a more comprehensive structure of the files generated. The key here is that all of the processes are run in the work directory and are given unique hash numbers.  If a process fails, you can fix the error and rerun the previous command using `-resume` and it will pick up where it failed.

**WARNING:** `-resume` only has one `-` not two as it is a nextflow command parameter and not a nextflow script parameter (as in --help above)

```
tree
.
|-- out_dir
|   |-- blast_output_combined.txt
|   |-- report.html
|   `-- timeline.html
`-- work
    |-- 70
    |   `-- 8c3ea2bbba3d4f31647d6410a4eced
    |       |-- blast.log
    |       |-- blastout
    |       `-- headtest.1.fasta -> /Users/severin/nextflow/nextflow-tutorial/blasttest/work/d3/4f9a6a126c4f0a3d6779f975ebcb7c/headtest.1.fasta
    |-- d3
    |   `-- 4f9a6a126c4f0a3d6779f975ebcb7c
    |       `-- headtest.1.fasta
    `-- fe
        `-- b6bdf808a6ebfdd7ef94e702821ac0
            `-- software_check.txt

8 directories, 8 files
```


Not all nextflow scripts will have a `--help` function or a `-profile test`.  If they are written well with the intention of sharing they should.

### Want to try it on your own dataset?

If you want to give it a go on your own dataset use the examples provided in the `--help` usage statement

#### 1. Create the blast database for me

Creates a blast database from the --genome file and then performs blast using `--query`

```
nextflow run parallelBLAST.nf --query QUERY.fasta --genome GENOME.fasta -profile local
```

#### 2. Use the blast database I provide

Uses the BLAST database specified by `--dbName` and `--dbDir` and then performs blast using `--query`

```
nextflow run parallelBLAST.nf -query QUERY.fasta --dbDir "blastDatabaseDirectory" --dbName "blastPrefixName" -profile local
```

## nextflow pull command

When you want the latest version of the nextflow script that you pulled from a github repository, you can execute the the `nextflow pull`

```
nextflow pull isugifNF/blast
```

## nextflow drop command

Sometimes the local copy of the nextflow github repo will cant `fast forward` and it may be easier to remove the local copy and start fresh.

```
nextflow drop isugifNF/blast
```
Once it is dropped you can start again using

```
nextflow run isugifNF/blast
```
and it will repull the latest version.

Note: if the `nextflow drop` command does not work you can manually remove this repository.  It is located: `~/.nextflow/assets/isugifNF/blast` or more generically `~/.nextflow/assets/OrganisationName/repoName`

## How to resume a nextflow pipeline

If your pipeline doesn't finish or errors at a specific process it can be restarted from where it left off using the `-resume` nextflow parameter.  Note: there is only one `-` (dash) as it is a nextflow parameter and not a pipeline parameter that your nextflow script has defined.


## Further Reading and Resources

* [NextFlow curated list of tutorials](https://nf-co.re/usage/nextflow)
* [NextFlow on YouTube](https://www.youtube.com/channel/UCB-5LCKLdTKVn2F4V4KlPbQ)
* [NextFlow Patterns](https://github.com/nextflow-io/patterns/blob/master/docs/mock-dependency.adoc) are examples of more complex functionality.
* [Searching for a biocontainer](https://quay.io/search?q=blast)
* [nf-core](https://github.com/nf-core) github repo of nextflow workflows.
* [nf-core modules](https://github.com/nf-core/modules)

### Why nextflow?

* [Using rapid prototyping to choose a bioinformatics workflow management system](https://www.biorxiv.org/content/10.1101/2020.08.04.236208v1.full.pdf+html)



## [Creating Your own workflow](dataAnalysis/nextflow/02_creatingAworkflow.md)
