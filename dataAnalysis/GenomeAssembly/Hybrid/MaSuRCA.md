# Genome Assembly with hybrid reads using MaSuRCA

### Background about the data Analysis

As with many examples in this workbook, we will explore publicly available data related to *Arabidopsis thaliana* which has a genome size of approximately 135 Mb.  

### Download the dataset

The data for this analysis is taken from the following NCBI project.  It contains both long and short reads.  We will only be using a fraction of these files so that the analysis will complete quickly but the full run could also performed with these data.

| SeqType              | Platform | ReadType   | BioProject                                                                                | Experiment                                                                                                                                                                 |
|----------------------|----------|------------|-------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Long Reads           | PacBio   | long-reads | [PRJNA314706](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA314706)                  | [Diploid Arabidopsis thaliana genome sequencing and assembly](https://www.ncbi.nlm.nih.gov/pubmed/27749838)              


The last column links directly to the article written about these data and experiment.



#### Let's start by organizing our data analysis folder,

* make a project directory called masurcatutorial
```
mkdir masurcatutorial
cd masurcatutorial
```

* make a directory for your raw data, let's call this step 00
```
mkdir 00_rawdata
cd 00_rawdata
```

#### Load the sra tool kit and download SRA file and convert to fastq file

*
```
module load sra-toolkit
```

* the --split-files parameter is used to separate the paired end data into two files rather than the default single file.
* [more info about SRA downloading found in this section of the workbook](https://github.com/ISUgenomics/bioinformatics-workbook/blob/3eae26f5ca6e7f16325372abdb8bbd9052b425d8/dataAcquisition/sra.md)
* Download the data into the masurcatutorial/00_rawdata folder  
```
# Illumina MiSeq data
fastq-dump --split-files SRR3703081
# Pacbio data
fastq-dump SRR3405330
```

While we are using MiSeq and PacBio data for this example, it is more common to have HiSeq and PacBio or HiSeq and Nanopore data.

### Background reading about the Assembling program MaSuRCA

* [MaSuRCA  home page](http://www.genome.umd.edu/masurca.html)
* [MaSuRCA  Paper](https://academic.oup.com/bioinformatics/article/29/21/2669/195975)
* [MaSuRCA github and Manual](https://github.com/alekseyzimin/masurca)


### Background on the required Configuration file

MaSuRCA requires a configuration file that states the type of data you have and the parameters you want to use in your run. See the github page for more information on the configuration file specific to your particular input files.  Here is a run down of some of the more common features.

* [MaSuRCA github and Manual](https://github.com/alekseyzimin/masurca)

In this example, we have three fastq files 2 files from a single Mi-Seq Illumina run and a single file from a PacBio run (SRR3405330.fastq  SRR3703081_1.fastq  SRR3703081_2.fastq).

* **sr_config.txt**

```
DATA
PE = pa 250 50 SRR3703081_1.fastq  SRR3703081_2.fastq
PACBIO = SRR3405330.fastq
CLOSE_GAPS=1
END

PARAMETERS
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 0
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.15
KMER_COUNT_THRESHOLD = 1
NUM_THREADS = 28
JF_SIZE = 20000000000
SOAP_ASSEMBLY=0
END
```

### Create the sr_config.txt file and place the above text into it.

```
# make a folder for the masurca run, we will call it 01_masurca
mkdir 01_masurca
cd 01_masurca
vim sr_config.txt

#Create softlinks to your raw data in the 01_masurca folder
ls /fullpath2analysisfolder/00_rawdata/ | xargs -I xx ln -s "/fullpath2analysisfolder/00_rawdata/"xx
```

If you have jump libraries or Nanopore you can enter them into the config file with lines that look similar to this.

```
JUMP = jb 15000 1000 Unicorn_R1_001.fastq Unicorn_R2_001.fastq
NANOPORE = Unicorn_NP.fastq
```

As stated in the [MaSuRCA Manual](https://github.com/alekseyzimin/masurca), here are some of the more common points in the manual to be aware of.

```
# DATA is specified as type {PE,JUMP,OTHER,PACBIO} and 5 fields:
1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads

#whether to attempt to close gaps in scaffolds with Illumina data
CLOSE_GAPS=1

PacBio/MinION data are supported. Note that you have to have 50x + coverage in Illumina Paired End reads to use PacBio of Oxford Nanopore MinION data. Supply PacBio or MinION reads (cannot use both at the same time) in a single fasta file as:

PACBIO=file.fa or NANOPORE=file.fa
```


### Running MaSuRCA using a Singularity container.

You don't have to run MaSuRCA using a singularity container as it may be already installed on your machine or you can install it locally using conda.  The benefit of a container is that it is fully reproducible by anyone else by providing them the container.  To make containers easier, we have created


#### Clone the repository

```
git clone https://github.com/ISUGIFsingularity/masurca.git
```
This repo contains the [singularity recipe](https://github.com/ISUGIFsingularity/masurca/blob/master/Singularity.1.0.0) if you want to build it from scratch.  However, You can also directly download the image from [singularity hub](https://www.singularity-hub.org/collections/1814).

```
cd masurca/SIMG
module load singularity
singularity pull shub://ISUGIFsingularity/masurca:1.0.0
# create a softlink if the singularity pull results in a file other than ISUGIFsingularity-masurca-master-1.0.0.simg as my wrappers are expecting the later and not the former
ln -s masurca_1.0.0.sif ISUGIFsingularity-masurca-master-1.0.0.simg
```

#### Set variables ```gitmasurca``` and ```PATH``` so we can find the container and the runscripts

modify and add the following lines to your .bashrc file.

* masurcagit is the path to the repo you just cloned.

```
export masurcagit="/home/severin/isugif/masurca/"
export PATH=$masurcagit/wrappers/:$PATH
```

- **runMaSuRCA.sh**

This script can be found in the repo you just cloned under bin

```
#!/bin/bash

module load singularity

MASURCA masurca sr_config.txt
MASURCA ./assemble.sh
```



### Example SLURM Job

* **masurca_AT.sub**

```
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 96:00:00
#SBATCH -J AT_0
#SBATCH -o AT_0.o%j
#SBATCH -e AT_0.e%j
#SBATCH --mail-user=youremail@email.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd $SLURM_SUBMIT_DIR
ulimit -s unlimited

module load singularity

$masurcagit/bin/runMaSuRCA.sh

scontrol show job $SLURM_JOB_ID

```

### Submit SLURM job

```
sbatch masurca_AT.sub
```
