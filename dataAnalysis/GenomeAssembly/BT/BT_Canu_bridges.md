---
title: Canu genome assembly of Bacillus thuringiensis on XSEDE bridges
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


  Bridges is part of the NSF funded XSEDE resources that anyone can apply for supercomputing (See [XSEDE tutorial](../../../Appendix/HPC/xsede/xsede.md) on how to get an allocation).


## Running on the large memory node

Due to some issues with the how SLURM is set up on bridges trying to run on the Regular Memory nodes (RM) will result in submission failure of the canu sub jobs. To overcome this issue, the following script will run on a single Large Memory node (XLM) rather than running on a grid. (Note: there is a workable solution on RM nodes below).  You will want to change the parameters and names accordingly in the canu line

```canu -s psc_lm.spec -p whiteabaloneBM -d $local_outdir genomeSize=1.8g correctedErrorRate=0.16  -nanopore-raw $local_disk/Abalone_nanopore.fastq```

This particular dataset had low coverage so I increased the correctedErrorRate to 0.16.  Your genome size, and prefix names will be different.  Also don't forget to update the rsync line.

```rsync -aP psc_lm.spec ../data/Abalone_nanopore.fastq $local_disk/```

#### canu.sub

```

#!/bin/bash
#SBATCH -p XLM
#SBATCH -t 14-00:00:00
#SBATCH --mem 12000GB
#SBATCH --mail-type=END
#SBATCH --mail-user=YOUREMAILADDRESS

module load gcc/7.3.0
module load perl/5.18.4-threads
module load java/jdk8u73
module load canu/1.7

jobid=$SLURM_JOBID
local_disk=/local/$jobid
local_outdir=$local_disk/Abalone_canu1_out
cd $SLURM_SUBMIT_DIR
rsync -aP psc_lm.spec ../data/Abalone_nanopore.fastq $local_disk/
cd $local_disk
canu -s psc_lm.spec -p whiteabaloneBM -d $local_outdir genomeSize=1.8g correctedErrorRate=0.16  -nanopore-raw $local_disk/Abalone_nanopore.fastq
```

In the same folder you will also need to create the following canu specification file.


#### psc_lm.spec

```
useGrid=0
#maxMemory=900
#maxThreads=100
ovsMethod=sequential
```


## Running on Regular Memory nodes

Download and install canu locally in your folder from (https://github.com/marbl/canu/releases).  I may create a singularity container for this and will post it here when I do.

Modify the following two files (Execution.pm and Grid_Slurm.pm), which can be found in the canu install directory (canu-1.7.1/Linux-amd64/lib/site_perl/canu).  These lines deal with how canu submits jobs to the SLURM scheduler and does not affect how canu assembles a genome.  I prefix my edits with something I can easily find later (GIFedit).

#### Execution.pm

```
#GIFedit $opts .= "$memOption "      if (defined($memOption));
#GIFedit $opts .= "$thrOption "      if (defined($thrOption));
```

#### Grid_Slurm.pm

```
#    setGlobalIfUndef("gridEngineMemoryOption",               "--mem-per-cpu=MEMORY");
```

#### Canu Submission Script on bridges for regular memory nodes

Be sure to modify the canu execution folder, prefix, error rate and input data.  This cannot just be cut and paste.

```
#!/bin/bash
#SBATCH -J canu4
#SBATCH -o canu4.o%j
#SBATCH -p RM
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -t 48:00:00

source /usr/share/Modules/init/bash


module load java
module load gcc
module load perl/5.18.4-threads

/pylon5/mc48o5p/severin/whiteAbaloneGenome/canu/canu-1.7.1/Linux-amd64/bin/canu \
 -p whiteabalone3 -d whiteAb-oxford4 \
 genomeSize=1.8g \
 gridOptions="-N 1 -p RM -n 28 -t 48:00:00" \
 correctedErrorRate=0.105 \
 -nanopore-raw Abalone_nanopore.fastq
```

#### Errors

When I first tried to run on bridges RM nodes I would get the following error.

```
CRASH:
CRASH: Canu 1.7
CRASH: Please panic, this is abnormal.
ABORT:
CRASH:   Failed to submit batch jobs.
CRASH:
CRASH: Failed at /opt/packages/canu/canu-1.7/Linux-amd64/bin/../lib/site_perl/canu/Execution.pm line 1101.
CRASH:  canu::Execution::submitOrRunParallelJob('whiteabalone', 'meryl', 'correction/0-mercounts', 'meryl', 1) called at /opt/packages/canu/canu-1.7/Linux-
amd64/bin/../lib/site_perl/canu/Meryl.pm line 538
CRASH:  canu::Meryl::merylCheck('whiteabalone', 'cor') called at /opt/packages/canu/canu-1.7/Linux-amd64/bin/canu line 622
CRASH:
CRASH: Last 50 lines of the relevant log file (correction/0-mercounts/meryl.jobSubmit-01.out):
CRASH:
CRASH:
```

When this error is followed up by looking in the ```meryl.jobSubmit-01.out``` you find that Slurm submission fails.

```

```

---

* [Bacillus thuringiensis data set Info](BT_background.md)
* [Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
