---
title: Introduction to MaSuRCA genome assembler
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---



## Background

* Main software website links
  * [MaSuRCA  home page](http://www.genome.umd.edu/masurca.html)
  * [MaSuRCA  Paper](https://academic.oup.com/bioinformatics/article/29/21/2669/195975)
  * [MaSuRCA github and Manual](https://github.com/alekseyzimin/masurca)


* #### Configuration file

   MaSuRCA requires a configuration file that states the type of data you have and the parameters you want to use in your run. See the github page for more information on the configuration file specific to your particular input files.  Here is a run down of some of the more common features.

   [MaSuRCA github and Manual](https://github.com/alekseyzimin/masurca)

   Here is an example of the config file (**sr_config.txt**)

     ```
     DATA
     PE = pa 250 50 SRR3166543_1.fastq  SRR3166543_2.fastq
     PE = pb 250 50 SRR3157034_1.fastq  SRR3157034_2.fastq
     JUMP = ja 8000 1600 SRR3156163_1.fastq  SRR3156163_2.fastq
     JUMP = jb 20000 4000 SRR3156596_1.fastq  SRR3156596_2.fastq
     PACBIO = SRR3405330.fastq
     END

     PARAMETERS
     GRAPH_KMER_SIZE = auto
     USE_LINKING_MATES = 0
     LIMIT_JUMP_COVERAGE = 300
     CA_PARAMETERS =  cgwErrorRate=0.15
     KMER_COUNT_THRESHOLD = 1
     NUM_THREADS = 28
     JF_SIZE = 200000000
     SOAP_ASSEMBLY=0
     END
     ```

  If you have jump libraries or Nanopore you can enter them into the config file with lines that look similar to this.

  ```bash
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

  More than one entry for each data type/set of files is allowed EXCEPT for PacBio/Nanopore data. That is if you have several pairs of PE fastq files, specify each pair on a separate line with a different two-letter prefix.
  ```


## How to run MaSuRCA

  * #### Using a Singularity container

  You don't have to run MaSuRCA using a singularity container as it may be already installed on your machine or you can install it locally using conda.  The benefit of a container is that it is fully reproducible by anyone else by providing them the container.  To make containers easier, we have created

  Clone the repository

  ```
  git clone https://github.com/ISUGIFsingularity/masurca.git
  ```
  This repo contains the [singularity recipe](https://github.com/ISUGIFsingularity/masurca/blob/master/Singularity.1.0.0) if you want to build it from scratch.  However, You can also directly download the image from [singularity hub](https://www.singularity-hub.org/).

  ```bash
  cd masurca/SIMG
  module load singularity
  singularity pull shub://ISUGIFsingularity/masurca:1.0.0
  # create a softlink if the singularity pull results in a file other than ISUGIFsingularity-masurca-master-1.0.0.simg as my wrappers are expecting the later and not the former
  ln -s masurca_1.0.0.sif ISUGIFsingularity-masurca-master-1.0.0.simg
  ```

  * #### Set variables ```gitmasurca``` and ```PATH``` so we can find the container and the runscripts

  modify and add the following lines to your .bashrc file.

  masurcagit is the path to the repo you just cloned.

  ```
  export masurcagit="/home/severin/isugif/masurca/"
  export PATH=$masurcagit/wrappers/:$PATH
  ```

  **runMaSuRCA.sh**

  This script can be found in the repo you just cloned under bin

  ```bash
  #!/bin/bash

  module load singularity

  MASURCA masurca sr_config.txt
  MASURCA ./assemble.sh
  ```



* #### Example SLURM Job

  **masurca_AT.sub**

  ```bash
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

## Expected files generated during assembly

|files|in|output |folder | |
|--|--|--|--|--|
|AT_0.e4290986|AT_0.o4290986|CA.mr.41.15.17.0.029|CA.mr.41.15.17.0.029.log|CA_dir.txt|
|--|--|--|--|--|
|ESTIMATED_GENOME_SIZE.txt|PLOIDY.txt|SRR3405330.fastq|SRR3703081_1.fastq|SRR3703081_2.fastq|
|assemble.sh|containees.txt|create_mega-reads.err|environment.sh|genome.uid|
|global_arrival_rate.txt|guillaumeKUnitigsAtLeast32bases_all.fasta|guillaumeKUnitigsAtLeast32bases_all.jump.fasta|guillaumeKUnitigsAtLeast32bases_all.41.fasta|k_u_hash_0|
|masurca_AT.sub|meanAndStdevByPrefix.pe.txt|mr.41.15.17.0.029.all.txt|mr.41.15.17.0.029.mr.txt|mr.41.15.17.0.029.single.txt|
|mr.41.15.17.0.029.sr.frg|mr.41.15.17.0.029.txt|mr.41.15.17.0.029.1.allowed|mr.41.15.17.0.029.1.fa|mr.41.15.17.0.029.1.frg|
|mr.41.15.17.0.029.1.mates.frg|mr.41.15.17.0.029.1.trims.txt|mr.41.15.17.0.029.1.unjoined.fa|mr.41.15.17.0.029.1.to_join.fa.tmp|mr.41.15.17.0.029.all_mr.fa|
|mr.41.15.17.0.029.all_mr.maximal.fa|mr.41.15.17.0.029.join_consensus.tmp|mr.41.15.17.0.029.maximal_mr.txt|pa.renamed.fastq|pacbio_nonredundant.fa|
|pe.cor.fa|pe.cor.tmp.log|pe_data.tmp|quorum.err|quorum_mer_db.jf|
|runCA.spec|sr_config.txt|super1.err|superReadSequences.named.fasta|tigStore.err|
|unitig_cov.txt|unitig_layout.txt|work1|work1_mr||

## SLURM standard output

  * This file provides time stamps of the steps that were run with MaSuRCA.  When you have more data this may come in handy to determine what step you are on.   This small dataset took about 33 minutes to run.  Your assembly stats may vary slightly if you rerun it multiple times.

  ```
  Verifying PATHS...
  jellyfish OK
  runCA OK
  createSuperReadsForDirectory.perl OK
  nucmer OK
  mega_reads_assemble_cluster.sh OK
  creating script file for the actions...done.
  execute assemble.sh to run assembly
  [Fri Oct 26 08:29:01 EDT 2018] Processing pe library reads
  [Fri Oct 26 08:30:02 EDT 2018] Average PE read length 249
  [Fri Oct 26 08:30:02 EDT 2018] Using kmer size of 127 for the graph
  [Fri Oct 26 08:30:02 EDT 2018] MIN_Q_CHAR: 33
  [Fri Oct 26 08:30:02 EDT 2018] Creating mer database for Quorum
  [Fri Oct 26 08:31:11 EDT 2018] Error correct PE
  [Fri Oct 26 08:36:18 EDT 2018] Estimating genome size
  [Fri Oct 26 08:37:06 EDT 2018] Estimated genome size: 147278815
  [Fri Oct 26 08:37:06 EDT 2018] Creating k-unitigs with k=127
  [Fri Oct 26 08:41:11 EDT 2018] Computing super reads from PE
  [Fri Oct 26 08:46:02 EDT 2018] Using CABOG from /opt/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/masurca-3.2.8-hgdshdatclkrv7kct6kqgfb3fa536iji/bin/../CA8/Linux-amd64/bin
  [Fri Oct 26 08:46:02 EDT 2018] Running mega-reads correction/assembly
  [Fri Oct 26 08:46:02 EDT 2018] Using mer size 15 for mapping, B=17, d=0.029
  [Fri Oct 26 08:46:02 EDT 2018] Estimated Genome Size 147278815
  [Fri Oct 26 08:46:02 EDT 2018] Estimated Ploidy 1
  [Fri Oct 26 08:46:02 EDT 2018] Using 28 threads
  [Fri Oct 26 08:46:02 EDT 2018] Output prefix mr.41.15.17.0.029
  [Fri Oct 26 08:46:02 EDT 2018] Pacbio coverage <25x, using the longest subreads
  [Fri Oct 26 08:46:02 EDT 2018] Reducing super-read k-mer size
  [Fri Oct 26 08:48:51 EDT 2018] Mega-reads pass 1
  [Fri Oct 26 08:48:51 EDT 2018] Running locally in 1 batch
  [Fri Oct 26 08:49:38 EDT 2018] Mega-reads pass 2
  [Fri Oct 26 08:49:38 EDT 2018] Running locally in 1 batch
  [Fri Oct 26 08:49:42 EDT 2018] Refining alignments
  [Fri Oct 26 08:49:43 EDT 2018] Joining
  [Fri Oct 26 08:49:43 EDT 2018] Gap consensus
  [Fri Oct 26 08:49:43 EDT 2018] Warning! Some or all gap consensus jobs failed, see files in mr.41.15.17.0.029.join_consensus.tmp, proceeding anyway, to rerun gap consensus erase mr.41.15.17.0.029.1.fa and re-run asse
  mble.sh
  [Fri Oct 26 08:49:43 EDT 2018] Generating assembly input files
  [Fri Oct 26 08:49:44 EDT 2018] Coverage of the mega-reads less than 5 -- using the super reads as well
  [Fri Oct 26 08:49:49 EDT 2018] Coverage threshold for splitting unitigs is 15 minimum ovl 250
  [Fri Oct 26 08:49:49 EDT 2018] Running assembly
  [Fri Oct 26 08:58:46 EDT 2018] Recomputing A-stat
  recomputing A-stat for super-reads
  [Fri Oct 26 09:01:15 EDT 2018] Mega-reads initial assembly complete
  [Fri Oct 26 09:01:15 EDT 2018] No gap closing possible
  [Fri Oct 26 09:01:15 EDT 2018] Removing redundant scaffolds
  [Fri Oct 26 09:02:07 EDT 2018] Assembly complete, final scaffold sequences are in CA.mr.41.15.17.0.029/final.genome.scf.fasta
  [Fri Oct 26 09:02:07 EDT 2018] All done
  [Fri Oct 26 09:02:07 EDT 2018] Final stats for CA.mr.41.15.17.0.029/final.genome.scf.fasta
  N50 32439
  Sequence 100640646
  Average 16714.9
  E-size 43482.6
  Count 6021
  JobId=4290986 JobName=AT_0
     UserId=severin(50922) GroupId=mc48o5p(15124) MCS_label=N/A
     Priority=2764 Nice=0 Account=mc48o5p QOS=rm
     JobState=RUNNING Reason=None Dependency=(null)
     Requeue=0 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
     RunTime=00:33:14 TimeLimit=2-00:00:00 TimeMin=N/A
     SubmitTime=2018-10-26T08:28:28 EligibleTime=2018-10-26T08:28:28
     StartTime=2018-10-26T08:28:54 EndTime=2018-10-28T08:28:55 Deadline=N/A
     PreemptTime=None SuspendTime=None SecsPreSuspend=0
     LastSchedEval=2018-10-26T08:28:54
     Partition=RM AllocNode:Sid=br005:26364
     ReqNodeList=(null) ExcNodeList=(null)
     NodeList=r203
     BatchHost=r203
     NumNodes=1 NumCPUs=28 NumTasks=16 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
     TRES=cpu=28,mem=123200M,node=1,billing/gpu=28
     Socks/Node=* NtasksPerN:B:S:C=16:0:*:* CoreSpec=*
     MinCPUsNode=16 MinMemoryNode=123200M MinTmpDiskNode=0
     Features=(null) DelayBoot=00:00:00

  ```



## Errors

* Consensus

    It appears the do_consensus.sh to perform gapfilling is currently having an issue that should be fixed in new releases.  This step can be performed manually.  [See here for more information](https://github.com/alekseyzimin/masurca/issues/53)

  The assembly is correct but it hasn't been gap filled.  This can occur on assemblies with lower levels of assembly coverage.


## Further Reading

  * Main software website links
    * [MaSuRCA  home page](http://www.genome.umd.edu/masurca.html)
    * [MaSuRCA  Paper](https://academic.oup.com/bioinformatics/article/29/21/2669/195975)
    * [MaSuRCA github and Manual](https://github.com/alekseyzimin/masurca)



## Tutorial examples with real datasets

* [Example with Arabidopsis](../Arabidopsis/AT_MaSuRCA.md)


---

* [Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
