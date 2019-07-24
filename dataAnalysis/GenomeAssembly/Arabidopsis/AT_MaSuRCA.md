---
title: MaSuRCA Genome Assembly of Arabidopsis
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


### Background about the data Analysis

As with many examples in this workbook, we will explore publicly available data related to *Arabidopsis thaliana* which has a genome size of approximately 135 Mb.  

###  How to download the data from SRA
The data for this analysis is taken from the following NCBI project.  It contains both long and short reads.  We will only be using a fraction of these files so that the analysis will complete quickly but the full run could also performed with these data.



#### Let's start by organizing our data analysis folder

* make a project directory called masurcatutorial

  ```bash
  mkdir masurcatutorial
  cd masurcatutorial
  ```

* make a directory for your raw data, let's call this step 00

  ```bash
  mkdir 00_rawdata
  cd 00_rawdata
  ```

### Download the data

* Load the sra tool kit and download SRA file and convert to fastq file

  ```bash
  module load sra-toolkit
  ```

  the --split-files parameter is used to separate the paired end data into two files rather than the default single file.

    [more info about SRA downloading found in this section of the workbook](https://isugenomics.github.io/bioinformatics-workbook/dataAcquisition/fileTransfer/sra.html)

  Download the data into the masurcatutorial/00_rawdata folder  

* Downloading from SRA will be performed using the [sra-toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc).

  create a file with the SRA ids and name it `srr.ids`

  ```
  SRR3156163
  SRR3156596
  SRR3157034
  SRR3166543
  ```

* Assuming the sra-toolkit is installed then load the module and run the following bash script on the command line.

  ```bash
  module load sra-toolkit
  while read line; do
    fastq-dump --split-files --origfmt ${line};
  done<srr.ids
  ```

* The pacbio data has to be downloaded separately

    ```
    fastq-dump --table SEQUENCE --origfmt SRR3156160
    ```
  **Note:** sra-toolkit will create a folder named ```ncbi`` in your home directory ```/home/userid/ncbi```  If you have a disk storage limit on your home directory (most supercomputers do), you will want to move that folder to a different location and then create a softlink in your home folder.

  **Error example:** 2019-04-16T19:43:49 fastq-dump.2.8.1 err: unknown while writing file within file system module - unknown system error errno='Disk quota exceeded(122)'

* Check the file sizes to see if they are all present and of reasonable size

  ```bash
  ls -lha | awk '{print $5, $NF}'
  125M SRR3156160_1.fastq
  8.1G SRR3156163_1.fastq
  8.1G SRR3156163_2.fastq
  218M SRR3156596_1.fastq
  218M SRR3156596_2.fastq
  1.2G SRR3157034_1.fastq
  1.2G SRR3157034_2.fastq
  71M SRR3166543_1.fastq
  71M SRR3166543_2.fastq
  ```

### MaSuRCA requires a Configuration file

In this example, we have 9 fastq files 8 files from a single Hi-Seq Illumina run and a single file from a PacBio run (SRR3156163_1.fastq  SRR3156163_2.fastq  SRR3156596_1.fastq  SRR3156596_2.fastq  SRR3157034_1.fastq  SRR3157034_2.fastq  SRR3166543_1.fastq  SRR3166543_2.fastq, SRR3156160_1.fastq).

* **sr_config.txt**

  ```bash
  DATA
  PE = pa 250 50 SRR3166543_1.fastq  SRR3166543_2.fastq
  PE = pb 250 50 SRR3157034_1.fastq  SRR3157034_2.fastq
  JUMP = ja 8000 1600 SRR3156163_1.fastq  SRR3156163_2.fastq
  JUMP = jb 20000 4000 SRR3156596_1.fastq  SRR3156596_2.fastq
  PACBIO = SRR3156160_1.fastq
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

### Create the sr_config.txt file and place the above text into it.

* make a folder for the masurca run, we will call it 01_masurca

  ```bash
  cd ..
  mkdir 01_masurca
  cd 01_masurca

  #add the text above into the sr_config.txt file
  vim sr_config.txt

  #Create softlinks to your raw data in the 01_masurca folder
  ls /fullpath2analysisfolder/00_rawdata/ | xargs -I xx ln -s "/fullpath2analysisfolder/00_rawdata/"xx
  ```


### Running MaSuRCA using a Singularity container.

In this example we use a MaSuRCA container.  See the main MaSuRCA tutorial for more information if you have not already downloaded the git repository and pulled the singularity container from singularity hub.  You can also use your locally installed MaSuRCA instead.

- **runMaSuRCA.sh**

  This script can be found in the repo you just cloned under bin in the git repo.  MASURCA in the script below is the wrapper to call the singularity container.  Delete that MASURCA from the beginning of each line if you are running using your local installation.

  ```bash
  #!/bin/bash

  module load singularity

  MASURCA masurca sr_config.txt
  MASURCA ./assemble.sh
  ```



### Example SLURM Job

* **masurca_AT.sub**

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

### Submit SLURM job

* batch script

  ```bash
  sbatch masurca_AT.sub
  ```

  **NOTE:** prior to running the following files should be in your directory

  ```
  SRR3156160_1.fastq
  SRR3156163_1.fastq
  SRR3156163_2.fastq
  SRR3156596_1.fastq
  SRR3156596_2.fastq
  SRR3157034_1.fastq
  SRR3157034_2.fastq
  SRR3166543_1.fastq
  SRR3166543_2.fastq
  masurca_AT.sub
  sr_config.txt
  ```

### This what to expect in your output folder from a successful run


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


### Assembly statistics

This script was originally part of the [assemblathon paper](https://gigascience.biomedcentral.com/articles/10.1186/2047-217X-2-10) that we keep locally in our shop to explore the quality of genomes.  You can download it to your machine from our [ISUGIFsingularity/utilities repo](https://github.com/ISUGIFsingularity/utilities/blob/master/utilities/new_Assemblathon.pl).  You will want to change the first line to #!/bin/perl.

* AT_0.o4290911

```bash
  module load perl/5.18.4-threads
  perl new_Assemblathon.pl final.genome.scf.fasta


  Information for assembly final.genome.scf.fasta


                                           Number of scaffolds       6037
                                       Total size of scaffolds  100660109
                                              Longest scaffold     296097
                                             Shortest scaffold        309
                                   Number of scaffolds > 1K nt       5543  91.8%
                                  Number of scaffolds > 10K nt       2909  48.2%
                                 Number of scaffolds > 100K nt         62   1.0%
                                   Number of scaffolds > 1M nt          0   0.0%
                                  Number of scaffolds > 10M nt          0   0.0%
                                            Mean scaffold size      16674
                                          Median scaffold size       9375
                                           N50 scaffold length      32181
                                            L50 scaffold count        878
                                                   scaffold %A      31.98
                                                   scaffold %C      18.05
                                                   scaffold %G      18.03
                                                   scaffold %T      31.94
                                                   scaffold %N       0.00
                                           scaffold %non-ACGTN       0.00
                               Number of scaffold non-ACGTN nt          0

                  Percentage of assembly in scaffolded contigs       0.0%
                Percentage of assembly in unscaffolded contigs     100.0%
                        Average number of contigs per scaffold        1.0
  Average length of break (>25 Ns) between contigs in scaffold          0

                                             Number of contigs       6037
                                Number of contigs in scaffolds          0
                            Number of contigs not in scaffolds       6037
                                         Total size of contigs  100660109
                                                Longest contig     296097
                                               Shortest contig        309
                                     Number of contigs > 1K nt       5543  91.8%
                                    Number of contigs > 10K nt       2909  48.2%
                                   Number of contigs > 100K nt         62   1.0%
                                     Number of contigs > 1M nt          0   0.0%
                                    Number of contigs > 10M nt          0   0.0%
                                              Mean contig size      16674
                                            Median contig size       9375
                                             N50 contig length      32181
                                              L50 contig count        878
                                                     contig %A      31.98
                                                     contig %C      18.05
                                                     contig %G      18.03
                                                     contig %T      31.94
                                                     contig %N       0.00
                                             contig %non-ACGTN       0.00
                                 Number of contig non-ACGTN nt          0
  ```



### SLURM output file



This file provides time stamps of the steps that were run with MaSuRCA.  When you have more data this may come in handy to determine what step you are on.   This small dataset took about 33 minutes to run.  Your assembly stats may vary slightly if you rerun it multiple times.

* AT_0.o4290911

  ```bash
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


### Errors

It appears the do_consensus.sh to perform gapfilling is currently having an issue that should be fixed in new releases.  This step can be performed manually.  [See here for more information](https://github.com/alekseyzimin/masurca/issues/53)

The assembly is correct but it hasn't been gap filled.

[Arabidopsis data set Info](Arabidopsis_background.md)
[Back to the Assembly and Annotation Index page](../../GenomeAnnotation/annotation_and_assembly_index.md)
