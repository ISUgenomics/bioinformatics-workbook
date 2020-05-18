---
title: "Slurm Basics"
layout: single
author: Siva Chudalayandi
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Introduction to SLURM: Simple Linux Utility for Resource Management

* Open source fault-tolerant, and highly scalable cluster management and job scheduling system for large and small Linux clusters.
* HPC systems admins use this system for smooth resource distribution among various users. A user can submit jobs with specific resources to the centralized manager.

## The three objectives of SLURM:

* Lets a user request a compute node to do an analysis (job)
* Provides a framework (commands) to start, cancel, and monitor a job
* Keeps track of all jobs to ensure everyone can efficiently use all computing resources without stepping on each others toes.

## SLURM functions:

The main SLURM user commands shown on the left give the user access to information pertaining to the super computing cluster and the ability to submit or cancel a job.  See table below for a description of the main SLURM user functions.

|command | Description |
| - | - |
|sbatch | Submit a batch script to SLURM |
|squeue| List all jobs currently running or in queue |
|scancel| Cancel a job you submitted |
|sinfo| Check the availability of nodes within all partitions|
|scontrol | See the configuration of a specific node |
|sacct| Displays accounting data for all jobs |



<!-- ![](assets/Slurm_components.png) Photo from [schedmd](https://slurm.schedmd.com/overview.html) -->

To check the availability of nodes within all partitions:
```
$ sinfo
PARTITION        AVAIL  TIMELIMIT  NODES  STATE NODELIST
debug               up    1:00:00      1  maint ceres19-compute-26
debug               up    1:00:00      1    mix ceres14-compute-4
debug               up    1:00:00      1   idle ceres19-compute-25
brief-low           up    2:00:00      2  maint ceres19-compute-[26,40]
brief-low           up    2:00:00      1  down* ceres19-compute-37
brief-low           up    2:00:00     59    mix ceres18-compute-[0-17,19-27],ceres19-compute-[0-5,7-9,12,21-24,35-36,38-39,41-42,44-45,47,55-63]
brief-low           up    2:00:00      4  alloc ceres18-compute-18,ceres19-compute-[6,28,43]
brief-low           up    2:00:00     26   idle ceres19-compute-[10-11,13-20,25,27,29-34,46,48-54]
mem768-low          up    2:00:00      3   idle ceres18-mem768-0,ceres19-mem768-[0-1]
mem-low             up    2:00:00      3    mix ceres18-mem-[0-1],ceres19-mem-1
```
To check availability within a specific partitio, in this case the partition called `short`
```

$ sinfo -p short
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
short*       up 2-00:00:00      4  maint ceres14-compute-[1-3],ceres19-compute-40
short*       up 2-00:00:00      1  down* ceres19-compute-37
short*       up 2-00:00:00     47    mix ceres14-compute-[4-6,26-28,34-39,48-50,52-56,58-59,62,65],ceres18-compute-[24-27],ceres19-compute-[27,35-36,38-39,41-42,44-45,47,55-63]
short*       up 2-00:00:00     16  alloc ceres14-compute-[7,29,32-33,44-47,51,60-61,63-64,66-67],ceres19-compute-28
short*       up 2-00:00:00     32   idle ceres14-compute-[8-24],ceres19-compute-[29-34,43,46,48-54]
      up 7-00:00:00      1   idle ceres19-mem-4

```
To see the configuration of a specific node for example `ceres14-compute-8`

```
$ scontrol show nodes ceres14-compute-8
NodeName=ceres14-compute-8 Arch=x86_64 CoresPerSocket=10
   CPUAlloc=0 CPUTot=40 CPULoad=0.01
   AvailableFeatures=AVX
   ActiveFeatures=AVX
   Gres=(null)
   NodeAddr=ceres14-compute-8 NodeHostName=ceres14-compute-8 Version=19.05.5
   OS=Linux 3.10.0-1062.12.1.el7.x86_64 #1 SMP Tue Feb 4 23:02:59 UTC 2020
   RealMemory=126000 AllocMem=0 FreeMem=85536 Sockets=2 Boards=1
   State=IDLE ThreadsPerCore=2 TmpDisk=975 Weight=1 Owner=N/A MCS_label=N/A
   Partitions=short,geneious
   BootTime=2020-02-17T17:14:55 SlurmdStartTime=2020-02-18T17:12:06
   CfgTRES=cpu=40,mem=126000M,billing=40
   AllocTRES=
   CapWatts=n/a
   CurrentWatts=0 AveWatts=0
   ExtSensorsJoules=n/s ExtSensorsWatts=0 ExtSensorsTemp=n/s
```


If you want to check your jobs after submission:
```
squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           2867457     short P3826e00 sivanand  R   21:50:29      1 ceres14-compute-53
           2867458     short P6370337 sivanand  R   21:50:29      1 ceres14-compute-53
           2867459     short Pa0567fb sivanand  R   21:50:29      1 ceres19-compute-38
           2867456      long   Falcon sivanand  R   21:50:45      1 ceres14-compute-55
           2867883     short       sh sivanand  R      48:03      1 ceres14-compute-64
```
## Slurm batch script: Guidelines

* Number of nodes
* Desired number of processors or jobs
* Type of partition/queue you want to use (optional)
* Memory requirement (Optional)
* Length of time you want to run the job (Each partition has a default)
* Where to write output and error files
* Name for your job while running on HPC
* Email ID to get job status (Optional)

****One of the most important takeaways in this tutorial is that a job is best run on `compute nodes` and not on the `login node`.**** We generally write a batch script where we can reserve the necessary resources and then write the commands or the actual job that you want to do. Obviously this example is trivial, however in reality most jobs run by users involve atleast some component of heavy computing or memory.

## A typical Job Script
* Every line in the script that begins with "##" is a comment and is a comment or an explanation*

* This script invokes the unix sleep and echo commands*

```bash
#!/bin/bash
##The shebang line or the absolute path to the bash interpreter

## All the lines below that start with a single `#SBATCH` is a slurm directive

#SBATCH -N 1
## Reserve a single node
#SBATCH -n 4
##The job steps will launch a max of 4 jobs
#SBATCH -p short
## Reserve in the short partition
#SBATCH -t 01:00:00
## Reserve for 1 hour
#SBATCH -J sleep
## the name of the job is "sleep"
#SBATCH -o sleep.o%j
## write any std output to a file named sleep.ojobid..>
#SBATCH -e sleep.e%j
## write any std output to a file named sleep.ejobid..>
#SBATCH --mail-user=user@domain.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
## Notify by email at the above address, specifically when the job starts and when it ends

## The following lines are the commands
cd $SLURM_SUBMIT_DIR
scontrol show job $SLURM_JOB_ID
## scontrol above is a slurm command to view the slurm configuration or state. It is useful to see how much of the resources you have used.
sleep 10 && echo "I slept for 10 seconds"
sleep 20 && ech "I slept for 20 seconds"
## Note in the above line, I deliberately mis spelt `ech`; this would cause a std error to be output
sleep 60 && echo "I slept for 1 min"
```
* save the job script as `slurm.batch.sh`

****Note****: Lines starting with `#SBATCH` are for `SLURM` resource manager to request resources for HPC.

This script can be submitted as follows:

```
sbatch slurm.batch.sh
```
This job will at least run for 1-2 mins, so soon after submitting you can actually issue commands to see the job run.

```
squeue -u $USER
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           2935316     short    sleep sivanand  R       0:04      1 ceres14-compute-34
```
****Notes****:                                              . We are using the `-u` option for `squeue` and supplying the variable `$USER`, which referes to your ****user name****. We notice that the job, ****sleep****, is running on the node `ceres14-compute-34` in the `short` partition and has a job ID `2935316`.

Once the job is completed the following files appear
```
sleep.o2935316 # this is the standard output
sleep.e2935316 # this is the std error
```
Let's take a look at the std output file
```
more sleep.o2935316

JobId=2935316 JobName=sleep
   UserId=sivanandan.chudalayandi(1727000561) GroupId=sivanandan.chudalayandi(1727000561) MCS_label=N/A
   Priority=213721 Nice=0 Account=scinet QOS=memlimit
   JobState=RUNNING Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:01 TimeLimit=01:00:00 TimeMin=N/A
   SubmitTime=2020-05-18T10:40:25 EligibleTime=2020-05-18T10:40:26
   AccrueTime=2020-05-18T10:40:26
   StartTime=2020-05-18T10:40:26 EndTime=2020-05-18T11:40:26 Deadline=N/A
   PreemptEligibleTime=2020-05-18T10:40:26 PreemptTime=None
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2020-05-18T10:40:26
   Partition=short AllocNode:Sid=ceres19-ipa-0:39699
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=ceres14-compute-34
   BatchHost=ceres14-compute-34
   NumNodes=1 NumCPUs=4 NumTasks=4 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=4,mem=12400M,node=1,billing=4
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryCPU=3100M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=OK Contiguous=0 Licenses=(null) Network=(null)
   Command=/project/isu_gif_vrsc/Siva/Service/Slurm/slurm.batch.sh
   WorkDir=/project/isu_gif_vrsc/Siva/Service/Slurm
   StdErr=/project/isu_gif_vrsc/Siva/Service/Slurm/sleep.e2935316
   StdIn=/dev/null
   StdOut=/project/isu_gif_vrsc/Siva/Service/Slurm/sleep.o2935316
   Power=

I slept for 10 seconds
I slept for 1 min

```
****Note****: the line starting with `JobID` through `Power=` is the slurm configuration and state (`scontrol`) and gives you an idea of how many resources you have used as mentioned before. The last two lines are directly from our `echo` command in the script.

Additionally, the error file `sleep.e2935316`:

```
more sleep.e2935316
/var/spool/slurmd/job2935316/slurm_script: line 16: ech: command not found
```
This tells us that the command `ech` (deliberately mis-spelt) is not found.


### Interactive Session

We could have also run the commands in the job script interactively by first reserving a node in the partion using `salloc`

```
# this command will give 1 Node with 4 cpu in the short partitio for a time of 00 hours: 30 minutes: 00 seconds

$ salloc -N 1 -n 4 -p short -t 00:30:00
```
```
salloc: Pending job allocation 2935626
salloc: job 2935626 queued and waiting for resources
salloc: job 2935626 has been allocated resources
salloc: Granted job allocation 2935626
salloc: Waiting for resource configuration
salloc: Nodes ceres14-compute-48 are ready for job
export TMPDIR=/local/bgfs//2935626
export TMOUT=5400

```

In an interactive session, we can primarily use it to run small test runs of a large job and/or run say, a bunch of file `compression` or `un-tarring`.

We can run the commands from out job script above directly in the interactive session.

```
sleep 10 && echo "I slept for 10 seconds"
I slept for 10 seconds
```
or

```
sleep 20 && ech "I slept for 20 seconds"
bash: ech: command not found

```

## Additional Resources

* Other Slurm tutorials
  * [schedmd.com](https://slurm.schedmd.com/tutorials.html)
  * [Fulton supercomputing](https://www.youtube.com/watch?v=U42qlYkzP9k&feature=player_embedded)


## References

This tutorial is a rehash of material found on [schedmd](https://slurm.schedmd.com/overview.html)
