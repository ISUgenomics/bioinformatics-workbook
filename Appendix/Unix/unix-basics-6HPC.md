---
title: "HPC Cluster Basics"
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# HPC cluster basics


### SLURM: Simple Linux Utility for Resource Management

* A simple text file with all the requirements for running your job
  * Memory requirement
  * Desired number of processors
  * Length of time you want to run the job
  * Type of queue you want to use (optional)
  * Where to write output and error files
  * Name foryour job while running on HPC



### Job Script Basics

A typical job script will look like this:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --mail-user=netid@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=JobName.%J.err
#SBATCH --output=JobName.%J.out

cd $SLURM_SUBMIT_DIR

module load modulename

your_commands_goes_here
```

Lines starting with `#SBATCH` are for `SLURM` resource manager to request resources for HPC. Some important options are as follows:


<table >
<thead><tr><th>Option</th><th>Examples</th><th>Description</th></tr></thead><tbody>
 <tr><td><blockcode>--nodes</blockcode></td><td><blockcode>#SBATCH --nodes=1</blockcode></td><td>Number of nodes</td></tr>
 <tr><td><blockcode>--cpus-per-task</blockcode></td><td><blockcode>#SBATCH --cpus-per-task=16</blockcode></td><td>Number of CPUs per node</td></tr>
 <tr><td><blockcode>--time</blockcode></td><td><blockcode>#SBATCH --time=HH:MM:SS</blockcode></td><td>Total time requested for your job</td></tr>
 <tr><td><blockcode>--output</blockcode></td><td><blockcode>#SBATCH -output filename</blockcode></td><td>STDOUT to a file</td></tr>
 <tr><td><blockcode>--error</blockcode></td><td><blockcode>#SBATCH --error filename</blockcode></td><td>STDERR to a file</td></tr>
 <tr><td><blockcode>--mail-user </blockcode></td><td><blockcode>#SBATCH --mail-user user@domain.edu</blockcode></td><td>Email address to send notifications</td></tr>
</tbody></table>

### Job Management Commands

<table>
<thead><tr><th>Job Status</th><th>Commands</th></tr></thead><tbody>
 <tr><td><blockcode>sinfo -a</blockcode></td><td>list all queues</td></tr>
 <tr><td><blockcode>squeue </blockcode></td><td>list all jobs</td></tr>
 <tr><td><blockcode>squeue -u userid</blockcode></td><td>list jobs for userid</td></tr>
 <tr><td><blockcode>squeue -t R</blockcode></td><td>list running jobs</td></tr>

</tbody></table>


Let's go ahead and give these job management commands a try.

```
sinfo -a
squeue
squeue -t R
#pick a name you saw when you typed squeue and specify all the jobs by that person with the following option
squeue -u first.lastname
```

There can be a lot of information using those two commands. I have created some useful alias' that change the output to something more informative.

```
alias sq='squeue -o "%8i %12j %4t %10u %20q %20a %10g %20P %10Q %5D %11l %11L %R"'
alias si='sinfo -o "%20P %5D %14F %8z %10m %10d %11l %16f %N"'
```

Where `(A/I/O/T)` = `available/idle/other/total`

You can place those alias' into your `~/.bashrc` file and it will automatically load every time you log in.

#### <span style="color:Green">Exercise:</span> Add these two alias' above to your `~/.bashrc` file
```
nano ~/.bashrc
```

### Job scheduling commands
<table>
<thead><tr><th>Commands</th><th>Function</th><th>Basic Usage</th><th>Example</th></tr></thead><tbody>
 <tr><td><blockcode>sbatch</blockcode></td><td>submit a slurm job</td><td>sbatch [script]</td><td>$ sbatch job.sub</td></tr>
 <tr><td><blockcode>scancel</blockcode></td><td>delete slurm batch job</td><td>scancel [job_id]</td><td>$ scancel 123456</td></tr>
</tbody></table>



### Interactive Session

To start a interactive session execute the following:

```
# this command will give 1 Node with 1 cpu in the brief-low queue for a time of 00 hours: 01 minutes: 00 seconds

salloc -N 1 -n 1 -p brief-low -t 00:01:00

# You can exit out of an interactive node by typing exit and hitting return
```

Interactive sessions are very helpful when you need more computing power than your laptop or desktop to wrangle the data or to test new software prior to submitting a full batch script.
