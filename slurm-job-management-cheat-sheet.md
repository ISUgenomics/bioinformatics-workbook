## Quick reference sheet for SLURM resource manager 

### Job scheduling commands

| Commands | Function | Basic Usage | Example |
| --- | --- | --- | --- |
| `sbatch` | submit a slurm job | sbatch [script] | $ sbatch job.sub |
| `scancel` | delete slurm batch job | scancel [job_id] | $ scancel 123456 |
| `scontrol hold` | hold slurm batch jobs | scontrol hold [job_id] | $ scontrol hold 123456 |
| `scontrol release` | release hold on slurm batch jobs | scontrol release [job_id] | $ scontrol release 123456 |

### Job management commands

| Job Status | Commands |
| --- | --- |
| `sinfo -a` | list all queues |
| `squeue` | list all jobs |
| `squeue -u userid` | list jobs for userid |
| `squeue -t R` | list running jobs |
| `smap` | show jobs, partitions and nodes in a graphical network topology |


### Job script basics

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

 
| Option | Examples | Description |
| --- | --- | --- |
| `--nodes` | `#SBATCH --nodes=1` | Number of nodes |
| `--cpus-per-task` | `#SBATCH --cpus-per-task=16` | Number of CPUs per node |
| `--time` | `#SBATCH --time=HH:MM:SS` | Total time requested for your job |
| `--output` | `#SBATCH -output filename` | STDOUT to a file |
| `--error` | `#SBATCH --error filename` | STDERR to a file |
| `--mail-user` | `#SBATCH --mail-user user@domain.edu` | Email address to send notifications |

### Interactive session 

To start a interactive session execute the following:

```
#this command will give 1 Node for a time of 4 hours
srun -N 1 -t 4:00:00 --pty /bin/bash 

```


### Aliases that provide useful information parsed from the Slurm commands

Bash
```
alias si="sinfo -o \"%20P %5D %14F %8z %10m %10d %11l %16f %N\""
alias sq="squeue -o \"%8i %12j %4t %10u %20q %20a %10g %20P %10Q %5D %11l %11L %R\""
```
More information about Slurm can be found here:

- [Slurm Cheat Sheat](https://www.chpc.utah.edu/presentations/SlurmCheatsheet.pdf)
- [Slurm Basics](http://researchit.las.iastate.edu/slurm-basics)
- [Slurm cluster commands](https://sites.google.com/a/case.edu/hpc-upgraded-cluster/slurm-cluster-commands)