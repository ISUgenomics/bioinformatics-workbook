## HPC cluster basics

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


<table >
<thead><tr><th>Option</th><th>Examples</th><th>Description</th></tr></thead><tbody>
 <tr><td><blockcode>--nodes</blockcode></td><td><blockcode>#SBATCH --nodes=1</blockcode></td><td>Number of nodes</td></tr>
 <tr><td><blockcode>--cpus-per-task</blockcode></td><td><blockcode>#SBATCH --cpus-per-task=16</blockcode></td><td>Number of CPUs per node</td></tr>
 <tr><td><blockcode>--time</blockcode></td><td><blockcode>#SBATCH --time=HH:MM:SS</blockcode></td><td>Total time requested for your job</td></tr>
 <tr><td><blockcode>--output</blockcode></td><td><blockcode>#SBATCH -output filename</blockcode></td><td>STDOUT to a file</td></tr>
 <tr><td><blockcode>--error</blockcode></td><td><blockcode>#SBATCH --error filename</blockcode></td><td>STDERR to a file</td></tr>
 <tr><td><blockcode>--mail-user </blockcode></td><td><blockcode>#SBATCH --mail-user user@domain.edu</blockcode></td><td>Email address to send notifications</td></tr>
</tbody></table>


### Job scheduling commands
<table>
<thead><tr><th>Commands</th><th>Function</th><th>Basic Usage</th><th>Example</th></tr></thead><tbody>
 <tr><td><blockcode>sbatch</blockcode></td><td>submit a slurm job</td><td>sbatch [script]</td><td>$ sbatch job.sub</td></tr>
 <tr><td><blockcode>scancel</blockcode></td><td>delete slurm batch job</td><td>scancel [job_id]</td><td>$ scancel 123456</td></tr>
 <tr><td><blockcode>scontrol hold</blockcode></td><td>hold slurm batch jobs</td><td>scontrol hold [job_id]</td><td>$ scontrol hold 123456</td></tr>
 <tr><td><blockcode>scontrol release </blockcode></td><td>release hold on slurm batch jobs</td><td>scontrol release  [job_id]</td><td>$ scontrol release  123456</td></tr>
</tbody></table>

### Job management commands

<table>
<thead><tr><th>Job Status</th><th>Commands</th></tr></thead><tbody>
 <tr><td><blockcode>sinfo -a</blockcode></td><td>list all queues</td></tr>
 <tr><td><blockcode>squeue </blockcode></td><td>list all jobs</td></tr>
 <tr><td><blockcode>squeue -u userid</blockcode></td><td>list jobs for userid</td></tr>
 <tr><td><blockcode>squeue -t R</blockcode></td><td>list running jobs</td></tr>
 <tr><td><blockcode>smap</blockcode></td><td> show jobs, partitions and nodes in a graphical network topology</td></tr>
</tbody></table>



### Interactive session

To start a interactive session execute the following:

```
#this command will give 1 Node for a time of 4 hours

srun -N 1 -t 4:00:00 --pty /bin/bash
```
