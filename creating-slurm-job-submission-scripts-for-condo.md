# Creating SLURM job submission scripts for Condo

For a quick list of frequently used SLURM commands/options, take a look at this [cheat sheet](https://gif.biotech.iastate.edu/slurm-slurm-job-management-cheat-sheet "SLURM job management cheat sheet"). This guide will show you how to generate SLURM job submission scripts for your commands, easily and efficiently!

## 1. Generate commands ##

We recommend using our run scripts. Run scripts are configured to run frequently used programs with settings that are tailored to fit our needs. The idea is to make it simple to run frequently used programs without having to read the entire manual. Some of the common run scripts are `runBWA.sh`, `runGSNAP.sh`, `runBLASTn` etc. Here we will take [runBLASTn.sh](https://github.com/ISUgenomics/common_analyses/blob/master/runBLASTn.sh "run script for blastn") as an example, but you can use this for almost any of the other run scripts as well as your own commands. First, lets generate the commands:

We have 5 fasta files: `file1.fsa`, `file2.fsa`, `file3.fsa`, `file4.fsa` and `file5.fsa` (each with 1000 sequences) for which we need to run BLAST against `nr` database, we will generate the commands as follows:

```
for file in file?.fsa; do
echo "./runBLASTn.sh $file;
done > blastn.cmds
```

## 2. Create submission files ##

The `blastn.cmds` now contains 5 lines, each set to run blastn on each of the file. To create slurm script for each one of these lines, we run the [makeSLURMs.sh](https://github.com/ISUgenomics/common_scripts/blob/master/makeSLURMs.py "creating slurm scripts") script as follows:
```
makeSLURMs.py 1 blastn.cmds
```
Here 1 is to tell that put one job per submission file, and `blastn.cmds` is the commands file that looks like this:
```
./runBLASTn.sh file1.fsa
./runBLASTn.sh file2.fsa
./runBLASTn.sh file3.fsa
./runBLASTn.sh file4.fsa
./runBLASTn.sh file5.fsa
```
Once the `makeSLURMs.py` command is run, you should see 5 `.slurm` files, each identified by numbers and the command file prefix. The content should look like this:
```
#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 96:00:00
#SBATCH -J blastn_1
#SBATCH -o blastn_1.o%j
#SBATCH -e blastn_1.e%j
#SBATCH --mail-user=${USER}@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
module use /work/GIF/software/modules
./unBLASTn.sh file1.fsa
scontrol show job $SLURM_JOB_ID
```
## 3. Submit the jobs ##

The next step is to submit the jobs to the cluster. Using the for loop
```
for f in blastn*.slurm; do
sbatch $f;
done
```
All commands will be submitted and will start to run as per the availability of nodes.