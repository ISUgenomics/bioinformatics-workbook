---
title: "BLAST"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Running NCBI-BLAST jobs in parallel
If there is a large file of sequences, then the traditional way for doing a BLAST search is to start with the first sequence and run them sequentially for all the sequences. This is very time consuming and waste of the processing power HPCs offer. There are several ways to speed up this process. Most of them split the input sequence file into multiple pieces (usually equal to the number of processors) and run the BLAST search simultaneously on each of the split file. So larger the number of processors, more faster the whole process.
## 1. Using GNU Parallel to split the large job to smaller chunks
One easy way is to use the ''parallel'' command from the [GNU Parallel](http://www.gnu.org/software/parallel).
To do this:
```bash
cat large.fasta | \
parallel -S :,server1,server2 \
--block 100k \
--recstart '>' \
--pipe blastp \
-evalue 0.01 \
-outfmt 6 \
-db db.fa
-query - > combined_results.txt
```

Here, the large file is divided into blocks of 100kb, making sure that each 'piece' starts with the symbol `>`. The total number of pieces is dependent on total number of processors (from all nodes combined) eg., (standard `fasta` format). If there are 64 processors and 2 nodes, there will be a total of 128 sequence file pieces. On each of these pieces, `blastp` program is called with options. Note that you need to follow your `query` with a dash `-` indicating that the input is coming from the stdout, rather than the file. The results are written to the final file `combined_results.txt`. The order of the results will vary depending on what sequence went through the `BLAST`pipe and doesn't match the input order.
## 2. Splitting input query in to smaller pieces

For this you have to use an external script (for splitting the input file). One such script can be found [here](https://github.com/ISUgenomics/common_scripts/blob/master/fasta-splitter.pl) on GitHub. The steps to run this are as follows. First set up a blast script (your favorite flavor against your chosen database. In this case `blastp` against `swissprot-db`). Name this file as ` runBLASTp.sh`
```bash
#!/bin/bash
# perfoms NR blast (blastx)
infile="$1"
outfile="$(basename "${infile%.*}").xml"
database="/data021/GIF/arnstrm/Baum/GenePrediction_Hg_20160115/05_databases/swissprot/uniprot_sprot"
module load ncbi-blast
blastp \
 -query "${infile}" \
 -db "${database}" \
 -out "${outfile}".xml \
 -evalue 1e-5 \
 -num_threads 8 \
 -max_target_seqs 50 \
 -outfmt 5
```
Split the input sequences into desired number of pieces (recommended is ~10000 sequences per file, so that the run finishes in < 24hrs wih 8 processors, but you can adjust it as per your needs). Here we do 10 splits as follows
```
fasta-splitter.pl --n-parts 10--measure count input_seq.fasta
```

Set up a commands to run these splits with the `runBLASTp.sh` script

```bash
for file in input.part-??.fasta; do echo "./runBLASTp.sh $file"; done >> blast.cmds
cat blast.cmds
/runBLASTp.sh input.part-01.fasta
/runBLASTp.sh input.part-02.fasta
/runBLASTp.sh input.part-03.fasta
/runBLASTp.sh input.part-04.fasta
/runBLASTp.sh input.part-05.fasta
/runBLASTp.sh input.part-06.fasta
/runBLASTp.sh input.part-07.fasta
/runBLASTp.sh input.part-08.fasta
/runBLASTp.sh input.part-09.fasta
/runBLASTp.sh input.part-10.fasta
```
You can request 5 nodes, with 16 procs for this job on condo to run them all in parallel. You can do this with a job file as shown below (note: if you have more than 10 splits, request the number of nodes accordingly):
```bash
#!/bin/bash
#PBS -l nodes=5:ppn=16
#PBS -l walltime=168:00:00
#PBS -N BLAST
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID} -e ${PBS_JOBNAME}.e${PBS_JOBID}
#PBS -m ae -M $USER@gmail.com
cd $PBS_O_WORKDIR
ulimit -s unlimited
chmod g+rw ${PBS_JOBNAME}.[eo]${PBS_JOBID}
parallel --env TMPDIR --jobs 2 --sshloginfile $PBS_NODEFILE --joblog blast_progress.log --workdir $PWD < blast.cmds
```

 Finally, you can merge all the blast output files but using `cat` command.

 Note that if you're using large number of splits to run the BLAST and all the jobs are trying the access the database files, there might be a huge bottleneck reading and writing files
