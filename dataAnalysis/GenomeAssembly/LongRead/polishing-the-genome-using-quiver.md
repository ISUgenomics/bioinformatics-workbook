---
title: Genome Assembly
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Running Quiver on the assembled genome (PacBio)

For running Quiver, you will need a finished assembly file (`fasta`) format and `bax.h5` reads from SMRTcells. You will need to run `quiver` for every single file of `bax.h5` separately. To make this more efficient, you can use `GNU parallel`.

First, set up a `runQuiver_prepare.sh` script as follows:

```
#!/bin/bash
if [ $# -lt 1 ] ; then
        echo ""
        echo "usage: runQuiver_prepare.sh <your.bas.h5> <genome.fasta>"
        echo "runs quiver for getting the genome consensus"
        echo ""
        exit 0
fi
module load smrtanalysis/64/2.3.0.140936
module load parallel
bash5="$1"
ref="$2"
out=$(basename ${bash5%%.**})
pbalign --forQuiver ${bash5} ${ref} aligned_reads.${out}.cmp.h5
```

Second, set up a `runQuiver_polish.sh` script as follows:

```
#!/bin/bash
if [ $# -lt 1 ] ; then
        echo ""
        echo "usage: runQuiver_polish.sh <aligned_reads_prefix> <genome.fasta>"
        echo "runs quiver for getting the genome consensus"
        echo ""
        exit 0
fi
module load smrtanalysis/64/2.3.0.140936
module load parallel
module load samtools
aligned=($1*)
ref="$2"
cpu=16
cmph5tools.py merge --outFile out_all.cmp.h5 ${aligned[@]}
cmph5tools.py sort --inPlace --deep out_all.cmp.h5
samtools faidx ${ref}
quiver out_all.cmp.h5 -j ${cpu} -r ${ref} -o ${ref%.*}_polished.fasta
```

Finally run these scripts in PBS/SLURM job script. For the first script:

* Run `runQuiver_prepare.sh` for every `bas.h5` file to create aligned reads to the genome.
* Run `runQuiver_polish.sh` to merge all aligned files, sort, and use it to polish indexed reference genome.

For running the `runQuiver_prepare.sh`, you need to find all `bas.h5` files in the path, and prepare slurm script for each one of them. For preparing Slurm script, you can use [makeSLURM_ceres.py](https://github.com/ISUgenomics/common_scripts/blob/master/makeSLURM_ceres.py)

```
cd directories_with_smrtcell_data
# create command file
for bash5 in $(find $(pwd) -name "*.bas.h5"); do
echo "./runQuiver_prepare.sh $bash5 genome.fasta";
done > quiver_prepare.cmds
# create slurm script
makeSLURMs.py 1 quiver_prepare.cmds
# submit the jobs
for sub in *.sub; do
sbatch $sub;
done
```

Once this completes, you will have many `aligned_reads.<file_prefix>.cmp.h5` files. Now, we will run `runQuiver_polish.sh` on all these files:

```
# create command file
echo "./runQuiver_polish.sh aligned_reads_prefix genome.fasta" > quiver_polish.cmds
# create slurm script
makeSLURMs.py 1 quiver_polish.cmds
# submit the job
sbatch quiver_polish_0.sub
```

You only have run this one as it will use all the alignment files you created in the first step. When the job completes sucessfully, you should be seeing a `genome_polished.fasta` file.
