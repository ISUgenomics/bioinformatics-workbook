# Running Quiver on the assembled genome (PacBio)

For running Quiver, you will need a finished assembly file (`fasta`) format and `bax.h5` reads from SMRTcells. You will need to run `quiver` for every single file of `bax.h5` separately. To make this more effecient, you can use `GNU parallel`.

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
cmph5tools.py merge --outFile out_all.cmp.h5 ${aligned[@]}
cmph5tools.py sort --inPlace --deep out_all.cmp.h5
samtools faidx ${ref} 
quiver out_all.cmp.h5 -j ${cpu} -r ${ref} -o ${ref%.*}_polished.fasta
```

Finally run these scripts in PBS/SLURM job script. For the first script:






