# Running Quiver on the assembled genome (PacBio)

For running Quiver, you will need a finished assembly file (`fasta`) format and `bax.h5` reads from SMRTcells. You will need to run `quiver` for every single file of `bax.h5` separately. To make this more effecient, you can use `GNU parallel`.

First, set up a `runQuiver.sh` script as follows:

```
#!/bin/bash
if [ $# -lt 1 ] ; then
        echo ""
        echo "usage: runQuiver.sh <your.bas.h5> <genome.fasta>"
        echo "runs quiver for getting the genome consensus"
        echo ""
        exit 0
fi
module load smrtanalysis/64/2.3.0.140936
module load parallel


```




