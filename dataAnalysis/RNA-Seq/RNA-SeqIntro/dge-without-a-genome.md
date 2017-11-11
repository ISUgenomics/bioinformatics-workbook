
Lets now assume that Arabidopsis doesnt have a sequenced genome. Then we will start with the reads and assemble them into transcripts. One such denovo assembler, that we will showcase here is ##Trinity


From Nathan Weeks: import the official Trinity Docker Image into a singularity image and run Trinity in a Singularity container:

```
salloc -p debug74
SINGULARITY_TMPDIR=$TMPDIR SINGULARITY_CACHEDIR=$TMPDIR singularity pull docker://trinityrnaseq/trinityrnaseq
....
....
....

Singularity container built: /scratch/sivanandan.chudalayandi/149863/trinityrnaseq.img

mv /scratch/sivanandan.chudalayandi/149863/trinityrnaseq.img .
```

## Script for Trinity assembly:
```
#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug74
#SBATCH --ntasks-per-node=16
#SBATCH -t 96:00:00
#SBATCH -J tri
#SBATCH -o tri.o%j
#SBATCH -e tri.e%j
#SBATCH --mail-user=csiva@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID
cat ../*1.*gz > left_1.gz
cat ../*2.*gz > right_2.gz

singularity exec --bind $PWD:/mnt --pwd /mnt /project/isu_gif_vrsc/Siva/transcriptomic_fastq/RNAseq/paired-end/trinity/trinityrnaseq.img Tri$
--seqType fq --max_memory 120G --CPU 16 --normalize_by_read_set --output TrinityOut --left left_1.gz --right right_2.gz --trimmomatic


```
Error: Unable to recognize read name formatting
```
Friday, November 10, 2017: 20:21:38     CMD: /usr/local/bin/trinityrnaseq/util/insilico_read_normalization.pl --seqType fq --JM 120G  --max_c
ov 50 --CPU 16 --output /mnt/TrinityOut/norm_for_read_set_1   --max_pct_stdev 10000  --left left_1.gz.PwU.qtrim.fq --right right_2.gz.PwU.qtr
im.fq --pairs_together --PARALLEL_STATS  
Converting input files. (both directions in parallel)CMD: seqtk-trinity seq -A /mnt/TrinityOut/left_1.gz.PwU.qtrim.fq >> left.fa
CMD: seqtk-trinity seq -A /mnt/TrinityOut/right_2.gz.PwU.qtrim.fq >> right.fa
Error, not recognizing read name formatting: [SRR4420293.1]

If your data come from SRA, be sure to dump the fastq file like so:

        SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra

Error, not recognizing read name formatting: [SRR4420293.1]

##If your data come from SRA, be sure to dump the fastq file like so:

        SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra

Thread 2 terminated abnormally: Error, cmd: seqtk-trinity seq -A /mnt/TrinityOut/right_2.gz.PwU.qtrim.fq >> right.fa died with ret 512 at /us
r/local/bin/trinityrnaseq/util/insilico_read_normalization.pl line 758.
Thread 1 terminated abnormally: Error, cmd: seqtk-trinity seq -A /mnt/TrinityOut/left_1.gz.PwU.qtrim.fq >> left.fa died with ret 512 at /usr/
local/bin/trinityrnaseq/util/insilico_read_normalization.pl line 758.
Error, conversion thread failed at /usr/local/bin/trinityrnaseq/util/insilico_read_normalization.pl line 329.
Error, cmd: /usr/local/bin/trinityrnaseq/util/insilico_read_normalization.pl --seqType fq --JM 120G  --max_cov 50 --CPU 16 --output /mnt/Trin
ityOut/norm_for_read_set_1   --max_pct_stdev 10000  --left left_1.gz.PwU.qtrim.fq --right right_2.gz.PwU.qtrim.fq --pairs_together --PARALLEL
_STATS   died with ret 6400 at /usr/local/bin/trinityrnaseq/Trinity line 2544.
        main::process_cmd("/usr/local/bin/trinityrnaseq/util/insilico_read_normalization"...) called at /usr/local/bin/trinityrnaseq/Trinity
line 3090
        main::normalize("/mnt/TrinityOut/norm_for_read_set_1", 50, ARRAY(0x22cc9b0), ARRAY(0x23139e8)) called at /usr/local/bin/trinityrnaseq
/Trinity line 3014
        main::run_normalization(50, ARRAY(0x1f89350), ARRAY(0x1f89398)) called at /usr/local/bin/trinityrnaseq/Trinity line 1297
```

So I decided to to use the -F option but still get the same error.

The error because Trinity expects left and right read names to be explicitly identified.
The following reads will result in error (note corresponding left and right reads have same names)

```

$zcat left_1.gz | head
@HWI-ST1136:361:HS250:1:1101:1130:2234
GCAACTTCTACAAGCAACTGGGTCGCTCTCTTCCATGACGACTCATCAGAACTGTCATACACGAAGACAGCTATGTCACATGCAGCTAAGGATTCCTTGG
+HWI-ST1136:361:HS250:1:1101:1130:2234
CCCFFFFFHHHHHJJJJJJJHIIJJJJJIJJJJIJJJJJJJJJJJJIJJJJJJJIIJJIJJJHHHFFFFEDEEEEEEEDDDDDDDDDDDDDDDDDEDDDD
@HWI-ST1136:361:HS250:1:1101:1409:2108
GGCTTTTCCTGCGCAGCTTAGGTGGAAGGCGAAGAAGGCCCCCTTCCGGGGGGCCCGAGCCATCAGTGAGATACTACTCTGGAAGAGCTAGAATTCTAAC
+HWI-ST1136:361:HS250:1:1101:1409:2108
CCCFFFFFHHHHHJJJJJJJJJFHJJJJJJJJJJJJJJJJJJJJJJJHHFDDDDDDDDDDDDDDDDEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDEEED

zcat right_2.gz | head
@HWI-ST1136:361:HS250:1:1101:1130:2234
AAGCGAGCTCCTGATCAGAGCTTTGAGTTGACAAATGCAGCAATCGATTTCCTGAAAGGAATGTATATGTTGTTTGATGACGACCAGGATAACAATCTGA
+HWI-ST1136:361:HS250:1:1101:1130:2234
CCCFFFFFHHHHHIJJJJJJIJJJJJJGIIJJJJIJJJIJJJJJJJJJJJJIJJJIJJJGJIIHJIIIJIGFFHHFFFFFFCCDDDDDDDDDDDDDDDED
@HWI-ST1136:361:HS250:1:1101:1409:2108
CGAGGAAACCTTTGCNNNNNTCCGTTACCTTTTGGGAGGCCTACGCCCCATAGAAACCGTCTACCTGAGACTGTCCCTTGGCCCGTAGGTCCTGACACAA
+HWI-ST1136:361:HS250:1:1101:1409:2108
<<<@@@@@@?@@@@@#####32@@@@@@???????????????????????????????==?????????=========<===<:::<<=========<<
```

We need to explicitly add /1 and /2 at the end of the left and right read names:

```
zcat left_1.gz | sed '/@HWI/ s_$_/1_' > correct_left_1.fq
zcat right_2.gz | sed '/@HWI/ s_$_/2_' > correct_right_2.fq
```

```
head correct_left_1.fq
@HWI-ST1136:361:HS250:1:1101:1130:2234/1
GCAACTTCTACAAGCAACTGGGTCGCTCTCTTCCATGACGACTCATCAGAACTGTCATACACGAAGACAGCTATGTCACATGCAGCTAAGGATTCCTTGG
+HWI-ST1136:361:HS250:1:1101:1130:2234
CCCFFFFFHHHHHJJJJJJJHIIJJJJJIJJJJIJJJJJJJJJJJJIJJJJJJJIIJJIJJJHHHFFFFEDEEEEEEEDDDDDDDDDDDDDDDDDEDDDD
@HWI-ST1136:361:HS250:1:1101:1409:2108/1
GGCTTTTCCTGCGCAGCTTAGGTGGAAGGCGAAGAAGGCCCCCTTCCGGGGGGCCCGAGCCATCAGTGAGATACTACTCTGGAAGAGCTAGAATTCTAAC
+HWI-ST1136:361:HS250:1:1101:1409:2108
CCCFFFFFHHHHHJJJJJJJJJFHJJJJJJJJJJJJJJJJJJJJJJJHHFDDDDDDDDDDDDDDDDEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDEEED

-$ head correct_right_2
@HWI-ST1136:361:HS250:1:1101:1130:2234/2
AAGCGAGCTCCTGATCAGAGCTTTGAGTTGACAAATGCAGCAATCGATTTCCTGAAAGGAATGTATATGTTGTTTGATGACGACCAGGATAACAATCTGA
+HWI-ST1136:361:HS250:1:1101:1130:2234
CCCFFFFFHHHHHIJJJJJJIJJJJJJGIIJJJJIJJJIJJJJJJJJJJJJIJJJIJJJGJIIHJIIIJIGFFHHFFFFFFCCDDDDDDDDDDDDDDDED
@HWI-ST1136:361:HS250:1:1101:1409:2108/2
CGAGGAAACCTTTGCNNNNNTCCGTTACCTTTTGGGAGGCCTACGCCCCATAGAAACCGTCTACCTGAGACTGTCCCTTGGCCCGTAGGTCCTGACACAA
+HWI-ST1136:361:HS250:1:1101:1409:2108
<<<@@@@@@?@@@@@#####32@@@@@@???????????????????????????????==?????????=========<===<:::<<=========<<
@HWI-ST1136:361:HS250:1:1101:1417:2131/2
CAAGGATACTGATATCTTGGCAGCATTCCGAGTAACTCCTCAACCTGGAGTTCCACCTGAAGAAGCAGGGGCTGCGGTAGCTGCTGAATCTTCTACTGGT
```
This solved the seqtk error, but raised a jellyfish error:
```
-------------------------------------------
----------- Jellyfish  --------------------
-- (building a k-mer catalog from reads) --
-------------------------------------------

CMD finished (0 seconds)
CMD: /usr/local/bin/trinityrnaseq/util/..//trinity-plugins/jellyfish/bin/jellyfish count -t 16 -m 25 -s 12947414142  --canonical  both.fa
Error, cmd: /usr/local/bin/trinityrnaseq/util/..//trinity-plugins/jellyfish/bin/jellyfish count -t 16 -m 25 -s 12947414142  --canonical  both.f
a died with ret 9 at /usr/local/bin/trinityrnaseq/util/insilico_read_normalization.pl line 758.
Error, cmd: /usr/local/bin/trinityrnaseq/util/insilico_read_normalization.pl --seqType fq --JM 120G  --max_cov 50 --CPU 16 --output /mnt/Trinit
yOut/norm_for_read_set_1   --max_pct_stdev 10000  --left correct_left_1.PwU.qtrim.fq --right correct_left_1.PwU.qtrim.fq --pairs_together --PAR
ALLEL_STATS   died with ret 512 at /usr/local/bin/trinityrnaseq/Trinity line 2544.
        main::process_cmd("/usr/local/bin/trinityrnaseq/util/insilico_read_normalization"...) called at /usr/local/bin/trinityrnaseq/Trinity li
ne 3090
        main::normalize("/mnt/TrinityOut/norm_for_read_set_1", 50, ARRAY(0x2145438), ARRAY(0x20fbc58)) called at /usr/local/bin/trinityrnaseq/T
rinity line 3014
        main::run_normalization(50, ARRAY(0x1db83f0), ARRAY(0x1db8438)) called at /usr/local/bin/trinityrnaseq/Trinity line 1297
gzip: correct_left_1.P.qtrim: No such file or directory
gzip: correct_left_1.U.qtrim: No such file or directory

```
