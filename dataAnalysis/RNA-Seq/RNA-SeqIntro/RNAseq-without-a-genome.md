---
title: RNA Sequence Analysis
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

Lets now assume that *Arabidopsis* doesn't have a sequenced genome. We then start with the RNAseq reads and assemble them *denovo* into transcripts. One such denovo assembler, that we will showcase here is [__Trinity__](https://github.com/trinityrnaseq/trinityrnaseq/wiki). It incorporates three software modules, Inchworm, Chrysalis and Butterfly in sequence. It divides the data to many smaller de bruijn graphs, each representing the transcriptional complexity at a locus.


There are many ways to run Trinity. The official download page is [here](https://github.com/trinityrnaseq/trinityrnaseq/releases).
Alternatively, you can import the official [Trinity Docker Image](https://hub.docker.com/r/trinityrnaseq/trinityrnaseq/) into a singularity image and run Trinity in a Singularity container. On Ceres singularity 2.4 is installed on nodes running linux Cent OS 7.4. This could be done from your project folder as follows.
*Note: The `pull` or import of the docker image could have also been done in the Trinity SBATCH script, but we demonstrate it separately just to show the process of building a singularity image. The pull command first saves a series of tar.gz files in the cache folder and then sequentially builds a squash image that we can use after mounting out working directory in the container, as shown in the SBATCH script.*

```
salloc -p debug74
singularity pull docker://trinityrnaseq/trinityrnaseq
WARNING: pull for Docker Hub is not guaranteed to produce the
WARNING: same image on repeated pull. Use Singularity Registry
WARNING: (shub://) to pull exactly equivalent images.
Docker image path: index.docker.io/trinityrnaseq/trinityrnaseq:latest
Cache folder set to /home/sivanandan.chudalayandi/.singularity/docker
Importing: base Singularity environment
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:ce640ed7800c29984dfd6c4cda9e5b6c759f630356e972c851125713628cdfd9.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:486cb8339a27635fa93dc47aa0c689326a0a7cce388966d16daf8d265436cf7f.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:dc6f0d824617ad8a5d1163a5b2084814665dd83156317ad06ccf14deb517a053.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:4f7a5649a30e3f318ce5d7e4dbcbbeb6c0938c4cbae4d4a641fe910562ff4978.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:672363445ad2c734e29221a6b47f4e614b5adc8a3cdca3364f62db2ed2bdff0c.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:fa5e03a7cea7fda9a438856574913bd4895091aa6db87feccabe56bab35de10d.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:60af5100289ac8a096793b2050c0dbb5040958661bd76bce44fe655d7800747b.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:71afd2414e9ba47b9b913096a234df60ec70cdff8e0d9dc136d849ade83ee1de.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:cfbe046dbb7a2d00a21e47136f0f65371cb44d471975a4010bec6503fce1fdde.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:ad7309631e7f9d0fde782912642154010ee230dcec0b31bda6262a6200860aaa.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:aa1d8e07c016851c7389bcefdeca671613fd5549cc8cb3716ec762bd7ecd6fb3.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:310877e7e48932846dd5fd88b594f1e7a0fa31a61401e10b31f71a09ee3841d5.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:229d926bd456b21601dcbb0b0148d092849fca1e589a149b141a13cb71f89667.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:dc3bf8f1f122e7c99d765a93aaec3db3011328e30d48f60c80f1e12d0794e223.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:4fe0a71b4fba93a388f18e50329f22dd078b4023640d5c0cddc17d873022870b.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:dd4a7aeb28caa9da4be497de7f7b65b61d790ae5ef4498b3fe7dd39768663beb.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:5aedeb9210d7db405b9e99bbefa2c64d9d05824ab067b831fdb95310193debe5.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:351c35ee3e8dcb6de90d3d2b81192e981660deee900baf94a41a144528222128.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:ae9af67d4f2aebce65b1ab4693adb4aee2e815146683052ee098f736b8c84b15.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:813e007a02f1d3b6d14180e90fd1de398cd453723474b566c6fc39d27e5e7182.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:38a8a24f9344eb0ded75c1005b17bd7163c88d1c61ef793122854bf9d6f5eae8.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:614634ba4025c2e0b2bf903367e353e4e75c68f60a5580513dcbe860c2ae73c9.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:245e88eb46fba5920fba64b8bc380e80ad844762ec8de67a408ac883cc79f2b4.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:56450ceb44fb47ca94a76e764eeb9c934712e0b15f8da1f54dd11c6e70be8941.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:923481697283cf0cb7b2bf4129428cd92ad6be9fb485701ef51d53c59d1352eb.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:68537cb83e308b416841a65c9061ef49615d005a800a7d647375de1cb65d191e.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:c9b6f189c784e7aa69ffd5f961d2f984210ccc160d3e8746b00ec0546bafeab4.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:19675e3d081b3482fcb1d982253203008212efeb32892edd8b1e79110d47c6e8.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:24549bf1ae8912a9d30c8efe69ffbc8979e67639a70bbfd7f5cba3f8b338f55b.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:16c5b4d662c99b94d81cf5902853bdbf98204bbb713d3df8254178151698defb.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:f2121024d426ab5040d03fa93fab380b1c296349f794767fc92a9d991143b4e7.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:d17c033ab3431aed136729f91794647c0d306c73991b9452f393ae1b1ee3f55d.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:717c8003e47e39f79d049eff0e82e25c106d959a4a9ebf7203932c9ca1a84456.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:6494a085bd622288605daaeab68ad49f1f00441358897520387fe1e2e81b724c.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/docker/sha256:76390ecbff1f342992c0b76eeab5bc68d90299cdf580ccff67fa9db4590d21c5.tar.gz
Importing: /home/sivanandan.chudalayandi/.singularity/metadata/sha256:caf5a16c57d3307da698a52ea556c652d0fd6a7f8fc292a52cbd7418e63c5aac.tar.gz
WARNING: Building container as an unprivileged user. If you run this container as root
WARNING: it may be missing some functionality.
Building Singularity image...
        Singularity container built: ./trinityrnaseq.img
Cleaning up...

exit # exit from the relevant node
```

## Script for Trinity assembly:
```
#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug74
#SBATCH --ntasks-per-node=40 # For Trinity always reserve the whole node
#SBATCH -t 96:00:00
#SBATCH -J tri
#SBATCH -o tri.o%j
#SBATCH -e tri.e%j
#SBATCH --mail-user=<user_email_address>
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
cd $SLURM_SUBMIT_DIR
ulimit -s unlimited
scontrol show job $SLURM_JOB_ID

# Trinity requires that the each set of reads are concatenated into a file each. make combined left and right reads.

cat ../*1.*gz > left_1.gz
cat ../*2.*gz > right_2.gz

# running Trinity after mounting our working directory inside the container using $PWD.

singularity exec --bind $PWD trinityrnaseq.img Trinity --seqType fq --max_memory 120G --CPU 16 --normalize_by_read_set --output TrinityOut --left left_1.gz --right right_2.gz --trimmomatic
```
After the complete run. We see the complete *de novo* assembly (trinity.fa) in the output directory, TrinityOut. We can perform a variety of downstream analyzes with this transcriptome assembly. For the purposes of this tutorial, we will demonstrate mapping the RNAseq reads back to the assembly using bowtie2, calculating transcript abundance, using FeatureCounts and then performing differential expression using DESeq2. We will use scripts that come packaged with Trinity.



We will now build a transcriptome index and align the RNAseq back to transcriptome and estimate the abundance of the transcripts. When using the singularity image, the absolute path to the folder containing the relevant scripts must be given as follows:


```
singularity exec --bind $PWD trinityrnaseq-2.5.0.img /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts path_to_assembly/TrinityOut/Trinity.fasta --seqType fq --left path_reads/93_1.gz --right transcriptomic_fastq/RNAseq/paired-end/Arabidopsis/93_2.gz --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir <RSEM_dir1> >& 93_AE_bt2.log &
```
The transcriptome index will be built in the folder containing the transcriptome. Also because we specify __--trinity_mode__, a gene to transcript map file is also prepared, which can be used to produce gene level counts in addition to transcript counts.

```
./Trinity.fasta.gene_trans_map
./Trinity.fasta.bowtie2.ok
./Trinity.fasta.bowtie2.4.bt2
./Trinity.fasta.bowtie2.3.bt2
./Trinity.fasta.bowtie2.1.bt2
./Trinity.fasta.bowtie2.rev.1.bt2
./Trinity.fasta.bowtie2.rev.2.bt2
./Trinity.fasta.bowtie2.2.bt2
```

We can run align_and_estimate_abundance.pl script for each set of reads. We store the bowtie2 alignment file and abundance estimates in separate directories for each set of reads.


Now we build a transcript expression matrix for all samples using the *abundance_estimates_to_matrix.pl* script. We will prefix all the output files with "all" and use the base name of the directory as sample names.
```
singularity exec --bind $PWD trinityrnaseq-2.5.0.img /usr/local/bin/trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method RSEM --name_sample_by_basedir --gene_trans_map transcriptomic_fastq/RNAseq/paired-end/trinity/TrinityOut/Trinity.fasta.gene_trans_map --out_prefix all RSEM_9*/*isoforms.results >& matrix1.log&
```
The matrix files generated are prefixed with "all" as under:
```
ls all*

all.gene.TMM.EXPR.matrix
all.gene.TPM.not_cross_norm.TMM_info.txt
all.isoform.TMM.EXPR.matrix
all.isoform.TPM.not_cross_norm.TMM_info.txt
all.gene.counts.matrix
all.gene.TPM.not_cross_norm
all.isoform.counts.matrix
all.isoform.TPM.not_cross_norm
```


```
head all.gene.counts.matrix
                       RSEM_93 RSEM_94 RSEM_95 RSEM_96 RSEM_97 RSEM_98
TRINITY_DN0_c0_g1       18.13   0.00    9.66    62.63   25.29   19.41
TRINITY_DN0_c0_g2       0.00    8.85    4.62    21.24   0.00    0.00
TRINITY_DN0_c0_g3       0.00    0.00    34.33   83.32   39.72   0.00
TRINITY_DN10000_c0_g1   37.74   3.54    13.85   88.67   31.02   32.56
TRINITY_DN10001_c0_g1   286.46  192.12  178.70  538.88  491.25  341.86
TRINITY_DN10001_c0_g2   505.17  151.43  275.19  743.57  380.79  553.27
TRINITY_DN10002_c0_g1   0.00    0.00    0.00    54.31   0.00    0.00
TRINITY_DN10003_c0_g1   0.00    43.16   0.00    22.91   0.00    0.00
TRINITY_DN10005_c0_g1   0.00    0.00    0.00    8.13    13.44   0.00
```
This matrix can be imported into R and differential expression analyses performed using DESeq2 as explained [here](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/Differential-Expression-Analysis).

---
[Table of contents](Differential-Expression-Analysis)
