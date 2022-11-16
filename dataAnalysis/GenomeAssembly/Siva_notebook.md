#### Genome Assembly:
```
[sivanandan.chudalayandi@ceres-dtn-1 mosquito]$ pwd
/90daydata/isu_gif_vrsc/Siva/mosquito
[sivanandan.chudalayandi@ceres-dtn-1 mosquito]$ more download.txt
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR121/085/SRR12121585/SRR12121585_subreads.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR121/086/SRR12121586/SRR12121586_subreads.fastq.gz
[sivanandan.chudalayandi@ceres-dtn-1 mosquito]$ more download.sh
parallel <download.txt
[sivanandan.chudalayandi@ceres-dtn-1 mosquito]$
```

* Install hifiasm

```
1040  cd /project/isu_gif_vrsc/Siva/software/
 1041  git init
 1042  git clone https://github.com/chhylp123/hifiasm
 1043  ls
 1044  cd hifiasm/
 1045  make
 1046  ls
 1047  pwd
 1048  nano ~/.bashrc
 1049  nano ~/.bash_profile
 1050  nano ~/.bashrc
 1051  pwd
 1052  /project/isu_gif_vrsc/Siva/software/hifiasm/hifiasm
 1053  nano ~/.bashrc
 1054  source ~/.bashrc
 1055  hifiasm
 1056  history

```

```
[sivanandan.chudalayandi@ceres19-compute-35 mosquito]$ pwd
/project/isu_gif_vrsc/Siva/Service/mosquito
[sivanandan.chudalayandi@ceres19-compute-35 mosquito]$ hifiasm -o anopheles.asm -t20 SRR12121585_subreads.fastq.gz SRR12121586_subreads.fastq.gz >& anopheles.log&
```

Files produced:

```
anopheles.log
anopheles.asm.a_ctg.lowQ.bed
anopheles.asm.a_ctg.gfa
anopheles.asm.a_ctg.noseq.gfa
anopheles.asm.p_ctg.lowQ.bed
anopheles.asm.p_ctg.gfa
anopheles.asm.p_ctg.noseq.gfa
anopheles.asm.p_utg.lowQ.bed
anopheles.asm.p_utg.noseq.gfa
anopheles.asm.p_utg.gfa
anopheles.asm.r_utg.lowQ.bed
anopheles.asm.r_utg.gfa
anopheles.asm.r_utg.noseq.gfa
anopheles.asm.ovlp.reverse.bin
anopheles.asm.ovlp.source.bin
anopheles.asm.ec.bin
```
* Get primary contigs in fasta:
`awk '/^S/{print ">"$2;print $3}' anopheles.asm.p_ctg.gfa > anopheles.asm.p_ctg.fa`

_weirdly, I am unable to see anything in /project/isu_gif_vrsc/Siva/software/_ though I had git cloned hifiasm there.

###### For Sylvie.Quiniou `https://usda-scinet.atlassian.net/browse/VRSC-3020`

```
6:05 AM; Sat Feb 20; 2021
# wonder why inputs are fasta?
# why using conda... when installation seems to be very simple
# test is working saw log files with kmer frequencies being written from 6:45 am

these are her conda commands and the batch job submision:

module load miniconda
source activate hifiasm
sbatch hifiasmdiploid.sh



[root@ceres19-compute-30 hifiasm_check]# pwd
/lustre/project/isu_gif_vrsc/Siva/Service/hifiasm_check
[root@ceres19-compute-30 hifiasm_check]# ls
hifiasm.man  SQ_hifiasm  test.log
[root@ceres19-compute-30 hifiasm_check]# ls -lth
total 4.0K
-rw-r--r--. 1 root proj-isu_gif_vrsc    0 Feb 20 06:05 test.log
-rw-r--r--. 1 root proj-isu_gif_vrsc 2.8K Feb 20 05:58 hifiasm.man
lrwxrwxrwx. 1 root proj-isu_gif_vrsc   39 Feb 20 05:16 SQ_hifiasm -> /project/tgbs_analysis/HIFIASM-CCSCCBL1
[root@ceres19-compute-30 hifiasm_check]#

# got kicked out of the job around 2 pm.

```

#### Nov 15/16 2022

Ecoli data for Genome Assembly:

Publication: https://www.tandfonline.com/doi/full/10.1080/19490976.2021.2014772
Bioproject: PRJNA725420
Details: https://ncbi.nlm.nih.gov/bioproject/?term=PRJNA725420
Data: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/043/SRR14417043/SRR14417043_subreads.fastq.gz
