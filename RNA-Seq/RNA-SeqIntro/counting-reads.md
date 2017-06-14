====== Generating counts data ======
 
Trinity package provides align_and_estimate_abundance.pl program to generate abundance estimates. \\

This program takes the reads file and the transcriptome to generate counts. So, we have to generate such commands for every sample. (Some samples have SE data and some have PE data). To make it easier, I have separate commands to SE and PE datasets.



```

## Paired end 
ls WBC-B1* | grep -v 'L005' | paste - - | tr '_' ' ' | awk '{print $1"_"$2"_"$3"_"$4"_"$5,$6"_"$7"_"$8"_"$9"_"$10,$1}' | awk '{print "align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --trinity_mode --transcripts transcripts.fa \\";print "--thread_count 16 \\";print "--seqType fq --left "$1,"--right "$2" \\";print "--output_prefix "$3,"--output_dir "$3"  &";print "wait";print " "}' > WBC-B1-PE-commands.RSEM

ls WBC-B2* | grep -v 'L005' | paste - - | tr '_' ' ' | awk '{print $1"_"$2"_"$3"_"$4"_"$5,$6"_"$7"_"$8"_"$9"_"$10,$1}' | awk '{print "align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --trinity_mode --transcripts transcripts.fa \\";print "--thread_count 16 \\";print "--seqType fq --left "$1,"--right "$2" \\";print "--output_prefix "$3,"--output_dir "$3"  &";print "wait";print " "}' > WBC-B2-PE-commands.RSEM

## Single end
ls *L005_R1* |  tr '_' ' ' | awk '{print $1"_"$2"_"$3"_"$4"_"$5,$6"_"$7"_"$8"_"$9"_"$10,$1"SE"}' | awk '{print "align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --trinity_mode --transcripts transcripts.fa \\";print "--thread_count 16 \\";print "--seqType fq --single "$1," \\";print "--output_prefix "$3,"--output_dir "$3"  &";print "wait";print " "}' > WBC-SE-commands.RSEM

##Bean vs corn
ls WBC-mg* | paste - - | tr '_' ' ' | awk '{print $1"_"$2"_"$3"_"$4"_"$5,$6"_"$7"_"$8"_"$9"_"$10,$1}' | awk '{print "align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --trinity_mode --transcripts transcripts.fa \\";print "--thread_count 64 \\";print "--seqType fq --left "$1,"--right "$2" \\";print "--output_prefix "$3,"--output_dir "$3"  &";print "wait";print " "}' > WBC-mgcommands.RSEM

ls WBC-mg* | paste - - | tr '_' ' ' | awk '{print $1"_"$2"_"$3"_"$4"_"$5,$6"_"$7"_"$8"_"$9"_"$10,$1}' | awk '{print "align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --trinity_mode --transcripts transcripts.fa \\";print "--thread_count 64 \\";print "--seqType fq --left "$1,"--right "$2" \\";print "--output_prefix "$3,"--output_dir "$3"  &";print "wait";print " "}' > WBC-mgcommands.RSEM

```








```
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=usha@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=alignestimateWBC-B2.%J.err
#SBATCH --output=alignestimateWBC-B2.%J.out

cd $SLURM_SUBMIT_DIR

module load trinityrnaseq/gcc/64/2.1.1



align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --trinity_mode --transcripts transcripts.fa \
--thread_count 16 \
--seqType fq --left WBC-B2-01_AGTTCC_L002_R1_001.fastq.gz --right WBC-B2-01_AGTTCC_L002_R2_001.fastq.gz \
--output_prefix WBC-B2-01 --output_dir WBC-B2-01  &
wait

align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --trinity_mode --transcripts transcripts.fa \
--thread_count 16 \
--seqType fq --left WBC-B2-02_ATGTCA_L002_R1_001.fastq.gz --right WBC-B2-02_ATGTCA_L002_R2_001.fastq.gz \
--output_prefix WBC-B2-02 --output_dir WBC-B2-02  &
wait

align_and_estimate_abundance.pl --est_method RSEM --aln_method bowtie --trinity_mode --transcripts transcripts.fa \
--thread_count 16 \
--seqType fq --left WBC-B2-03_CCGTCC_L002_R1_001.fastq.gz --right WBC-B2-03_CCGTCC_L002_R2_001.fastq.gz \
--output_prefix WBC-B2-03 --output_dir WBC-B2-03  &
wait
```






Once the align_and_estimate_abundance jobs are done, each folder has *.isoforms.results file. Concatenate these files to generate counts matrix.
<sxh bash>
module load trinityrnaseq/gcc/64/2.1.1
/software/apps/trinityrnaseq/gcc/64/2.1.1/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix all */*isoform*results

## Above command generates all.counts.matrix and some other files
</sxh>

Look at the ERROR and OUTPUT files from RSEM to see whether the reads aligned well to the transcriptome. Pay attention to what percentage of reads in each sample aligned. For samples that have very low mapping rate, exclude them from downstream analyses.