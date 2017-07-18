#Transcriptome Assembly
In this tutorial, we will learn about how to assemble the transcriptome from RNA-Seq reads..

**Software required:** Trinity https://github.com/trinityrnaseq/trinityrnaseq/wiki

**Files required:** Illumina RNA sequencing reads either in Fastq or Fasta format

**Step1:** In case of SE reads, concatenate all the fastq files (no need to unzip) into a single file. This can be achieved by


```
cat file1.fq file2.fq ... > R1.fq
```

In case of PE reads, concatenate R1 and R2 files separately.
In case of mixed SE and PE data, concatenate all single end data along to PE R1 file.

**Step 2:** Decide (1) how many CPUs and (2) how much memory will be used to assemble the transcriptome. Let's say you have decided on 16 CPUS and 32G memory. The command to assemble is as follows:

```
## Single-end Illumina
Trinity --seqType fq --max_memory 32G --CPU 16 --trimmomatic --normalize_by_read_set --output TrinityOut --single R1.gz &

##Paired-end Illumina or a combination of SE and PE
Trinity --seqType fq --max_memory 32G --CPU 16 --trimmomatic --normalize_by_read_set --output TrinityOut --left all_left_R1.gz --right all_right_R2.gz &

## To show all available options for running Trinity
Trinity --show_full_usage_info

```

Here is an explanation of the parameters used.

--seqType - Specify fq for fastq files or fa for fasta files
--max_memory - specify available memory
--trimmomatic - include this option to trim reads
--normalize_by_read_set - include this option to perform normalization
--output - output_folder_name
--left - all R1 fastq/fasta files concatenated
--right - all R2 fastq/fasta files concatenated

and calculate differential expression.

 When there is no reference genome or when you are interested in isoform-level differential expression, you may want to assemble the transcriptome before proceeding with DE analysis.

 We will use Trinity software package for transcriptome assembly. XX version of trinity includes reads trimming and normalization options. SO, there is no need to pre-trim or normalize the reads before starting Trinity.

Command to start Trinity


```
## Single-end Illumina
Trinity --seqType fq --max_memory 32G --CPU 16 --trimmomatic --normalize_by_read_set --output TrinityOut --single R1.gz &

##Paired-end Illumina or a combination of SE and PE
Trinity --seqType fq --max_memory 32G --CPU 16 --trimmomatic --normalize_by_read_set --output TrinityOut --left all_left_R1.gz --right all_right_R2.gz &

## To show all available options for running Trinity
Trinity --show_full_usage_info
```


--seqType Specify fq for fastq files or fa for fasta files
--max_memory specify available memory
--trimmomatic include this option to trim reads
--normalize_by_read_set include this option to perform normalization
--output output_folder_name
--left all R1 fastq/fasta files concatenated
--right all R2 fastq/fasta files concatenated

## FAQ
**
What input is required for transcriptome assembly?
**Illumina sequencing reads. The reads can be single-end, or paired-end or a combination of both.
**Where is the output of the final assembly?
**Trinity.fasta in the TrinityOut folder is a fasta file of transcripts.
**What if I have a combination of single-end and paired-end sequencing data?
**No problem, just concatenate all single-end data to R1 of paired-end data files.
cat *R1_001.fastq.gz > all_left_R1.gz &
cat *R2_001.fastq.gz > all_right_R2.gz &
**Where can I get more information on Trinity?**
[https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)