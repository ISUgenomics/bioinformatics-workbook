---
title: "Introduction to GNU parallel"
layout: single
author: Siva Chudalayandi
author1: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


Some analyses take a long time because it is running on a single processor and there is a lot of data that needs processing.  A problem is consider **trivially parallelizable** if the data can be chunked into pieces and processed separately.  

### Examples of data that can be **trivially parallelized** include:

* When each line of a file can be processed individually
* Each chromosome of a genome can be processed individually
* Each scaffold of an assembly can be processed individually

### Examples of problems that are **trivially parallelizable**  

* Zipping or unzipping 10s to 100s of files
* Counting the number of lines in a large file
* Aligning raw sequencing data files of many samples to a genome


### Examples of problems that are not trivially parallelizable

* Genome assembly is **not** trivially parallelizable because the first step requires alignment of each read to each other read in order to find which ones are similar and should be joined (assembled).  Take a subset of the reads would result in a bunch of small poor assemblies.


## GNU parallel can be used on trivially parallelizable problems

The program that we use to parallelize a bioinformatics problem is **GNU parallel**. It is "a shell tool for executing jobs in parallel using one or more compute nodes". GNU parallel helps you run jobs that you would have otherwise run sequentially one by one or in a loop. You can check the [GNU parallel website](https://www.gnu.org/software/parallel/) to determine how to install parallel on your cluster and/or learn how to use it. On ceres, we have parallel version 20181222.

Here, we load the module and look at the version

```
module load parallel
parallel --version
```

```
GNU parallel 20181222
Copyright (C) 2007-2018 Ole Tange and Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
This is free software: you are free to change and redistribute it.
GNU parallel comes with no warranty.

Web site: http://www.gnu.org/software/parallel
```


## Example dataset

We will be using COVID-19 data collated by [New York Times github repository](https://github.com/nytimes/covid-19-data)

```
mkdir GNU-parallel
cd GNU-parallel
wget https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv
```

This is a comma separated file so let's convert that to a tab delimitated file

```
more us-counties.csv  | tr ',' '\t' > us-counties.tab
```

As you can see, this data contains the county and state information about the pandemic over time.

```
head us-counties.tab

date    county  state   fips    cases   deaths
2020-01-21      Snohomish       Washington      53061   1       0
2020-01-22      Snohomish       Washington      53061   1       0
2020-01-23      Snohomish       Washington      53061   1       0
2020-01-24      Cook    Illinois        17031   1       0
2020-01-24      Snohomish       Washington      53061   1       0
2020-01-25      Orange  California      06059   1       0
2020-01-25      Cook    Illinois        17031   1       0
2020-01-25      Snohomish       Washington      53061   1       0
2020-01-26      Maricopa        Arizona 04013   1       0
```

Instead of one large file let's separate this data by county-state

Using sort and awk we can first sort the file by county/state and then using awk to print each line ($0) to a file named county-state.tab.

```
sort -k 2,3 us-counties.tab | awk '{print $0 > $2"-"$3".tab"}'
```

This will generate 2578 files + the original 2 files we downloaded

```
ls | wc
   2580    2580   50550
```

## GNU examples

### Example 1) Gzipping 2580 text files

#### Learning objectives

  * How to use parallel with gzip to zip many files simultaneously
  * The anatomy of a parallel command

Let's make a copy of the data and compare how long it takes to run gzip using a for loop vs using parallel

```
mkdir -p gzip/parallel
mkdir -p gzip/forloop
cp *.tab gzip/parallel
cp *.tab gzip/forloop
```

#### Gzip using a for loop
We can do this using a for loop as follows.

```
time for f in *.tab; do gzip $f; done
real    0m15.801s
user    0m1.414s
sys     0m5.045s
```

#### Gzip using parallel

However, we can make better use of all the available CPUs by using GNU parallel. The anatomy of the function is

* `parallel` command
* `-j10` number of jobs or cpus to use for processing
* `"command"` in this case `gzip {}` where `{}` is a place holder for substituting a list of files defined after the delimiter
* ':::' the delimiter
* `*.tab` the list of files using the `*` operator for any file that ends in tab

```
time parallel -j10 "gzip {}" ::: *.tab

real    0m6.405s
user    0m6.872s
sys     0m11.623s
```

As you can see this sped up the gziping command by a factor of 2.3.

Here is how you can unzip all the files in parallel It takes about the same amount of time.
```
parallel -j10 "gunzip {}" ::: *.tab.gz

real    0m5.519s
user    0m0.376s
sys     0m1.367s

```

### Example 2) listing all the files using `ls`

#### Learning objectives
  * parallel has a time overhead when it starts and as it runs

Interestingly, this is about the same amount of time it takes to run a simple ls command on each file using parallel.  `ls` is a simple command and shouldn't take very long to run in parallel.

```
time parallel -j10 "ls {}" ::: *.tab > /dev/null

real    0m5.955s
user    0m6.713s
sys     0m10.994s
```

Above, I am redirecting the output to the null device so It doesn't print it to standard output

However, if we add the `-X` parameter it finishes in a fraction of the time.

```
time parallel -X -j10 "ls {}" ::: *.gz > /dev/null

real    0m0.333s
user    0m0.210s
sys     0m0.214s
```

This is because there is some overhead from the parallel command starting a new shell and opening files for buffering.  The `-X` will tell parallel not to do this as it is all on the same node.  

This is actually slower than just typing `ls *.gz` So parallel is great when you have lots of data and lots of files in the data or if a task is taking a really long time using a single processor but overhead of the command will make it slower if the dataset is too small.

```
time ls *.gz > /dev/null

real    0m0.277s
user    0m0.017s
sys     0m0.079s

```

Unfortunately, Repeating the gzipping example from above using this parameter does not give as significant of an improvement as the `ls` example.

```
time parallel -X -j10 "gzip {}" ::: *.tab

real    0m5.438s
user    0m0.492s
sys     0m1.398s
```



### Example 3) Using parallel for programs that require more than one input file.

#### Learning objectives

  * How to use parallel with programs that take more than one input
  * `:::+`

A good example of this in bioinformatics is aligning paired-end reads to a genome. Every sample has a forward (Left) and reverse (right) read pair.  For this example, we will be aligning several small sample files to a genome using bowtie2.

Let's download a toy example from Arabidopsis.  I grabbed the first 250 sequences from four Arabidposis samples taken from NCBI's SRA database.

* need to find a place that people can easily download this from to their local or remote machines.  Github has some issues with this. Perhaps a gitrepo of just super simple raw data set examples then I could just git clone the entire repo and hopefully it will remain small


```
mkdir rawdata
cd rawdata/
git init
git remote add origin -f https://github.com/ISUgenomics/GNU_parallel.git
git config core.sparseCheckout true
echo fastqfiles >> .git/info/sparse-checkout
git pull origin master

```


```
https://raw.githubusercontent.com/ISUgenomics/GNU_parallel/master/fastqfiles/SRR4420293_1.fastq?token=ABQPC5HDNOZWZFTT72GKZZC6TIOLM

```

```
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
module load bowtie2
bowtie2-build  TAIR10_chr_all.fas tair
parallel -j2 "bowtie2 --threads 4 -x tair -k1 -q -1 {1} -2 {2} -S {1/.}.sam >& {1/.}.log" ::: *_1.fastq :::+ *_2.fastq
```




The reads are in a directory named `reads_dir`


****Left Reads:****
```
first_R1.fq.gz
second_R1.fq.gz
third_R1.fq.gz
fourth_R1.fq.gz
fifth_R1.fq.gz
```
****Right reads****
```
first_R2.fq.gz
second_R2.fq.gz
third_R2.fq.gz
fourth_R2.fq.gz
fifth_R2.fq.gz

```
the genome index is located in a directory called `bwt_index`

We decide to write the output SAM and log files to <out_dir>.

We can design a bash script with the bowtie2 command for each pair separately . Obviously if dealing with hundreds of files this could become cumbersome.

```
module load bowtie2
bowtie2 --threads 4 -x <bwt_index> -k1 -q -1 first_R1.fq.gz -2 first_R2.fq.gz -S <out_dir>/first_R1.sam >& <out_dir>/first.log
..........................
..........................
bowtie2 --threads 4 -x <bwt_index> -k1 -q -1 fifth_R1.fq.gz -2 fifth_R2.fq.gz -S <out_dir>/fifth_R1.sam >& <out_dir>/fifth.log

```

GNU parallel let's us automate this task by using a combination of `substitution` and `separators` notably `:::` and `:::+`. We can also make optimum use of the available threads.

```
module load bowtie2
parallel -j2 "bowtie2 --threads 4 -x <bwt_index> -k1 -q -1 {1} -2 {2} -S <out_dir>/{1/.}.sam >& <out_dir?/{1/.}.log" ::: reads_dir/*R1*gz :::+ reads_dir/*R2*gz
```







****III) Splitting a big job to make use of all the available cpus****

Lets assume we have a large file `test.fa`. Our aim is simply to count the lines in the fasta file.
```
head test.fa

>TRINITY_DN22368_c0_g1_i1 len=236 path=[427:0-235] [-1, 427, -2]
ATTGGTTTTCTACGGAGAGAGAGAAAATGGAGACGGCGAGTGTCTAAAGCTAGAGCTTGT
GTTGGAGAAGGAAACGGAGATTTGCGTAGTAGTGGAAGCTTTAGGTATTTGTTGTGGTTA
CTCACGGCGGCGATATTTGACGGCGGGAGGAGGAAAAGAGAGAGGAAAGAACAGAGGAAG
AAGATGAGAGGAAACATTGAGAGAGAGTGAGAAGGGTTTTGTGATTTTTGTGTCTG
>TRINITY_DN22345_c0_g1_i1 len=615 path=[593:0-614] [-1, 593, -2]
GCCGGATTCAGATACGCAAGGAGAATCTGAGCAGGTCGAATGTTGATGGTATGCTTTCAT
CGGCACTTCCAGGTGGTCAGGAGAAGATCCCCATACGACTGCACTCTCTTTGCTATATGA
TGAAGCAGGAACTGTCACAAGAGGCAGAGAAGTACTGGACTCTGCCATTTGCTCATTTGT
AGCATGATTTCCTTCCCCATTCTCAGTTCCGGGAGTGCAGTGAAAGCAACAATCATTATT
```

We have reserved a node and are making use of 10 cpus on Ceres and we decide to run the `wc -l` command to find the number of lines. We can also use the unix command `time` to see how much time it takes to run the command.

```
time wc -l test.fa
1082567 test.fa

real	0m1.237s
user	0m0.025s
sys	0m0.057s

```

Now using parallel we can take advantage of the 10 cpus and spread this job over all the cpus;

```
parallel -a test.fa --pipepart --block -1  time wc -l
```
****Note:
i) `-a` option means input is read from a file            
ii)`--pipepart`: pipe parts of a file****


```
111748
real	0m0.021s
user	0m0.005s
sys	0m0.002s

104450

real	0m0.023s
user	0m0.002s
sys	0m0.005s
111050

real	0m0.026s
user	0m0.002s
sys	0m0.005s
108797

real	0m0.023s
user	0m0.005s
sys	0m0.002s
108114

real	0m0.021s
user	0m0.005s
sys	0m0.002s
109001

real	0m0.021s
user	0m0.003s
sys	0m0.003s
106528

real	0m0.022s
user	0m0.002s
sys	0m0.004s
109069

real	0m0.019s
user	0m0.002s
sys	0m0.005s
104735

real	0m0.020s
user	0m0.003s
sys	0m0.003s
109075

real	0m0.021s
user	0m0.003s
sys	0m0.003s

```

Notice we use less time with each block being counted by each cpu. The longest time in this case is `0.026 seconds` compared to `1.237 seconds` when not making use of parallel.

****A relevant Bioinformatics example****             
We have seen above how we can use GNU parallel to take advantage of available computing power. We can use this to our advantage in performing blast searches of e.g, `test.fa` that we have seen above. We can separate the input fasta file into various chunks (say `size of 1024K`) based on any defined record separators/starter (`--recstart`) in this case `>`, i.e the fasta header then use `--pipe` to read from std input and spread the job to the available cpus.

```
cat test.fa | parallel --block 1024k --recstart '>' --pipe blastp -evalue 0.01 -outfmt 6 -db <nr> -query - > blast_results
```

 The alternative of this command by using the faster `--pipepart` is as follows (*note: we can use the -a argument to directly read from the file*):

 ```
 parallel -a test.fa --pipepart --block 1024k --recstart '>' blastp -evalue 0.01 -outfmt 6 -db <nr> -query - > blast_results
 ```
