## Genome guided assembly of transcripts

This section covers different methods for genome guided transcript assembly. You will need a `bam` file after aligning RNAseq reads to the genome using any splice aware alignment program (HiSat2, STAR, GSNAP etc) and should be coordinate sorted.

Following four methods will be covered in this section:

1. [StringTie](https://www.nature.com/articles/nbt.3122)
2. [Class2](https://www.ncbi.nlm.nih.gov/pubmed/26975657)
3. [Trinity](https://www.ncbi.nlm.nih.gov/pubmed/21572440)
4. [Cufflinks](https://www.ncbi.nlm.nih.gov/pubmed/22383036)


### StringTie

The input file is the BAM format file that is coordinate sorted (`alignments_sorted.bam`). It also requires that there is `xs` tag included for the splice aligned reads.

You can run the transcript assembly as follows:

```
stringtie -v -o alignments_stringtie.gtf -p 16 alignments_sorted.bam
```
There are many more options that can be fine tuned for your dataset. Type `stringtie -h` for the help page.

### Class2

Again the input is BAM file. The command is as follows:

```
run_class.pl -a alignments_sorted.bam -o alignments_class2.gtf -p 16 --verbose
```

### Trinity

For Trinity, you can run the following command:

```
Trinity --genome_guided_bam alignments_sorted.bam --max_memory 100G --genome_guided_max_intron 10000 --CPU 16
```

### Cufflinks

```
cufflinks -o alignments_cufflinks.gtf -p 16 --multi-read-correct --frag-bias-correct genome.fasta alignments_sorted.bam
```


### FAQs:

**1. How do I convert my sam file to bam format?**

There are many programs that lets you convert sam to bam. One option is,  use `samtools` for doing the conversion as follows:
```
samtools view --threads 16 -b -o output.bam input.sam
```

**2. How to sort my BAM file?**

Again, there are many alternatives, but using `samtools` you can sort them as follows:*

```
samtools sort--threads 16 -o output_sorted.bam  input.bam
```

**3. How to include XS tag for spliced alignments?**

For this, you need to make sure that you run the alginment program with proper options. If you are using HiSat2, ensure that you use `--dta` option while aligning. If you are using STAR, then add ` --outSAMattributes XS`. with GSNAP, you `XS` tag is always included for splice alignments.

**4. My alignment file doesn't have XS tag for spliced alignments, how can I fix this?**

Fortunately, `class2` program comes with an accessory executable that lets you add the `XS` flag post alignment! You can run this as follows:

```
samtools view -h input.bam | addXS ref_genome.fasta | samtools view -bS - > output_with_xs_tag.bam
```


**5. How much time does the program take for running assembly?**

Each program differs in time and the data size greatly affect the performance of each program. If your datasize is small, `cufflinks` might be fastest. If you are running on HPC, request plenty of walltime before running the jobs.

**6. I have several RNAseq libraries, should I run the assembly separately?**

Normally, you should merge the files (all R1 and all R2 seprately before aligning OR merge bam files after alignment step), before running the transcript assembly. This will provide better coverage for the transcript and chances of assembling the complete transcript is higher. If you want the results quickly, then you can run them separately and then merge the gtf file/transcript file obtained. However, you will have to do additional step of removing duplicates in your final file.
