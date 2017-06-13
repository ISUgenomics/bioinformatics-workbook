#Downloading all SRA files related to a BioProject/study

NCBI Sequence Read Archive (SRA) stores sequence and quality data (fastq files) in aligned or unaligned formats from NextGen sequencing platforms. A BioProject is a collection of biological data related to a single initiative, originating from a single organization or from a consortium. A BioProject record provides users a single place to find links to the diverse data types generated for that project. Often times, once single BioProject will hold a considerable number of experiments and it gets tedious to download them all individually. Here is the guide to show how to do this in a effecient way:

First load the modules that are needed:
<pre>
module load edirect
module load sratoolkit/2.5.4-1
module load parallel
</pre>

To get the SRR numbers associated with the project:
<pre>
esearch -db sra -query PRJNA301162 | efetch --format runinfo |cut -d "," -f 1 > SRR.numbers
</pre>

To download them all in parallel (limit the number to 3 concurrent downloads)
<pre>
parallel --jobs 3 "fastq-dump --split-files --origfmt --gzip {}" ::: SRR.numbers
</pre>
Make sure you do this on Condoddtn node or as a PBS job



