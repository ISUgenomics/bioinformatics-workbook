# Downloading files from SRA

## Using SRAtoolkit
SRA toolkit has been configured to connect to NCBI SRA and download via FTP. The simple command to fetch a SRA file you can use this command:

```
module load sratoolkit
fastq-dump SRR1234567
```
This will download the SRA file (in `sra` format) and then convert them to `fastq` file for you. 
If your SRA file is paired, you will still end up with a single `fastq` file, since, `fastq-dump`, by default writes them as interleaved file. To change this, you can provide `--split-files` argument.
```
module load sratoolkit
fastq-dump --split-files SRR1234567
```
The downloaded fastq files will have `sra` number suffixed on all header lines of `fastq` file
```
head -n 12 SRR447882_1.fastq
@SRR447882.1.1 HWI-EAS313_0001:7:1:6:844 length=84
ATTGATCATCGACCAGAGNCTCATACACCTCACCCCACATATGTTTCCTTGCCATAGATCACATTCTTGNNNNNNNGGTGGANA
+SRR447882.1.1 HWI-EAS313_0001:7:1:6:844 length=84
BBBBBB;BB?;>7;?<?B#AA3@CBAA?@BAA@)=6ABBBBB?ACA;0A=257?A7+;;&########################
@SRR447882.2.1 HWI-EAS313_0001:7:1:6:730 length=84
AGTTGATTGTGATATAGGNGTCTATCGACATTGATGCATAGGTCCTCTATTAAACTTGTTTTGTGATGTNNNNNNNTTTTTTNA
+SRR447882.2.1 HWI-EAS313_0001:7:1:6:730 length=84
A?@B:@CA:=?BCBC:2C#7>BACB??@4@B@<=>;'>@>3:86>=6@=B@B<;)@@###########################
@SRR447882.3.1 HWI-EAS313_0001:7:1:6:1343 length=84
CATCAATGCAAGGATTGTNCCATTGGTAACAATTCCACTCCTAACTTGTCAATTGATTTTCATATAACTNNNNNNNCCAAAANT
+SRR447882.3.1 HWI-EAS313_0001:7:1:6:1343 length=84
BCB@BBC+5BCA>BABBA#@4BCCA>?CBBB4CB(*ABB?ABBAACCB8ABBB?(<<B?:########################
```
Although, this normally does not affect any programs, some programs might throw an error saying that it can't process these `fastq` files. To avoid this, you an request the file to be in the orignal format (`--origfmt`). Also, note that if you're downloading files in bulk, you can save a lot of space by compressing them in gzip format (`--gzip`).

```
module load sratoolkit
fastq-dump --split-files --origfmt --gzip SRR1234567
```
The `fastq-dump` is also capable of doing:
 * *Additional filtering or clipping of the downloaded reads*: to remove reads with poor quality or to trim adapters. Although, this will work for the single end reads, for paired-end reads it may cause differential treatment for each pairs and might not be usable for mapping programs that needs strict pairs. 
 * *Compressed format:* either as gzipped or bzipped files using `--gzip` or `--bzip2` options.
  
## Using Linux commands:

In cases were you cannot run the SRA toolkit or any other programs to download the file, you can still use the inbuilt commands of Linux such as `wget` and `curl`. The standard web link for downloading the SRA files is:

```
http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRRNNNNNN&format=fastq
```
You need to replace the `SRRNNNNNN` with the actual SRR number for it to work. 

You can either use `wget`
```
wget http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRRNNNNNN&format=fastq
```
or `curl`
```
curl -O http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRRNNNNNN&format=fastq
```
If you have a large list of ids, you can simply loop it over using a `while` loop

```
while read line; do
wget http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=${line}&format=fastq;
done<list_of_ids
```

## Using 




##Downloading all SRA files related to a BioProject/study

NCBI Sequence Read Archive (SRA) stores sequence and quality data (fastq files) in aligned or unaligned formats from NextGen sequencing platforms. A BioProject is a collection of biological data related to a single initiative, originating from a single organization or from a consortium. A BioProject record provides users a single place to find links to the diverse data types generated for that project. Often times, once single BioProject will hold a considerable number of experiments and it gets tedious to download them all individually. Here is the guide to show how to do this in a effecient way:

First load the modules that are needed:
```
module load edirect
module load sratoolkit/2.5.4-1
module load parallel
```

To get the SRR numbers associated with the project:
```
esearch -db sra -query PRJNA301162 | efetch --format runinfo |cut -d "," -f 1 > SRR.numbers
```

To download them all in parallel (limit the number to 3 concurrent downloads)
```
parallel --jobs 3 "fastq-dump --split-files --origfmt --gzip {}" ::: SRR.numbers
```
Make sure you do this on Condoddtn node or as a PBS job



