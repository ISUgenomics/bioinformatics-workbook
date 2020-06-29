---
title: "Introduction to Data Acquisition"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Downloading files from SRA

Introduced in 2007, the [SRA file format](https://en.wikipedia.org/wiki/Sequence_Read_Archive) helps standardize and link sequencing data  within the NCBI databases ([Wheeler et al., 2008](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2238880/)). There are multiple methods for downloading SRA files including the [SRA toolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/), [wget/curl commands](https://www.baeldung.com/linux/curl-wget), and [Aspera Connect](https://downloads.asperasoft.com/connect2/). There is also a method for downloading all SRA files relating to a given NCBI BioProject. 

## Using SRA toolkit
The [SRA toolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/) has been configured to connect to [NCBI SRA Database](https://www.ncbi.nlm.nih.gov/sra) and download files via [FTP](https://en.wikipedia.org/wiki/File_Transfer_Protocol). Assuming you are on the HPC cluster, the basic command to fetch a SRA file will look similar to:

```
module load sratoolkit    # <= may need to search for module name
fastq-dump SRR1234567
```
This will download a SRA file (in `sra` format) and then convert it to a single `fastq` file in [fastq format](https://en.wikipedia.org/wiki/FASTQ_format).
If your SRA file is paired, you will still end up with a single `fastq` file, since `fastq-dump` writes them as interleaved file by default. To get two separate `fastq` files (left- and right-reads), provide the `--split-files` argument:

```
module load sratoolkit
fastq-dump --split-files SRR1234567
```
The downloaded fastq files will have a `sra` number suffixed on all header lines of `fastq` file. Click the `head` command below, to see the first 12 lines of an example fetched file for `SRR447882`.

<!-- todo: Is there a smaller SRR file to use in these examples?-->

<details><summary> `head -n 12 SRR447882.fastq` </summary>

```
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

</details>

Although, this normally does not affect any programs, some programs might throw an error saying that it can't process these `fastq` files. To avoid this, you an request the file to be in the orignal format (`--origfmt`). Also, note that if you're downloading files in bulk, you can save a lot of space by compressing them in gzip format (`--gzip`).

```
module load sratoolkit
fastq-dump --split-files --origfmt --gzip SRR1234567
```
The `fastq-dump` is also capable of doing:
 * *Additional filtering or clipping of the downloaded reads*: to remove reads with poor quality or to trim adapters. Although, this will work for the single end reads, for paired-end reads it may cause differential treatment for each pairs and might not be usable for mapping programs that needs strict pairs.
 * *Compressed format:* either as gzipped or bzipped files using `--gzip` or `--bzip2` options.
 * *fasta format*: by using the `--fasta` option

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
#! /usr/bin/env bash
# USAGE: bash thisscript.sh list_of_ids.txt

list_of_ids=$1

while read line; do
  wget http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=${line}&format=fastq;
done < list_of_ids
```
The datasets can also be downloaded from DDBJ or EMBL using the FTP links, but the transfer speeds might be affected if you're not near their servers.

## Using Aspera Connect (ascp)

[Aspera]() uses high-speed file transfer to rapidly transfer large files and data sets over an existing WAN infrastructure.

To get the `sra` files:

```
prefetch --max-size 100G \
  --transport ascp \
  --ascp-path "/path/to/aspera/3.6.2/bin/ascp_or_etc/asperaweb_id_dsa.openssh" \
  SRRNNNNNN
```

This usually prefetches the SRA file to your home directory in folder named ncbi. If your home directory does not contain enough space to store all data, you may want to create another directory and softlink to the home. To do this:

```
cd ~
mkdir -p /project/storage/your_dir/ncbi
ln -s /project/storage/your_dir/ncbi .
```

when you run this, you will have a directory named `ncbi` in your home, but the data is actually stored in `/project/storage/your_dir/ncbi`

```
ncbi
└── public
    └── sra
        ├── SRR006189.sra
        └── SRR006190.sra
```

Then you can convert the SRA files back to fastq format using `fastq-dump` command.

```
for sra file in ~/ncbi/public/sra/*; do
  fastq-dump --split-files --origfmt --gzip ${sra}
done
```

##Downloading all SRA files related to a BioProject/study

[NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) stores sequence and quality data (fastq files) in aligned or unaligned formats from NextGen sequencing platforms. A [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/) is a collection of biological data related to a single initiative, originating from a single organization or from a consortium. A BioProject record provides users a single place to find links to the diverse data types generated for that project. Often times, once single BioProject will hold a considerable number of experiments and it gets tedious to download them all individually. Here is the guide to show how to do this in a effecient way:

First load the modules that are needed:

<!-- scrap: module load sratoolkit/2.5.4-1 (why do we need version number?) -->
<!-- todo: do we need a link about searching for module names every time we mention loading modules? -->

```
module load edirect
module load sratoolkit
module load parallel
```

To get the SRR numbers associated with the [BioProject PRJNA301162](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA301162):

```
esearch -db sra -query PRJNA301162 |\
  efetch --format runinfo |\
  cut -d "," -f 1 > SRR.numbers
```

To download them all in parallel (limit the number to 3 concurrent downloads)

```
parallel --jobs 3 "fastq-dump --split-files --origfmt --gzip {}" ::: SRR.numbers
```

Make sure you do this on Condoddtn node or as a PBS job
