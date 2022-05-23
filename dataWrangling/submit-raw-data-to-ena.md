---
title: "Submitting raw datasets to ENA"
layout: single
author: Margaret Woodhouse
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Submitting raw datasets to ENA

This is part 2 of the ENA data submission. Here we demonstrate how to submit raw reads and associated data and connect it to your existing BioProject, BioExperiments and SRA datasets. This assumes that you have already finished the part 1 (BioProject creation) step.

An overview of the reads submission process to ENA can be found here:

https://ena-docs.readthedocs.io/en/latest/submit/reads.html

This section will address two different types of submission:

1. Raw read and other data submission for genomes;
2. Raw read and other data submission for functional data such as RNA-seq, bisulfite sequencing, and other data types not related to genome sequencing or assembly.

Before beginning, it should be noted that ENA frequently changes its submission protocols,  so that while we will make every attempt to keep this information up to date, there may be times when some of the steps have become modified or obsolete. Be sure to check the link above to see if protocols have changed from what is described here.

## Uploading Files Using Command Line FTP Client

This section explains how to upload files to ENA using a command line FTP client in Linux or Mac. The built-in FTP tool for Windows command line does not support FTPS so Windows users are recommended to use an alternative:

    Using Webin File Uploader
    Using FileZilla On Windows

The below instructions describe how you may upload your files to us through a command line FTP client in Linux or Mac.

    Open a terminal and type

```
lftp webin2.ebi.ac.uk -u Webin-xxxxx
```

filling in your Webin username that had been created when you created your Bioproject in Section 1.

    -Enter your password when prompted.
    -Type ls command to check the content of your drop box.
    -Use mput <filename> command to upload files.
    -Use bye command to exit the ftp client.

Note that in your Webin username, the ‘W’ should be upper case. In this section, we will use the fictional Webin-00000 account as our example, and a fictional password PaSSW0rd.

```
lftp webin2.ebi.ac.uk -u Webin-00000
new password: PaSSW0rd #Enter your password here
```

## Using Aspera ascp Command Line Program

Aspera is a commercial file transfer protocol that may provide better transfer speeds than FTP.

Please select the correct operating system. The ascp command line client is distributed as part of the Aspera Cli in the cli/bin folder.

Your command should look similar to this:

```
ascp -QT -l300M -L- <file(s)> <Webin-N>@webin.ebi.ac.uk:.
```

The `-l300M` option sets the upload speed limit to 300MB/s. You may wish to lower this value to increase the reliability of the transfer.

The `-L-` option is for printing logs out while transferring,

The `<file(s)>` can be a file mask (e.g. `*.cram`), a list of files or a single file.

`<Webin-N>` is your Webin submission account name.

```
ascp -QT -l300M -L- <file(s)> Webin-00000@webin.ebi.ac.uk:.
```

## 1. Submitting raw reads data for genomes


In this example, we will use the submission of PacBio subread bam files to demonstrate the raw read genomic submission process. We are submitting information related to the maize cultivar B104. We are using Aspera ascp (described above) as the mode of submission via the command line.

To submit all PacBio bam files associated with the B104 genome submission, we used the following command:

```
ascp -QT -l300M -L- B104.CELL01.bam B104.CELL03.bam B104.CELL05.bam B104.CELL07.bam B104.CELL09.bam B104.CELL11.bam B104.CELL13.bam B104.CELL02.bam B104.CELL04.bam B104.CELL06.bam B104.CELL08.bam B104.CELL10.bam B104.CELL12.bam B104.CELL14.bam Webin-00000@webin.ebi.ac.uk:.
```

Make an md5 of all files:

```
md5sum *.bam > B104_pacbio_subreads_md5
```

You will need this md5 data for the next step.

Next, you will return to the ENA website and go to the Read Submission page:

https://www.ebi.ac.uk/ena/submit/webin/read-submission

There, you will select your file type (in this case, "BAM" file format), and then download the pre-formatted spreadsheet for BAM file submission on the next page; you will have the option of de-selecting certain metadata types if you wish. Fill out the spreadsheet as indicated. Column 1 under "sample" will be the sample name of the genome that you had registered in Section 1, followed by the study ID, then instrument model and other information if you chose to include those in the spreadsheet. Lastly, you will include the filename of each bam file which you submitted via Aspera, followed by its md5 which you generated above. You will then upload your completed spreadsheet to the ENA reads submission portal and check to make sure all the data is correct.


For other genomic data types, such as 10Xgenomics reads data, a similar submission process is followed:

```
ascp -QT -l300M -L- B104_S100_L001_R1_001.fastq.gz B104_S201_L001_R1_001.fastq.gz B104_S301_L001_R1_001.fastq.gz B104_S401_L001_R1_001.fastq.gz B104_S100_L001_R2_001.fastq.gz B104_S201_L001_R2_001.fastq.gz Webin-00000@webin.ebi.ac.uk:.
```
(Note that the list of files is shortened for brevity.)

Once again, an md5 of all the 10xgenomics files were made:

```
md5sum *.fastq.gz > B104_10x_fastq_md5
```

However, this time, in the ENA reads submission portal https://www.ebi.ac.uk/ena/submit/webin/read-submission, you will select "Fastq" as the file type instead of bam files.

To complete the submission of your genome(s), please see Section 3, "Submission Guidelines for Chromosomal level Assemblies to ENA":

https://github.com/ISUgenomics/bioinformatics-workbook/blob/master/dataWrangling/ena-genome-submission.md
