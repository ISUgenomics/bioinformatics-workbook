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

This is part 2 of the ENA data submission. Here we demonstrate how to submit raw reads and associated data and connect it to your existing BioProject, BioExperiments and SRA datasets. This assumes that you have already finished Part 1, the BioProject creation step.

An overview of the reads submission process to ENA [can be found here](https://ena-docs.readthedocs.io/en/latest/submit/reads.html).

This section will address two different types of raw reads submission:

1. Raw read and other data submission for *genomes*;
2. Raw read and other data submission for *functional data* such as RNA-seq, bisulfite sequencing, and other data types not related to genome sequencing or assembly.

Before beginning, it should be noted that ENA frequently changes its submission protocols,  so that while we will make every attempt to keep this information up to date, there may be times when some of the steps have become modified or obsolete. Be sure to check the link above to see if protocols have changed from what is described here.

## Uploading Files Using Command Line FTP Client

This section explains how to upload files to ENA using a command line FTP client in Linux or Mac.

Open a terminal and type

```
lftp webin2.ebi.ac.uk -u Webin-xxxxx
```

filling in your Webin username that had been created when you created your Bioproject in Section 1.

Enter your password when prompted. Type ls command to check the content of your drop box. Use `mput <filename>` command to upload files. Use `bye` command to exit the ftp client.

Note that in your Webin username, the ‘W’ should be upper case. In this section, we will use the fictional `Webin-00000` account as our example, and a fictional password `PaSSW0rd`.

```
lftp webin2.ebi.ac.uk -u Webin-00000
new password: PaSSW0rd #Enter your password here
```

## Using Aspera ascp Command Line Program

Aspera is a commercial file transfer protocol that may provide better transfer speeds than FTP. When uploading many large files, such as bam, fastq, or other types of read files, Aspera is much more efficient. Aspera can be downloaded here:

http://downloads.asperasoft.com/en/downloads/62

The basic Aspera command is structured like so:

```
ascp -QT -l300M -L- <file(s)> <Webin-N>@webin.ebi.ac.uk:.
```

The `-l300M` option sets the upload speed limit to 300MB/s. You may wish to lower this value to increase the reliability of the transfer.

The `-L-` option is for printing logs out while transferring.

The `<file(s)>` can be a file mask (e.g. *.cram), a list of files, or a single file.

`<Webin-N>` is the Webin submission account name that was generated when you created your account.


## 1. Submitting raw reads data for genomes


In this example, we will use the submission of PacBio subread bam files to demonstrate the raw read genomic submission process. We are submitting reads related to the maize cultivar genome B104. We are using Aspera ascp (described above) as the mode of submission via the command line.

To submit all PacBio bam files associated with the B104 genome submission, we used the following command:

```
ascp -QT -l300M -L- B104.CELL01.bam B104.CELL02.bam B104.CELL03.bam Webin-00000@webin.ebi.ac.uk:.
```
(Note that the list of files is shortened for brevity.)

Next, make an md5 of all files. You will need this md5 data for the next step.

```
md5sum *.bam > B104_pacbio_subreads_md5
```


Next, you will return to the ENA website and go to the [ENA reads submission portal](https://www.ebi.ac.uk/ena/submit/webin/read-submission) (you will need to be logged into your Webin account to see the reads portal page). There, you will select your file type (in this case, "BAM" file format), and then download the pre-formatted spreadsheet for BAM file submission on the next page; you will have the option of de-selecting certain metadata types if you wish. Fill out the downloaded spreadsheet as indicated. Column 1 under "sample" will be the sample name of the genome that you had registered in Section 1, followed by the study ID or BioProject ID, then instrument model and other information if you chose to include those in the spreadsheet. Lastly, you will include the filename of each bam file which you submitted via Aspera, followed by its md5 which you generated above. You will then upload your completed spreadsheet to the ENA reads submission portal and check to make sure all the data is correct.


For other genomic data types, such as 10Xgenomics reads data, a similar submission process is followed:

```
ascp -QT -l300M -L- B104_S100_L001_R1_001.fastq.gz B104_S201_L001_R1_001.fastq.gz B104_S301_L001_R1_001.fastq.gz B104_S401_L001_R1_001.fastq.gz B104_S100_L001_R2_001.fastq.gz B104_S201_L001_R2_001.fastq.gz Webin-00000@webin.ebi.ac.uk:.
```
(Note that the list of files is shortened for brevity.)

Once again, an md5 of all the 10xgenomics files were made:

```
md5sum *.fastq.gz > B104_10x_fastq_md5
```

However, this time, in the [ENA reads submission portal](https://www.ebi.ac.uk/ena/submit/webin/read-submission), you will select "Fastq" as the file type instead of bam files to generate the metadata download spreadsheet.

To complete the submission of your genome(s), please see Section 3, [Submission Guidelines for Chromosomal level Assemblies to ENA](https://bioinformaticsworkbook.org/dataWrangling/ena-genome-submission.html#gsc.tab=0).



## 2. Submitting raw reads for functional (non-genomic) data to EBI:

For data such as RNA-seq, bisulfite, and other data types, you will first need to create an [Annotare account](https://www.ebi.ac.uk/fg/annotare/).

An overview of Annotare file transfer can be found [here](https://www.ebi.ac.uk/fg/annotare/help/ftp_upload.html#TerminalFTP).

To transfer files using Aspera command line client:

```
ascp -P 33001 <yourfile.txt> aexpress@fasp.ebi.ac.uk:<your upload directory>/

aexpress = username
```

When prompted for a password enter `aexpress1`.

Proceed to import the uploaded files to your submission, as described [here](https://www.ebi.ac.uk/fg/annotare/help/sample_attributes.html#annotate) and [here](https://www.ebi.ac.uk/arrayexpress/help/UHTS_submissions.html#HowToSubmit).

Please note that there is a private directory with a unique name for each submission where files must be copied to. Before starting to transfer files, click on `FTP/Aspera Upload` on your Annotare page. The dialogue will show you the directory name for your submission (e.g. `ftp-private-2.ebi.ac.uk/abcfefdh-o12e456789sup/`). Copy the data files to this directory following the FTP/Aspera transfer instructions below. Please transfer the compressed files one by one (and not bundling multiple fastq.gz files in one tar.gz archive) to avoid time-out issues and to allow ENA to process your files promptly.

To associate the transferred files with your experiment submission and ensure the file integrity, follow the on-screen instructions of the dialogue, and fill in file names and their corresponding md5 checksums (remember to use the checksum of the actual compressed file that is sent to ENA). Annotare will then verify the presence of the files on their FTP site and the md5 checksums. If verification passes, you will be able to assign data files to each of your samples.

To begin uploading files via Aspera on the command line:


1. Compress fastq files;
2. get private directory name from Annotare "FTP/Aspera download" button:

```
ftp://aexpress:aexpress1@ftp-private-2.ebi.ac.uk/abcfefdh-o12e456789sup/
```

Next, get checksums of all files, which you enter into the popup that contains the directory name from the Annotare "FTP/Aspera download" button:

```
md5 *.fastq.gz > fastq_checksums
```

Since the password "aexpress1" is required for all submissions, and you can only submit one file at a time, you can set up a script to recursively submit your reads.

First, set your computer to automatically use the aexpress1 password when submitting all files. This can be done by typing in the command line

```
export ASPERA_SCP_PASS="aexpress1"
```

then make a bash script with a submission command for each file with the following structure:

```
ascp -P 33001 file.fastq.gz aexpress@fasp.ebi.ac.uk:abcfefdh-o12e456789sup/
```

You can do this by typing:

```
ls *fastq.gz | awk '{print "ascp -P 33001 "$1" aexpress@fasp.ebi.ac.uk:abcfefdh-o12e456789sup/"}' > arrayexpress.sh
```

Then run the bash script in your command line:

```
sh arrayexpress.sh
```

and this will recursively submit your reads.


Once the files are uploaded, return to Annotare and submit the md5 information for your reads in the pop-up; the reads will then be processed. Enter all metadata in the tabs on the right. When that is finished, click on the "Validate" button. If validated and you are happy with it, then hit the button at the top "Submit to ArrayExpress".
