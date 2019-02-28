---
title: "Introduction to Data Acquisition"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

SRA toolkit has been configured to connect to NCBI SRA and download via FTP. The simple command to fetch a SRA file and to split it to forward/reverse reads, use this command:

<pre>
module load sratoolkit
fastq-dump --split-files --origfmt --gzip SRR1234567
</pre>
You will see 2 files <blockcode>SRR1234567_1.fastq</blockcode> and <blockcode>SRR1234567_2.fastq</blockcode> downloaded directly from NCBI. If the file size is more than 1Gb, submit this within a PBS script.

If this doesn't work for you (or too slow, because of FTP) then, you can try <blockcode>aspera</blockcode>  which will be fast (very useful if you have large number of files to download)
<pre>
ascp -i $KEY/asperaweb_id_dsa.openssh -k 1 -QT -l 200m \
anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR390/SRR390728/SRR390728.sra \ # url for the file
./ # save location
# the above command should be in single line
</pre>
If you have large number of files to download (usually organized as <blockcode>project</blockcode>  with a SRR project number), you can save all the IDs in a file and loop through the lines.

Finally, if you want to download using <blockcode>wget</blockcode> (which will be very slow), you can use this template:
<pre>
wget http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list=SRRXXXXXX&format=fastq
</pre>

<p>If you have a previously downloaded FASTQ file, without using the <blockcode>--origfmt</blockcode>, option, then the first field (before the instrument name) is SRA id number. Although. most of the programs won't mind, some of them might throw an error example example [Khmer] (http://khmer.readthedocs.org/en/v1.0). </p>

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


<p>To remove this using the below script:</p>

<pre>
for f in SRR447882_[12]_paired.fq; do\
awk '$1 ~ /@SRR447882*/ {$1="@"}{print}' $f | sed 's/^@ /@/g' | \
sed 's/^+SRR447882.\+/+/g' > $f.cleaned; \
done
</pre>

<p>Now the sequences should look like these:</p>

    @HWI-EAS313_0001:7:1:6:844 length=84
    ATTGATCATCGACCAGAGNCTCATACACCTCACCCCACATATGTTTCCTTGCCATAGATCACATTCTTGNNNNNNNGGTGGANA
    +
    BBBBBB;BB?;>7;?<?B#AA3@CBAA?@BAA@)=6ABBBBB?ACA;0A=257?A7+;;&########################
    @HWI-EAS313_0001:7:1:6:730 length=84
    AGTTGATTGTGATATAGGNGTCTATCGACATTGATGCATAGGTCCTCTATTAAACTTGTTTTGTGATGTNNNNNNNTTTTTTNA
    +
    A?@B:@CA:=?BCBC:2C#7>BACB??@4@B@<=>;'>@>3:86>=6@=B@B<;)@@###########################
    @HWI-EAS313_0001:7:1:6:1343 length=84
    CATCAATGCAAGGATTGTNCCATTGGTAACAATTCCACTCCTAACTTGTCAATTGATTTTCATATAACTNNNNNNNCCAAAANT
    +
    BCB@BBC+5BCA>BABBA#@4BCCA>?CBBB4CB(*ABB?ABBAACCB8ABBB?(<<B?:########################


Instead of this, you can also redonwload the original SRA file using <blockcode>--origfmt</blockcode> option, if it saves time.

<h2>Download all SRR files related to a project </h2>

If you have large number of SRR files to donwload, see if they belong to a specific project. Eg., project [[http://www.ncbi.nlm.nih.gov/Traces/sra/?study=SRP011907 | SRP011907 ]] has 283 SRR files. You can use aspera to download all 283 files at once.

Here are the steps:

  - Get the FTP link: go to [[http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=download_reads | SRA download reads]] and look for the project id (eg., SRP011907, is located under <blockcode> reads --> ByStudy --> sra --> srp --> SRP011 --> SRP011907 </blockcode>. Once you reach the link, clicking the "save icon" (next to the total size) will give the ftp link.
  - Now, download it using Aspera as expalined before.

<pre>
module load aspera/3.3.3.81344
ascp -i $KEY/asperaweb_id_dsa.openssh -k 1 -QT -l 200m \
  anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra/SRP/SRP011/SRP011907 \
  ./Destination_dir
</pre>
