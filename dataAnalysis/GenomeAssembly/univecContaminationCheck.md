---
title: Checking your Assembly for contaminants using VecScreen
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# **Learning Objectives**

Upon completion of this section on Checking your Assembly for contaminants using VecScreen, you will understand the following:

* What the UniVec sequence database is
* How to create a test dataset with a known contaminant.
* How to check for genome contamination locally using python's jcvi module
* How to check for genome contamination locally using BLAST and the UniVec database


# **Resources**

This tutorial is based on the following webpages.  

* [About vecscreen and BLAST parameters](https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/)
* [UniVec database ftp download site](ftp://ftp.ncbi.nih.gov/pub/UniVec/)
* [Biostars: running VecScreen locally](https://www.biostars.org/p/69003/)
* [jcvi python module](https://github.com/tanghaibao/jcvi/tree/297644ecac0d660b42f450447c0b9085e64bcf5d)
  * [jcvi vecscreen](https://github.com/tanghaibao/jcvi/blob/297644ecac0d660b42f450447c0b9085e64bcf5d/jcvi/apps/vecscreen.py)

# **What is the UniVec sequence database?**

Here is a description taken directly from the README.uv file.

```
UniVec is a non-redundant database of sequences commonly attached to
cDNA or genomic DNA during the cloning process.  It was developed by
staff at the National Center for Biotechnology Information, part of
the National Library of Medicine at the National Institutes of
Health. UniVec_Core is a subset of the full UniVec database.
```


# **How to create a test dataset with a known contaminant.**

If you assembled a genome, you will have a genome to test.  If you need a genome for this exercise. You can download one.  This genome is from a fish (_Seriola dorsalis_).


```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/814/215/GCF_002814215.1_Sedor1/GCF_002814215.1_Sedor1_genomic.fna.gz
gunzip GCF_002814215.1_Sedor1_genomic.fna.gz
```

In this example we are going to artificially add the PhiX phage genome so that we can find it in later steps.  In this way you can see the difference between a strong hit and random noise from short alignments that blast will also find.  If you have a genome that you want to see if it has PhiX you should skip this step.


```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna

seqlen.awk NC_001422.fna
gi|9626372|ref|NC_001422.1| 5386

```
* [seqlen.awk script](https://github.com/ISUGIFsingularity/utilities/blob/master/utilities/seqlen.awk)

Add the PhiX scaffold to the end of the genome file.

```
cat NC_001422.fna >> GCF_002814215.1_Sedor1_genomic.fna
```

# **How to check for genome contamination locally using python's jcvi module**

## Load a version of python 2 and install jcvi module locally

The jcvi python module only works with python2. To check which python you are using using the following commands.
```
python --version
which python
```
Load module that will give you a version of python 2 then install jcvi as a local user.

```
module avail
module load python_2/2.7.14
pip install --user jcvi
```
If the install fails and you are sure you are using python2 and not python3 then you can run univec without jcvi.  See the next section.  If it all installed successfully.  Continue below.

## Run the jcvi vecscreen on the univec_core database

jcvi requires the use of blast so load blast.  The latest blast you have on your system is fine.
```

module load blast/2.7.1
python -m jcvi.apps.vecscreen mask --db=UniVec GCF_002814215.1_Sedor1_genomic.fna &
```

Skip down to the BLAST Output section below for an example of what the output will look like.

# **How to check for genome contamination locally using BLAST and the UniVec database**

## Download UniVec sequences to blast genome to.


The file is a fasta file.
```
wget ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core
```

## Make a Blast database for UniVec

```
module load blast
makeblastdb -in UniVec_Core -input_type fasta -dbtype nucl -parse_seqids -out UniVecDB

```

## Perform blast to find potential contamination in the genome

```
blastn -db UniVecDB -query GCF_002814215.1_Sedor1_genomic.fna -num_threads 16 -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue .01 -searchsp 1750000000000 -outfmt 6 > univec_blastOut &

```

### Header information for blast output
```
 1.	 qseqid	 query (e.g., gene) sequence id
 2.	 sseqid	 subject (e.g., reference genome) sequence id
 3.	 pident	 percentage of identical matches
 4.	 length	 alignment length
 5.	 mismatch	 number of mismatches
 6.	 gapopen	 number of gap openings
 7.	 qstart	 start of alignment in query
 8.	 qend	 end of alignment in query
 9.	 sstart	 start of alignment in subject
 10.	 send	 end of alignment in subject
 11.	 evalue	 expect value
 12.	 bitscore	 bit score
 ```

## BLAST Output

There are 67 blast hits of these the top 3 are to the PhiX construct we deliberately placed in the file as a positive control.  The next ten have evalus between 1e-3 and 1e-6.  All blast hits have a score greater than 30 and would be considered a strong match to a vector.  But this along is not sufficient to say it is contamination.  Each one will need to be looked at in the context of the sequencing technology used and the sequence data in that region should be examined more closely to rule each hit as contamination or not.

```
more univec_blastOut | sort -gk 11
gi|9626372|ref|NC_001422.1|     uv:J02482.1:1-5386-49   100.000 5386    0       0       1       5386    1       5386    0.0     10801
NW_019525453.1  uv:AY390769.1:1-1000    93.069  101     7       0       1862465 1862565 136     236     3.16e-24        118
gi|9626372|ref|NC_001422.1|     uv:J02482.1:1-5386-49   100.000 49      0       0       1       49      5387    5435    3.44e-18        98.7
NW_019524969.1  uv:DQ666273.1:28-245    100.000 29      0       0       52281   52309   117     145     4.07e-06        58.6
NW_019524890.1  uv:U89960.1:4738-7717   97.059  34      1       0       1115904 1115937 2757    2724    1.64e-05        56.6
NW_019524969.1  uv:DQ666273.1:28-245    100.000 28      0       0       52283   52310   113     140     1.64e-05        56.6
NW_019525353.1  uv:JX025643.1:1-2941    97.059  34      1       0       54161   54194   4       37      1.64e-05        56.6
NW_019525353.1  uv:JX025643.1:5331-5389-49      97.059  34      1       0       54161   54194   63      96      1.64e-05        56.6
NW_019525353.1  uv:JX025643.1:5331-5389-49      97.059  34      1       0       54161   54194   92      59      1.64e-05        56.6
NW_019525453.1  uv:AY390769.1:1-1000    97.059  34      1       0       1862398 1862431 72      105     1.64e-05        56.6
NW_019525353.1  uv:JX025643.1:1-2941    96.970  33      1       0       54161   54193   33      1       6.57e-05        54.6
NW_019524828.1  uv:JN835257.1:94-190    100.000 26      0       0       115169  115194  71      46      2.64e-04        52.6
NW_019524929.1  uv:NGB01023.1:946-1024  96.875  32      1       0       311919  311950  73      42      2.64e-04        52.6
NW_019520528.1  uv:JN835257.1:94-190    100.000 25      0       0       5604    5628    47      71      0.001   50.6
NW_019521395.1  uv:DQ666273.1:28-245    100.000 25      0       0       2662    2686    140     116     0.001   50.6
NW_019524181.1  uv:JN835257.1:94-190    100.000 25      0       0       21183   21207   73      49      0.001   50.6
NW_019524181.1  uv:NGB00118.1:1-154     100.000 25      0       0       21183   21207   25      1       0.001   50.6
NW_019524181.1  uv:NGB00119.1:1-81      100.000 25      0       0       21183   21207   25      1       0.001   50.6
NW_019524181.1  uv:NGB00120.1:1-86      100.000 25      0       0       21183   21207   25      1       0.001   50.6
NW_019524181.1  uv:NGB00122.1:1-104     100.000 25      0       0       21183   21207   25      1       0.001   50.6
NW_019524181.1  uv:NGB00399.1:2658-2791 100.000 25      0       0       21183   21207   61      85      0.001   50.6
NW_019524228.1  uv:JN835257.1:94-190    100.000 25      0       0       297938  297962  43      67      0.001   50.6
NW_019525019.1  uv:NGB00123.1:1691-1849 100.000 25      0       0       230801  230825  159     135     0.001   50.6
NW_019525137.1  uv:AY330212.1:119-333   100.000 25      0       0       1080501 1080525 177     153     0.001   50.6
NW_019525137.1  uv:U29899.1:1847-4463   96.875  32      0       1       1930702 1930733 2045    2075    0.001   50.6
NW_019525235.1  uv:DQ666273.1:28-245    100.000 25      0       0       62534   62558   140     116     0.001   50.6
NW_019525243.1  uv:JN835257.1:94-190    100.000 25      0       0       336403  336427  44      68      0.001   50.6
NW_019525301.1  uv:JN835257.1:94-190    100.000 25      0       0       5990505 5990529 46      70      0.001   50.6
NW_019525348.1  uv:JN835257.1:94-190    100.000 25      0       0       869897  869921  71      47      0.001   50.6
NW_019525429.1  uv:U29899.1:1847-4463   100.000 25      0       0       17326   17350   571     595     0.001   50.6
NW_019456168.1  uv:J02400.1:1-5243-49   100.000 24      0       0       299     322     2905    2928    0.004   48.5
NW_019521395.1  uv:DQ666273.1:28-245    100.000 24      0       0       2663    2686    145     122     0.004   48.5
NW_019521753.1  uv:U09128.1:3334-10218  100.000 24      0       0       136722  136745  2438    2415    0.004   48.5
NW_019524159.1  uv:JN835257.1:94-190    100.000 24      0       0       10688   10711   69      46      0.004   48.5
NW_019524181.1  uv:NGB00118.1:1811-1938 100.000 24      0       0       21184   21207   105     128     0.004   48.5
NW_019524181.1  uv:NGB00120.1:1895-2001 100.000 24      0       0       21184   21207   84      107     0.004   48.5
NW_019524181.1  uv:NGB00121.1:1125-1248 100.000 24      0       0       21184   21207   101     124     0.004   48.5
NW_019524181.1  uv:NGB00123.1:1691-1849 100.000 24      0       0       1161    1184    159     136     0.004   48.5
NW_019524781.1  uv:NGB00399.1:2658-2791 100.000 24      0       0       61029   61052   86      63      0.004   48.5
NW_019524828.1  uv:NGB00843.1:1-38      100.000 24      0       0       119755  119778  38      15      0.004   48.5
NW_019524838.1  uv:NGB00399.1:2658-2791 100.000 24      0       0       932657  932680  90      67      0.004   48.5
NW_019524881.1  uv:U89960.1:3-1162      96.667  30      1       0       184581  184610  1133    1104    0.004   48.5
NW_019524971.1  uv:JN835257.1:94-190    100.000 24      0       0       1722    1745    69      46      0.004   48.5
NW_019525017.1  uv:JN835257.1:94-190    100.000 24      0       0       1394355 1394378 46      69      0.004   48.5
NW_019525072.1  uv:JN835257.1:94-190    100.000 24      0       0       1183803 1183826 71      48      0.004   48.5
NW_019525102.1  uv:KC702166.1:1-523     100.000 24      0       0       184131  184154  367     344     0.004   48.5
NW_019525111.1  uv:JN835257.1:94-190    100.000 24      0       0       124086  124109  44      67      0.004   48.5
NW_019525125.1  uv:JN835257.1:94-190    100.000 24      0       0       62628   62651   46      69      0.004   48.5
NW_019525138.1  uv:EU545996.1:1461-1586 96.667  30      1       0       1352776 1352805 86      57      0.004   48.5
NW_019525199.1  uv:JN835257.1:94-190    100.000 24      0       0       377772  377795  69      46      0.004   48.5
NW_019525235.1  uv:DQ666273.1:28-245    100.000 24      0       0       62535   62558   145     122     0.004   48.5
NW_019525242.1  uv:NGB00272.1:1-57      96.667  30      1       0       51556   51585   28      57      0.004   48.5
NW_019525243.1  uv:JN835257.1:94-190    100.000 24      0       0       1717558 1717581 68      45      0.004   48.5
NW_019525390.1  uv:JN835257.1:94-190    100.000 24      0       0       104486  104509  46      69      0.004   48.5
NW_019525446.1  uv:NGB00399.1:2658-2791 100.000 24      0       0       242023  242046  90      67      0.004   48.5
NW_019525459.1  uv:JN835257.1:94-190    100.000 24      0       0       232931  232954  46      69      0.004   48.5
NW_019525479.1  uv:NGB00225.1:1-59      96.667  30      1       0       1944011 1944040 27      56      0.004   48.5
NW_019525483.1  uv:NGB00842.1:1-38      100.000 24      0       0       373303  373326  38      15      0.004   48.5
NW_019525485.1  uv:JN835257.1:94-190    100.000 24      0       0       3562625 3562648 69      46      0.004   48.5
NW_019525485.1  uv:NGB00399.1:2658-2791 100.000 24      0       0       3793263 3793286 67      90      0.004   48.5
NW_019525487.1  uv:JN835257.1:94-190    100.000 24      0       0       301263  301286  69      46      0.004   48.5
NW_019525496.1  uv:NGB00842.1:1-38      100.000 24      0       0       2454921 2454944 38      15      0.004   48.5
NW_019525520.1  uv:JN835257.1:94-190    100.000 24      0       0       826124  826147  46      69      0.004   48.5
NW_019525522.1  uv:KF680544.1:323-942   100.000 24      0       0       285040  285063  178     155     0.004   48.5
NW_019525533.1  uv:JN835257.1:94-190    100.000 24      0       0       189206  189229  46      69      0.004   48.5
NW_019525549.1  uv:JN835257.1:94-190    100.000 24      0       0       1907391 1907414 70      47      0.004   48.5
NW_019525551.1  uv:J02400.1:1-5243-49   100.000 24      0       0       415782  415805  2898    2921    0.004   48.5
```

```
VecScreen Match Categories
Vector contamination usually occurs at the beginning or end of a sequence; therefore, different criteria are applied for terminal and internal matches. VecScreen considers a match to be terminal if it starts within 25 bases of the beginning of the query sequence or stops within 25 bases of the end of the sequence. Matches that start or stop within 25 bases of another match are also treated like terminal matches. Matches are categorized according to the expected frequency of an alignment with the same score occurring between random sequences.

Strong Match to Vector
(Expect 1 random match in 1,000,000 queries of length 350 kb.)
Terminal match with Score ≥ 24.
Internal match with Score ≥ 30.
Moderate Match to Vector
(Expect 1 random match in 1,000 queries of length 350 kb.)
Terminal match with Score 19 to 23.
Internal match with Score 25 to 29.
Weak Match to Vector
(Expect 1 random match in 40 queries of length 350 kb.)
Terminal match with Score 16 to 18.
Internal match with Score 23 to 24.
Segment of Suspect Origin
Any segment of fewer than 50 bases between two vector matches or between a match and an end.
```

# Further Reading

* [Large-scale contamination of microbial isolate genomes by Illumina PhiX control](https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/1944-3277-10-18)
* [About vecscreen and BLAST parameters](https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/)
* [UniVec database ftp download site](ftp://ftp.ncbi.nih.gov/pub/UniVec/)
* [Biostars: running VecScreen locally](https://www.biostars.org/p/69003/)
* [jcvi python module](https://github.com/tanghaibao/jcvi/tree/297644ecac0d660b42f450447c0b9085e64bcf5d)
  * [jcvi vecscreen](https://github.com/tanghaibao/jcvi/blob/297644ecac0d660b42f450447c0b9085e64bcf5d/jcvi/apps/vecscreen.py)
