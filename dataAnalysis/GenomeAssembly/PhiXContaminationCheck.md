---
title: Checking your Assembly for PhiX contamination
layout: single
author: Andrew Severin
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Learning Objectives

Upon completion of this section on Checking your Assembly for PhiX contamination you will understand the following:

* Where can I find the sequence for PhiX that is a common contaminent?
* How do I check my genome for PhiX contamination using BLAST?

PhiX is a very common contaminant that can be misassembled into genomes

  - Read [Large-scale contamination of microbial isolate genomes by Illumina PhiX control](https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/1944-3277-10-18)for more details on this issue.


# Download PhiX sequence from NCBI

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna

seqlen.awk NC_001422.fna
gi|9626372|ref|NC_001422.1| 5386

```
# Download a Genome

If you assembled a genome, you will have a genome to test.  If you need a genome for this exercise. You can download one.  This genome is a fish (_Seriola dorsalis_).

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/814/215/GCF_002814215.1_Sedor1/GCF_002814215.1_Sedor1_genomic.fna.gz
gunzip GCF_002814215.1_Sedor1_genomic.fna.gz
```

# Prep example genome

In this example we are going to artificially add the PhiX phage genome so that we can find it in later steps.  In this way you can see the difference between a strong hit and random noise from short alignments that blast will also find.  If you have a genome that you want to see if it has PhiX you should skip this step.

```
cat NC_001422.fna >> GCF_002814215.1_Sedor1_genomic.fna
```

# Make a Blast database for the _Seriola_ genome

```
module load blast
makeblastdb -in GCF_002814215.1_Sedor1_genomic.fna -input_type fasta -dbtype nucl -parse_seqids -out seriola_blastDB
makeblastdb -in GCF_002814215.1_Sedor1_genomic.fna -input_type fasta -dbtype prot -parse_seqids -out seriola_blastDB

```

# Perform blast to find potential contamination in the genome

```
blastn  -db seriola_blastDB -query NC_001422.fna -outfmt 6 | sort -k 7n > phiX_2_seriola.blastnout

tblastx -db seriola_blastDB -query NC_001422.fna -outfmt 6 | sort -k 7n > phiX_2_seriola.tblastxout
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

# BLAST Output

As you can see the only match in the blastn is the phiX genome we added and it aligned perfectly with 100% match and the full 5386 bp length.  You will not likely get a perfect match that is full length.  If you do not have PhiX contamination, this file will be empty.

```
phiX_2_seriola.blastnout
::::::::::::::
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 5386    0       0       1       5386    1       5386    0.0     9947
```

Similarly, if we look at the tblastx output.  We see a lot of perfect matches 100% identity and longish lengths >50 bp to the known PhiX genome we placed in there (NC_001422.1) and some that have very poor sequence identity 32%, 34%, 50% and so on with very short length of 43 bp, 26 bp and 28 bp, respectively.  The poor sequence identity combined with the very short alignment lengths indicate that those are false positives.


```
more phiX_2_seriola.tblastxout
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 295     0       0       1       885     1       885     0.0     626
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 152     0       0       2       457     2       457     0.0     352
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 19      0       0       3       59      3       59      0.0     42.3
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 315     0       0       108     1052    108     1052    0.0     658
gi|9626372|ref|NC_001422.1|     NW_019525106.1  32.558  43      29      0       313     441     889397  889525  3.6     37.7
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 92      0       0       497     772     497     772     0.0     194
gi|9626372|ref|NC_001422.1|     NW_019453341.1  34.615  26      17      0       569     646     365     442     0.32    28.1
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 198     0       0       594     1       594     1       0.0     437
gi|9626372|ref|NC_001422.1|     NW_019525231.1  50.000  28      14      0       599     682     173316  173233  1.0     39.6
gi|9626372|ref|NC_001422.1|     NW_019432001.1  50.000  16      8       0       602     649     402     355     4.4     25.4
gi|9626372|ref|NC_001422.1|     NW_019430130.1  42.857  14      8       0       605     646     604     563     1.2     23.1
gi|9626372|ref|NC_001422.1|     NW_019431419.1  54.545  11      5       0       605     637     333     301     0.089   24.0
gi|9626372|ref|NC_001422.1|     NW_019437257.1  42.857  14      8       0       605     646     597     556     2.3     22.6
gi|9626372|ref|NC_001422.1|     NW_019437257.1  50.000  14      7       0       605     646     402     361     8.1     25.8
gi|9626372|ref|NC_001422.1|     NW_019437257.1  50.000  14      7       0       605     646     432     391     0.26    26.7
gi|9626372|ref|NC_001422.1|     NW_019437257.1  50.000  14      7       0       605     646     567     526     0.048   25.8
gi|9626372|ref|NC_001422.1|     NW_019439991.1  50.000  14      7       0       605     646     225     184     1.6     26.3
gi|9626372|ref|NC_001422.1|     NW_019467615.1  46.667  15      8       0       605     649     226     182     0.049   23.1
gi|9626372|ref|NC_001422.1|     NW_019503050.1  50.000  14      7       0       605     646     50      91      0.011   26.3
gi|9626372|ref|NC_001422.1|     NW_019503426.1  42.857  14      8       0       605     646     265     224     4.7     21.7
gi|9626372|ref|NC_001422.1|     NW_019478677.1  45.000  20      11      0       623     682     304     363     3.9     25.8
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 539     0       0       815     2431    815     2431    0.0     1301
gi|9626372|ref|NC_001422.1|     NW_019434829.1  35.135  37      24      0       833     943     1000    890     9.5     36.3
gi|9626372|ref|NC_001422.1|     NW_019525440.1  35.135  37      24      0       833     943     153264  153154  9.5     36.3
gi|9626372|ref|NC_001422.1|     NC_001422.1     100.000 283     0       0       850     2       850     2       0.0     578
```

To look at this more closely, we can find the boundary where we only find the PhiX genome alignments

Here I sort the output by the subject (target) sequence in the assembly then require that the sequence identity (column 3) is greater than 50 and the alignment length (column 4) is greater than 30 bp.  With that filter on, all the false positives fall away and we are left with only the known PhiX scaffold.

```
more phiX_2_seriola.tblastxout | sort -k 2n | awk '$3>50 && $4>30' | awk '{print $2}' | sort | uniq -c
     45 NC_001422.1
```


# Further Reading

* [Large-scale contamination of microbial isolate genomes by Illumina PhiX control](https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/1944-3277-10-18)
