---
title: "Introduction to Data Wrangling"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Calculate Sequence Length (fasta)
Sometimes it is essential to know the length distribution of your sequences. It may be your newly assembled scaffolds or it might be a genome, that you wish to know the size of chromosomes, or it could just be any multi fasta sequence file.
##  1. Using biopython ##
Save this as a script, make it an executable and run on a fasta file:
```
#!/usr/bin/python
from Bio import SeqIO
import sys
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
 output_line = '%s\t%i' % \
(seq_record.id, len(seq_record))
 print(output_line)
```
To run:
```
chmod +x seq_length.py
seq_length.py input_file.fasta
```
This will print length for all the sequences in that file.


##  2. Using bioawk

Bioawk is an extension of the <blockcode>awk</blockcode> written by [Heng Li](https://github.com/lh3).  It is available to donwload from this [link](https://github.com/lh3/bioawk/releases). Installation is easy too. To get sequence length, run it as:

```
bioawk -c fastx '{print $name length($seq)}' input.fasta
```
Output will be similar to the above script and can be redicrected to any file if you want.

# More information
*  [Introduction to Bioawk](../../Appendix/Unix/bioawk-basics.md)

[Table of contents](https://isugenomics.github.io/bioinformatics-workbook/)
