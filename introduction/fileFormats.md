# File Formats
---
## Learning Objective
Upon completion of this section you will have a better understanding of the following file formats, how to read them and interpret the information they contain.

* FASTA – plain sequences
* FASTQ – sequencing reads
* GFF – gene models
* VCF – sequence variants
* SAM – sequence alignments
* BAM – alignments in binary


[Table of contents](/index.md)

## FASTA

Text file format for storing sequences
for nucleotide & amino acid data.  For a given sequence, a single line description and ID is supplied followed by one or more lines of sequence.  Multiple sequences can be placed in a single file and empty lines are typically ignored by programs.  The recommended number of sequence characters per line is 60 – 80.

```Line 1:``` starts with “>” followed by ID  
```Line 2:``` Sequence data  

##### Examples
```
>gi|296581|emb|Z22600.1| D.tigrina homeodomain mRNA
ttcgcggttcataactacctgacgaggttgagacggtacgagctggcggtggccctcaatcttaacgaaa
gacagataaaagtttgg
```
```
>gi|425153|gb|L26238.1|MUSHOME Mus domesticus (lbx) homeodomain mRNA, partial cds
CCATTTCAACAAGTACCTGACCAGGGCTCGGCGAGTGGAAGTTGCCGCTATTCTCGAGCTCAACGAAACT
CAAGTGAAAATT
```
```
>gi|2695691|gb|AAC36493.1| BRCA1 [Rattus norvegicus]
MDLSAVRIQEVQNVLHAMQKILECPICLELIKEPVSTQCDHIFCKFCMLKLLNQKKGPSQCPLCKNEITK
RSLQGSARFSQLVEELLKIIDAFELDTGMQCANGFSFSKKKNSSSELLNEDASIIQSVGYRNRVKKLQQI
ESGSATLKDSLSVQLSNLGIVRSMKKNRQTQPQNKSVYIALESDSSEERVNAPDGCSVRDQELFQIAPGG
AGDEGKLNSAKKAACDFSEGIRNIEHHQCSDKDLNPTENHATERHPEKCPRISVANVHVEPCGTDARASS
```

##### Common Errors that occur with this file type

* Program requires the sequences to all be on a single but the fasta file is on multiple lines
* Program requies the sequences to be on multiple lines with a string length per line less than 80 characters but the sequences are written on a single lines

---

## FASTQ File Format

FASTQ files are similar to FASTA but contain the quality score of the sequence data (only nucleotide sequences). The format contains two additional lines beyond FASTA format.

```Line 1:``` starts with “@” followed by ID  
```Line 2:``` Sequence data  
```Line 3:``` Starts with “+”       rest of the description is optional  
```Line 4:``` Quality score for each base in the sequence


Text file format for storing sequences and its quality scores (only nucleotide sequences)
File extension is .fastq or .fq
4 lines per entry

```
@HISEQ:402:H147CADXX:1:1101:1250:2208 1:N:0:CGATGT
TGATGCTGCNAATTTTATTCAGTCAGCGGAGGGGGCTTACGTGTATTTTCTGCAACCTTT
+
CCCFFFFFH#4AFIJJJJJJJJIJJJJJJJJJJJJJJJJJJHHHHHHFFFFFFFEEEEED
```

##### Quality score

* Probability of an error in base calling
* Higher score means low probability of error
* Different types of encoding available (Sanger, Phred, etc..)
* Different sequencing technologies use different encoding

For more information on Quality Score encoding see [Fastq Quality score Encoding](introduction/fastqquality-Score-encoding.md)


---

## GFF3: General Feature Format

This is a nine column tab separated text file that stores information about gene annotation.

```Column 1```  seqID (e.g. chromosome/scaffold, genome id, etc..)  
```Column 2```  Source (program used to generate or location of download)  
```Column 3```  Feature type (gene, mRNA, CDS, exon, etc.)  
```Column 4```  Start position of feature  
```Column 5```  End position of feature  
```Column 6```  Score (some program outputs will have a score of confidence for feature)  
```Column 7```  Strand (+,-,.)  
```Column 8```  Phase  
```Column 9```  List of attributes in the format tag=value. Multiple attributes are separated by “;”

Undefined fields are replaced with “.” character


```
##gff-version 3
##date Thu Nov  7 15:29:10 2013
##source gbrowse gbgff gff3 dumper
Chr1	TAIR9	gene	813471	816749	.	+	.	ID=AT1G03310;Note=protein_coding_gene;Name=AT1G03310
Chr1	TAIR9	mRNA	813471	816749	.	+	.	ID=AT1G03310.1;Parent=AT1G03310;Name=AT1G03310.1;Index=1
Chr1	TAIR9	protein	813975	816623	.	+	.	ID=AT1G03310.1-Protein;Name=AT1G03310.1;Derives_from=AT1G03310.1
Chr1	TAIR9	exon	813471	813581	.	+	.	Parent=AT1G03310.1
Chr1	TAIR9	five_prime_UTR	813471	813581	.	+	.	Parent=AT1G03310.1
Chr1	TAIR9	exon	813929	816749	.	+	.	Parent=AT1G03310.1
Chr1	TAIR9	five_prime_UTR	813929	813974	.	+	.	Parent=AT1G03310.1
Chr1	TAIR9	CDS	813975	816623	.	+	0	Parent=AT1G03310.1,AT1G03310.1-Protein;
Chr1	TAIR9	three_prime_UTR	816624	816749	.	+	.	Parent=AT1G03310.1
Chr1	TAIR9	mRNA	813486	816749	.	+	.	ID=AT1G03310.2;Parent=AT1G03310;Name=AT1G03310.2;Index=1
Chr1	TAIR9	protein	813975	816623	.	+	.	ID=AT1G03310.2-Protein;Name=AT1G03310.2;Derives_from=AT1G03310.2
Chr1	TAIR9	exon	813486	813604	.	+	.	Parent=AT1G03310.2
Chr1	TAIR9	five_prime_UTR	813486	813604	.	+	.	Parent=AT1G03310.2
Chr1	TAIR9	exon	813929	816749	.	+	.	Parent=AT1G03310.2
Chr1	TAIR9	five_prime_UTR	813929	813974	.	+	.	Parent=AT1G03310.2
Chr1	TAIR9	CDS	813975	816623	.	+	0	Parent=AT1G03310.2,AT1G03310.2-Protein;
Chr1	TAIR9	three_prime_UTR	816624	816749	.	+	.	Parent=AT1G03310.2
Chr1	TAIR9	gene	795532	798463	.	-	.	ID=AT1G03260;Note=protein_coding_gene;Name=AT1G03260
Chr1	TAIR9	mRNA	795532	798463	.	-	.	ID=AT1G03260.1;Parent=AT1G03260;Name=AT1G03260.1;Index=1
Chr1	TAIR9	protein	795678	798102	.	-	.	ID=AT1G03260.1-Protein;Name=AT1G03260.1;Derives_from=AT1G03260.1
Chr1	TAIR9	five_prime_UTR	798103	798463	.	-	.	Parent=AT1G03260.1
Chr1	TAIR9	CDS	798001	798102	.	-	0	Parent=AT1G03260.1,AT1G03260.1-Protein;
Chr1	TAIR9	exon	798001	798463	.	-	.	Parent=AT1G03260.1
Chr1	TAIR9	gene	799191	802436	.	+	.	ID=AT1G03270;Note=protein_coding_gene;Name=AT1G03270
Chr1	TAIR9	mRNA	799191	802436	.	+	.	ID=AT1G03270.1;Parent=AT1G03270;Name=AT1G03270.1;Index=1
Chr1	TAIR9	protein	799191	802436	.	+	.	ID=AT1G03270.1-Protein;Name=AT1G03270.1;Derives_from=AT1G03270.1
Chr1	TAIR9	exon	799191	799431	.	+	.	Parent=AT1G03270.1
```

Example from GFF3 definition by Lincoln Stein

![](https://github.com/The-Sequence-Ontology/Specifications/blob/master/img/figure1.png)

#####  More information
* [GFF3 definition](https://raw.githubusercontent.com/The-Sequence-Ontology/Specifications/master/img/figure1.png)
