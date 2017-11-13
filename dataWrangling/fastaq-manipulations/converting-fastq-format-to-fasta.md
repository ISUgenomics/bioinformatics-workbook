# Converting FASTQ format to FASTA

There are several ways you can convert fastq to fasta sequences. Some methods are listed below.

### Using SED
`sed` can be used to selectively print the desired lines from a file, so if you print the first and 2rd line of every 4 lines, you get the sequence header and sequence needed for fasta format.
```
sed -n '1~4s/^@/>/p;2~4p' INFILE.fastq > OUTFILE.fasta
```

###  Using PASTE
You can linerize every 4 lines in a tabular format and print first and second field using `paste`
```
cat INFILE.fastq | paste - - - - |cut -f 1, 2| sed 's/@/>/'g | tr -s "/t" "/n" > OUTFILE.fasta
```

###  EMBOSS:seqret
Standard script that can be used for many purposes. One such use is fastq-fasta conversion
```
seqret -sequence reads.fastq -outseq reads.fasta
```
`awk` can be used for conversion as follows:
###  Using AWK

```
cat infile.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > file.fa
```


###  FASTX-toolkit
`fastq_to_fasta` is available in the FASTX-toolkit that scales really well with the huge datasets
```
fastq_to_fasta -h
usage: fastq_to_fasta [-h] [-r] [-n] [-v] [-z] [-i INFILE] [-o OUTFILE]
# Remember to use -Q33 for illumina reads!
version 0.0.6
	   [-h]         = This helpful help screen.
	   [-r]         = Rename sequence identifiers to numbers.
	   [-n]         = keep sequences with unknown (N) nucleotides.
			       Default is to discard such sequences.
	   [-v]         = Verbose - report number of sequences.
			       If [-o] is specified,  report will be printed to STDOUT.
			       If [-o] is not specified (and output goes to STDOUT),
			       report will be printed to STDERR.
	   [-z]         = Compress output with GZIP.
	   [-i INFILE]  = FASTA/Q input file. default is STDIN.
	   [-o OUTFILE] = FASTA output file. default is STDOUT.
```

###  Bioawk
Another option to convert fastq to fasta format using `bioawk`
```
bioawk -c fastx '{print ">"$name"\n"$seq}' input.fastq > output.fasta
```

###  Seqtk
From the same developer, there is another option using a tool called `seqtk`
```
seqtk seq -a input.fastq > output.fasta
```
Note that you can use either compressed or uncompressed files for this tool

[Table of contents](index.md)
