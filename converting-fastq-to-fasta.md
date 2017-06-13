<p>There are several ways you can convert fastq to fasta sequences. Some methods are listed below. </p>

<h3>Using SED </h3>
<p><blockcode>sed</blockcode> can be used to selectively print the desired lines from a file, so if you print the first and 2rd line of every 4 lines, you get the sequence header and sequence needed for fasta format. </p>
<pre>
sed -n '1~4s/^@/>/p;2~4p' INFILE.fastq > OUTFILE.fasta
</pre>

<h3> Using PASTE</h3>
<p>You can linerize every 4 lines in a tabular format and print first and second field using <blockcode>paste</blockcode></p>
<pre>
cat INFILE.fastq | paste - - - - |cut -f 1, 2| sed 's/@/>/'g | tr -s "/t" "/n" > OUTFILE.fasta
</pre>

<h3> EMBOSS:seqret </h3>
<p>Standard script that can be used for many purposes. One such use is fastq-fasta conversion</p>
<pre>
seqret -sequence reads.fastq -outseq reads.fasta
</pre>
<p><blockcode>awk</blockcode> can be used for conversion as follows:</p>
<h3> Using AWK</h3>

<pre>
cat infile.fq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > file.fa
</pre>


<h3> FASTX-toolkit</h3>
<p><blockcode>fastq_to_fasta</blockcode> is available in the FASTX-toolkit that scales really well with the huge datasets</p>
<pre>
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
</pre>

<h3> Bioawk</h3>
<p>Another option to convert fastq to fasta format using <blockcode>bioawk</blockcode> </p>
<pre>
bioawk -c fastx '{print ">"$name"\n"$seq}' input.fastq > output.fasta
</pre>

<h3> Seqtk</h3>
<p>From the same developer, there is another option using a tool called <blockcode>seqtk</blockcode> </p>
<pre>
seqtk seq -a input.fastq > output.fasta
</pre>
Note that you can use either compressed or uncompressed files for this tool
