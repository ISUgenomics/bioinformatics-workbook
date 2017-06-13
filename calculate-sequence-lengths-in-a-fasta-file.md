Sometimes it is essential to know the length distribution of your sequences. It may be your newly assembled scaffolds or it might be a genome, that you wish to know the size of chromosomes, or it could just be any multi fasta sequence file.
<h3> 1. Using biopython </h3>
<p>Save this as a script, make it an executable and run on a fasta file:</p>
<python>
#!/usr/bin/python
from Bio import SeqIO
import sys
cmdargs = str(sys.argv)
for seq_record in SeqIO.parse(str(sys.argv[1]), "fasta"):
 output_line = '%s\t%i' % \
(seq_record.id, len(seq_record))
 print(output_line)
</python>
<p>To run:</p>
<code>
chmod +x seq_length.py
seq_length.py input_file.fasta
</code>
<p>This will print length for all the sequences in that file.</p>
