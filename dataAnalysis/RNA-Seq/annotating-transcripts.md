# Running Trinotate for annotating the transcripts:

[Trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki) can be used to annotate the transcripts. The files used in this example are as follows:

1. Input fasta file `trinity.fasta`
2. Databases:
  * uniprot_sprot.pep.gz
  * Pfam-A.hmm.gz
3. TransDecoder output: `longest_orfs.pep`


## Database downloads
```

wget https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/uniprot_sprot.pep.gz
wget https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Pfam-A.hmm.gz
gunzip uniprot_sprot.pep.gz
gunzip Pfam-A.hmm.gz
makeblastdb -in uniprot_sprot.pep -dbtype prot
hmmpress Pfam-A.hmm
```


Next, database searches and predictions were carried out:
If you haven't run the TransDecoder on your `trinity.fasta`, you can run it as follows:

## TransDecoder

```
module load transdecoder
TransDecoder.LongOrfs -m 10 -t trinity.fa
# you will need `longest_orfs.pep` for next steps
```


## Searches

```
blastx -query trinity.fasta \
  -db uniprot_sprot.pep \
  -num_threads 8 \
  -max_target_seqs 1 \
  -outfmt 6 > blastx.outfmt6
  
blastp -query longest_orfs.pep \
  -db uniprot_sprot.pep \
  -num_threads 8 \
  -max_target_seqs 1 \
  -outfmt 6 > blastp.outfmt6
  
hmmscan --cpu 8 \
  --domtblout TrinotatePFAM.out \
  Pfam-A.hmm longest_orfs.pep > pfam.log
  
signalp -f short \
  -n signalp.out longest_orfs.pep
  
tmhmm --short < longest_orfs.pep > tmhmm.out

RnammerTranscriptome.pl --transcriptome ttrinity.fasta \
  --path_to_rnammer /usr/bin/software/rnammer_v1.2/rnammer
```

## Loading results

Trinotate SQLite Database was updated with the new predictions:

```
wget "https://data.broadinstitute.org/Trinity/Trinotate_v3_RESOURCES/Trinotate_v3.sqlite.gz" -O Trinotate.sqlite.gz

gunzip Trinotate.sqlite.gz

get_Trinity_gene_to_trans_map.pl trinity.fasta >  Trinity.fasta.gene_trans_map

Trinotate Trinotate.sqlite init \
  --gene_trans_map Trinity.fasta.gene_trans_map \
  --transcript_fasta trinity.fasta \
  --transdecoder_pep longest_orfs.pep
  
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
Trinotate Trinotate.sqlite LOAD_signalp signalp.out
```

and finally, report was generated as follows:

## Report:

```
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
```


## GO Annotation (optional)

Using the above report, you can assign GO for your sequences as follows:

```
${TRINOTATE_HOME}/util/extract_GO_assignments_from_Trinotate_xls.pl  \
  --Trinotate_xls trinotate_annotation_report.xls \
  -G --include_ancestral_terms \
  > go_annotations.txt
```

## More information

Trinotate does not yet have its own paper but is recommended to cite the following:

* Bryant, D.M., Johnson, K., DiTommaso, T., Tickle, T., Couger, M.B., Payzin-Dogru, D., Lee, T.J., Leigh, N.D., Kuo, T.H., Davis, F.G. and Bateman, J., 2017. [A tissue-mapped axolotl de novo transcriptome enables identification of limb regeneration factors](https://pubmed.ncbi.nlm.nih.gov/28099853/). Cell reports, 18(3), pp.762-776. 
