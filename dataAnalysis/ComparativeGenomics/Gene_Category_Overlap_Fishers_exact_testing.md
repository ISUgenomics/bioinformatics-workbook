# GeneOverlap R package / Fisher's exact testing
Identifying genes that are significantly associated with a feature or characteristic of a transcriptome/genome is a common bioinformatics test. <br/>
The [GeneOverlap]("https://www.rdocumentation.org/packages/GeneOverlap/versions/1.8.0") package offers a quick and easy way to accomplish this task.<br/>
<br/>
Perhaps you want to see if a mutation has LTR retrotransposons

Grab relevant files
```
#get your genome's gff
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_genomic.gff.gz

#grab your differentially expressed gene list. (I took this from the RNA-seq tutorial on DESEQ)
less DEGArabidopsisSiva.csv

#get your LTR retrotransposon coordinates. (I took this from the LTR_finder tutorial)
ln -s ../08_LTRFinder/LTRfinder.bed
```

###Rename the scaffolds so that the DESEQ analysis genes  will match the genes identified in LTRfinder analysis
Since the gene overlap will essentially be found using bedtools, the scaffold names need to match.
```
#the easiest and fastest way to do this is to remove all other columns except scaffold name, so sed isnt searching the entire gff
less GCF_000001735.3_TAIR10_genomic.gff |awk '$3=="gene" {print $1}' >GeneScaffNames

Then grab the unique scaffold names and renamed them with sed.
less GCF_000001735.3_TAIR10_genomic.gff |awk '$3=="gene" {print $1}' |sort|uniq|less

sed -i 's/NC_000932.1/Plastid/g' GeneScaffNames
sed -i 's/NC_001284.2/Mitochondria/g' GeneScaffNames
sed -i 's/NC_003070.9/Chr1/g' GeneScaffNames
sed -i 's/NC_003071.7/Chr2/g' GeneScaffNames
sed -i 's/NC_003074.8/Chr3/g' GeneScaffNames
sed -i 's/NC_003075.7/Chr4/g' GeneScaffNames
sed -i 's/NC_003076.8/Chr5/g' GeneScaffNames
paste GeneScaffNames <(awk '$3=="gene"' GCF_000001735.3_TAIR10_genomic.gff) |awk '{print $1,$3,$4,$5,$6,$7,$8,$9,$10}' |tr " " "\t" >MOD.gff
```


### Create gene lists for GeneOverlap
Now we need to create some genes lists: genes in ltr retrotransposons, genes in the genome, and a set of genes that are differentially upregulated (say top 10 percentile of genes in log fold change).
```
 bedtools intersect -wo -a <(sed 's/^/Chr/g' LTRfinder.bed) -b MOD.gff  |awk '{print $13}' |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq >LTRGene.list
#total count of all genes in the genome
  less MOD.gff  |awk '{print $9}' |sed 's/ID=//g' |sed 's/;/\t/g' |cut -f 1 |sort|uniq|wc -l >All.genes

  #figure out the top 10 percentile
  less DEGArabidopsisSiva.csv |wc
    38159   38159  868807
#so find the top 3816 genes with the high log fold change
less DEGArabidopsisSiva.csv |sed 's/,/\t/g' |sort -k2,2nr |awk 'NR<3817' >Top10Percentile
#what the heck, lets find the lowest also
less DEGArabidopsisSiva.csv |sed 's/,/\t/g' |sort -k2,2nr |awk 'NR>34343' >Bottom10Percentile
```

### Time to load it all into R

```
module load R
library(GeneOverlap)
LTRGene.list <- read.table("LTRGene.list")
Top10Percentile <- read.table("Top10Percentile")
Bottom10Percentile <- read.table("Bottom10Percentile")
All.genes <- read.table("All.genes")
```
Perform the comparison
```
go.obj <- newGeneOverlap(LTRGene.list$V1, Top10Percentile$V1, genome.size=All.genes)
go.obj1 <- testGeneOverlap(go.obj)
go.obj <- newGeneOverlap(LTRGene.list$V1, Bottom10Percentile$V1, genome.size=All.genes)
go.obj2 <- testGeneOverlap(go.obj)
```
Concatenate the output and push into a file.
```
AllCombine <- c(go.obj1, go.obj2)

library(RJSONIO)
 exportJSON <- toJSON(AllCombine)
 write(exportJSON,"LTRExpression.out")
#exit R
cat <(ls -1 Top10Percentile ) <(ls -1 Bottom10Percentile ) |paste - <(less LTRExpression.out|sed 's/\[/\t/g' |awk 'NF<15' |sed 's/{//g' |sed 's/}//g' |sed 's/]//g' |sed 's/"//g' |sed 's/,/\t/g' |sed '/^$/d' | tr "\n" "\t" |sed 's/true/true\n/g' |sed 's/false/false\n/g') |less
 ```
 ### Results
 ```
 Top10Percentile                    notA: 29431            inA: 68                          notA: 3805             inA: 11               odds.ratio: 1.2512      pval: 0.29241   Jaccard: 0.0028321      is.tested: true
Bottom10Percentile                         notA: 29427            inA: 72                          notA: 3809             inA: 7                odds.ratio: 0.75111     pval: 0.81441   Jaccard: 0.0018004      is.tested: true

 ```
 Well, it looks like there are 11 LTR retrotransposon associated genes that are upregulated in expression and 7 that are downregulated.  As a whole, LTR retrotransposon gene expression does not seem to be associated with the changes that have occurred between the two arabidopsis lines analyzed here.   
