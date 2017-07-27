# Gene Ontology enrichment analysis

The [Gene Ontology](http://www.geneontology.org/) (GO) project provides structured description for genes known biological information at different levels. The controlled vocabularies (fixed words) describe gene function at 3 levels:
- biological processes
- molecular functions
- cellular components

These terms are species independent and hierarchically structured (see more information [here](http://geneontology.org/page/ontology-documentation)).

Enrichment analysis is a test to see a small subset of genes when sampled from large set of genes (reference set), what is the probability that small subset of genes (or statistically large proportion of subset genes) belong to a functional category as opposed to a randomly sampled subset of genes. This is normally done using either hypergeometric test (test without replacement) or using binomial test (with replacement).

## GO enrichment using Ontologizer
Before you begin, you need to have 2 important files. 
  1. GO ontology file, where you describe every possible GO term numbers with what they are. You can easily get this from here: http://purl.obolibrary.org/obo/go.obo 
  2. Mapping file, this will describe the association of genes of the organism with GO terms (this can be generated using Blast2Go- by loading swissprot blast xml files and interproscan and exporting it in .anno file). 

Once you have these 2 files, you can run Ontologizer. Ontologizer is a command line tool which very effecient when you have large number of list to analyze. You need to download the `jar` file to run this
The command:

```
java -jar Ontologizer.jar \
   -a association.anno \
   -g gene_ontology.obo \
   -s your_input_list.txt \
   -p population.txt \
   -c Parent-Child-Union \
   -m Westfall-Young-Single-Step \
   -d 0.05 \
   -r 1000
```

Here the `population.txt` is basically the full list of genes that you have in anno file. The only time you need to change this is when you have a different background set to test. You can also play around the other settings like `--mode`, `--resamplingsteps` to optimize. With the `-d` it will also generate a `dot` file, that can be used with `GraphViz` to show the pathway where these genes are enriched.

You will see a `dot` file after the run is complete. You can convert this to `png` format to view them easily:

```
dot -Tpng input.dot -o output.png
```
