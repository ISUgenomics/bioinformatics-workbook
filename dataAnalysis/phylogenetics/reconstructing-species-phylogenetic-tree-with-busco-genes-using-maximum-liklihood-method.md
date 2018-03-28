# Building maximum likelihood phylogenetic tree using BUSCO genes

In this guide, we will explain how to get the species phylogenetic tree using just the draft genomes (genomes that still don't have the annotations yet). For data, we will use the sequenced legume genomes. You will need the following programs installed and access to a large cluster (we use Condo in this exercise).



## Data

Like mentioned before, we will use the sequenced legume genomes that are available on NCBI (both RefSeq and GenBank). Since it is most likely that the new genome you're working won't have good gene predictions, we will use the BUSCO genes to construct phylogenetic tree. [BUSCO](https://doi.org/10.1093/bioinformatics/btv351) (Benchmarking Universal Single-Copy Orthologs) are ideal for this because, these genes are found in a genome in single-copy are evolutionarily sound. There are 1,440 BUSCO's for `embryophyta` profile, which will be a good sample size for building phylogenetic tree in reasonable amount of time.

**Dataset** : At the time of writing this tutorial, there were 38 genomes in NCBI belonging to Legumes group ([Fabaceae](https://www.ncbi.nlm.nih.gov/assembly/?term=Fabaceae)). Out of these 38, 17 genomes have chromosome level assemblies, 11 genomes as scaffolds, and 10 genomes as contigs. See table 1 for more information on these genomes.

**Table 1:** Legume genomes used in this study.

| Organism                                            | SpeciesTaxid | RefSeq_category       | Genbank         | ScaffoldN50 |
|-----------------------------------------------------|--------------|-----------------------|-----------------|-------------|
| _Arachis duranensis_ (eudicots)                     | 130453       | representative genome | GCA_000817695.2 | 1,068,268   |
| _Arachis duranensis_ (eudicots)                     | 130453       | na                    | GCA_001687015.1 | 149,039     |
| _Arachis ipaensis_ (eudicots)                       | 130454       | representative genome | GCA_000816755.2 | 6,169,993   |
| _Cajanus cajan_ (pigeon pea)                        | 3821         | representative genome | GCA_000340665.1 | 555,764     |
| _Cajanus cajan_ (pigeon pea)                        | 3821         | na                    | GCA_000230855.2 | 5,341       |
| _Cicer arietinum_ (chickpea)                        | 3827         | representative genome | GCA_000331145.1 | 697,963     |
| _Cicer arietinum_ (chickpea)                        | 3827         | na                    | GCA_000347275.3 | 74,024      |
| _Cicer arietinum_ (chickpea)                        | 3827         | na                    | GCA_002896005.1 | 232,601     |
| _Cicer echinospermum_ (eudicots)                    | 90897        | representative genome | GCA_002896215.1 | 206,896     |
| _Cicer reticulatum_ (eudicots)                      | 90898        | representative genome | GCA_002896235.1 | 109,263     |
| _Glycine max_ (soybean)                             | 3847         | representative genome | GCA_000004515.3 | 48,577,505  |
| _Glycine max_ (soybean)                             | 3847         | na                    | GCA_002905335.1 | 15,020,773  |
| _Glycine max_ (soybean)                             | 3847         | na                    | GCA_001269945.2 | 42,937      |
| _Glycine soja_ (wild soybean)                       | 3848         | representative genome | GCA_002907465.1 | 4,430,511   |
| _Glycine soja_ (wild soybean)                       | 3848         | na                    | GCA_000722935.2 | 404,776     |
| _Lotus japonicus_ (eudicots)                        | 34305        | representative genome | GCA_000181115.2 | 25,054      |
| _Lupinus angustifolius_ (narrow-leaved blue lupine) | 3871         | representative genome | GCA_001865875.1 | 702,610     |
| _Lupinus angustifolius_ (narrow-leaved blue lupine) | 3871         | na                    | GCA_002285895.2 | 30,409,451  |
| _Lupinus angustifolius_ (narrow-leaved blue lupine) | 3871         | na                    | GCA_000338175.1 | 15,485      |
| _Medicago truncatula_ (barrel medic)                | 3880         | representative genome | GCA_000219495.2 | 49,172,423  |
| _Medicago truncatula_ (barrel medic)                | 3880         | na                    | GCA_002024945.1 | 12,848,239  |
| _Medicago truncatula_ (barrel medic)                | 3880         | na                    | GCA_002251925.1 | 575,111     |
| _Medicago truncatula_ (barrel medic)                | 3880         | na                    | GCA_002251935.1 | 859,387     |
| _Medicago truncatula_ (barrel medic)                | 3880         | na                    | GCA_002251955.1 | 503,143     |
| _Phaseolus vulgaris_ (eudicots)                     | 3885         | representative genome | GCA_000499845.1 | 50,367,376  |
| _Phaseolus vulgaris_ (eudicots)                     | 3885         | na                    | GCA_001517995.1 | 433,759     |
| _Trifolium pratense_ (eudicots)                     | 57577        | representative genome | GCA_900079335.1 | 22,682,783  |
| _Trifolium pratense_ (eudicots)                     | 57577        | na                    | GCA_000583005.2 | 2,432       |
| _Trifolium subterraneum_ (eudicots)                 | 3900         | representative genome | GCA_001742945.1 | 287,605     |
| _Trifolium subterraneum_ (eudicots)                 | 3900         | na                    | GCA_002003065.1 | 1,510       |
| _Vicia faba_ (fava bean)                            | 3906         | representative genome | GCA_001375635.1 | 1,723       |
| _Vigna angularis_ (adzuki bean)                     | 3914         | representative genome | GCA_001190045.1 | 1,292,063   |
| _Vigna angularis var. angularis_ (adzuki bean)      | 3914         | na                    | GCA_001723775.1 | 8,174,047   |
| _Vigna angularis var. angularis_ (adzuki bean)      | 3914         | na                    | GCA_000465365.1 | 21,697      |
| _Vigna radiata var. radiata_ (mung bean)            | 157791       | representative genome | GCA_000741045.2 | 25,360,630  |
| _Vigna radiata var. radiata_ (mung bean)            | 157791       | na                    | GCA_001584445.1 | 683,756     |
| _Vigna radiata_ (mung bean)                         | 157791       | na                    | GCA_000180895.1 | 228         |
| _Vigna unguiculata_ subsp. unguiculata (cowpea)     | 3917         | representative genome | GCA_001687525.1 | 7,412       |

The data can be downloaded from using the commands:


### Download the Data

For downloading, you can do it directly from the NCBI website. To do so, go to [NCBI genomes](https://www.ncbi.nlm.nih.gov/assembly/?term=Fabaceae) page for Fabaceae, click on `Download Assemblies`, select the `GenBank` as source database and `Genomic FASTA` as the file type, followed by `Download`

![](assets\fig1.png)

To easily indetify these files, we will rename them.


Now we are all set to begin the analysis!

## Run BUSCO to predict conserved genes

Running BUSCO is straight forward. The BUSCO for plants have 1,440 total genes and running BUSCO will predict those genes from the genomes. Of course, the number will vary depending on how complete the genomes are but we can use the common genes that were found in 3 or more genomes to build the gene tree and then use it for inferring the species tree.

For running BUSCO we will use `runBUSCO.sh`. This script will make it easier to run on several genomes by limiting the complexity and standardizing analyses to best parameters. The script can be found [here](https://github.com/ISUgenomics/common_analyses/blob/master/runBUSCO_genome.sh)

To run this in parallel we will use a for loop to generate commands for each species.

```
for file in *.fasta; do
echo "./runBUSCO.sh $file";
done >> busco.cmds
```


Next, generate multiple PBS scripts to run them on a cluster. For Condo the  `makeSLURMs.py` script can be used (script can be obtained from our repo [here](https://github.com/ISUgenomics/common_scripts/blob/master/makeSLURMs.py).

After this step, you should have 21 BUSCO submission scripts. Submit them to the cluster using for loop:

```
for file in busco_?.sub; do
sbatch $file;
done
```

### Optional step: Get the BUSCO plot for your genomes

Copy the short summary files from BUSCo run directories to another directory for getting the BUSCO plots

```
mkdir busco_summaries
for file in $(find . -type f -name "short_summary_*"); do
cp $file ./busco_summaries/;
done
```

Generate the BUSCO plot using the script in the BUSCO package


```
python3 ${BUSCO_HOME}/BUSCO_plot.py -wd busco_summaries
```
This should generate `busco_figure.R` and `busco_figure.png` files.

## Extract complete BUSCO genes

Next, we need to find all the "complete" BUSCO genes predicted from all of the genomes. Since 1,440 genes have unique ids, we can quickly get a list of all gene names that are complete (from BUSCO logs) and use them for downstream analysis.

```
for file in $(find . -name "full_table_*.tsv"); do
grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >> complete_busco_ids.txt;
done
```

To filter out the complete BUSCOs that are only present less than 3 genomes, we can use `uniq` and `awk` commands

```
sort complete_busco_ids.txt |uniq -c > complete_busco_ids_with_counts.txt
awk '$NF > 2 {print $1}' complete_busco_ids_with_counts.txt > final_busco_ids.txt
```

Now let's copy the genes (both nucleotides and amino acids) to different directory. Note that the gene names are just indetified as BUSCO id with either faa or fna extenstions. If you try to cooy them directly to a single directory, you will likely overwrite all the files and end up with only the last set of files. You can  give them a unique name and then merge them all together, writing them to a single busco id file. before you do this, you will also need to edit the sequence name (fasta header) so that it includes organism identifier in it.

```
mkdir -p busco_aa
mkdir -p busco_nt
for dir in $(find . -type d -name "single_copy_busco_sequences"); do
  sppname=$(basename $(dirname $dir)|cut -f 2-3 -d "_" | sed 's/_/ /g');
  abbrv=$(echo $sppname | sed 's/\(\w\w\w\)\w*\( \|$\)/\1/g')
  for file in ${dir}/*.faa; do
    cp $file busco_aa/${abbrv}_${file}
    sed -i 's/^>/>'${abbrv}'|/g' busco_aa/${abbrv}_${file}
	cut -f 1 -d ":" busco_aa/${abbrv}_${file} | tr '[:lower:]' '[:upper:]' > busco_aa/${abbrv}_${file}.1
	mv busco_aa/${abbrv}_${file}.1 busco_aa/${abbrv}_${file}
  done
  for file in ${dir}/*.fna; do
    cp $file busco_nt/${abbrv}_${file}
    sed -i 's/^>/>'${abbrv}'|/g' busco_nt/${abbrv}_${file}
	cut -f 1 -d ":" busco_nt/${abbrv}_${file} | tr '[:lower:]' '[:upper:]' > busco_nt/${abbrv}_${file}.1
	mv busco_nt/${abbrv}_${file}.1 busco_nt/${abbrv}_${file}  done
done
```

This will convert the spp name to an abbreviation (from `Vicia_faba` to `Vicfab`) and edit the fasta header to include it in the name (eg. >EOG09360J1D to `>VICFAB|EOG09360J1D`). Next step is to generate single file for each of the BUSCO id. We will do this with a simple bash loop. Remember the file `final_busco_ids.txt` we created previously? we need that for this purpose:

```
while read line; do
cat ??????_${line}.faa >> ${line}_aa.fasta;
cat ??????_${line}.fna >> ${line}_nt.fasta;
done<final_busco_ids.txt
```
As we won't be needing the individual fasta files (because we now how a single file for each BUSCO gene that contains complete sequence from all the species), let's go ahead and delete them (alternatively, you can move them to another directory as backup copy, just in case)

```
rm ??????_EOG09360???.f?a
```
## Alignment

NOTE: we will only use the peptide sequence for building tree here, but you can easily change this to do the same with nucleotide sequence as well. We will use the `makeSLURMs.py` found [here](https://github.com/ISUgenomics/common_scripts) and `runGuidenceAA.sh` found [here](https://github.com/ISUgenomics/common_scripts/blob/master/runGuidenceAA.sh) for running alignments.

```
for file in *_aa.fasta; do
echo "./runGuidenceAA.sh $file";
done > guidence.cmds
makePBSs.py 200 guidence.cmds
for sub in guidence_*.sub; do
sbatch $sub;
done
```

Once the alignment completes, you will see a folder (same name as the input sequence but without fasta extension) and in each directory will have many files. we will only be needing `Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names` file for the next step. You also need to make sure that the alignment completed successfully by checking if there is any error message in the `log` file.

Next, let's copy the final alignment file to a different directory where we will run the RAxML.

```
mkdir raxml_aa
# from the directory where you ran guidence from
for aln in $(find . -name "Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names"); do
busco=$(dirname $(basename $aln));
cp $aln ./raxml_aa/${busco}.aln;
done
```

## Gene tree reconstruction

Gene phylogenetic tree will be reconstructed with maximum likelihood (ML) approach using RAxML  program. We will use a rapid bootstrap analysis and search for the best-scoring ML tree with 1000 bootstrap replicates. Best protein model with respect to the likelihood on a fixed, reasonable tree will be automatically determined using the PROTGAMMAAUTO option of RAxML (other option is to test which evolutionary model fits your data by running `prottest` or other similar programs). _Arabidopsis thaliana_ will be the outgroup.

Again, we will use the `runRAXML.sh` (found [here](https://github.com/ISUgenomics/common_scripts/blob/master/runRAxML.sh)) and `makeSLURMs.py` scripts for running the jobs:

```
cd raxml_aa
for aln in *.aln; do
echo "./runRAXML.sh $aln AA"
done > raxml.cmds
makePBSp.py 200 raxml.cmds
for sub in raxml_*.sub; do
sbatch $sub;
done
```

Once complete, you will have `RAxML_bestTree.EOG09360???` file that contains the best scoring tree (not the bootstrap consensus tree). The branch scores represent the likelihood score. If you want to get the bootstrap support trees you should be using `RAxML_bootstrap.EOG09360???` file that will have 1000 bootstrap replicates to generate a consensus tree. Here, we will use the ML tree.

To generate a single file with ML trees for all BUSCO genes:

```
cat RAxML_bestTree.EOG09360??? >> busco_genes_ML.tree
```

## Species tree reconstruction

Before we proceed with the species tree, we need to fix the title (or the taxa name) in the gene trees. As even single tax in these gene trees are identified as `genspp|buscoid`, they all have unique ids and it will be considered as seperate taxa. We don't want this as we only have 22 tax here. For this, we need to remove the `buscoid` from the name (i.e, convert names that look like this `LOTJAP_03I2` to just `LOTJAP`). We will use the sed command for that purpose.

```
# create a file with all the abbreviation's
# pick any one RAXML best tree from teh raxml_aa folder
grep -oE "......_...." RAxML_bestTree.EOG093603IM |cut -f 1 -d "_" > capabrv
# edit the merged ML tree file to replace it
while read line; do
sed -i 's/'${line}'_....:/'${line}':/g' busco_genes_ML.tree;
done<capabrv
```

As a last step, we will generate a species tree using the [ASTRAL](https://www.ncbi.nlm.nih.gov/pubmed/25161245 "ASTRAL: genome-scale coalescent-based species tree estimation") coalescent-based species tree estimation program. ASTRAL is the fast method for estimating species trees from multiple genes which is statistically consistent, and can run on multiple gene trees efficiently.


```
java -jar ${ASTRAL_HOME}/astral.jar --input busco_genes_ML.tree --output busco_genes_ASTRAL_spp.tree
```

You will get a newick format tree file that can be uploaded to `EvolView` to generate publication quality graphics (tree). It can be accessed from this [link](http://www.evolgenius.info/evolview "online phylogenetic tree viewer and customization tool").
