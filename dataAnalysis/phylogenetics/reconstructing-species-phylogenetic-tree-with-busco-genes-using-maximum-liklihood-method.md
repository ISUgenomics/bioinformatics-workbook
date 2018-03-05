# Building maximum liklihood phylogenetic tree using BUSCO genes

In this guide, we will explain how to get the species phylogenetic tree using just the draft genomes (genomes that still don't have the annotations yet). For data, we will use the sequenced legume genomes. You will need the following programs installed and access to a large cluster (we use Condo in this exercise).



## Data

Like mentioned before, we will use the sequenced legume genomes for this purpose. Not all genomes had gene predictions hence we will predict the BUSCO genes before we construct phylogenetic tree. Following are the genomes used:

| Common Name  | Scientific Name                             | Genomes | Source (website)                                                           |
|--------------|---------------------------------------------|---------|----------------------------------------------------------------------------|
| Chickpea     | _Cicer arietinum_                           | 2       | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000331145.1)              |
| Wild peanut  | _Arachis duranensis_ and _Arachis ipaensis_ | 2       | [PeanutBase](https://peanutbase.org/peanut_genome)                         |
| peanut       | _Arachis hypogaea_ cv. Tifrunner            | 1       | [PeanutBase](https://peanutbase.org/peanut_genome)                         |
| Soybean      | _Glycine max_                               | 1       | [SoyBase](https://soybase.org)                                             |
| Wild soybean | _Glycine soja_                              | 1       | [NCBI](https://www.ncbi.nlm.nih.gov/genome/13239)                          |
| White Clover | _Trifolium pratense_                        | 2       | [CloverGARDEN](http://clovergarden.jp/)                                    |
| Red Clover   | _Trifolium repens_                          | 2       | [CloverGARDEN](http://clovergarden.jp/)                                    |
| Alfa-alfa    | _Medicago truncatula_                       | 1       | [EnsemblPlants](https://plants.ensembl.org/Medicago_truncatula/Info/Index) |
| Adzuki bean  | _Vigna angularis_                           | 1       | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_001190045.1)              |
| Mung bean    | _Vigna radiata_                             | 1       | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000741045.1)              |
| Lotus        | _Lotus japonicus_                           | 1       | [kazusa](http://www.kazusa.or.jp/lotus/)                                   |
| Lupin        | _Lupinus angustifolius_                     | 1       | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_001865875.1)              |
| Cowpea       | _Cajanus cajan_                             | 1       | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCA_001687525.1)              |
| Broad bean   | _Vicia faba_                                | 1       | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCA_001375635.1)              |
| Common bean  | _Phaseolus vulgaris_                        | 1       | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000499845.1)              |
| Pigeon bean  | _Cajanus cajan_                             | 1       | [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCA_000230855.2)              |
| Tale cress (outgrp)   | _Arabidopsis thaliana_                      | 1       | [araport11](https://www.araport.org/data/araport11)                        |


The data can be downloaded from using the commands:




```
# chickpea
wget ftp://ftp.bioinfo.wsu.edu/species/Cicer_arietinum/C.arietinum_CDCFrontier_v1.0/assembly/Cicer_arietinum_GA_v1.0_kabuli.chr_scf.fa
wget https://www.coolseasonfoodlegume.org/sites/default/files/files/C_arietinum_ICC4958_v2_0.fasta
# wild peanut
wget http://peanutbase.org/files/genomes/Arachis_ipaensis/assembly/Araip_v1.0_20140908.fa.gz
wget http://peanutbase.org/files/genomes/Arachis_duranensis/assembly/Aradu_v1.0_20140908.fa.gz
# soybean
wget https://soybase.org/GlycineBlastPages/archives/Gmax_275_v2.0.softmasked.20140305.fasta.zip
# wild soybean
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/722/935/GCA_000722935.2_W05v1.0/GCA_000722935.2_W05v1.0_genomic.fna.gz
# clover
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/fasta/trifolium_pratense/dna/Trifolium_pratense.Trpr.dna.toplevel.fa.gz
wget http://legumeinfo.org/data/public/Trifolium_pratense/MilvusB.gnm2.1/tripr.MilvusB.gnm2.1.gNmT.pchr_plus_unanchored.fna.gz
# red clover
wget ftp://ftp.kazusa.or.jp/pub/clover/Woogenellup/TSUw_r1.0.fasta.gz
wget ftp://ftp.kazusa.or.jp/pub/clover/Daliak/TSUd_r1.1.fasta.gz
# alfa-alfa
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-34/fasta/medicago_truncatula/dna/Medicago_truncatula.MedtrA17_4.0.dna.toplevel.fa.gz
# Adzuki bean
wget http://legumeinfo.org/data/public/Vigna_angularis/Gyeongwon.gnm3/vigan.Gyeongwon.gnm3.JyYC.pchr_plus_unanchored.fna.gz
# Mung bean
wget http://legumeinfo.org/data/public/Vigna_radiata/VC1973A.gnm6/vigra.VC1973A.gnm6.3nL8.pchr_plus_unassigned.fna.gz
# lotus
wget ftp://ftp.kazusa.or.jp/pub/lotus/lotus_r3.0/Lj3.0_pseudomol.fna.gz
# lupin
wget http://legumeinfo.org/data/public/Lupinus_angustifolius/Tanjil.gnm1/lupan.Tanjil.gnm1.Qq0N.scaf.fna.gz
# cowpea
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/687/525/GCA_001687525.1_Cowpea_0.03/GCA_001687525.1_Cowpea_0.03_genomic.fna.gz
# Broad bean
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/375/635/GCA_001375635.1_VfEP_Reference-Unigene/GCA_001375635.1_VfEP_Reference-Unigene_genomic.fna.gz
# Common Bean
wget http://legumeinfo.org/data/public/Phaseolus_vulgaris/G19833.gnm1/phavu.G19833.gnm1.zBnF.fa.gz
# Pigeon pea
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/340/665/GCA_000340665.1_C.cajan_V1.0/GCA_000340665.1_C.cajan_V1.0_genomic.fna.gz
# outgroup, tale cress
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
```

Next we will rename them to identify them easily.

```
mv GCA_000340665.1_C.cajan_V1.0_genomic.fna Cajanus_cajan.fasta
mv GCA_000722935.2_W05v1.0_genomic.fna Glycine_soja.fasta
mv Lj3.0_pseudomol.fna Lotus_japonicus.fasta
mv lupan.Tanjil.gnm1.Qq0N.scaf.fna Lupinus_angustifolius.fasta
mv tripr.MilvusB.gnm2.1.gNmT.pchr_plus_unanchored.fna Trifolium_pratense.fasta
mv GCA_001687525.1_Cowpea_0.03_genomic.fna Vicia_faba.fasta
mv vigan.Gyeongwon.gnm3.JyYC.pchr_plus_unanchored.fna Vigna_angularis.fasta
mv vigra.VC1973A.gnm6.3nL8.pchr_plus_unassigned.fna Vigna_radiata.fasta
mv GCA_001375635.1_VfEP_Reference-Unigene_genomic.fna Vigna_unguiculata.fasta
mv Arabidopsis_thaliana_TAIR10.fa Arabidopsis_thaliana.fasta
mv Aradu_v1.0_20140908.fa Arachis_duranensis.fasta
mv Araip_v1.0_20140908.fa Arachis_ipaensis.fasta
mv Cicer_arietinum_GA_v1.0_kabuli.chr_scf.fa Cicer_arietinum_kabuli.fasta
mv Medicago_truncatula.MedtrA17_4.0.dna.toplevel.fa Medicago_truncatula.fasta
mv Pvulgaris_442_v2.0.fa Phaseolus_vulgaris_jgi.fasta
mv phavu.G19833.gnm1.zBnF.fa Phaseolus_vulgaris.fasta
mv Trifolium_pratense.Trpr.dna.toplevel.fa Trifolium_pratense_ensembl.fasta
mv C_arietinum_ICC4958_v2_0.fasta Cicer_arietinum_ICC4958.fasta
mv Gmax_275_v2.0.softmasked.20140305.fasta Glycine_max.fasta
mv TSUd_r1.1.fasta Trifolium_subterraneum_r1.1.fasta
mv TSUw_r1.0.fasta Trifolium_subterraneum_r1.0.fasta
```

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
