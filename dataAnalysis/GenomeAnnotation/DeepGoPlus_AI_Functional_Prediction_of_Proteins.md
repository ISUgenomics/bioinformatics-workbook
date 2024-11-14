---
title: DeepGoPlus: Using AI to associate GO terms with novel proteins 
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

DeepGOPlus is a bioinformatics tool designed for protein function prediction using deep learning. It combines convolutional neural networks with sequence-based features to predict Gene Ontology (GO) terms for proteins. This tool is particularly useful for large-scale genome annotation and proteome analysis, though a webserver is currently availble for smaller requests: https://deepgo.cbrc.kaust.edu.sa/ . Frequently genome annotation projects will produce tens of thousands of genes, and depending upon the model status of your chosen species, many of these genes may lack significant homology to anything availabe in current databases. This software is not a replacement for identifying ontology terms using interproscan, but an alternative for novel gene candidates. 


**Software used in this tutorial**
* Kulmanov, M., & Hoehndorf, R. (2020). DeepGOPlus: Improved protein function prediction from sequence. *Bioinformatics*, 36(2), 422–429. https://doi.org/10.1093/bioinformatics/btz595
* Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. *Nature Methods*, 12(1), 59–60. https://doi.org/10.1038/nmeth.3176




### Installation

**Installation of DeepGOPlus**
```
#create conda environment in python 3.7.9
ml miniconda3
conda create -n deepgoplus_env python=3.7.9

#if you havent done this previously you'll need to run this to use conda
conda init
source ~/.bashrc

#activate the environment
conda activate deepgoplus_env

#create and actiavte a python virtual environment
python3 -m venv DeepGoPlus_pyenv
source DeepGoPlus_pyenv/bin/activate

#upgrade pip and setuptools
pip install --upgrade pip setuptools

#download the repository 
git clone https://github.com/bio-ontology-research-group/deepgoplus.git
cd deepgoplus/
pip install -r requirements.txt
pip install deepgoplus

#download the metadata and training
wget http://deepgoplus.bio2vec.net/data/data.tar.gz

#automatically unpacks into data folder if within deepgoplus
tar -zxvf data.tar.gz


cp -rf metadata/ data/.
```

**Install Diamond BLAST to path**
```
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

#add to .bashrc and source or run in terminal for temporary access
export PATH="/work/gif3/masonbrink/Software/:$PATH"
```

**I had to modify the update.py script to get around the check for the uniprot header. This header must not exist anymore.** 
```
# Release information
Current version is 1.0.21. The model in the current release was trained using the Gene Ontology
released on 2024-09-08 and the SwissProt data with version 2024_05.

python update.py

#allow check for header to pass
if last_release_version == new_release_version:
#to 
if last_release_version != new_release_version:

#Updated uniprot address
#this line
response = requests.get('https://www.uniprot.org')
#to
response = requests.get('https://rest.uniprot.org/uniprotkb/search')

# this allows the update and retraining of the databases 
python update.py 
```


### Proceed with DeepGoPlus Run 
```
#get some fun data from a plant parasitic nematode
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ditylenchus_dipsaci/PRJNA498219/ditylenchus_dipsaci.PRJNA498219.WBPS19.protein.fa.gz

# this will save formatting steps later, reducing the fasta header to just the basic name (first column).
less ditylenchus_dipsaci.PRJNA498219.WBPS19.protein.fa |awk '{print $1}' >Cleanditylenchus_dipsaci.PRJNA498219.WBPS19.protein.fa

# use the same threshold as the online database
ml miniconda3; conda activate deepgoplus_env; source bin/activate ; deepgoplus -if Cleanditylenchus_dipsaci.PRJNA498219.WBPS19.protein.fa -dr data/ -of results_DDipsaci -t 0.3


awk 'NF>3' results_DDipsaci |sed 's/\t/#/1' |sed 's/\t/,/g' |sed 's/#/\t/g' >DeepGoPlusResults.tab
```



[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)