---
title: DeepGoPlus | Using AI to associate GO terms with novel proteins 
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

DeepGOPlus is a bioinformatics tool designed for protein function prediction using deep learning. It combines convolutional neural networks with sequence-based features to predict Gene Ontology (GO) terms for proteins. This tool is particularly useful for large-scale genome annotation and proteome analysis, though a [webserver](https://deepgo.cbrc.kaust.edu.sa/) is currently available for smaller requests. Frequently, genome annotation projects will produce tens of thousands of genes, and depending upon the model status of your chosen species, many of these genes may lack significant homology to anything available in current databases. This software is not a replacement for identifying ontology terms using interproscan, but an alternative for novel gene candidates. 


**Software used in this tutorial**
* Kulmanov, M., & Hoehndorf, R. (2020). DeepGOPlus: Improved protein function prediction from sequence. *Bioinformatics*, 36(2), 422–429. https://doi.org/10.1093/bioinformatics/btz595
* Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. *Nature Methods*, 12(1), 59–60. https://doi.org/10.1038/nmeth.3176
* [Github Repository](https://github.com/bio-ontology-research-group/deepgoplus)


## Installation  of DeepGOPlus and prerequisite software

#### Create conda environment in python 3.7.9

* create conda environment
* activate conda environment
* create python virtual environment
* upgrade pip and setuptools

```bash
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
```



#### Download the repository 

* Download the repository
* install the requirements

```bash
git clone https://github.com/bio-ontology-research-group/deepgoplus.git
cd deepgoplus/
pip install -r requirements.txt
pip install deepgoplus
```


#### Download the metadata and training datasets

* download datasets
* copy metadata to data folder

```bash
wget http://deepgoplus.bio2vec.net/data/data.tar.gz

#automatically unpacks into data folder if within deepgoplus
tar -zxvf data.tar.gz

cp -rf metadata/ data/.
```

#### Install Diamond BLAST

* download tar ball
* add to PATH variable

```bash
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

#add to .bashrc and source or run in terminal for temporary access
export PATH="/work/gif3/masonbrink/Software/:$PATH"
```

## Run the database update.py script

In order for this to run properly, I had to modify the update.py script to get around the check for the uniprot header. This header is not required and is not found in the dataset we downloaded.


* Release information
  * Current version is 1.0.21
    * Gene Ontology released on 2024-09-08
    * SwissProt data with version 2024_05.



The following lines need to be modified in the `update.py` script. 

```python
#allow check for header to pass
if last_release_version == new_release_version:
```

to 

```python
#allow check for header to pass
if last_release_version != new_release_version:
```

And Updated uniprot address

```python
response = requests.get('https://www.uniprot.org')
```

to

```python
response = requests.get('https://rest.uniprot.org/uniprotkb/search')
```

Once these have been updated we can run 

```python
python update.py
```




## Run DeepGoPlus analysis

* Download some fun data from a plant parasitic nematode
* Format the data
* Run DeepGoPlus
* Format the output

```

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ditylenchus_dipsaci/PRJNA498219/ditylenchus_dipsaci.PRJNA498219.WBPS19.protein.fa.gz

# this will save formatting steps later, reducing the fasta header to just the basic name (first column).
less ditylenchus_dipsaci.PRJNA498219.WBPS19.protein.fa |awk '{print $1}' >Cleanditylenchus_dipsaci.PRJNA498219.WBPS19.protein.fa

# use the same threshold as the online database
ml miniconda3; conda activate deepgoplus_env; source bin/activate ; deepgoplus -if Cleanditylenchus_dipsaci.PRJNA498219.WBPS19.protein.fa -dr data/ -of results_DDipsaci -t 0.3

# format the output
awk 'NF>3' results_DDipsaci |sed 's/\t/#/1' |sed 's/\t/,/g' |sed 's/#/\t/g' >DeepGoPlusResults.tab
```



[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)