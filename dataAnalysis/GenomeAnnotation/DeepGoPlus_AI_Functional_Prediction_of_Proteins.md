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

Once this has been updated we can run 

```python
python update.py
```


## Run DeepGoPlus analysis

* Download some fun data from a plant parasitic nematode (Ditylenchus dipsaci)
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


| Gene ID   | GO Terms                                                                                           |
|-----------|----------------------------------------------------------------------------------------------------|
| jg20583   | GO:0110165-0.355, GO:0005575-0.361                                                                 |
| jg20588   | GO:0110165-0.349, GO:0005575-0.357                                                                 |
| jg16733   | GO:0003674-0.567, GO:0005488-0.500, GO:0005515-0.377, GO:0005575-0.725, GO:0005622-0.613, GO:0005737-0.575, GO:0005886-0.331, GO:0007154-0.330, GO:0008150-0.727, GO:0009987-0.583, GO:0016020-0.499, GO:0023052-0.332, GO:0043226-0.400, GO:0043227-0.329, GO:0043229-0.359, GO:0043231-0.311, GO:0050789-0.481, GO:0050794-0.450, GO:0050896-0.393, GO:0051716-0.344, GO:0065007-0.494, GO:0071944-0.352, GO:0110165-0.720 |
| jg20471   | GO:0110165-0.387, GO:0005575-0.396                                                                 |
| jg2136    | GO:0003674-0.503, GO:0005488-0.369, GO:0005515-0.349, GO:0005575-0.788, GO:0005622-0.688, GO:0005737-0.635, GO:0005886-0.434, GO:0008104-0.443, GO:0008150-0.742, GO:0009987-0.687, GO:0016020-0.654, GO:0016043-0.382, GO:0032991-0.303, GO:0033036-0.444, GO:0033365-0.308, GO:0043226-0.566, GO:0043227-0.528, GO:0043229-0.555, GO:0043231-0.503, GO:0045184-0.435, GO:0051179-0.510, GO:0051234-0.500, GO:0051641-0.465, GO:0051668-0.362, GO:0070727-0.443, GO:0071840-0.387, GO:0071944-0.435, GO:0072657-0.359, GO:0090150-0.353, GO:0110165-0.783 |
| jg12535   | GO:0110165-0.345, GO:0005575-0.345                                                                 |
| jg2135    | GO:0110165-0.395, GO:0005575-0.397, GO:0005622-0.365, GO:0016020-0.309, GO:0043229-0.329, GO:0043226-0.337 |
| jg2134    | GO:0110165-0.369, GO:0005575-0.370, GO:0005622-0.343                                               |
| jg12533   | GO:0110165-0.360, GO:0005575-0.364, GO:0016020-0.329                                               |
| jg14715   | GO:0110165-0.370, GO:0005575-0.383, GO:0005622-0.305                                               |
| jg21581   | GO:0110165-0.380, GO:0005575-0.380                                                                 |
| jg2133    | GO:0003674-0.732, GO:0005215-0.637, GO:0005575-0.637, GO:0005622-0.327, GO:0005737-0.327, GO:0005886-0.602, GO:0006810-0.737, GO:0006811-0.616, GO:0006812-0.598, GO:0008150-0.948, GO:0008324-0.354, GO:0009987-0.771, GO:0015075-0.362, GO:0015095-0.319, GO:0015318-0.575, GO:0015693-0.555, GO:0016020-0.603, GO:0022857-0.637, GO:0022890-0.566, GO:0030001-0.589, GO:0034220-0.616, GO:0046873-0.340, GO:0050896-0.334, GO:0051179-0.737, GO:0051234-0.737, GO:0051716-0.307, GO:0055085-0.679, GO:0071944-0.602, GO:0098655-0.596, GO:0098660-0.607, GO:0098662-0.595, GO:0110165-0.637, GO:1903830-0.445 |
| jg12534   | GO:0110165-0.394, GO:0005575-0.396, GO:0005622-0.314                                               |


[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)