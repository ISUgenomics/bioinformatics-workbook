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
|Gene_name       GoTerms|
|:----|
|jg20583 GO:0110165#0.355        GO:0005575#0.361|
|jg20588 GO:0110165#0.349        GO:0005575#0.357|
|jg16733 GO:0003674#0.567        GO:0005488#0.500        GO:0005515#0.377        GO:0005575#0.725        GO:0005622#0.613        GO:0005737#0.575        GO:0005886#0.331        GO:0007154#0.330        GO:0008150#0.727 GO:0009987#0.583        GO:0016020#0.499        GO:0023052#0.332        GO:0043226#0.400        GO:0043227#0.329        GO:0043229#0.359        GO:0043231#0.311        GO:0050789#0.481        GO:0050794#0.450 GO:0050896#0.393        GO:0051716#0.344        GO:0065007#0.494        GO:0071944#0.352        GO:0110165#0.720|
|jg20471 GO:0110165#0.387        GO:0005575#0.396|
|jg2136  GO:0003674#0.503        GO:0005488#0.369        GO:0005515#0.349        GO:0005575#0.788        GO:0005622#0.688        GO:0005737#0.635        GO:0005886#0.434        GO:0008104#0.443        GO:0008150#0.742 GO:0009987#0.687        GO:0016020#0.654        GO:0016043#0.382        GO:0032991#0.303        GO:0033036#0.444        GO:0033365#0.308        GO:0043226#0.566        GO:0043227#0.528        GO:0043229#0.555 GO:0043231#0.503        GO:0045184#0.435        GO:0051179#0.510        GO:0051234#0.500        GO:0051641#0.465        GO:0051668#0.362        GO:0070727#0.443        GO:0071840#0.387        GO:0071944#0.435 GO:0072657#0.359        GO:0090150#0.353        GO:0110165#0.783|
|jg12535 GO:0110165#0.345        GO:0005575#0.345|
|jg2135  GO:0110165#0.395        GO:0005575#0.397        GO:0005622#0.365        GO:0016020#0.309        GO:0043229#0.329        GO:0043226#0.337|
|jg2134  GO:0110165#0.369        GO:0005575#0.370        GO:0005622#0.343|
|jg12533 GO:0110165#0.360        GO:0005575#0.364        GO:0016020#0.329|
|jg14715 GO:0110165#0.370        GO:0005575#0.383        GO:0005622#0.305|
|jg21581 GO:0110165#0.380        GO:0005575#0.380|
|jg2133  GO:0003674#0.732        GO:0005215#0.637        GO:0005575#0.637        GO:0005622#0.327        GO:0005737#0.327        GO:0005886#0.602        GO:0006810#0.737        GO:0006811#0.616        GO:0006812#0.598 GO:0008150#0.948        GO:0008324#0.354        GO:0009987#0.771        GO:0015075#0.362        GO:0015095#0.319        GO:0015318#0.575        GO:0015693#0.555        GO:0016020#0.603        GO:0022857#0.637 GO:0022890#0.566        GO:0030001#0.589        GO:0034220#0.616        GO:0046873#0.340        GO:0050896#0.334        GO:0051179#0.737        GO:0051234#0.737        GO:0051716#0.307        GO:0055085#0.679 GO:0071944#0.602        GO:0098655#0.596        GO:0098660#0.607        GO:0098662#0.595        GO:0110165#0.637        GO:1903830#0.445|
|jg12534 GO:0110165#0.394        GO:0005575#0.396        GO:0005622#0.314|
|jg24909 GO:0000902#0.352        GO:0003674#0.510        GO:0005488#0.356        GO:0005515#0.326        GO:0005575#0.791        GO:0005622#0.389        GO:0005886#0.490        GO:0007154#0.433        GO:0007165#0.395 GO:0007166#0.343        GO:0007275#0.622        GO:0007399#0.531        GO:0007409#0.345        GO:0007411#0.325        GO:0008150#0.946        GO:0009653#0.424        GO:0009987#0.768        GO:0016020#0.579 GO:0016043#0.472        GO:0022008#0.405        GO:0023052#0.432        GO:0030030#0.360        GO:0030154#0.441        GO:0030182#0.384        GO:0031175#0.360        GO:0032501#0.665        GO:0032502#0.645 GO:0048468#0.419        GO:0048518#0.435        GO:0048522#0.377        GO:0048583#0.318        GO:0048666#0.372        GO:0048667#0.347        GO:0048699#0.391        GO:0048731#0.606        GO:0048812#0.348 GO:0048856#0.645        GO:0048858#0.348        GO:0048869#0.441        GO:0050789#0.629        GO:0050794#0.563        GO:0050896#0.579        GO:0051716#0.428        GO:0061564#0.345        GO:0065007#0.646 GO:0071840#0.472        GO:0071944#0.562        GO:0097485#0.341        GO:0110165#0.791        GO:0120036#0.360        GO:0120039#0.348|
|jg14716 GO:0110165#0.373        GO:0005575#0.379|
|jg21582 GO:0110165#0.382        GO:0005575#0.389        GO:0005622#0.307|
|jg2132  GO:0003674#0.662        GO:0005215#0.547        GO:0005575#0.643        GO:0005886#0.483        GO:0006810#0.656        GO:0006811#0.580        GO:0006812#0.566        GO:0008150#0.934        GO:0008324#0.327 GO:0009987#0.734        GO:0015075#0.328        GO:0015095#0.319        GO:0015318#0.541        GO:0015693#0.550        GO:0016020#0.534        GO:0022857#0.547        GO:0022890#0.539        GO:0030001#0.561 GO:0034220#0.575        GO:0046873#0.319        GO:0050896#0.411        GO:0051179#0.661        GO:0051234#0.656        GO:0051716#0.341        GO:0055085#0.587        GO:0071944#0.490        GO:0098655#0.565 GO:0098660#0.568        GO:0098662#0.565        GO:0110165#0.640        GO:1903830#0.448|
|jg12531 GO:0110165#0.364        GO:0005575#0.374|
|jg14717 GO:0110165#0.370        GO:0005575#0.371        GO:0005622#0.305|
|jg2131  GO:0003674#0.693        GO:0003824#0.533        GO:0004596#0.302        GO:0005575#0.794        GO:0005622#0.767        GO:0005634#0.345        GO:0005737#0.700        GO:0005829#0.433        GO:0008080#0.452 GO:0008150#0.728        GO:0008152#0.422        GO:0009058#0.373        GO:0009059#0.340        GO:0009987#0.509        GO:0010467#0.331        GO:0016020#0.550        GO:0016407#0.452        GO:0016410#0.452 GO:0016740#0.486        GO:0016746#0.464        GO:0016747#0.456        GO:0019538#0.347        GO:0031248#0.356        GO:0031414#0.356        GO:0032991#0.465        GO:0034212#0.302        GO:0043170#0.374 GO:0043226#0.548        GO:0043227#0.501        GO:0043229#0.538        GO:0043231#0.478        GO:0044237#0.378        GO:0044238#0.405        GO:0044249#0.349        GO:0110165#0.791        GO:0140535#0.376 GO:1902493#0.356        GO:1902494#0.387        GO:1990234#0.366|
|jg12532 GO:0003674#0.493        GO:0005488#0.322        GO:0005575#0.874        GO:0005622#0.810        GO:0005737#0.716        GO:0005773#0.331        GO:0005886#0.503        GO:0006810#0.393        GO:0006811#0.343 GO:0008150#0.874        GO:0009987#0.794        GO:0016020#0.788        GO:0030154#0.353        GO:0032502#0.434        GO:0034220#0.314        GO:0043226#0.772        GO:0043227#0.719        GO:0043229#0.728 GO:0043231#0.677        GO:0048856#0.412        GO:0048869#0.353        GO:0050789#0.322        GO:0050896#0.322        GO:0051179#0.440        GO:0051234#0.397        GO:0055085#0.316        GO:0065007#0.537 GO:0071944#0.534        GO:0098590#0.305        GO:0110165#0.872|
|jg14718 GO:0110165#0.401        GO:0005575#0.405        GO:0005622#0.352        GO:0016020#0.338        GO:0043226#0.308|
|jg21580 GO:0110165#0.351        GO:0005575#0.367|
|jg2130  GO:0003674#0.357        GO:0005575#0.859        GO:0005622#0.811        GO:0005737#0.680        GO:0005773#0.408        GO:0005794#0.341        GO:0005886#0.440        GO:0008150#0.580        GO:0009987#0.440 GO:0012505#0.384        GO:0016020#0.702        GO:0032991#0.306        GO:0043226#0.693        GO:0043227#0.667        GO:0043229#0.669        GO:0043231#0.652        GO:0051179#0.309        GO:0071944#0.457 GO:0110165#0.851|
|jg20474 GO:0110165#0.349        GO:0005575#0.353|
|jg14711 GO:0003674#0.429        GO:0005488#0.396        GO:0005515#0.376        GO:0005575#0.836        GO:0005622#0.775        GO:0005737#0.724        GO:0005886#0.353        GO:0006810#0.537        GO:0008104#0.304 GO:0008150#0.843        GO:0009987#0.787        GO:0016020#0.694        GO:0016192#0.307        GO:0031090#0.352        GO:0032501#0.486        GO:0032502#0.308        GO:0032991#0.384        GO:0033036#0.307 GO:0043226#0.649        GO:0043227#0.625        GO:0043229#0.623        GO:0043231#0.597        GO:0046907#0.423        GO:0048856#0.308        GO:0051179#0.648        GO:0051234#0.642        GO:0051641#0.557 GO:0051649#0.438        GO:0070727#0.305        GO:0071944#0.370        GO:0110165#0.830|
|jg12530 GO:0000785#0.355        GO:0000976#0.566        GO:0000977#0.496        GO:0000978#0.441        GO:0000981#0.346        GO:0000987#0.449        GO:0001067#0.566        GO:0003674#0.785        GO:0003676#0.663 GO:0003677#0.625        GO:0003690#0.615        GO:0003700#0.453        GO:0005488#0.785        GO:0005515#0.674        GO:0005575#0.820        GO:0005622#0.682        GO:0005634#0.515        GO:0005667#0.380 GO:0005694#0.356        GO:0005737#0.324        GO:0006139#0.749        GO:0006351#0.725        GO:0006355#0.720        GO:0006357#0.494        GO:0006366#0.502        GO:0007275#0.460        GO:0007399#0.345 GO:0008134#0.430        GO:0008150#0.945        GO:0008152#0.775        GO:0009058#0.770        GO:0009059#0.770        GO:0009889#0.754        GO:0009891#0.666        GO:0009893#0.679        GO:0009987#0.899 GO:0010467#0.770        GO:0010468#0.748        GO:0010556#0.754        GO:0010557#0.654        GO:0010604#0.674        GO:0016020#0.643        GO:0016043#0.347        GO:0016070#0.749        GO:0019219#0.730 GO:0019222#0.775        GO:0022008#0.336        GO:0030154#0.428        GO:0030182#0.329        GO:0031323#0.770        GO:0031325#0.676        GO:0031326#0.754        GO:0031328#0.662        GO:0032501#0.485 GO:0032502#0.533        GO:0032774#0.749        GO:0032991#0.464        GO:0034654#0.749        GO:0043170#0.770        GO:0043226#0.652        GO:0043227#0.565        GO:0043228#0.450        GO:0043229#0.646 GO:0043231#0.550        GO:0043232#0.448        GO:0043425#0.330        GO:0043565#0.583        GO:0044237#0.770        GO:0044238#0.765        GO:0044249#0.770        GO:0045893#0.640        GO:0045935#0.653 GO:0045944#0.473        GO:0048468#0.316        GO:0048518#0.719        GO:0048522#0.711        GO:0048699#0.331        GO:0048731#0.440        GO:0048856#0.533        GO:0048869#0.430        GO:0050789#0.871 GO:0050794#0.859        GO:0050896#0.346        GO:0051252#0.728        GO:0051254#0.645        GO:0060255#0.763        GO:0065007#0.872        GO:0071840#0.352        GO:0080090#0.765        GO:0090304#0.749 GO:0097159#0.663        GO:0110165#0.738        GO:0140110#0.467        GO:0140297#0.387        GO:0141187#0.749        GO:1902680#0.640        GO:1990837#0.574        GO:2001141#0.721|
|jg14712 GO:0003674#0.487        GO:0005488#0.446        GO:0005575#0.906        GO:0005622#0.864        GO:0005737#0.808        GO:0005739#0.572        GO:0006091#0.344        GO:0006119#0.340        GO:0006120#0.340 GO:0006139#0.390        GO:0006163#0.343        GO:0006164#0.343        GO:0006753#0.345        GO:0006754#0.340        GO:0006793#0.358        GO:0006796#0.358        GO:0008150#0.714        GO:0008152#0.557 GO:0009058#0.540        GO:0009060#0.343        GO:0009117#0.345        GO:0009141#0.340        GO:0009142#0.340        GO:0009144#0.340        GO:0009145#0.340        GO:0009150#0.340        GO:0009152#0.340 GO:0009165#0.344        GO:0009199#0.340        GO:0009201#0.340        GO:0009205#0.340        GO:0009206#0.340        GO:0009259#0.340        GO:0009260#0.340        GO:0009987#0.658        GO:0015980#0.344 GO:0015986#0.340        GO:0016020#0.822        GO:0019637#0.351        GO:0019646#0.340        GO:0019693#0.340        GO:0022900#0.344        GO:0022904#0.340        GO:0031974#0.306        GO:0032991#0.390 GO:0034654#0.382        GO:0042773#0.340        GO:0042775#0.340        GO:0043226#0.809        GO:0043227#0.785        GO:0043229#0.805        GO:0043231#0.780        GO:0043233#0.306        GO:0044237#0.543 GO:0044238#0.542        GO:0044249#0.404        GO:0044281#0.356        GO:0045333#0.343        GO:0046034#0.340        GO:0046390#0.340        GO:0055086#0.346        GO:0070013#0.304        GO:0072521#0.343 GO:0072522#0.343        GO:0090407#0.349        GO:0098796#0.316        GO:0110165#0.903        GO:1901135#0.356        GO:1901137#0.348        GO:1901293#0.344        GO:1902494#0.321|

[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)