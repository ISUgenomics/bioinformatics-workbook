---
title: Identify secreted proteins and predict subcellular localization
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# How to identify the secreted protein from an annotation and predict subcellular localization

Here we will be using a set of predicted proteins from a plant parasitic nematode genome to predict secretion, transmembrane domains, and subcellular localization.

**Software used in this tutorial**
1.  Signalp 6.0
2.  Tmhmm 2.0c
3.  Samtools 1.16.1
4.  Localizer 1.0.5
5.  DeepLoc 2.0



### Installation of SignalP 6.0
```
# create and activate a python virtual environment
python -m venv signalp6_env
source signalp6_env/bin/activate

# register for signalp use here to obtain the software
https://services.healthtech.dtu.dk/services/SignalP-6.0/

#install SignalP
tar -zxvf signalp-6.0h.fast.tar.gz
cd signalp6_fast/
pip install signalp-6-package/

#it installs with numpy 2, which is incorrect, so had to run this. 
pip install numpy==1.24.3

#This code stores the installation path of the signalp module in the variable SIGNALP_DIR, so that you can use this path in further shell commands.
SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
```

### Run Signalp 6.0
```
#softlink your proteins
ln -s /work/gif3/masonbrink/Baum/01_ReannotateAllSCNGenomes/42_Fix67CDSs/TN10FinalManualAnnotation_proteins.fasta

#run signalp 6 in a compute node
echo "source signalp6_env/bin/activate; signalp6 --fastafile TN10FinalManualAnnotation_proteins.fasta --organism euk --output_dir Signalp6_out  --format txt --mode fast --model_dir signalp6_fast/signalp-6-package/models/" >signalp.sh
```

Excerpt of results from Signalp 6.0. We want to get the fasta names that are secreted and the last position of the cleavage site.  
```
# SignalP-6.0   Organism: Eukarya       Timestamp: 20240910112251
# ID    Prediction      OTHER   SP(Sec/SPI)     CS Position
Hg_chrom1_TN10mRNA_1 gene=Hg_chrom1_TN10gene_1  OTHER   1.000000        0.000000
Hg_chrom1_TN10mRNA_2 gene=Hg_chrom1_TN10gene_2  OTHER   1.000000        0.000001
Hg_chrom1_TN10mRNA_3 gene=Hg_chrom1_TN10gene_3  OTHER   1.000000        0.000000
Hg_chrom1_TN10mRNA_4 gene=Hg_chrom1_TN10gene_4  OTHER   1.000000        0.000003
Hg_chrom1_TN10mRNA_5 gene=Hg_chrom1_TN10gene_5  OTHER   1.000000        0.000023
Hg_chrom1_TN10mRNA_6 gene=Hg_chrom1_TN10gene_6  OTHER   1.000000        0.000006
Hg_chrom1_TN10mRNA_7 gene=Hg_chrom1_TN10gene_7  OTHER   1.000000        0.000000
Hg_chrom1_TN10mRNA_8 gene=Hg_chrom1_TN10gene_8  OTHER   0.999587        0.000437
Hg_chrom1_TN10mRNA_9 gene=Hg_chrom1_TN10gene_9  OTHER   1.000000        0.000006
Hg_chrom1_TN10mRNA_10 gene=Hg_chrom1_TN10gene_10        OTHER   1.000000        0.000000
Hg_chrom1_TN10mRNA_11 gene=Hg_chrom1_TN10gene_11        SP      0.000243        0.999714        CS pos: 19-20. Pr: 0.9780
Hg_chrom1_TN10mRNA_12 gene=Hg_chrom1_TN10gene_12        OTHER   1.000000        0.000000
Hg_chrom1_TN10mRNA_13 gene=Hg_chrom1_TN10gene_13        SP      0.000200        0.999755        CS pos: 19-20. Pr: 0.9861
Hg_chrom1_TN10mRNA_14 gene=Hg_chrom1_TN10gene_14        OTHER   0.976381        0.023633
Hg_chrom1_TN10mRNA_15 gene=Hg_chrom1_TN10gene_15        OTHER   1.000000        0.000039
Hg_chrom1_TN10mRNA_16 gene=Hg_chrom1_TN10gene_16        OTHER   1.000000        0.000000
Hg_chrom1_TN10mRNA_17 gene=Hg_chrom1_TN10gene_17        OTHER   1.000000        0.000000
Hg_chrom1_TN10mRNA_18 gene=Hg_chrom1_TN10gene_18        OTHER   1.000000        0.000000
Hg_chrom1_TN10mRNA_19 gene=Hg_chrom1_TN10gene_19        OTHER   0.844543        0.155449
Hg_chrom1_TN10mRNA_20 gene=Hg_chrom1_TN10gene_20        OTHER   1.000000        0.000005
Hg_chrom1_TN10mRNA_21 gene=Hg_chrom1_TN10gene_21        SP      0.000452        0.999497        CS pos: 25-26. Pr: 0.9702
continued ...
```


### Identify secreted proteins in signalp6 results and remove signal peptide for transmembrane domain prediction
```
# create a faidx index of our proteins to subtract the signal peptide from the fasta of secreted proteins
ml samtools;samtools faidx TN10FinalManualAnnotation_proteins.fasta

# this gets the signalp results file and the protein lengths in the same order 
awk '$2=="SP"{print $1"\t"$7}'  Signalp6_out/prediction_results.txt|sed -e 's/\.//g' -e 's/-/\t/g' |awk -F"\t" '{print "samtools faidx TN10FinalManualAnnotation_proteins.fasta "$1":"$3"- >>SignalPeptidesSubtracted6.fasta"}' >ExtractSignalPContainingProts6.sh

sh ExtractSignalPContainingProts6.sh
```

**Excerpt from SignalPeptidesSubtracted6.fasta**
```
>Hg_chrom1_TN10mRNA_11:20-
HDQQLRSENLGDVAAKDHRATGGDTLPKYSKDQLHQRAQKIGEAKTAEHQQQNFQNAVAD
QQQHQAAAPSAFKVENAVGHDNQEQQQKAAEHLENKFVASGTKSKLMNTDARKNVHEQQQ
QQAVGQQAVQQRTDELRQHLEQPSTSAIEGDQQRHARVYVVKQVAAEPIFNSNAEQPEAE
EFNAGRKLNKDMGNELLRPVYLNDNNRFATAARIVSKVFGGHYYEQQQQQRQFQYQRGYG
SNKMGQPEQGEQFENAETRQQRVQFSGKNGQTELVLDNLNQHQQQQYNNNMGKDSQSFKS
ERLDNSRGMQRVQEHQQQQQNKLGQDSRVILGDLQQQKLEQDNARTFQGEQQQQQKQKLE
QGNSRTFQGEQQQKLEQDSMRTFQGEQQQKLEQDSMRTFQGEQQHQQKEQHEIQQFEQQE
NKQQEPSLALERARLEELRIKERARVAKLEAERQARLVEIQARERQAKLEAEMRVRNVEN
RVRQAKLEAEARANQAILEAKVRQAERRARLGHGGGFADASEELQKDEQQAQQRQQQQLL
EQPQINTPQQLRENVAQPIVAAPVVAVAGVH.
>Hg_chrom1_TN10mRNA_13:20-
AHQQQLGTPVIGGICKLDTPDVHIGGKQTQFFLRCEPNTDSAQGEGVWVVKSRHAGSTAA
KIPAAAVVSQQEIEQGEQNSQKIAATKQISSYVCDLVSDAVEHGYCSTSENCLQPVYTDR
TAFLQCDPTSRRWTKKHCQSGFNFDFERQACIAHTAKMSHHFRQFPTQVLPHSQRGGGIV
CTFAQCSVADPCNVGSCNNGYCCTAAVGSMPKAVLPISTNKTISSGNSNNIAKPLSPPSS
TIAVSASETVPPIPSILVIELNNNNNHNDDVIEYDQLRLNDIPSTLIGHISSNSQPMPTI
FDHCSSGFHSPIRCGGTSGGECPVGLICELATRFCCPFPPKKLKRKSTSTNSGGIGGVVG
TTATGNWAQSLVPNGAKVVRRRTALRRHGRFERSFFSSFCCVPAAEQLLIDANDNETPEQ
QLNMFICPTGAPALGPCSPGYCVAPLQCAGALCCAPPAPVMAPPPPIAYACPGGLPALGP
CIGGRCAAPAICASPSNVCCAQPSVPTAIPATVCPDGTQAAGACVNGQCGAGFTCNQGLC
CTNSSQTPRCLDGSQSIGACIQGRCGTGYTCTTGNICCPSQLNLCPPGQTSVGQAVNGRC
PAGYTNFNGQCCGPQSAQSQVSCSVEDSFGRCDSNQQCAEPGYACDVANNWCCPQVIGDA
IGPCIQGEGGNRLCPEGYACVGEGEGQCFRLDTGTCAPEEQSGPCAPDGTCPPGYECIEG
FCCQTDAGGANPAAAALRLFHQRRRRRRASASSSAFLKMINYDRGSRQPTQNRNSKK.
>Hg_chrom1_TN10mRNA_21:26-
TISTAFSQFLRNYFNDEVEKNIARRDLGTDGSFGGGKGRIKAKKRPVIFVHGLTNLAGEL
DYVRRLFREKGGYRDSELYATTYGYGMKGWKWLRDAMKCDHVKMIRIMIRAITQFTHSSQ
AILGGVCVDTEEQLGPPLTHLVHTFIGVGGANREAVHLCRLFEWAMPCNPVNGMRCHSQF
LYDINSNVGYEATTRIFVIRSMDDGTVGTRDCEGRSVSAIDGQNDEIVLRNYSHQMVIFG
TGEQQLKLLTF.
```

# Install and run Tmhmm to identify transmembrane domains
```
#Download and extract tmhmm from here. 
https://services.healthtech.dtu.dk/services/TMHMM-2.0/

add this to your ~/.bashrc 
#export PATH="/path/to/your/tmhmm-2.0c/bin/:$PATH"

#run tmhmm
tmhmm SignalPeptidesSubtracted6.fasta >SignalPeptidesSubtracted6.tmhmmout
```

**Excerpt of the results from Tmhmm**
```
# Hg_chrom1_TN10mRNA_11:20- Length: 571
# Hg_chrom1_TN10mRNA_11:20- Number of predicted TMHs:  0
# Hg_chrom1_TN10mRNA_11:20- Exp number of AAs in TMHs: 0.00222
# Hg_chrom1_TN10mRNA_11:20- Exp number, first 60 AAs:  0
# Hg_chrom1_TN10mRNA_11:20- Total prob of N-in:        0.00112
Hg_chrom1_TN10mRNA_11:20-       TMHMM2.0        outside      1   571
# Hg_chrom1_TN10mRNA_13:20- Length: 777
# Hg_chrom1_TN10mRNA_13:20- Number of predicted TMHs:  0
# Hg_chrom1_TN10mRNA_13:20- Exp number of AAs in TMHs: 0.01067
# Hg_chrom1_TN10mRNA_13:20- Exp number, first 60 AAs:  0.00038
# Hg_chrom1_TN10mRNA_13:20- Total prob of N-in:        0.00040
Hg_chrom1_TN10mRNA_13:20-       TMHMM2.0        outside      1   777
# Hg_chrom1_TN10mRNA_21:26- Length: 251
# Hg_chrom1_TN10mRNA_21:26- Number of predicted TMHs:  0
# Hg_chrom1_TN10mRNA_21:26- Exp number of AAs in TMHs: 0.02999
# Hg_chrom1_TN10mRNA_21:26- Exp number, first 60 AAs:  0.02192
# Hg_chrom1_TN10mRNA_21:26- Total prob of N-in:        0.01958
Hg_chrom1_TN10mRNA_21:26-       TMHMM2.0        outside      1   251
# Hg_chrom1_TN10mRNA_25:23- Length: 2168
# Hg_chrom1_TN10mRNA_25:23- Number of predicted TMHs:  0
# Hg_chrom1_TN10mRNA_25:23- Exp number of AAs in TMHs: 30.69038
# Hg_chrom1_TN10mRNA_25:23- Exp number, first 60 AAs:  0.01395
# Hg_chrom1_TN10mRNA_25:23- Total prob of N-in:        0.00072
continued ...
```

# Subcellular Localization
Here we are using two distinct subcellular localization predictors.  Each of these two software's uses a different approach to identifying the cellular compartment localization. 


| Feature                  | Localizer                                  | DeepLoc                                  |
|--------------------------|--------------------------------------------|------------------------------------------|
| **Focus**                 | Plant-specific, focused on chloroplast, mitochondrion, and nucleus. | General eukaryotic protein localization across 10 compartments. |
| **Subcellular compartments** | Chloroplast, mitochondrion, nucleus.     | Nucleus, cytoplasm, extracellular, mitochondrion, ER, Golgi, lysosome/vacuole, peroxisome, plasma membrane, chloroplast. |
| **Methodology**           | Rule-based predictions using known motifs/signals. | Deep learning model (RNN) based on sequence data. |
| **Organism specificity**  | Plant proteins.                            | General eukaryotes, including plants. |
| **Prediction accuracy**   | Focuses on high accuracy for three compartments. | Predicts across a broad range of compartments using neural networks. |
| **Ease of interpretation**| Results are interpretable based on known localization signals. | Predictions from a neural network may be harder to interpret. |
| **Input format**          | Protein sequences (FASTA).                 | Protein sequences (FASTA). |
| **Publication**           | Sperschneider et al. (2017).               | Almagro Armenteros et al. (2017). |


### Localizer
Localizer is a highly specialized software for plants and only predicts chlorplast, mitochondrial, and nuclear localization using known motifs/signals.

**Install Localizer**
```
git clone https://github.com/JanaSperschneider/LOCALIZER.git
tar -xvf LOCALIZER-1.0.5.tar.gz
cd LOCALIZER-1.0.5/Scripts


python -m venv localizer
source localizer/bin/activate

tar xvf emboss-latest.tar.gz
cd EMBOSS-6.5.7/
./configure
make
cd ../ 
unzip weka-3-6-12.zip

# added this to my ~/.bashrc
#export PATH="/work/gif3/masonbrink/Baum/01_ReannotateAllSCNGenomes/39_localizer/LOCALIZER/Scripts/:$PATH"
```
### Run Localizer 1.0.5
```
echo "source localizer/bin/activate; LOCALIZER.py -e -M -o SP6Out -i SignalPeptidesSubtracted6.fasta ">localizerSP6.sh
```

**Excerpt of results from Localizer**
```
# -----------------
# LOCALIZER 1.0.5 Predictions (-e mode)
# -----------------
Identifier                              Chloroplast             Mitochondria            Nucleus
Hg_chrom1_TN10mRNA_1003:25-587          Y (0.964 | 62-86)       -                       Y (RKRNRRLSASSVHRRRRHS,RRNCNEKLCNFTQETKRWR)
Hg_chrom1_TN10mRNA_1006:26-1143         -                       -                       Y (KWRRRNNGNKKGERSKKEKRRP)
Hg_chrom1_TN10mRNA_1028:26-246          -                       -                       Y (APAPHSQRKEKQGKKEKAKSSSSSSSAEEEEKKHRKKGNGRQIVKHFGKGKKKD)
Hg_chrom1_TN10mRNA_102:19-109           -                       -                       -
Hg_chrom1_TN10mRNA_1045:36-353          -                       -                       -
Hg_chrom1_TN10mRNA_1046:24-361          -                       -                       -
Hg_chrom1_TN10mRNA_105:30-452           -                       -                       -
Hg_chrom1_TN10mRNA_1069:30-355          -                       -                       -
Hg_chrom1_TN10mRNA_1090:22-475          -                       -                       Y (KKGKIGGFLIQKIQRQRRQ)
Hg_chrom1_TN10mRNA_1092:44-1224         -                       -                       -
Hg_chrom1_TN10mRNA_1097:25-161          -                       -                       -
Hg_chrom1_TN10mRNA_1108:25-397          -                       -                       -
Hg_chrom1_TN10mRNA_1115:30-4744         -                       -                       Y (PKPS)
Hg_chrom1_TN10mRNA_1116:23-494          -                       -                       -
Hg_chrom1_TN10mRNA_112:24-315           -                       -                       Y (RKIMQQTNAKKAFKK,KKGDKLKPKKPPAKNAPEAPVTGKSKAEEKGQEKEEKKGTTDKKEAGAEKKKKEL)
Hg_chrom1_TN10mRNA_113:26-479           -                       -                       -
Hg_chrom1_TN10mRNA_1150:30-127          Y (1.0 | 5-45)          -                       Y (RRIPKMSRTEKRKLL,KRHRATAQSDNRARRI,KGRKRKRTDGRQKRHRA)
Hg_chrom1_TN10mRNA_1154:27-652          -                       -                       -
Hg_chrom1_TN10mRNA_1157:25-598          -                       -                       Y (RKSK,PANKRRC)
Continued ...
```

### DEEPLOC 2.0
DeepLoc is a general tool for eukaryotic protein localization, covering a much wider range of subcellular compartments using advanced deep learning techniques. It is more flexible but less focused on specific plant biology needs compared to Localizer. The software will provide different predictions based on whether you've removed the signal peptide in your protein, so here I am going to predict on my truncated proteins and proteins that were unmodified.


**Install Deeploc**
```
You will need to register for this software and they will send you a download link 
https://services.healthtech.dtu.dk/services/DeepLoc-2.0/

tar -zxvf deeploc-2.0.All.tar.gz
cd deeploc2_package
pip install . --user

```
### Run deeploc2 -- 2 runs: one on subtracted SP proteins, and one on unmodified proteins.
```
deeploc2 -f SignalPeptidesSubtracted6.fasta -o Signap6Accurate -m Accurate -d cuda
deeploc2 -f TN10FinalManualAnnotation_proteins.fasta -o AllProteins -m Accurate -d cuda
```
**Excerpt of results from Deeploc 2.0**
```
Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,Lysosome/Vacuole,Golgi apparatus,Peroxisome
Hg_chrom1_TN10mRNA_1,Nucleus,Nuclear localization signal,0.4494999945163727,0.6108999848365784,0.1216999962925911,0.07580000162124634,0.3386000096797943,0.008999999612569809,0.13439999520778656,0.048700001090765,0.042399998754262924,0.005799999926239252
Hg_chrom1_TN10mRNA_2,Cytoplasm,Nuclear localization signal,0.7009999752044678,0.4546999931335449,0.09380000084638596,0.39579999446868896,0.09120000153779984,0.017500000074505806,0.05480000004172325,0.024700000882148743,0.02889999933540821,0.011900000274181366
Hg_chrom1_TN10mRNA_3,Cytoplasm|Nucleus,Nuclear localization signal,0.5720000267028809,0.7562000155448914,0.03790000081062317,0.12630000710487366,0.052400000393390656,0.016200000420212746,0.029400000348687172,0.016100000590085983,0.03009999915957451,0.009100000374019146
Hg_chrom1_TN10mRNA_4,Nucleus,,0.4250999987125397,0.451200008392334,0.24070000648498535,0.14259999990463257,0.20229999721050262,0.03359999880194664,0.11720000207424164,0.06080000102519989,0.04450000077486038,0.05620000138878822
Hg_chrom1_TN10mRNA_5,Cytoplasm,Nuclear localization signal,0.5914000272750854,0.37220001220703125,0.13279999792575836,0.15309999883174896,0.16269999742507935,0.054999999701976776,0.09549999982118607,0.11569999903440475,0.07930000126361847,0.01510000042617321
Continued ...
```

### Create feature lists for each mRNA
I always create a excel chart for each feature of a gene, so it is nice to have a tabular list of gene name "\t" feature.  In this case I have the Signalp 6 secretion score and the number of transmembrane domains after the signal peptide is cleaved from the protein. 
```
#Signalp scores for those that are secreted
 less Signalp6_out/prediction_results.txt |awk '$2=="SP" {print $2"\t"$4}' >signalp6Scores.tab

#Number of transmembrane domains in each secreted protein 
grep "Number of predicted" SignalPeptidesSubtracted6.tmhmmout |sed 's/:/\t/g' |awk '{print $2"\t"$8}' >Signalp6Tmhmm.tab

#subcellular localization Localizer for secreted proteins
less SP6Out/Results.txt |awk 'NR>4' |awk -F"\t" '{if(substr($2,1,1)=="Y") {print $1"\tChloroplast",$2} else if(substr($3,1,1)=="Y" ) {print $1"\tMitochondria",$3} else if(substr($4,1,1)=="Y") {print $1"\tNucleus",$4} else {next;}}' |sed 's/:/\t/g' |awk '{print $1"\t"$3}' >LocalizerSP6.tab

#subcellular localization Deeploc secreted proteins
cat Signap6Accurate/results_20240910-133321.csv Signap5Accurate/results_20240910-133354.csv  |sed 's/,/\t/g'  |cut -f 1,2,3 |sed 's/:/\t/g' |cut -f 1,3,4,5 |sed 's/\t/#/1' |sed 's/\t/ /g' |sed 's/#/\t/g' >SecretedProteinsDeepLoc.tab

#subcellular localization Deeploc non-secreted proteins
cat AllProteins/results_20240910-142450.csv Signap5Accurate/results_20240910-133354.csv  |sed 's/,/\t/g'  |cut -f 1,2,3 |sed 's/:/\t/g' |cut -f 1,3,4,5 |sed 's/\t/#/1' |sed 's/\t/ /g' |sed 's/#/\t/g' >AllOtherProteinsDeepLoc.tab
```

| mRNA_Names              | Gene_Names              | SignalPv6.0_Scores | Transmembrane_Domain_Counts| Localizer 1.0.5  | DeepLoc2.0 Secreted                     | DeepLoc2.0 All Proteins                              |
|-------------------------|-------------------------|-------------------|-------------------------------------------------------|------------------|------------------------------------------|-------------------------------------------------------|
| Hg_chrom1_TN10mRNA_1     | Hg_chrom1_TN10gene_1     | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_10    | Hg_chrom1_TN10gene_10    | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| Hg_chrom1_TN10mRNA_100   | Hg_chrom1_TN10gene_96    | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| Hg_chrom1_TN10mRNA_1000  | Hg_chrom1_TN10gene_956   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Transmembrane domain                                  |
| Hg_chrom1_TN10mRNA_1001  | Hg_chrom1_TN10gene_957   | 0.499             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide                                        |
| Hg_chrom1_TN10mRNA_1002  | Hg_chrom1_TN10gene_958   | 0.006             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| Hg_chrom1_TN10mRNA_1003  | Hg_chrom1_TN10gene_959   | 0.999             | 0                                                     | Chloroplast      | Extracellular Signal peptide            | Signal peptide                                        |
| Hg_chrom1_TN10mRNA_1004  | Hg_chrom1_TN10gene_960   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| Hg_chrom1_TN10mRNA_1005  | Hg_chrom1_TN10gene_960   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear export signal                                 |
| Hg_chrom1_TN10mRNA_1006  | Hg_chrom1_TN10gene_961   | 1.000             | 1                                                     | Nucleus          | Golgi apparatus Signal peptide|Transmembrane domain | Signal peptide|Transmembrane domain                   |
| Hg_chrom1_TN10mRNA_1007  | Hg_chrom1_TN10gene_962   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1008  | Hg_chrom1_TN10gene_963   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal|Nuclear export signal     |
| Hg_chrom1_TN10mRNA_1009  | Hg_chrom1_TN10gene_964   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_101   | Hg_chrom1_TN10gene_97    | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| Hg_chrom1_TN10mRNA_1010  | Hg_chrom1_TN10gene_963   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1011  | Hg_chrom1_TN10gene_965   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| Hg_chrom1_TN10mRNA_1012  | Hg_chrom1_TN10gene_966   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide                                        |
| Hg_chrom1_TN10mRNA_1013  | Hg_chrom1_TN10gene_967   | 0.427             | #N/A                                                  | #N/A             | Cytoplasm Nuclear localization signal    | Signal peptide                                        |
| Hg_chrom1_TN10mRNA_1014  | Hg_chrom1_TN10gene_968   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| Hg_chrom1_TN10mRNA_1015  | Hg_chrom1_TN10gene_969   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| Hg_chrom1_TN10mRNA_1016  | Hg_chrom1_TN10gene_970   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1017  | Hg_chrom1_TN10gene_971   | 0.162             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| Hg_chrom1_TN10mRNA_1018  | Hg_chrom1_TN10gene_972   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1019  | Hg_chrom1_TN10gene_973   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| Hg_chrom1_TN10mRNA_102   | Hg_chrom1_TN10gene_98    | 0.997             | 0                                                     | #N/A             | Nucleus                                 | Signal peptide                                        |
| Hg_chrom1_TN10mRNA_1020  | Hg_chrom1_TN10gene_974   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Transmembrane domain                                  |
| Hg_chrom1_TN10mRNA_1021  | Hg_chrom1_TN10gene_975   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| Hg_chrom1_TN10mRNA_1022  | Hg_chrom1_TN10gene_976   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1023  | Hg_chrom1_TN10gene_977   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1024  | Hg_chrom1_TN10gene_978   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear export signal                                 |
| Hg_chrom1_TN10mRNA_1025  | Hg_chrom1_TN10gene_979   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1026  | Hg_chrom1_TN10gene_979   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1027  | Hg_chrom1_TN10gene_980   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| Hg_chrom1_TN10mRNA_1028  | Hg_chrom1_TN10gene_981   | 0.981             | 0                                                     | Nucleus          | Cytoplasm|Nucleus|Cell membrane      | Nuclear localization signal                           |
| Hg_chrom1_TN10mRNA_1029  | Hg_chrom1_TN10gene_982   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| Hg_chrom1_TN10mRNA_103   | Hg_chrom1_TN10gene_99    | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear export signal                                 |


[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)