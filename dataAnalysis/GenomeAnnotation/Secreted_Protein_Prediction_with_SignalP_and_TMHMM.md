---
title: Identifying secreted proteins and predicting their subcellular localization
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

Here we will be using a set of predicted proteins from a plant parasitic nematode genome to predict protein secretion, transmembrane domains, and subcellular localization.

**Software used in this tutorial**
- SignalP 6.0 [Teufel et al., 2022](https://www.nature.com/articles/s41587-021-01156-3)
- TMHMM 2.0c [Krogh et al., 2001](https://pubmed.ncbi.nlm.nih.gov/11152613/)
- Samtools 1.16.1 [Li et al., 2009](https://pubmed.ncbi.nlm.nih.gov/19505943/)
- Localizer 1.0.5 [Sperschneider et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5353544/)
- DeepLoc 2.0 [Thumuluri et al., 2022](https://academic.oup.com/nar/article/50/W1/W228/6576357)
- EMBOSS [Rice et al., 2000](http://emboss.open-bio.org/html/use/pr02s04.html)

# Secretion
### Signalp 6.0
SignalP 6.0 leverages deep neural networks to predict the presence, location, and cleavage sites of signal peptides in protein sequences. 

**Installation of SignalP 6.0**
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

#On my system it installed with numpy 2, which is incorrect, so I had to run this. 
pip install numpy==1.24.3

#This code stores the installation path of the signalp module in the variable SIGNALP_DIR, so that you can use this path in further shell commands.
SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
```

### Run Signalp 6.0
```
#run signalp 6 in a compute node
echo "source signalp6_env/bin/activate; signalp6 --fastafile YourSpecies_proteins.fasta --organism euk --output_dir Signalp6_out  --format txt --mode fast --model_dir signalp6_fast/signalp-6-package/models/" >signalp.sh
```

**Excerpt of results from Signalp 6.0** 
We want to get the fasta names that are secreted and the last position of the cleavage site.  
```
# SignalP-6.0   Organism: Eukarya       Timestamp: 20241007151909
# ID    Prediction      OTHER   SP(Sec/SPI)     CS Position
mRNA_1    OTHER   1.000000        0.000000
mRNA_2    OTHER   1.000000        0.000001
mRNA_3    OTHER   1.000000        0.000000
mRNA_4    OTHER   1.000000        0.000003
mRNA_5    OTHER   1.000000        0.000023
mRNA_6    OTHER   1.000000        0.000006
mRNA_7    OTHER   1.000000        0.000000
mRNA_8    OTHER   0.999587        0.000437
mRNA_9    OTHER   1.000000        0.000006
mRNA_10   OTHER   1.000000        0.000000
mRNA_11   SP      0.000243        0.999714        CS pos: 19-20. Pr: 0.9780
mRNA_12   OTHER   1.000000        0.000000
mRNA_13   SP      0.000200        0.999755        CS pos: 19-20. Pr: 0.9861
mRNA_14   OTHER   0.976381        0.023633
mRNA_15   OTHER   1.000000        0.000039
mRNA_16   OTHER   1.000000        0.000000
mRNA_17   OTHER   1.000000        0.000000
mRNA_18   OTHER   1.000000        0.000000
mRNA_19   OTHER   0.844543        0.155449
mRNA_20   OTHER   1.000000        0.000005
mRNA_21   SP      0.000452        0.999497        CS pos: 25-26. Pr: 0.9702
mRNA_22   OTHER   1.000000        0.000004
mRNA_23   OTHER   0.999822        0.000190
mRNA_24   OTHER   0.999777        0.000249
mRNA_25   SP      0.005658        0.994305        CS pos: 22-23. Pr: 0.5986
mRNA_26   SP      0.000441        0.999528        CS pos: 22-23. Pr: 0.7915less re
continued ...
```


### Identify secreted proteins in signalp6 results and remove the signal peptide for transmembrane domain prediction
```
# create a faidx index of our proteins to subtract the signal peptide from the fasta of secreted proteins
ml samtools;samtools faidx YourSpecies_proteins.fasta

#creates a small script to extract each secreted protein with the signal peptide truncated from the sequnce.  
awk '$2=="SP"{print $1"\t"$7}'  Signalp6_out/prediction_results.txt|sed -e 's/\.//g' -e 's/-/\t/g' |awk -F"\t" '{print "samtools faidx YourSpecies_proteins.fasta "$1":"$3"- >>SignalPeptidesSubtracted6.fasta"}' >ExtractSignalPContainingProts6.sh

sh ExtractSignalPContainingProts6.sh
```

**Excerpt from SignalPeptidesSubtracted6.fasta**
```
>mRNA_11:20-
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
>mRNA_13:20-
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
>mRNA_21:26-
TISTAFSQFLRNYFNDEVEKNIARRDLGTDGSFGGGKGRIKAKKRPVIFVHGLTNLAGEL
DYVRRLFREKGGYRDSELYATTYGYGMKGWKWLRDAMKCDHVKMIRIMIRAITQFTHSSQ
AILGGVCVDTEEQLGPPLTHLVHTFIGVGGANREAVHLCRLFEWAMPCNPVNGMRCHSQF
LYDINSNVGYEATTRIFVIRSMDDGTVGTRDCEGRSVSAIDGQNDEIVLRNYSHQMVIFG
TGEQQLKLLTF.
```
# Transmembrane domains 
### TMHMM 2.0
TMHMM 2.0 uses a hidden Markov model (HMM) to predict transmembrane helices in protein sequences. 

**Install and run TMHMM** 
```
#Download and extract tmhmm from here. 
https://services.healthtech.dtu.dk/services/TMHMM-2.0/

add this to your ~/.bashrc 
#export PATH="/path/to/your/tmhmm-2.0c/bin/:$PATH"

#run tmhmm
tmhmm SignalPeptidesSubtracted6.fasta >SignalPeptidesSubtracted6.tmhmmout
```

**Excerpt of the results from TMHMM**
```
# mRNA_11:20- Length: 571
# mRNA_11:20- Number of predicted TMHs:  0
# mRNA_11:20- Exp number of AAs in TMHs: 0.00222
# mRNA_11:20- Exp number, first 60 AAs:  0
# mRNA_11:20- Total prob of N-in:        0.00112
mRNA_11:20-       TMHMM2.0        outside      1   571
# mRNA_13:20- Length: 777
# mRNA_13:20- Number of predicted TMHs:  0
# mRNA_13:20- Exp number of AAs in TMHs: 0.01067
# mRNA_13:20- Exp number, first 60 AAs:  0.00038
# mRNA_13:20- Total prob of N-in:        0.00040
mRNA_13:20-       TMHMM2.0        outside      1   777
# mRNA_21:26- Length: 251
# mRNA_21:26- Number of predicted TMHs:  0
# mRNA_21:26- Exp number of AAs in TMHs: 0.02999
# mRNA_21:26- Exp number, first 60 AAs:  0.02192
# mRNA_21:26- Total prob of N-in:        0.01958
mRNA_21:26-       TMHMM2.0        outside      1   251
# mRNA_25:23- Length: 2168
# mRNA_25:23- Number of predicted TMHs:  0
# mRNA_25:23- Exp number of AAs in TMHs: 30.69038
# mRNA_25:23- Exp number, first 60 AAs:  0.01395
# mRNA_25:23- Total prob of N-in:        0.00072
continued ...
```

# Subcellular Localization
We are using two distinct subcellular localization prediction tools, each of which employs a different method to determine the cellular compartment.


| Feature                  | Localizer                                  | DeepLoc                                  |
|--------------------------|--------------------------------------------|------------------------------------------|
| **Focus**                 | Plant-specific, focused on chloroplast, mitochondrion, and nucleus. | General eukaryotic protein localization across 10 compartments. |
| **Subcellular compartments** | Chloroplast, mitochondrion, nucleus.     | Nucleus, cytoplasm, extracellular, mitochondrion, ER, Golgi, lysosome/vacuole, peroxisome, plasma membrane, chloroplast. |
| **Methodology**           | Rule-based predictions using known motifs/signals. | Deep learning model (RNN) based on sequence data. |
| **Organism specificity**  | Plant proteins.                            | General eukaryotes, including plants. |
| **Prediction accuracy**   | Focuses on high accuracy for three compartments. | Predicts across a broad range of compartments using neural networks. |
| **Ease of interpretation**| Results are interpretable based on known localization signals. | Predictions from a neural network may be harder to interpret. |
| **Input format**          | Protein sequences (FASTA).                 | Protein sequences (FASTA). |


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
#export PATH="/path/to/your/software/LOCALIZER/Scripts/:$PATH"
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
mRNA_1003:25-587          Y (0.964 | 62-86)       -                       Y (RKRNRRLSASSVHRRRRHS,RRNCNEKLCNFTQETKRWR)
mRNA_1006:26-1143         -                       -                       Y (KWRRRNNGNKKGERSKKEKRRP)
mRNA_1028:26-246          -                       -                       Y (APAPHSQRKEKQGKKEKAKSSSSSSSAEEEEKKHRKKGNGRQIVKHFGKGKKKD)
mRNA_102:19-109           -                       -                       -
mRNA_1045:36-353          -                       -                       -
mRNA_1046:24-361          -                       -                       -
mRNA_105:30-452           -                       -                       -
mRNA_1069:30-355          -                       -                       -
mRNA_1090:22-475          -                       -                       Y (KKGKIGGFLIQKIQRQRRQ)
mRNA_1092:44-1224         -                       -                       -
mRNA_1097:25-161          -                       -                       -
mRNA_1108:25-397          -                       -                       -
mRNA_1115:30-4744         -                       -                       Y (PKPS)
mRNA_1116:23-494          -                       -                       -
mRNA_112:24-315           -                       -                       Y (RKIMQQTNAKKAFKK,KKGDKLKPKKPPAKNAPEAPVTGKSKAEEKGQEKEEKKGTTDKKEAGAEKKKKEL)
mRNA_113:26-479           -                       -                       -
mRNA_1150:30-127          Y (1.0 | 5-45)          -                       Y (RRIPKMSRTEKRKLL,KRHRATAQSDNRARRI,KGRKRKRTDGRQKRHRA)
mRNA_1154:27-652          -                       -                       -
mRNA_1157:25-598          -                       -                       Y (RKSK,PANKRRC)
Continued ...
```

### DEEPLOC 2.0
DeepLoc is a general tool for eukaryotic protein localization, covering a much wider range of subcellular compartments using advanced deep learning techniques. It is more flexible but less focused on specific plant biology needs compared to Localizer. The software will provide different predictions based on whether you've removed the signal peptide in your protein, so here I am going to predict on my truncated secreted proteins but also all proteins together.

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
deeploc2 -f SignalPeptidesSubtracted6.fasta -o Signap6Accurate -m Fast
deeploc2 -f YourSpecies_proteins.fasta -o AllProteins -m Fast


```
**Excerpt of results from Deeploc 2.0**
```
Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,Lysosome/Vacuole,Golgi apparatus,Peroxisome
mRNA_1,Nucleus,Nuclear localization signal,0.4494999945163727,0.6108999848365784,0.1216999962925911,0.07580000162124634,0.3386000096797943,0.008999999612569809,0.13439999520778656,0.048700001090765,0.042399998754262924,0.005799999926239252
mRNA_2,Cytoplasm,Nuclear localization signal,0.7009999752044678,0.4546999931335449,0.09380000084638596,0.39579999446868896,0.09120000153779984,0.017500000074505806,0.05480000004172325,0.024700000882148743,0.02889999933540821,0.011900000274181366
mRNA_3,Cytoplasm|Nucleus,Nuclear localization signal,0.5720000267028809,0.7562000155448914,0.03790000081062317,0.12630000710487366,0.052400000393390656,0.016200000420212746,0.029400000348687172,0.016100000590085983,0.03009999915957451,0.009100000374019146
mRNA_4,Nucleus,,0.4250999987125397,0.451200008392334,0.24070000648498535,0.14259999990463257,0.20229999721050262,0.03359999880194664,0.11720000207424164,0.06080000102519989,0.04450000077486038,0.05620000138878822
mRNA_5,Cytoplasm,Nuclear localization signal,0.5914000272750854,0.37220001220703125,0.13279999792575836,0.15309999883174896,0.16269999742507935,0.054999999701976776,0.09549999982118607,0.11569999903440475,0.07930000126361847,0.01510000042617321
Continued ...
```

### Create feature lists for each mRNA
I always create a excel chart for each feature of a gene, so it is nice to have a tabular list of gene name "\t" feature.  
```
#Signalp scores for proteins that are secreted
 less Signalp6_out/prediction_results.txt |awk '$2=="SP" {print $2"\t"$4}' >signalp6Scores.tab

#Number of transmembrane domains in each secreted protein 
grep "Number of predicted" SignalPeptidesSubtracted6.tmhmmout |sed 's/:/\t/g' |awk '{print $2"\t"$8}' >Signalp6Tmhmm.tab

#subcellular localization Localizer for secreted proteins
less SP6Out/Results.txt |awk 'NR>4' |awk -F"\t" '{if(substr($2,1,1)=="Y") {print $1"\tChloroplast",$2} else if(substr($3,1,1)=="Y" ) {print $1"\tMitochondria",$3} else if(substr($4,1,1)=="Y") {print $1"\tNucleus",$4} else {next;}}' |sed 's/:/\t/g' |awk '{print $1"\t"$3}' >LocalizerSP6.tab

#subcellular localization Deeploc using secreted proteins only
cat Signap6Accurate/results_20240910-133321.csv Signap5Accurate/results_20240910-133354.csv  |sed 's/,/\t/g'  |cut -f 1,2,3 |sed 's/:/\t/g' |cut -f 1,3,4,5 |sed 's/\t/#/1' |sed 's/\t/ /g' |sed 's/#/\t/g' >SecretedProteinsDeepLoc.tab

#subcellular localization Deeploc all proteins
cat AllProteins/results_20240910-142450.csv Signap5Accurate/results_20240910-133354.csv  |sed 's/,/\t/g'  |cut -f 1,2,3 |sed 's/:/\t/g' |cut -f 1,3,4,5 |sed 's/\t/#/1' |sed 's/\t/ /g' |sed 's/#/\t/g' >AllOtherProteinsDeepLoc.tab
```

| mRNA_Names              | Gene_Names              | SignalPv6.0_Scores | Transmembrane_Domain_Counts| Localizer 1.0.5  | DeepLoc2.0 Secreted                     | DeepLoc2.0 All Proteins                              |
|-------------------------|-------------------------|-------------------|-------------------------------------------------------|------------------|------------------------------------------|-------------------------------------------------------|
| mRNA_1     | gene_1     | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_10    | gene_10    | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| mRNA_100   | gene_96    | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| mRNA_1000  | gene_956   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Transmembrane domain                                  |
| mRNA_1001  | gene_957   | 0.499             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide                                        |
| mRNA_1002  | gene_958   | 0.006             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| mRNA_1003  | gene_959   | 0.999             | 0                                                     | Chloroplast      | Extracellular Signal peptide            | Signal peptide                                        |
| mRNA_1004  | gene_960   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| mRNA_1005  | gene_960   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear export signal                                 |
| mRNA_1006  | gene_961   | 1.000             | 1                                                     | Nucleus          | Golgi apparatus Signal peptide|Transmembrane domain | Signal peptide|Transmembrane domain                   |
| mRNA_1007  | gene_962   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_1008  | gene_963   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal|Nuclear export signal     |
| mRNA_1009  | gene_964   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_101   | gene_97    | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| mRNA_1010  | gene_963   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_1011  | gene_965   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| mRNA_1012  | gene_966   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide                                        |
| mRNA_1013  | gene_967   | 0.427             | #N/A                                                  | #N/A             | Cytoplasm Nuclear localization signal    | Signal peptide                                        |
| mRNA_1014  | gene_968   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| mRNA_1015  | gene_969   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| mRNA_1016  | gene_970   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_1017  | gene_971   | 0.162             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| mRNA_1018  | gene_972   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_1019  | gene_973   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| mRNA_102   | gene_98    | 0.997             | 0                                                     | #N/A             | Nucleus                                 | Signal peptide                                        |
| mRNA_1020  | gene_974   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Transmembrane domain                                  |
| mRNA_1021  | gene_975   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| mRNA_1022  | gene_976   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_1023  | gene_977   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_1024  | gene_978   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear export signal                                 |
| mRNA_1025  | gene_979   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_1026  | gene_979   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear localization signal                           |
| mRNA_1027  | gene_980   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | #N/A                                                  |
| mRNA_1028  | gene_981   | 0.981             | 0                                                     | Nucleus          | Cytoplasm|Nucleus|Cell membrane      | Nuclear localization signal                           |
| mRNA_1029  | gene_982   | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Signal peptide|Transmembrane domain                   |
| mRNA_103   | gene_99    | 0.000             | #N/A                                                  | #N/A             | #N/A                                     | Nuclear export signal                                 |


[Back to the Assembly and Annotation Index page](annotation_and_assembly_index.md)