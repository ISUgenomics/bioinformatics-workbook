# Finding orthologs using OrthoMCL program

This is the first section of Phylogenomics chapter, where we predict orthologs from the gene predictions for carrying out downstream analysis. We will use OrthoMCL program for this purpose, using the Singularity container. We will use many different species of plasmodium as an example dataset for this exercise.

Following is the outline for this chapter:

1. Obtaining the data
2. Cleaning and formating sequences
3. Running orthoMCL program
4. Formatting and processing results (for downstream analysis)


### Downloading data

At the time of writing this tutorial, there were 13 RefSeq assemblies available on [NCBI](https://www.ncbi.nlm.nih.gov/assembly/?term=txid5820%5BOrganism%3Aexp%5D), with chromosomal level assembly. We will use them all for finding orthologs. Although we will only need the protein sequences for this section, we will also download the CDS sequences that is necessary for other sections following this chapter. This step may vary if you want to use another dataset or custom dataset.

![Downloading Proteins](assets/Fig2.png)

Download it and extract the archive using `tar` command, and explore the contents.

```
tar xf genome_assemblies.tar
cd
grep -c ">" *.faa
```

Following is the summary:

| Filename                                               | Species                         | Sequences |
|:-------------------------------------------------------|---------------------------------|---------:|
| GCF_000002415.2_ASM241v2_protein.faa                   | *Plasmodium vivax* Sal-1        | 5,392    |
| GCF_000002765.4_ASM276v2_protein.faa                   | *Plasmodium falciparum* 3D7     | 5,392    |
| GCF_000006355.1_ASM635v1_protein.faa                   | *Plasmodium knowlesi* strain H  | 5,101    |
| GCF_000321355.1_PcynB_1.0_protein.faa                  | *Plasmodium cynomolgi* strain B | 5,716    |
| GCF_000524495.1_Plas_inui_San_Antonio_1_V1_protein.faa | *Plasmodium inui* San Antonio 1 | 5,832    |
| GCF_000709005.1_Plas_vinc_vinckei_V1_protein.faa       | *Plasmodium vinckei vinckei*    | 4,954    |
| GCF_000956335.1_Plas_frag_nilgiri_V1_protein.faa       | *Plasmodium fragile*            | 5,672    |
| GCF_001601855.1_ASM160185v1_protein.faa                | *Plasmodium reichenowi*         | 5,279    |
| GCF_001602025.1_ASM160202v1_protein.faa                | *Plasmodium gaboni*             | 5,354    |
| GCF_001680005.1_ASM168000v1_protein.faa                | *Plasmodium coatneyi*           | 5,516    |
| GCF_900002335.2_PCHAS01_protein.faa                    | *Plasmodium chabaudi chabaudi*  | 5,158    |
| GCF_900002375.1_PBANKA01_protein.faa                   | *Plasmodium berghei* ANKA       | 4,896    |
| GCF_900002385.1_PY17X01_protein.faa                    | *Plasmodium yoelii*             | 5,926    |

### Cleaning and formatting the data to desired specifications

Next step is to clean up this data so that it is easier to handle for both program and for the user. Since these files are identified by accessing id, we will rename them to some easily readable name. We will use the first 4 letters of genus name and first 4 letters of spp name (eg: `plaviva` for *Plasmodium vivax* Sal-1).

#### Files before:
```
$ ls -1
GCF_000002415.2_ASM241v2_protein.faa
GCF_000002765.4_ASM276v2_protein.faa
GCF_000006355.1_ASM635v1_protein.faa
GCF_000321355.1_PcynB_1.0_protein.faa
GCF_000524495.1_Plas_inui_San_Antonio_1_V1_protein.faa
GCF_000709005.1_Plas_vinc_vinckei_V1_protein.faa
GCF_000956335.1_Plas_frag_nilgiri_V1_protein.faa
GCF_001601855.1_ASM160185v1_protein.faa
GCF_001602025.1_ASM160202v1_protein.faa
GCF_001680005.1_ASM168000v1_protein.faa
GCF_900002335.2_PCHAS01_protein.faa
GCF_900002375.1_PBANKA01_protein.faa
GCF_900002385.1_PY17X01_protein.faa
```
##### script to rename files

```
for faa in *.faa; do
new=$(head -n 1 ${faa} |grep -wEo "Plasmodium .*\b" |cut -f 1-2 -d " "| sed 's/\(\w\w\w\w\)\w*\( \|$\)/\1/g');
mv ${faa} ${new}.fasta;
done
```
#### Files after:

```
Plasberg.fasta
Plaschab.fasta
Plascoat.fasta
Plascyno.fasta
Plasfalc.fasta
Plasfrag.fasta
Plasgabo.fasta
Plasinui.fasta
Plasknow.fasta
Plasreic.fasta
Plasvinc.fasta
Plasviva.fasta
Plasyoel.fasta
```

Next, we will clean the sequence headers. Although this is not absolutely required for these files, when using custom files it is necessary to keep consistency among all the files being used for the analysis. Here our file fasta header looks like this (all 13 are from NCBI and they are consistent)

#### Headers before:

```
$ head -n 1 *.fasta |grep "^>"
>XP_022711843.1 BIR protein, partial [Plasmodium berghei ANKA]
>XP_016652856.1 CIR protein [Plasmodium chabaudi chabaudi]
>XP_019912384.1 Dual specificity phosphatase [Plasmodium coatneyi]
>XP_004220499.1 hypothetical protein PCYB_021010, partial [Plasmodium cynomolgi strain B]
>XP_001347288.1 erythrocyte membrane protein 1, PfEMP1 [Plasmodium falciparum 3D7]
>XP_012333075.1 hypothetical protein AK88_00001, partial [Plasmodium fragile]
>XP_018638642.1 putative EMP1-like protein, partial [Plasmodium gaboni]
>XP_008813840.1 hypothetical protein C922_00001 [Plasmodium inui San Antonio 1]
>XP_002257498.1 SICA antigen (fragment), partial [Plasmodium knowlesi strain H]
>XP_012760262.2 rifin [Plasmodium reichenowi]
>XP_008621913.1 hypothetical protein YYE_00001, partial [Plasmodium vinckei vinckei]
>XP_001608309.1 hypothetical protein [Plasmodium vivax Sal-1]
>XP_022810663.1 YIR protein [Plasmodium yoelii]
```
#### Script to modify headers:

```
for fasta in *.fasta; do
cut -f 1 -d " " $fasta > ${fasta%.*}.temp;
mv ${fasta%.*}.temp $fasta
done
```

#### Headers after:

```
$ head -n 1 *.fasta |grep "^>"
>XP_022711843.1
>XP_016652856.1
>XP_019912384.1
>XP_004220499.1
>XP_001347288.1
>XP_012333075.1
>XP_018638642.1
>XP_008813840.1
>XP_002257498.1
>XP_012760262.2
>XP_008621913.1
>XP_001608309.1
```

### Running OrthoMCL program

![The overview of running OrthoMCL](assets/Fig4.png)


Since we don't want to setup/configure a MySQL server we will use a `docker` container (using the `singularity` program) for running.

```
module load singularity
singularity pull docker://docker pull granek/orthomcl
```
 Now this should create a `docker.img` file.  You can rename it if you want, but for this tutorial, we will simply use `docker.img` for simplicity.
