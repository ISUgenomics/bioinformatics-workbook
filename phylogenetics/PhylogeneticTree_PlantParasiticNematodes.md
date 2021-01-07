---
title: Call orthologues and create phylogenetic tree
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---


### Collect data
```
#/work/GIF/remkv6/USDA/17_OrthofinderTreeSecretome/01_GFFandProteins

#Collect protein fastas and gffs for each plant parasitic nematode

#From wormbase
bursaphelenchus_xylophilus.PRJEA64437.WBPS10.protein.fa
ditylenchus_destructor.PRJNA312427.WBPS10.protein.fa
globodera_pallida.PRJEB123.WBPS10.protein.fa
globodera_rostochiensis.PRJEB13504.WBPS10.protein.fa
meloidogyne_floridensis.PRJEB6016.WBPS10.protein.fa
meloidogyne_hapla.PRJNA29083.WBPS10.protein.fa
meloidogyne_incognita.PRJEA28837.WBPS10.protein.fa
bursaphelenchus_xylophilus.PRJEA64437.WBPS11.annotations.gff3
meloidogyne_javanica.PRJEB8714.WBPS11.annotations.gff3
meloidogyne_javanica.PRJEB8714.WBPS11.protein.fa
meloidogyne_incognita.PRJEB8714.WBPS11.annotations.gff3
meloidogyne_hapla.PRJNA29083.WBPS11.annotations.gff3
meloidogyne_floridensis.PRJEB6016.WBPS11.annotations.gff3
meloidogyne_arenaria.PRJEB8714.WBPS11.annotations.gff3
meloidogyne_arenaria.PRJEB8714.WBPS11.protein.fa
globodera_rostochiensis.PRJEB13504.WBPS11.annotations.gff3
globodera_pallida.PRJEB123.WBPS11.annotations.gff3
ditylenchus_destructor.PRJNA312427.WBPS11.annotations.gff3

#This is one I created from G. ellingtonae genome using RNAseq with braker
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/augustus.gff3 G.ellingtonae.gff3
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/augustus.aa G.ellingtonae.protein.fa

Published in BIOXIV at this time.
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/augustus.aa H.glycines.protein.fa
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/58_Renamatorium/1_genomeNgff/fixed.augustus.gff3 H.glycines.gff3
```

### Rename genes and gff so that they are reasonably easy to work with
```
#had to make modifications to names of proteins to get the gff and protein.fa to agree
sed 's/\.t/\./g' H.glycines.protein.fa |sed 's/Hetgly\.G/Hetgly\.T/g'  > H.glycinesFixed.protein.fa

module load maker
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Hgly --justify 8 H.glycines.gff3 >H.glycines.map
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids H.glycines.map H.glycinesFixed.protein.fa
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids H.glycines.map H.glycines.gff3

##Run this for each gff and protein fasta

#B.xylophilus
sed 's/transcript://g' B.xylophilus.gff3 |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' >B.xylophilusFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Bxyl. --justify 8 B.xylophilusFixed.gff3 >B.xylophilus.map
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids B.xylophilus.map B.xylophilus.protein.fa

#D. destructor
sed 's/transcript://g' D.destructor.gff3 |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' >D.destructorFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Ddes. --justify 8 D.destructorFixed.gff3 >D.destructor.map
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids D.destructor.map D.destructor.protein.fa
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids D.destructor.map D.destructorFixed.gff3

#G. pallida
sed 's/transcript://g' G.pallida.gff3 |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' >G.pallidaFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Gpal. --justify 8 G.pallidaFixed.gff3 >G.pallida.map
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids G.pallida.map G.pallida.protein.fa
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids G.pallida.map G.pallidaFixed.gff3

#G. rostochiensis
sed 's/transcript://g' G.rostochiensis.gff3 |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' >G.rostochiensisFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Gros. --justify 8 G.rostochiensisFixed.gff3 >G.rostochiensis.map
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids G.rostochiensis.map G.rostochiensis.protein.fa
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids G.rostochiensis.map G.rostochiensis.gff3

#G.ellingtonae
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/augustus.gff3 G.ellingtonae.gff3
cp /work/GIF/remkv6/Baum/CamTechGenomeComparison/26_GloboderaSynteny/otherGloboderaGenomes/ellingtonae/braker/GCA_001723225.1_ASM172322v1_genomic/augustus.aa G.ellingtonae.protein.fa


less G.ellingtonae.gff3 |sed 's/\.exon.*;/;/g' |sed 's/\.gene.*;/;/g' |sed 's/\.CDS.*;/;/g' | sed 's/\.transcript.*;/;/g' > G.ellingtonaeFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Gell. --justify 8 G.ellingtonaeFixed.gff3 > G.ellingtonae.map
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids G.ellingtonae.map G.ellingtonae.protein.fa
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids G.ellingtonae.map G.ellingtonaeFixed.gff3

#M. arenaria
sed 's/transcript://g' M.arenaria.gff3 |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' >M.arenariaFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Mare --justify 8 M.arenariaFixed.gff3 >M.arenaria.map

/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids  M.arenaria.map M.arenariaFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids  M.arenaria.map M.arenaria.protein.fa

#M. hapla
sed 's/transcript://g' M.hapla.gff3 |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' >M.haplaFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Mhap --justify 8 M.haplaFixed.gff3 >M.hapla.map
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids  M.hapla.map M.hapla.protein.fa
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids  M.hapla.map M.haplaFixed.gff3

#M.floridensis
less M.floridensis.gff3 |sed 's/transcript://g' |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' > M.floridensisFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Mflo --justify 8 M.floridensisFixed.gff3 >M.floridensis.map

/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids M.floridensis.map M.floridensis.protein.fa
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids M.floridensis.map M.floridensisFixed.gff3

#M.incognita
less M.incognita.gff3 |sed 's/transcript://g' |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' > M.incognitaFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Minc --justify 8 M.incognitaFixed.gff3 >M.incognita.map

##The M.incognita gff had scaffold names attached to gene names, so I had to make a separate map.
less M.incognita.map |sed 's/g/\tMinc/g' |cut -f 2,3 >M.incognitaFasta.map

#While this was successful for primary transcripts, it did not change the alternative isoforms.  
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids M.incognitaFasta.map M.incognita.protein.fa
##removing alternative isoforms that were left in the protein fasta file
module load cdbfasta
cdbfasta M.incognita.protein.fa
#only getting primary isoforms, all of which end with "RA"
grep ">" M.incognita.protein.fa |grep "RA" |sed 's/>//g' |cdbyank M.incognita.protein.fa.cidx >M.incognitaFixed.protein.fa

##This map and gff worked fine without  modification to the map, but still had cds, exon, mrna, and gene attached to gene names
less M.incognita.gff3 |sed 's/transcript://g' |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' > M.incognitaFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids M.incognita.map M.incognitaFixed.gff3

# M.javanica
less M.javanica.gff3  |sed 's/transcript://g' |sed 's/exon://g' |sed 's/gene://g' |sed 's/cds://g' > M.javanicaFixed.gff3
/shared/software/GIF/programs/maker/2.31.8b/bin/maker_map_ids --prefix Mjav --justify 8 M.javanicaFixed.gff3 >M.javanica.map
/shared/software/GIF/programs/maker/2.31.8b/bin/map_fasta_ids  M.javanica.map M.javanica.protein.fa
/shared/software/GIF/programs/maker/2.31.8b/bin/map_gff_ids  M.javanica.map M.javanicaFixed.gff3
```

### Gather primary protein isoforms and run orthofinder
```
#/work/GIF/remkv6/USDA/17_OrthofinderTreeSecretome/02_orthofinder
for f in ../01_GFFandProteins/*protein.fa; do ln -s $f; done

#this one has the original names, so we dont want it
unlink M.incognita.protein.fa
unlink H.glycines.protein.fa

#gettting ready to extract only primary isoforms
module load cdbfasta
for f in *fa; do cdbfasta $f;done

#grab the fasta ids, RA is on every primary isoform, get rid of ">", because cdbyank needs it removed
for f in *fa; do grep ">" $f |grep "RA" |sed 's/>//g' |cdbyank $f.cidx >Primary$f; done

#move all primary proteins to a new folder and run orthofinder
#just a note, the reason we want primary is that orthofinder will call alternative spliced isoforms as duplicated, and so flaw the analysis
mkdir 01_orthofinderRun
cd 01_orthofinderRun/
for f in ../Primary*fa; do ln -s $f;done
printf "module load orthofinder\northofinder -t 16 -a 16  -f 01_orthofinderRun/\n" >orthofinder.sh
################################################################
module load orthofinder
orthofinder -t 16 -a 16 -f 01_orthofinderRun/
################################################################
```

### Run signalp to get proteins with a signal peptide and lacking a transmembrane domain (secreted proteins)
```
#/work/GIF/remkv6/USDA/17_OrthofinderTreeSecretome/03_SignalP
#softlink all of the protein files
for f in ../02_orthofinder/01_orthofinderRun/*fa; do ln -s $f;done

#need to split up the protein fasta so that there are fewer than 10,000 proteins for each run, otherwise signalp complains
fasta-splitter.pl --n-parts 10 PrimaryM.hapla.protein.fa

#run signalp on all of the split fasta files
module load signalp

## Im using the notm database, outputting a signal peptide file, and feeding the output from each process into the same file. This only takes a couple minutes for 20k proteins
for f in PrimaryM.hapla.protein.part*; do echo "signalp -s notm -m -V "$f" >>"${f%%.*}".SignalP.out &"; done >M.haplaSignalp.sh
sh M.haplaSignalp.sh

for f in *fa; do fasta-splitter.pl --n-parts 20 $f ;done

#run signalp for each split fasta file, all the extra commands are making sure that the output files are the same for each species
for f in *part*fa; do echo "signalp -s notm -m -V "$f" >>"${f%.*}".signalp.out &";done |sed 's/\.part/\t/2' |sed 's/\.signalp.out/\t\.signalp.out/g' |sed '0~15 s/$/\nwait/g' >signalp.sh

```

### Run orthofinder on secreted proteins
```
# softlink the secreted protein fastas to subfolder below
#/work/GIF/remkv6/USDA/17_OrthofinderTreeSecretome/04_OrthofinderSignalP/01_orthofinderRun
for f in ../../03_SignalP/*secreted.fa; do ln -s $f;done

#/work/GIF/remkv6/USDA/17_OrthofinderTreeSecretome/04_OrthofinderSignalP
module load orthofinder
orthofinder -t 16 -a 16 -f 01_orthofinderRun/


```
