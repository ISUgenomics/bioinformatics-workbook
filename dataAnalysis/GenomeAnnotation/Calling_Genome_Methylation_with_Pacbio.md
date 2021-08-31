# Identify methylation status in the genome using existing pacbio reads



### Example run of smrttools ipdSummary
```
#explanation of use found here
https://www.pacb.com/wp-content/uploads/SMRT-Tools-Reference-Guide-v8.0.pdf

# this is their example command
ipdSummary aligned.bam --reference ref.fasta  m6A,m4C --gff basemods.gff --csv_h5 kinetics.h5

```


### smrttools setup, attempt to use previously aligned data *FAIL
```
/work/gif/remkv6/USDA/02_PacbioMethylation

ml smrtlink/10.1.0
ln -s  /work/gif/remkv6/Baum/04_DovetailSCNGenome/49_RenameChromosomes/01_Transfer2Box/SCNgenome.fasta

#Create a list subread bam files available for methylation calling
for f in  /work/GIF/archive1/Baum/091615_SCN/RawReads/*/Analysis_Results/subreads.bam; do ls $f;done >bamFile.fofn

vi bamFile.fofn
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_A01_1/Analysis_Results/m150811_212317_00125_c100874372550000001823188403031620_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_B01_1/Analysis_Results/m150812_014206_00125_c100874372550000001823188403031621_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_C01_1/Analysis_Results/m150812_060219_00125_c100874372550000001823188403031622_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_D01_1/Analysis_Results/m150812_102027_00125_c100874372550000001823188403031623_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_E01_1/Analysis_Results/m150812_144206_00125_c100874372550000001823188403031624_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10_1_E01_1/Analysis_Results/m150827_192450_00125_c100873012550000001823192503031684_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_E02_1/Analysis_Results/m150810_043842_00125_c100866742550000001823186502121680_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_F01_1/Analysis_Results/m150812_190110_00125_c100874372550000001823188403031625_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10_1_F01_1/Analysis_Results/m150827_234403_00125_c100873012550000001823192503031685_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_F02_1/Analysis_Results/m150810_085742_00125_c100866742550000001823186502121681_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_G01_1/Analysis_Results/m150828_040822_00125_c100873012550000001823192503031686_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_G02_1/Analysis_Results/m150810_132337_00125_c100866742550000001823186502121682_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_H01_1/Analysis_Results/m150828_082752_00125_c100873012550000001823192503031687_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10-1_H02_1/Analysis_Results/m150810_174223_00125_c100866742550000001823186502121683_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10_A01_1/Analysis_Results/m150829_012521_00125_c100873002550000001823192503031690_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10_B01_1/Analysis_Results/m150829_054404_00125_c100873002550000001823192503031691_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10_C01_1/Analysis_Results/m150829_100646_00125_c100873002550000001823192503031692_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10_D01_1/Analysis_Results/m150829_142812_00125_c100873002550000001823192503031693_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10_E01_1/Analysis_Results/m150829_185026_00125_c100873002550000001823192503031694_s1_p0.subreads.bam
/work/GIF/archive1/Baum/091615_SCN/RawReads/TN10_F01_1/Analysis_Results/m150829_231139_00125_c100873002550000001823192503031695_s1_p0.subreads.bam

# align the reads using blasr, it will not accept any other aligner
ml dafoam;ml blasr/5.1; blasr bamFile.fofn SCNgenome.fasta  --nproc 36 --useQuality --bam --out subreadsPBbam.bam;
ml smrttools/4.0; pbindex subreadsPBbam.bam;
ipdSummary subreadsPBbam.bam --reference SCNgenome.fasta  --gff basemods.gff


```
