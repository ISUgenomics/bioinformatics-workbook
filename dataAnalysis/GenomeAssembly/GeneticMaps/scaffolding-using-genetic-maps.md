# Applications of Genetic Map

In the previous section, we developed a genetic map for maize NAM founder line CML247 using the GBS dataset from [_Wallace et. al._](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004845). The original paper ([_McMullen et. al._](http://science.sciencemag.org/content/325/5941/737)) that created this NAM population also had created a genetic map, albeit with less number of markers. So in this exercise, we will use both these genetic maps to scaffold CML247 genome.

## Dataset

1. From the previous section, we need the `results_onemap.txt`, genetic map that we created for CML247.
2. MaizeGDB hosts the original markers developed by [_McMullen et. al._](http://science.sciencemag.org/content/325/5941/737), that are accessible on [MaizeGDB](https://www.maizegdb.org/data_center/pos?id=1167939). The markers can be downloaded as separate files for each chromosome. We refer them as GoldenGate markers (based on technology used to develop this).
3. A genome to scaffold. We need to have the same genome for which the map was developed and preferably as scaffolds. [MaizeGDB](https://www.maizegdb.org/genome/genome_assembly/Zm-CML247-DRAFT-PANZEA-1.0) hosts a version of CML247 that can be used for this exercise.

### Obtaining datasets

GBS markers from pervious section is available [here](assets/results_onemap.txt).

GoldenGate Markers can be downloaded using `wget` commands

```bash
wget -O CML247_1203657.map https://www.maizegdb.org/map_text?id=1203657
wget -O CML247_1203882.map https://www.maizegdb.org/map_text?id=1203882
wget -O CML247_1203682.map https://www.maizegdb.org/map_text?id=1203682
wget -O CML247_1203707.map https://www.maizegdb.org/map_text?id=1203707
wget -O CML247_1203732.map https://www.maizegdb.org/map_text?id=1203732
wget -O CML247_1203757.map https://www.maizegdb.org/map_text?id=1203757
wget -O CML247_1203782.map https://www.maizegdb.org/map_text?id=1203782
wget -O CML247_1203807.map https://www.maizegdb.org/map_text?id=1203807
wget -O CML247_1203832.map https://www.maizegdb.org/map_text?id=1203832
wget -O CML247_1203857.map https://www.maizegdb.org/map_text?id=1203857
```

Scaffolds for `CML247` can also be downloaded from MaizeGDB using `wget` command

```
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-CML247-DRAFT-PANZEA-1.0/Zm-CML247-DRAFT-PANZEA-1.0.fasta.gz
gunzip Zm-CML247-DRAFT-PANZEA-1.0/Zm-CML247-DRAFT-PANZEA-1.0.fasta.gz
mv Zm-CML247-DRAFT-PANZEA-1.0/Zm-CML247-DRAFT-PANZEA-1.0.fasta CML247-scaf.fasta
```

Apart from these primary datasets, we will also need additional data for this tutorial. First, the GoldenGate markers were developed for version 3, And GBS map was developed using the version 4, so we need both these genomes to extract marker sequences.

```bash
wget https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-22/fasta/zea_mays/dna/Zea_mays.AGPv3.22.dna.toplevel.fa.gz
gunzip *.gz
```

Second, for getting the positions for the markers that we generated, we need the co-ordinates for them. This information is located in the VCF file (cleaned), that we used in the previous section (found [here](/assets/B73xCML247_cleaned.vcf)).

Third, we need to singularity container for running ALLMAPS

```bash
module load singularity
singularity pull docker://tanghaibao/jcvi
```

## Steps

1. Prepare the genome: index it with HiSat2 aligner that is suitable to align short reads such as markers.
2. Extract the markers for both maps.
3. Map the markers to the genome (separately)
4. Merge these markers and create a weights file (GoldenGate markers need to be weighed higher than GBS markers)
5. Run ALLMAPS and create AGP file as well as Pseduomolecules for CML247.

### Preparing genome

Create a script for indexing genome `00_build_index_hisat2.sh`

```bash
#!/bin/bash
if [ "$#" -ne 2 ] ; then
echo "please provide:"
echo -e "\t\t(1) name for index, preferably the NAM line name"
echo -e "\t\t(2) genome sequences, use only scaffolds"
echo "";
echo "./00_build_index_hisat2 <NAM name> <FASTA-scaffolds-file>" ;
echo "";
exit 0;
fi

NAM=$1
file=$2
module load hisat2
hisat2-build ${file} $NAM
```

run it as:

```bash
./00_build_index_hisat2.sh CML247 CML247-scaf.fasta
```

### Processing OneMap genome map:

Split the OneMap genetic map to different chromosomes and check if they are consistent. Genetic map programs identify chromosome names based on the size (largest to smallest). While it works for most genomes, in maize the sizes are not in descending order, so we will have to rename them to correct chromosomes. We will create a series of scripts and run them to achieve this:

`om1_splitOnemap.sh`

```bash
#!/bin/bash
# run this onemap result file
# splits them based on linkage groups to check consistency
if [ "$#" -ne 1 ] ; then
echo "please provide the onemap results file";
fi

file="$1"
base=$(basename ${file%.*} |cut -f 1 -d "_")
min=$(cut -f 1 -d " " ${file} | awk 'NF>0' | sort -n | uniq | head -n 1)
max=$(cut -f 1 -d " " ${file} | awk 'NF>0' | sort -n | uniq | tail -n 1)
for i in $(seq ${min} ${max}); do
awk -v x=$i '$1==x' $file | sed 's/ /,/g' > ${base}_onemap_chr${i}.csv;
done
```

`om2_CheckSplitsConsistency.sh`
```bash
#!/bin/bash
# checks the split onemap results for consistency
if [ "$#" -ne 1 ] ; then
echo "please provide split onemap result file";
exit 0;
fi
file="$1"
lg=$(cut -f 1 -d "," $file | sort | uniq)
markers=$(cut -f 2 -d "," $file | cut -f 1 -d "_" | sort | uniq -c | awk 'BEGIN{ORS=";"}{print $2"("$1")"}');
echo -e "${file%.*}\t${lg}\t${markers}";
# count number of markers in each linkage group and warn if they have less than 10
count=$(cat $file |wc -l)
if [ "${count}" -lt "10" ] ; then
echo "$(tput setaf 1)$file has less than 10 markers$(tput sgr0)";
fi
```

`om3_FixLinkageGroupName.sh`
```bash
#!/bin/bash
# correct the linkage group to the actual chromosome name based on B73 naming
if [ "$#" -ne 1 ] ; then
echo "please provide a NAM name for the marker file";
exit 0;
fi
file="$1"
correct=$(cut -f 2 -d "," $file | cut -f 1 -d "_" | sed 's/S//g' | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k 2 -n | tail -n 1 |cut -f 1);
current=$(cut -f 1 -d "," $file | sort | uniq )
awk -v x=${correct} 'BEGIN{OFS=FS=","} $1=x'  $file > ${file%.*}_right.csv
echo "the linkage group was listed as $current but was changed to $correct"
```
`om4_MakeMarkerFile.sh`
```bash
#!/bin/bash
# run this on right files
# merges them to one file
if [ "$#" -ne 1 ] ; then
echo "please provide a NAM name for the marker file";
exit 0;
fi
base=$1
rm ${base}_onemap-cM.csv 2>/dev/null;
cat *_right.csv | grep -v "^chr" >> ${base}_onemap-cM.csv;
```

```bash
om5_ExtractAndMapMarkers.sh
#!/bin/bash
if [ "$#" -ne 4 ] ; then
echo "please provide:"
echo -e "\t\t(1) path for hisat2 indexed NAM genome (without ht2 extension)"
echo -e "\t\t(2) csv marker file checked, processed and merged"
echo -e "\t\t(3) VCF file to obtain co-ordinates for the markers"
echo -e "\t\t(4) NAM line name for naming the output file"
echo "";
echo "./04_onemap_ExtractAndMapMarkers.sh <index_path> <csv_file> <vcf_file> <NAM-name>" ;
echo "";
exit 0;
fi
module load bedtools2
module load hisat2
module load samtools
B73v4="/work/LAS/mhufford-lab/arnstrm/Canu_1.8/required-files/B73.fasta"
# full path for ht2 files
index="$1"
# full path for markers csv file
csv="$2"
# NAM genome
#genome="$3"
# vcf file for locations
vcf="$3"
# full nam name eg Ki3
full="$4"
# convert the genetic map to tab and only keep SNP name and genetic distance
sed 's/,/\t/g' $csv | grep -v "chr" | cut -f 2,3 > marker_distnace.txt
# get the B73.V4 locations for SNPs from the subsetted VCF file
grep -Fw -f <(cut -f 1 marker_distnace.txt) ${vcf} | awk '{print $3"\t"$1"\t"$2}' > marker_locations.txt
# merge locations with distance and create a bed file for extracting marker sequence
# uses 50bp up and down from the SNP location
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' marker_distnace.txt marker_locations.txt |\
  sed 's/ //g' |\
  awk 'BEGIN{OFS=FS="\t"}{print $2,$3-50,$3+50,$1,$4,"+"}' > markers_locations_distance.bed
# use bedtools to extract marker sequence
bedtools getfasta -fi ${B73v4} -fo ${full}_markers.fasta -bed markers_locations_distance.bed -name
# save marker sequnce as variable
markers=${full}_markers.fasta
# generate index for the genome
#hisat2-build ${genome} ${base}_index
# map marker sequence to the indexed genome
hisat2 -p 12 --mp 1,1 --no-softclip -f -x $index -U ${markers}  -S ${markers%.*}_mapped.sam &> mapping_stats.txt
# convert sam to bam and then sort
samtools view -b -o ${markers%.*}_mapped.bam ${markers%.*}_mapped.sam
samtools sort -o ${markers%.*}_mapped_sorted.bam ${markers%.*}_mapped.bam
# save the marker name and location they map
bedtools bamtobed -i ${markers%.*}_mapped_sorted.bam |  awk '{print $4"\t"$1"\t"$2}' > ${markers%.*}_part1.txt
# get the linkage group from the genetic map
awk 'BEGIN{OFS=FS=","} {print $2,$1,$3}' $csv | sed 's/,/\t/g' > ${markers%.*}_part2.txt
# write a csv file merging the part1 and part2 file
echo "Scaffold ID,scaffold position,LG,genetic position" > ${full}_mapped_onemap-cM.csv
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' ${markers%.*}_part2.txt ${markers%.*}_part1.txt | \
  sed 's/ //g' | \
  cut -f 2- | \
  sed 's/\t/,/g' >> ${full}_mapped_onemap-cM.csv
```


Now run all these scripts as follows:

```
./om1_splitOnemap.sh results_onemap.txt
for f in *.csv; do
   ./om2_CheckSplitsConsistency.sh $f;
done
for f in *.csv; do
   ./om3_FixLinkageGroupName.sh $f;
done
./om4_MakeMarkerFile.sh CML247
om5_ExtractAndMapMarkers.sh CML247 CML247_onemap-cM.csv B73xCML247_cleaned.vcf CML247
```

The final result will be the mapping position of the markers on the draft genome, in a CSV format (`CML247_mapped_onemap-cM.csv`). The scripts will also output other useful information such as how many markers were found per chromosomes, how many was extractable (sequence), chromosome names identified by OneMap and what was it changed to, number of mapping markers etc., Please double check these numbers to make sure the script is doing what it is supposed to do.

 
