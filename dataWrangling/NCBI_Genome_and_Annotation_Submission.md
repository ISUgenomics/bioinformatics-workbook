---
title: "How to submit a genome and gene annotation NCBI"
layout: single
author: Rick Masonbrink
author_profile: true
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

Depositing a genome and its gene annotation into NCBI ensures that it is accessible to the scientific community, is properly archived, and links the relative metadata. Submission typically involves preparing the sequence data, genome, gene annotation, and assembled organelle genomes.

# The most efficient order for submitting a genome to NCBI
1. Create a Bioproject - This serves as an umbrella for your study, linking multiple datasets together.
2. Create Biosamples - This describes each of the biological samples that were acquired for assembling or annotating your genome.
3. Submit sequencing to SRA - Each independent sequencing sample that is not a technical replicate will need to be deposited as a separate file. 
4. Submit genome and annotation - Once the information for the previous submissions is complete, you should have all of the information you need to deposit your genome, annotation,assembled organellar genomes, and plasmids (if any).

# Outline of NCBI submissions described here
1. Bioproject creation
2. Biosample creation 
3. Genome-only submission (no gene annotation)
4. Genome and gene annotation submission

# Bioproject submission

If you do not have an account with NCBI yet, you will need to create one for the submission. Once that is accomplished you can go to "My submissions" and click on Bioproject at the top left of the window. Then click the 'New submission' at the top left.

The first tab to fill out is the 'SUBMITTER' tab

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/SubmitterNoName.png)

Since I created this genome and am submitting it, my submitter information is simple. This may become more complex if submitting for a group. 

Then you will need to choose the type of data you will be submitting. Here we assembled a genome, so checking this box will give us the required options needed for later genome submission steps. Then you need to choose the scope of your sample, which is often a monoisolate for genome submissions. We will also be needing locus tag prefixes to submit our annotation in a later step, so check this box.

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BioprojectProjectType.png)

This takes you to the 'TARGET' tab, where you will list the scientific name of your organism and fill out one of the five options in the second line of boxes. The description is entirely optional, but including as much information as possible will help future users utilize your data in subsequent analyses. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BioprojectTarget.png)

Then you will be taken to the 'GENERAL INFO' tab, where you will specify the release date of your data, add a public description, and select the relevance of your genome. You'll have to answer a few questions about external links to your data, grants that funded your genome assembly/annotation, consortium associations and if you are using a data provider. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BioprojectGeneralInfo.png)

The next tab is BIOSAMPLE, where you will list the previously created biosample accessions that are associated with the sequencing data used to assemble the genome. In this case, I have not created my Biosample yet, so I leave the area blank. If you do have previously submitted Biosamples, then add them here. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BioprojectBiosample.png)

The next tab is PUBLICATIONS, where you would list any published composition associated with this genome and annotation.

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BioprojectPublications.png)

The last tab is 'REVIEW & SUBMIT'. This is your last opportunity to change any of the information in your submission, as after it is accepted you will have to email NCBI to make changes. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BioprojectReviewAndSubmit.png)


# Biosample

Creating a Biosample is the standard next step, essentially to describe the organism/organisms that were sequenced. The Biosample exists to link this description with your unique Bioproject, sequence read archive(SRA), genome, etc.

Starting at 'My submissions' you should click Biosample in the top left box. Then select 'New submission' at the top left. Once there you will have to fill out the 'SUBMITTER' information if it is different from your bioproject submission. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/SubmitterNoName.png)

Then, in the 'GENERAL INFO' tab, you will need to choose a release date and whether you want to submit a single biosample or multiple.
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BiosampleGeneralInformation.png)

Then you get to the SAMPLE TYPE tab, which gives you lots of options for packages to submit your genome. Choosing a package that best represents your sample will require reading about each package, but here I will be using the Invertebrate package as an example. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BiosampleSampleType.png)

Pressing continue will move you to the 'ATTRIBUTES' tab, where you will be given the option to use a built-in editor to enter your data one at a time, or download a TSV or Excel file for editing multiple samples simultaneously. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BiosampleAttributes.png)

I submitted multiple biosamples for this project using the Excel file download. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/BiosampleExcel.png)

Then you will be given one final time to review the biosample related information before it cannot be changed, short of an email to NCBI. In my experience, if the input fits the expected format, your submission should be accepted in a few minutes. 

# Genome Only

A genome-only submission is commonly used when a gene annotation is not yet available or is to remain private. Uploading your genome-only sequence is easier than uploading a genome and gene annotation. The benefits of a genome-only submission is that if your assembly is in pseudomolecules/chromosomes then NCBI to annotate your assembly, with chromosome-level assemblies being prioritized.  

To deposit a genome you will need to go back to 'My Submissions' and then click Genome in the box at the top left corner. After that, you will see a blue 'New Submission' that allow you to start your Genome submission.   
You start with the familiar Submitter tab that you've seen previously. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/SubmitterNoName.png)

Next is the GENERAL INFO tab, where you will be asked a lot of information about your assembly and how it associates with your Bioproject and Biosamples. Since we have already submitted our Bioproject and Biosample, we can move forward by adding these accessions.

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeGeneralInfo.png)

PRJNA1153013 is the assigned Bioproject accession and SAMN43384812 is the assigned Biosample accession.


Next you have to choose the release date, which I always put as the furthest out I can (four years). This way my data will stay private long enough for me publish the research it is used in. Once the bioproject is published in a paper, all of the associated data is released.
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeReleaseDate.png)

Then you move to the 'GENOME INFO' tab, where you will need to input the various details of your data and assembly methods. I have an approximate date that the assembly was done, put in the software title, and version for my assembly method. I named my genome assembly, as this is the third update to this genome. 
Then you will need to compute the coverage of your reads to your genome size. I used nanopore and pacbio reads to assemble the genome, so I will need to get the entire length of those reads used for assembly.

**Compute genome coverage**
One of the first questions is to estimate the coverage of your reads to your genome assembly. This can be calculated by summing the total length of your input sequence reads and dividing by total genome size.
```bash
#I ran this multiple times, as I had multiple read sets

zcat *.fastq.gz | seqtk seq -A - | awk '{if(NR%2==0) total+=length($0)} END {print total}'

5,754,629,444 bp is the total length of my nanopore reads.
18,212,104,559 bp is the total length of my pacbio RSII reads

23,966,734,003 bp is the total length of all reads used in assembly
My genome size is 114,996,807, so reads/genome size = 208.4
```
I then fill out the sequencing technology of the reads used for assembly of the genome, which are GridION and RSII in my case. Because I used PacBio RSII reads for assembly, I was asked if I wanted to submit the motif_summary.csv. In this case it does not apply, since these reads are already deposited and released (SRX2692203-SRX2692222).

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeGenomeInfo.png)

There are now some follow-up questions that apply to a select number of genomes:

* Did you attempt to include the entire genome in your submission? Or is it only partial coverage of the genome?
* Is this the final version of your genome? So, will you be updating this genome any time in the near future? I will be updating this genome again, so I selected "no" here.
* Is it a de novo assembly? Did you assemble the genome from reads, or were these reference guided assemblies, etc? 
* Is this genome submission an update of an existing genome submission? There are older versions of this genome that this genome will replace, so I selected yes here. 

Then you will need to answer whether you are the assembler of the genome or if you are someone depositing the genome for someone else. If YOU assembled the genome then, you would select original. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeSubmissionCategory.png)


The next step is to create a title for your submission and to write up any details that may be pertinent to the NCBI curators. With this submission, I am planning on submitting a further improved genome from the same line.

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeSubmissionTitle.png)

The next step is filling out the 'FILES' tab.  In thi scase I was lucky that my genome was assembled into 9 complete chromosomal scaffolds without extraneous contigs. 
If your genome does not have all contigs integrated into pseudomolecule scaffolds then you will need to 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeFiles.png)

Now you will need to decide how you want to upload the files for your genome submission. I typically use Aspera over FTP when uploading due to speed, reliability, security, efficiency, and scalability.  All we need now is to find our correct genome to upload and set up the folder structure that aspera requests. 

```bash
#create this folder 
/work/gif3/masonbrink/Baum/01_ReannotateAllSCNGenomes/14_SCNBase/01_GenomesAndAnnotations/01_TN10Deposit
#copy the contents from the key_file link and paste into a file on my server. 
vi key.file
#Create another folder with just my genome
/work/gif3/masonbrink/Baum/01_ReannotateAllSCNGenomes/14_SCNBase/01_GenomesAndAnnotations/01_TN10Deposit/DepositFolder
#softlink my genome
ln -s ../../TN10Genome.fasta
#go back to the folder with the key
Then modify the command that they give us. 
ascp -i key.file -QT -l100m -k1 -d  /work/gif3/masonbrink/Baum/01_ReannotateAllSCNGenomes/14_SCNBase/01_GenomesAndAnnotations/01_TN10Deposit/DepositFolder subasp@upload.ncbi.nlm.nih.gov:uploads/arnstrm_iastate.edu_jCiP2Yhq
```
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeAspera.png)

Depending on file size the genome should upload in minutes, though the folder will not appear in the "Select preload folder' window for approximately ten minutes. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeSelectPreloadFolder.png)
Once the folder appears, then select "Use selected folder" and continue to the bottom of the page. 

The next step is the GAPS tab, where you will be propositioned to define the gaps in your genome.
* First you are prompted with "Did you randomly merge sequences together into scaffolds?"
Since it is never appropriate to merge contigs together randomly, you should select "no".

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeGapsAnalysis.png)

* Describe your gaps: are they an estimated size or are they all the same length, representing an unknown size?
In my case all of my gaps were 500bp.   
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeGapsAnalysis2.png)

If they are not etimated size then put no and proceed.  If your assembly was an all-in-one software from reads to pseudomolecule, then you will likely have estimated gap sizes. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeGapsAnalysis3.png)

What kind of technology was used to scaffold your genome assembly, if scaffolding was performed? I selected proximity-ligation, as I used HiC to scaffold my assembly. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeTaxonomy.png)

The next step is the ASSIGNMENT tab, where you will name the scaffolds that represent pseudomolecules/chromosomes, organellar genomes, and/or plasmids. 
* Do you have any sequences in your genome that represent mitochondria or chloroplast DNA? If so, you get to fill in the information about your organellar scaffolds.
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeAssignment.png)

Do you have any sequences in your genome that are representations of a plasmid? If so, you need the names of the scaffolds that are plasmids, if they are complete and if they are circular. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomePlasmidAssignment.png)

My submission is simpler than most, only needing the nine chromosomes to be labeled.  Others may find the csv more efficient.
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeChromosomeAssignments.png)

For the REFERENCES you will need to list all of the people in your group who contributed to the development of your genome assembly. Below this there is the opportunity to link any publications associated with this genome, though most times during submission a publication is not ready.
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/GenomeReferences.png)

You should have one last step to review the information before you submit. If the genome lacks errors then your submission should say processing. Depending on the size, NCBI will check your  for contamination within a few hours. If contamination is found, then you will have to modify your genomic scaffolds to remove it.  

# Genome and Annotation

To deposit an annotation with a genome the process is more complex and requires modification of your genome fasta file and genome gff3 file. 

There is software specifically designed to convert your genome and gff3 into the file format NCBI require for submitting an annotation. The next step is converting your gff3 and annotation into the desired format (.sqn), which involves using the NCBI software table2asn. Table2asn is an executable file that you can download from NCBI. There is some explanation on how it works and the download link here: https://www.ncbi.nlm.nih.gov/genbank/table2asn/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=table2asn-updated-20230706

In order to use this software you will need to make a template with the metadata for your genome here:  https://submit.ncbi.nlm.nih.gov/genbank/template/submission/

The next step is to create a 5 column feature table for your annotation. 
```bash
ml micromamba; micromamba activate cufflinks
gffread SortedTN10FinalManualAnnotation.gff3 -T -o SortedTN10FinalManualAnnotation.gtf
vi gtf_to_tbl.py
python gtf_to_tbl.py

gtf_to_tbl.py
###################################################################################################################

###################################################################################################################
```

As of March 1, 2025, an .agp (a golden path) file is no longer needed to submit a genome, though they can be included. There is a script in Juicebox software that can create this .agp file, though it does not make it perfectly. The Juicebox script uses 0 based numbering, while NCBI expects 1 based numbering.

Once you've made your template and the received the locus tag from your Bioproject, you can start using the table2asn script to create your .sqn files (ASN.1 format). 
```bash
./table2asn -t template.sbt -i TN10GenomeSNPFixedFinal.fasta -f SortedTN10FinalManualAnnotation.gff3 -gaps-unknown 500 -M n -o annotation.sqn -locus-tag-prefix ACNJWR -Z
```

This script will produce four files. 
* annotation.stats -- overview of errors in annotation
* annotation.dr -- discrepancy report has more specific error information
* annotation.val -- Different types of specific error information
* annotation.sqn -- the .sqn output you need to upload to NCBI to get your annotation deposited. 

This information is enough to deposit an annotation and genome to NCBI, however if you'd like your genes to have functional information associated with them, you will need to add "product=" and "Dbxref=" to column 9 of the mRNA and CDS features in the gff3.

Once you've made a .sqn file for your annotation, you will need to create one for your genome as well. 
```bash
./table2asn -t template.sbt -i TN10GenomeSNPFixedFinal.fasta -j "[organism=Heterodera_glycines]" -M n -o genome.sqn
```

Once you've created the .tbl, the .agp, and the two .sqn files, you are ready to begin the submission. It is very important to select the correct options here, or NCBI will not give you the opportunity to upload an annotation. 


To deposit these files you will need to go back to 'My Submissions' and then click 'Genome' in the top left corner. After that, you will see a blue 'New Submission' that allow you to start your submission.   

First, you will need to decide if you are submitting a single or multiple genomes. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/AnnotationDataType.png)

Then you will start with the familiar Submitter tab that you've previously seen. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/SubmitterNoName.png)

Then you will be forwarded into the GENERAL INFO tab where you can list your previously submitted Bioproject and Biosample accessions to be associated with your genome. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/AnnotationGeneral.png)

Then you will have to set the genome and annotation release date, which can be pushed back as far as four years. Setting four years is likely the best option for most people, as the presence of the Bioproject accession number in a publication will auto-release all associated data. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/AnnotationGenomeAndAnnotationReleaseDate.png)

The next part works in two possible ways. You can try to get all of this metadata into your .sqn file and check the box at the top, or fill it out here. Filling it out here is the easier option in my opinion, so leave the box at the top unchecked and fill out the respective information.

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/AnnotationGenomeInfo.png)

Then some more basic questions that I will elaborate on.
* Did you assemble the whole genome? 
* Are you planning to reassemble this data when you obtain more, or is this the final version of the genome?
* Did you assemble the genome using assembly software with reads or was it a reference-guided assembly?
* Is this an update of existing submission?  If so, then you can list the existing genome accession. 


![](assets/NCBI_Genome_and_Annotation_Submission_Assets/AnnotationGenomeInfo2.png)

Then check Original if you are the person that assembled the genome. Then you can make a descriptive title and if you have extra information that may be important for NCBI curators add it into the Private 'comments'.

The next step is the tricky part, at the 'FILES tab' you can miss your opportunity to deposit your annotation if you select the wrong options. Even if your chromosomes are in a single sequence and there are no extra sequences, to be able to deposit an annotation you need to choose number 2 here. 

![](assets/NCBI_Genome_and_Annotation_Submission_Assets/AnnotationFilesForSubmission.png)

Then you get to choose how you want to upload your data and, at least currently, you need to answer yes to the bottom question. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/AnnotationDataUpload.png)


Then you can answer if you have unplaced scaffolds (scaffolds not attributed to a chromosome, organelle, or plasmid).  We produced the .agp files earlier and we can use that file here. The following question is related, and asks if you have unplaced scaffolds, but you know their chromosomal assignment (just not order). The third question gives you th eoption to deposit an annotation with your genome, so select 'yes'. Then there is an odd question at the bottom that I am not entirely sure of the meaning on, but you need to select yes to this question as well so you can deposit the .sqn and .tbl files. 
![](assets/NCBI_Genome_and_Annotation_Submission_Assets/AnnotationAGPFileQs.png)
