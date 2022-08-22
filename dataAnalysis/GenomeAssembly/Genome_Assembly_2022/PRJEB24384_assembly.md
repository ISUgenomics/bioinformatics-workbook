* File download
Details: https://www.ebi.ac.uk/ena/browser/view/PRJEB24384?show=reads

Working Directory: `/project/isu_gif_vrsc/Siva_2022/USDA/PRJEB24384`
```
more report.txt
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/002/ERR2244302/ERR2244302_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/002/ERR2244302/ERR2244302_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/003/ERR2244303/ERR2244303_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/003/ERR2244303/ERR2244303_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/004/ERR2244304/ERR2244304_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/004/ERR2244304/ERR2244304_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/005/ERR2244305/ERR2244305_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR224/005/ERR2244305/ERR2244305_2.fastq.gz

parallel -j3 < report.txt >& download.log&
```
