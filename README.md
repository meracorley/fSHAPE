#Download vitro/vivo icSHAPE struture probing reads from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149767
#Process and map reads from any or all cell type sets (K562,HepG2,HeLa,293T) with /bin/bash processReadsPipelineMaster.sh
#First specify path to Hg38 STAR reference with "STARREF=/path/to/star/ref/" in first line of processReadsPipeline.sh
#Download .fa file with all hg28 chromosomes: 
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
#Outputs (.coverage files) to /coverage folders within folder for each sample.
#Generate fSHAPE reactivity profiles for any region or transcript (expressed as BED regions) with python getReactivities.py -h

#DEPENDENCIES:
Bedtools
Samtools
Bamtools
UMItools
STAR aligner
Cutadapt
Python libraries: numpy, hmmlearn, pybedtools
