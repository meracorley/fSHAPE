## fSHAPE processes reads from fSHAPE (footprinting with selective 2'-hydroxyl acylation followed by primer extension) libraries.
In regions with sufficient read coverage, this pipeline calculates a reactivity value at each nucleotide, where high reactivities
indicate bases that hydrogen bond with protein.

## DEPENDENCIES:
Bedtools

Samtools

Bamtools

UMItools

STAR

Cutadapt

Python3

Python libraries: numpy, hmmlearn, joblib, sklearn, pybedtools

## Download vitro/vivo icSHAPE struture probing reads from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149767
Or use example reads included: K562_vitro_N1.ftl.fastq, K562_vitro_N2.ftl.fastq, K562_vivo_N1.ftl.fastq,K562_vivo_N2.ftl.fastq


## Download .fa file with all hg38 chromosomes: 
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz

$ gunzip hg38.fa.gz

## Build a STAR reference for hg38 (requires > 32 Gb of memory):

$ STAR --runThreadN 8 --runMode genomeGenerate --genomeDir STAR_hg38 --genomeFastaFiles hg38.fa

## Generate coverage data files for each of the fastq samples with processReadsPipeline.sh:

$ /bin/bash processReadsPipeline.sh -i sample1.fastq,sample2.fast2,…,sampleN.fastq -r STAR_reference -a adapterSequence –u UMIpattern -p threads -l logfile.txt

See UMItools documentation for how to specify the UMI specific to your reads.

With the example fastq files the pipeline would be run as:

$ /bin/bash processReadsPipeline.sh -i K562_vitro_N1.ftl.fastq,K562_vitro_N2.ftl.fastq,K562_vivo_N2.ftl.fastq,K562_vivo_N1.ftl.fastq STAR_hg38/ -p 8 -l LOGFILE 

This calculates coverages across each chromosome for each input fastq file, which are output to sampleN/coverage/ . For example, for K562_vitro_N1.ftl its coverage files are output to ./K562_vitro_N1.ftl/coverage/
For very large input fastq files, running a separate instance of processReadsPipeline.sh per sample will complete in a more reasonable amount of time.

## Generate reactivities for regions or full transcripts with bedReactivities.py:

$ python bedReactivities.py -i input_bed_regions.bed -g reference_genome.fa -v invivo_sample1/coverage/,invivo_sample2/coverage/,...,invivo_sampleNcoverage/ -t invitro_sample1/coverage/,...invitro_sampleN/coverage/ 

The unless the "-s" option is used, the input bed regions are grouped together by name (transcript identifier in column 4 of the bed region)
such that reactivities are grouped, calculated, and normalized by transcript.

With the example files, generating reactivities for the mRNA FTL would be done as follows:

$ python bedReactivities.py -i ftl.bed -g hg38.fa -v K562_vivo_N1.ftl/coverage/,K562_vivo_N2.ftl/coverage/ -t K562_vitro_N1.ftl/coverage/,K562_vitro_N2.ftl/coverage/

This outputs .map and .rx reactivity files for each transcript with regions in the input bed file. 
If a transcript has no coverage, however, no files are output for it.

Map file format is a tab-delimited file with four columns and one row per base. Columns are:

base_number(1-based) Average_reactivity Standard_error Base(A,C,T,G)

Rx file format (or .shape file format) is a tab-delimited file with at least two columns containing the reactivity value for each replicate for each base in the region/transcript.

Traditionally, reactivity values of "-999.0" are understood to indicate "no data." High fSHAPE reactivities (typically > 2.0) indicate bases that hydrogen bond with protein. Obviously this is most certain for bases with higher coverage (>200 reads)  and good agreement between replicates.
