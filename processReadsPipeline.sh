STARREF=""
wd="./"
cd $wd
module load star
module load umitools
module load samtools
module load bamtools
module load cutadapt
module load bedtools

trim() {
zcat "$file".fastq.gz | umi_tools extract --bc-pattern=NNNNCCCCNNNNN -L extract.log > "$file".umi.fastq
cutadapt -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG --minimum-length=23 --overlap=5 -o "$file".trim.fastq "$file".umi.fastq
}

runStar() {
mkdir "$file"
STAR --runThreadN 8 --genomeDir $STARREF --readFilesIn "$file".trim.fastq --outFileNamePrefix "$file"/"$file" --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate
samtools index "$file"/"$file"Aligned.sortedByCoord.out.bam "$file"/"$file"Aligned.sortedByCoord.out.bam.bai
}

dedup(){
#dedup on full aligned file is too big
#umitools dedup -I "$file"/"$file"Aligned.sortedByCoord.out.bam --spliced-is-unique -S "$file"/"$file".dedup.bam -L "$file"/dedup.log --output-stats="$file"/dedup.stats
cd "$wd"/"$file"

##take out unmapped reads
samtools view -b -F 4 "$file"Aligned.sortedByCoord.out.bam -o "$file"mapped.bam
samtools sort "$file"mapped.bam -o "$file"mapped.sort.bam
samtools index "$file"mapped.sort.bam "$file"mapped.sort.bam.bai

##split the bam file up by chromosome and dedup each separately
bamtools split -in "$file"mapped.sort.bam -reference
mkdir splitBam
mv *REF*.bam splitBam
cd splitBam
mkdir dedup

##for each split bam, index it and dedup it
for bam in *.bam; do
samtools index "$bam" "$bam".bai
name=$(echo $bam | sed 's/.bam//')
loc=$(pwd)
umi_tools dedup -I "$bam" --spliced-is-unique -S dedup/"$bam"

cd dedup ##make sorted deduped bam:
samtools view -q 3 "$name".bam -o temp.bam #Filter to leave out reads with poor mapping quality.
samtools sort temp.bam -o "$name".sort.bam
samtools index "$name".sort.bam "$name".sort.bam.bai
cd ../
done

##re-combine split, deduped bams
cd ../
ls splitBam/dedup/*.sort.bam > deduped_bams.txt
samtools merge -O BAM -b deduped_bams.txt "$file".dedup.bam
cd $wd
}

mappingStats(){
start=$(wc -l "$file".umi.fastq | awk '{print $1/4}')
trim=$(wc -l "$file".trim.fastq | awk '{print $1/4}')
#-F 256 to avoid counting secondary alignments
mapped=$(samtools view -c -F 256 "$file"/"$file"mapped.bam)
uniq=$(samtools view -c -F 256 "$file"/"$file".dedup.bam)
yield=$(echo "scale=2; $uniq/$start*100" | bc)
duprate=$(echo "scale=2; (1-$uniq/$mapped)*100" | bc)
echo "$file Starting:$start Trimmed:$trim Mapped:$mapped Unique:$uniq %Yield:$yield %Duplicates:$duprate"
}

genomeCov(){
mkdir "$wd"/"$file"/coverage/
cd "$wd"/"$file"/splitBam/dedup

for bam in *.sort.bam; do
chrom=$(echo $bam | sed s/"$file"singleMapped.sort.REF_// | sed s/.sort.bam//) #the chromosome name of this bam file
bedtools genomecov -split -strand + -dz -ibam "$bam" -g hg38.chrom.sizes > "$wd"/"$file"/coverage/"$chrom"+.coverage
bedtools genomecov -split -strand - -dz -ibam "$bam" -g hg38.chrom.sizes > "$wd"/"$file"/coverage/"$chrom"-.coverage
bedtools genomecov -5 -strand + -dz -ibam "$bam" -g hg38.chrom.sizes > "$wd"/"$file"/coverage/"$chrom".5prime+.coverage
bedtools genomecov -5 -strand - -dz -ibam "$bam" -g hg38.chrom.sizes > "$wd"/"$file"/coverage/"$chrom".5prime-.coverage
done

}

cell="$1"
file="$2" #if using processReadsPipelineMaster.sh to run all samples at once, sample name is provided as argument to this script

echo "TRIMMING..."
trim

echo "RUNNING STAR..."
runStar

echo "DEDUPING..."
dedup

echo "CALC MAPPING STATS..."
mappingStats

echo "GENOMECOV"
genomeCov

