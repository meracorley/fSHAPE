wd=$(pwd)

#module load star
#module load umitools
#module load samtools
#module load bamtools
#module load cutadapt
#module load bedtools

trim() {
if [ ${S: -3} == ".gz" ]; then
zcat "$S" | umi_tools extract --bc-pattern="$bcPattern" -L extract.log > "$file".umi.fastq
else
cat "$S" | umi_tools extract --bc-pattern="$bcPattern" -L extract.log > "$file".umi.fastq
fi
cutadapt -a "$adapter" --minimum-length=23 --overlap=5 -o "$file".trim.fastq "$file".umi.fastq
}

runStar() {
mkdir "$file"
STAR --runThreadN "$nthreads" --genomeDir "$STARREF" --readFilesIn "$file".trim.fastq --outFileNamePrefix "$file"/"$file" --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate
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
mkdir splitBam
mv "$file"mapped.sort.bam splitBam
cd splitBam
bamtools split -in "$file"mapped.sort.bam -reference -refPrefix REF_
mv "$file"mapped.sort.bam ../
mkdir dedup

##for each split bam, index it and dedup it
for bam in *.bam; do
samtools index "$bam" "$bam".bai
name=$(echo $bam | sed 's/.bam//')
loc=$(pwd)
umi_tools dedup -I "$bam" --spliced-is-unique --mapping-quality=3 -S dedup/"$name".sort.bam

cd dedup ##make sorted deduped bam:
#samtools view -q 3 "$name".bam -o temp.bam #Filter to leave out reads with poor mapping quality.
#samtools sort temp.bam -o "$name".sort.bam
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

if [[ -f "$LOGFILE" ]];	then
echo "$file Starting:$start Trimmed:$trim Mapped:$mapped Unique:$uniq %Yield:$yield %Duplicates:$duprate" >> "$wd"/"$LOGFILE"
fi
}

genomeCov(){
mkdir "$wd"/"$file"/coverage/
cd "$wd"/"$file"/splitBam/dedup

for bam in *.sort.bam; do
chrom=$(echo $bam | sed s/"$file"mapped.sort.REF_// | sed s/.sort.bam//) #the chromosome name of this bam file
bedtools genomecov -split -strand + -dz -ibam "$bam" > "$wd"/"$file"/coverage/"$chrom"+.coverage
bedtools genomecov -split -strand - -dz -ibam "$bam" > "$wd"/"$file"/coverage/"$chrom"-.coverage
bedtools genomecov -5 -strand + -dz -ibam "$bam" > "$wd"/"$file"/coverage/"$chrom".5prime+.coverage
bedtools genomecov -5 -strand - -dz -ibam "$bam" > "$wd"/"$file"/coverage/"$chrom".5prime-.coverage
done
cd $wd
}


#########Begin command line options#########
if [ $# -lt 1 ] ; then
   echo "Usage: processReadsPipeline.sh -i file1.fastq,...,filen.fastq -r STAR_reference [...]"
   echo "processReadsPipeline.sh -h for more information"
   exit

fi
usage()
{
cat << EOF

/bin/bash processReadsPipeline.sh -i file1.fastq,...filen.fastq -r STAR_reference [...]

This pipeline processes input .fastq files from fSHAPE libraries to produce coverage bedgraph files for each chromosome-strand.
Fastq files can also be in .gz format.
If .fastq files are very large it is recommended to run multiple processReadPipeline.sh instances in parallel, rather than providing
all samples as a input in a list. Very large files may take 48 hours to be processed.

This pipeline assumes an RT-stop SHAPE method was used (as opposed to mutational profiling (MaP) method), such that both total coverage 
and 5' read end coverage is calculated for each chromosome and strand. (for ex: chr1+.cov, chr1+.5cov, chr1-.cov, chr1-.5cov).

Mutational profiling methods are also possible with fSHAPE, but require different analysis steps than implemented in this pipeline. 

OPTIONS:

REQUIRED
   -i, --inputs           A .fastq file of fSHAPE sequencing reads or comma separated list of .fastq files: file1.fastq,...filen.fastq
   -r, --reference	  Path to the STAR reference to align reads to.

OPTIONAL(has defaults)
   -h, --help      Print this help message.
   -a, --adapter   The adapter sequence that will be trimmed by CutAdapt. Default is AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG.
   -u, --UMI       The UMI sequence pattern to be recognized by UMItools. Default is NNNNCCCCNNNNN.
   -p, --threads   The number of threads to use during read aligning by STAR.
   -l, --log       Name of logfile that records completion status of major steps of pipeline. Path assumed to be current working directory.
EOF
}

file="" #use processReadsPipelineMaster.sh to run all samples at once
STARREF=""
bcPattern=NNNNCCCCNNNNN
adapter=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
nthreads=1
SAMPLES=""
LOGFILE=""
while [[ -n $1 ]]; do
  case $1 in
    -h | --help) usage; exit 1;;
    -i | --input) shift; SAMPLES="$1";;
    -r | --reference) shift; STARREF="$1";;
    -a | --adapter) shift; adapter="$1" ;;
    -u | --UMI) shift; bcPattern="$1";;
    -p | --threads) shift; nthreads="$1";;
    -l | --log) shift; LOGFILE="$1";;
    -* | --*) usage; exit 1;;
    *) echo "Argument $1 needs switch..."; exit ;;
   esac
   shift
done
###############End command line options###########

##echo "$STARREF $SAMPLES $nthreads $bcPattern $adapter"

if [[ "$SAMPLES" == "" ]]; then echo "No input files provided."
exit
fi

if [[ -f "$LOGFILE" ]]; then
rm $LOGFILE
fi

for S in $(echo $SAMPLES |  sed "s/,/ /g"); do
if [[ ! -f "$S" ]]; then echo "Cannot find file $S. Please check filename/filepath."; exit ; fi

file=$(echo $S | sed "s/.fastq/ /" | awk '{print$1}' | tr '/' '\t' | rev | cut -f1 | rev)
cd $wd
trim
if [[ "$LOGFILE" != "" ]]; then
echo "TRIMMED $file" >> "$wd"/"$LOGFILE"
echo "" >> "$wd"/"$LOGFILE"
fi

cd $wd
runStar
if [[ -f "$LOGFILE" ]];	then
echo "MAPPED with STAR on $file" >> "$wd"/"$LOGFILE"
echo "STAR mapping stats are in $wd/$file/$fileLog.final.out" >> "$wd"/"$LOGFILE"
echo "" >> "$wd"/"$LOGFILE"
fi

cd $wd
dedup
if [[ -f "$LOGFILE" ]];	then
echo "DEDUPED $file with Umi_tools" >> "$wd"/"$LOGFILE"
echo "" >> "$wd"/"$LOGFILE"
fi

cd $wd
mappingStats

cd $wd
genomeCov
if [[ -f "$LOGFILE" ]];	then
echo "GENOMECOV on $file" >> "$wd"/"$LOGFILE"
echo "" >> "$wd"/"$LOGFILE"
fi
done




