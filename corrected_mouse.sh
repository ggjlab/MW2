#BSUB -q normal
#BSUB -J col03
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=18]"
#BSUB -n 54

cd `pwd`
str1="_1.fq.gz"
str2="_2.fq.gz"
sample_name=$(basename `pwd`)
# sample_name="E100016613_L01_$(basename `pwd`)"

## merge bam file without filtering linkers
java -jar /share/home/guoguoji/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar FastqToSam F1=${sample_name}$str1 F2=${sample_name}$str2  O=H.bam QUALITY_FORMAT=Standard SAMPLE_NAME=sample_name 

## Example Cell Barcode:
/share/home/guoguoji/tools/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended INPUT=H.bam OUTPUT=unaligned_tagged_Cell.bam  SUMMARY=unaligned_tagged_Cellular.bam_summary.txt BASE_RANGE=1-18,25-34 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1

# read counts of H.bam
read_counts_total=`samtools view -c H.bam`
echo "H.bam: ${read_counts_total}"


## Example Molecular Barcode:
/share/home/guoguoji/tools/Drop-seq_tools-1.12/TagBamWithReadSequenceExtended INPUT=unaligned_tagged_Cell.bam OUTPUT=unaligned_tagged_CellMolecular.bam SUMMARY=unaligned_tagged_Molecular.bam_summary.txt BASE_RANGE=19-24 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1

rm unaligned_tagged_Cell.bam

## FilterBAM:
/share/home/guoguoji/tools/Drop-seq_tools-1.12/FilterBAM TAG_REJECT=XQ INPUT=unaligned_tagged_CellMolecular.bam OUTPUT=unaligned_tagged_filtered.bam

read_counts=`samtools view -c unaligned_tagged_CellMolecular.bam`
echo "unaligned_tagged_CellMolecular.bam: ${read_counts}"

rm unaligned_tagged_CellMolecular.bam

## PolyATrimmer
/share/home/guoguoji/tools/Drop-seq_tools-1.12/PolyATrimmer INPUT=unaligned_tagged_filtered.bam OUTPUT=unaligned_mc_tagged_polyA_filtered.bam OUTPUT_SUMMARY=polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=10

read_counts=`samtools view -c unaligned_tagged_filtered.bam`
echo "unaligned_tagged_filtered.bam: ${read_counts}"

rm unaligned_tagged_filtered.bam
## bam to sam file
bamfile="unaligned_mc_tagged_polyA_filtered.bam"
samfile="unaligned_mc_tagged_polyA_filtered.sam"
samtools view -h $bamfile > $samfile
echo "bam to sam file done"


read_counts=`samtools view -c unaligned_mc_tagged_polyA_filtered.bam`
echo "unaligned_mc_tagged_polyA_filtered.bam: ${read_counts}"
rm unaligned_mc_tagged_polyA_filtered.bam

## corrected bam for one mismatch
python_path="/share/home/guoguoji/Desktop/callsnptools/anaconda2/bin"
script_path="/share/home/guoguoji/tools/Microwellseq_barcode"
path=`pwd`
input_folder=$path
output_folder=$path
#input_folder=$(pwd)
#output_folder=$(pwd)
barcodepath="/share/home/guoguoji/tools/Microwellseq_barcode"
script=$script_path/Microwellseq2_correctBC_mismatch.py

$python_path/python2.7 $script $input_folder $output_folder $barcodepath $samfile
echo "correct sam files done"

rm unaligned_mc_tagged_polyA_filtered.sam
## sam to bam file
samtools view -b unaligned_tagged_filtered_corrected.sam > filtered.bam

rm unaligned_tagged_filtered_corrected.sam

## read counts of filtered.bam
read_counts=`samtools view -c filtered.bam`
echo "filtered.bam: ${read_counts}"

ratio=`echo "scale=2; $read_counts/$read_counts_total" | bc`
echo -e "Total reads count: ${read_counts_total}\nClean reads count: ${read_counts}\nRatio: ${ratio}"


#BAM_FILE="unaligned_tagged_filtered_corrected.bam"
#filter_file="/share/home/guoguoji/tools/Microwellseq_barcode/barcode_88K_tag.txt"
# Save the header lines
#samtools view -H $BAM_FILE > SAM_header
# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
#samtools view $BAM_FILE | LC_ALL=C grep -F -f $filter_file > filtered_SAM_body
# Combine header and body
#cat SAM_header filtered_SAM_body > filtered.sam
# Convert filtered.sam to BAM format
#samtools view -b filtered.sam > filtered.bam



java -Xmx100g -jar /share/home/guoguoji/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar SamToFastq INPUT=filtered.bam FASTQ=unaligned_mc_tagged_polyA_filtered.fastq


## Alignment STAR
/share/home/guoguoji/tools/STAR-2.5.2a/source/STAR --genomeDir /share/home/guoguoji/tools/STAR_Reference_Mouse/genomeDir --readFilesIn unaligned_mc_tagged_polyA_filtered.fastq --outFileNamePrefix star

#rm unaligned_mc_tagged_polyA_filtered.fastq
## SortSam
java -Xmx100g -jar /share/home/guoguoji/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar SortSam I=starAligned.out.sam O=aligned.sorted.bam SO=queryname

## MergeBamAlignment
java -Xmx100g -jar /share/home/guoguoji/tools/Drop-seq_tools-1.12/3rdParty/picard/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=/share/home/guoguoji/tools/STAR_Reference_Mouse/genomeDir/Mus_musculus.GRCm38.88.fasta UNMAPPED_BAM=filtered.bam ALIGNED_BAM=aligned.sorted.bam OUTPUT=merged.bam INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false

rm filtered.bam
rm aligned.sorted.bam
## TagReadWithGeneExon
/share/home/guoguoji/tools/Drop-seq_tools-1.12/TagReadWithGeneExon I=merged.bam O=star_gene_exon_tagged.bam ANNOTATIONS_FILE=/share/home/guoguoji/tools/STAR_Reference_Mouse/genomeDir/Mus_musculus.GRCm38.88.gtf TAG=GE

rm merged.bam

# Digital Gene Expression
/share/home/guoguoji/tools/Drop-seq_tools-1.12/DigitalExpression I=star_gene_exon_tagged.bam O=_dge.txt.gz SUMMARY=_dge.summary.txt NUM_CORE_BARCODES=50000

# Digital Gene Expression
/share/home/guoguoji/tools/Drop-seq_tools-1.12/BAMTagHistogram I=star_gene_exon_tagged.bam O=_out_cell_readcounts.txt.gz TAG=XC

files_to_delete="unaligned_tagged_Cell.bam unaligned_tagged_CellMolecular.bam unaligned_tagged_filtered.bam unaligned_mc_tagged_polyA_filtered.fastq merged.bam star.Aligned.out.sam aligned.sorted.bam unaligned_mc_tagged_polyA_filtered.fastq unaligned_mc_tagged_polyA_filtered.sam unaligned_mc_tagged_polyA_filtered.bam filtered.bam aligned.sorted.bam starAligned.out.sam unaligned_tagged_filtered_corrected.sam"

rm $files_to_delete