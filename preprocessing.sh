#!/bin/bash

# Substitute the command to run your code below:
# Input: 
# [1] /u/INPUT_DIRECTORY (do not include FILENAME) 
# [2] SAMPLE_NAME (do not include _R1 or _R2, use ONLY sample name)
# [3] /u/OUTPUT_DIRECTORY
# [4] Enter 1 if single read, 2 if paired read

# Example command: 
# ex: srun preprocessing.sh /home/user/preprocessing_input/ sample_name /home/user/preprocessing_output/ 2


mkdir $3$2"_processing"
echo "___________"$2" sample folder made___________"


if [ "$4" -eq "2" ]; then
    ~/utils/fastp -i /home/al018041/utils/sratoolkit.3.0.7-ubuntu64/bin/$2_1.fastq.gz -I /home/al018041/utils/sratoolkit.3.0.7-ubuntu64/bin/$2_2.fastq.gz -o $3$2"_processing"/$2_1.trim.fastq -O $3$2"_processing"/$2_2.trim.fastq

    echo "___________fastp___________"

    hisat2 -p 8 --dta -x $1"hg38/hisat2index/hg38/genome" -1 $3$2"_processing/"$2"_1".trim.fastq -2 $3$2"_processing/"$2"_2".trim.fastq  -S $3$2"_processing/"$2.sam 2> $3$2"_processing/"$2.summary

    echo "___________HISAT2___________"

elif [ "$4" -eq "1" ]; then
    ~/utils/fastp -i /home/al018041/utils/sratoolkit.3.0.7-ubuntu64/bin/"$2".fastq.gz -o "$3""$2"_processing/"$2.trim.fastq"

    echo "___________fastp___________"

    hisat2 --phred33 --dta -x "$1"hg38/hisat2index/hg38/genome -U "$3""$2"_processing/"$2.trim.fastq" -S "$3""$2"_processing/"$2.sam" 2> "$3""$2"_processing/"$2.summary"

else
    echo "Invalid argument: 4th argument should be either '1' or '2'"
    echo "[1] should be for single reads."
    echo "[2] should be for paired reads."
fi

# Convert SAM file to BAM file 
module purge 
module load samtools 
samtools view -bS $3$2"_processing/"$2.sam > $3$2"_processing/"$2.bam

echo "___________SAMtoBAM___________"

module purge
module load java/17
module load picard3

picard  ValidateSamFile \
I=$3$2"_processing/"$2.bam \
MODE=SUMMARY

echo "___________ValidateSamFile___________"


# Mark Duplicates 
# Load Modules
module purge
module load java/8
module /share/apps/gatk-4.0.2.1

/share/apps/gatk-4.0.2.1/gatk MarkDuplicatesSpark \
-I $3$2"_processing/"$2.bam \
-O $3$2"_processing/"$2"_MarkDuplicatesPicard".bam \
-M $3$2"_processing/"$2"_MarkedDuplicates_Metrics".txt

echo "___________MarkDuplicates___________"
module purge
module load java/17
module load picard3

picard  AddOrReplaceReadGroups \
I=$3$2"_processing/"$2"_MarkDuplicatesPicard".bam \
O=$3$2"_processing/"$2"_AddReplaceRG".bam \
RGID=XXXX \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=NULL \
RGSM=$2 

echo "___________AddOrReplaceReadGroups___________"


#samtools sort $3$2"_processing/"$2"_MarkDuplicatesPicard".bam > $3$2"_processing/"$2"_sorted".bam

samtools sort $3$2"_processing/"$2"_AddReplaceRG".bam > $3$2"_processing/"$2"_sorted".bam
samtools index $3$2"_processing/"$2"_sorted".bam

echo "___________Sorted_and_Indexed___________"

module purge
module load java/8
module load /home/al018041/utils/gatk-4.3.0.0

/home/al018041/utils/gatk-4.3.0.0/gatk Mutect2 \
-R $1"hg38/hg38".fa \
-I $3$2"_processing/"$2"_sorted".bam \
-O $3$2"_processing/"$2"_Mutect2".vcf.gz

echo "___________Mutect2___________"


# Generate SNP and Indel Files 

# Load Modules
module purge
module load java/8
module load /share/apps/gatk-4.0.2.1

# 2/8/24 does not need modding outside gatk step and removing INDEL version
# This step does not take a lot of time (Less than 1 minute)
/share/apps/gatk-4.0.2.1/gatk SelectVariants \
-R $1/hg38/hg38.fa \
-V $3$2"_processing/"$2"_Mutect2".vcf.gz \
--select-type-to-include SNP \
-O $3$2"_processing/"$2"_raw_snps".vcf.gz


echo "___________SelectVariants___________"

# Filter SNP Variants 
# These are the standard values

# Load Modules
module purge
module load java/8
module load /share/apps/gatk-4.0.2.1

/share/apps/gatk-4.0.2.1/gatk VariantFiltration \
-R $1/hg38/hg38.fa \
-V $3$2"_processing/"$2"_raw_snps".vcf.gz \
-O $3$2"_processing/"$2"_filtered_snps_final".vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \

echo "___________Variant Filtration ___________"


module purge
module load samtools 

bcftools annotate -a $1"hg38/hg38".knownGeneAddGene_name.bed.gz -c CHROM,-,INFO/FEATURE,FROM,TO,-,-,INFO/FRAME,-,INFO/ATTRIBUTES -h $1"extractVariants"/hdr.txt $3$2"_processing/"$2"_filtered_snps_final".vcf -o $3$2"_processing/"$2"_annotatedSNPsfinalbcf".vcf

echo "___________bcftools ___________"

python find_cSNP.py $2"_processing/"$2"_annotatedSNPsfinalbcf.vcf" $2"_processing/"$2"_cSNPs.txt"

echo "___________python_cSNP ___________"
