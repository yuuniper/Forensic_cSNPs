# Recommended cores

#!/bin/bash
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=8
#SBATCH --error=preprocessing_part1-%J.err
#SBATCH --output=preprocessing_part1-%J.out
#SBATCH --job-name=preprocessing_part1

# -------------------------------------------
# Preprocessing Part 1: 
# Hisat2 mapping 
# Convert SAM to BAM 
# Mark Duplicates 

# Only if needed: Add Read Groups 
# -------------------------------------------

# Input: 
# use 3 arguments: 
# [1] /u/INPUT_DIRECTORY (do not include FILENAME) 
# [2] FILENAME (do not include _R1 or _R2, use ONLY sample name)
# [3] /u/OUTPUT_DIRECTORY
# ex: srun preprocessing_part1.sh /home/nisreen/NewRNAseq/ 14D-C3 /home/aliceyu/hisat2test/

# Load Modules
module purge
module load hisat2

hisat2 -p 8 --dta -x 
$1/hg38/hg38/genome 
-1 $1/trimmomatic/$2"_R1".trim.fastq 
-2 $1/trimmomatic/$2"_R2".trim.fastq 
-S $3$2.sam 2> $3$2.summary

# Convert SAM file to BAM file 
module load samtools 
samtools view -bS $3$2.sam > $3$2.bam

# Mark Duplicates 
# Load Modules
module purge
module /utils/gatk-4.0.2.1
 
/utils/gatk-4.0.2.1/gatk MarkDuplicatesSpark \
-I $3$2.bam \
-O $3$2"MarkDuplicatesPicard".bam \
-M $3"marked_duplicates_metrics".txt

# Add Read Groups if needed 
# Load Modules
module purge
module load java/17
module load picard3

# Modify as needed  
picard  AddOrReplaceReadGroups \
I=$2"MarkDuplicatesPicard".bam \
O=$3"14D-C3-AddReplaceRG".bam \
RGID=XXX-XXXX\
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=00000






