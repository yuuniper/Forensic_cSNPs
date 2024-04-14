# Forensic_cSNPs 

## Introduction 

This bioinformatics pipeline is designed to detect forensically relevant coding region single nucleotide polymorphisms (cSNPs) from RNA-seq data.

The cSNPs in this pipeline are a set 35 cSNPs sourced from [Dorum et. al](10.1016/j.fsigen.2022.102685) in body fluid-specific mRNA transcripts that represent a direct link between body fluids and their donors. 

The pipeline is adaptable for use on compute clusters with job submission engines like SLURM, as well as on standalone machines.

The pipeline offers seamless execution from start to finish, allowing users to initiate the analysis from raw FASTQ files and proceed through somatic variant calling with a single submission command.

The pipeline supports both single-end and paired-end data. 

To read more about this pipeline's methodology, please click here for undergraduate thesis link: (Will be updated by 5/10/24 once published by UCF STARS.)


## Install all the required packages, softwares, and tools: 
> NOTE:
> You only need to do this once. If you are working on a bioinformatics compute cluster, it may have most dependencies already preinstalled. 
> Follow the steps below to install all the ones your cluster or machine does not already have.


**Java 17**

This pipeline uses Java version 17 to be compatible with Picard 3.0.0: https://www.oracle.com/java/technologies/javase/jdk17-archive-downloads.html


**Human Genome**

Please have hg38 downloaded from the NCBI website: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

**GATK**

The GATK version of this pipeline is gatk-4.0.2.1: https://github.com/broadinstitute/gatk/releases/tag/4.0.2.1


**picard**

Method 1: 
1. Download *picard.jar* from this GitHub.
2. Drag and drop *picard.jar* into your home directory.

Method 2: 
1. Download the 3.0.0 Version on Picard  here: https://github.com/broadinstitute/picard/releases/tag/3.0.0

**HISAT2**
1. Download the correct binary version of HISAT2 here: https://daehwankimlab.github.io/hisat2/download/

2. Download the Index for GRCh38 here: 	https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz 


**samtools**
Plase download samtools from here: https://www.htslib.org/download/


## How to use the preprocessing script 

> NOTE:
> Please make sure all the packages are downloaded before running the scripts!

**Scripts use 3 arguments:**
1. Directory of your input file (DO NOT INCLUDE FILENAME!)
```
/u/INPUT_DIRECTORY/
```
2. Name of the sample file (DO NOT INCLUDE DIRECTORY!)

> For example, if the name is *SAMPLE_R1.fastq.gz* then your sample name is SAMPLE
```
SAMPLE
```
3. Directory of your output file (DO NOT INCLUDE FILENAME!)
```
/u/OUTPUT_DIRECTORY/
```

4. Single End or Paired End Processing 

Enter 1: For Single End 
```
1 
```
Enter 2: For Paired End
```
2
```

For example, you can run the following command to submit as a job on the University of Central Florida (UCF) Coombs cluster: 
```
srun preprocessing.sh /u/home/user/ SRRXXXXXX /u/home/user/output/ 1
```
This means I have a single end file */u/home/user/SRRXXXXXX.fastq.gz I want to process.
I want to output my processed files into the directory */u/home/user/output/*.


