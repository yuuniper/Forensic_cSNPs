# Forensic_cSNPs

An RNA-seq pipeline for detecting forensically relevant coding-region SNPs (cSNPs) in body-fluid-specific mRNA transcripts.

## Research question

mRNA profiling is used in forensics to identify which body fluid (blood, saliva, semen, etc.) a stain came from. Can the same RNA-seq reads also be used to genotype coding SNPs inside those transcripts, linking the body fluid to its donor without a separate DNA assay?

The pipeline targets a panel of 35 cSNPs in body-fluid-specific transcripts, taken from Dorum et al., *Forensic Sci Int Genet* 2022 ([doi:10.1016/j.fsigen.2022.102685](https://doi.org/10.1016/j.fsigen.2022.102685)).

## Data

- **Input:** public bulk RNA-seq runs (single-end or paired-end FASTQ, e.g. SRA accessions such as `SRRXXXXXX`), downloaded with the SRA Toolkit.
- **Reference:** human genome hg38 (NCBI GCF_000001405.26) with a HISAT2 GRCh38 index, and a gene annotation BED (`hg38.knownGeneAddGene_name.bed.gz`) used to annotate variants.
- **Not in this repo:** no sequencing data, VCFs, or reference files are committed. Raw reads from human subjects can be identifiable or access-restricted, so obtain them from the original source under its terms of use. `.gitignore` blocks common data file types.

## Methods

`preprocessing.sh` runs, per sample:

1. **Trim** reads with fastp.
2. **Align** to hg38 with HISAT2 (`--dta`).
3. **Convert** SAM to BAM with samtools; validate with Picard `ValidateSamFile`.
4. **Mark duplicates** with GATK `MarkDuplicatesSpark`; add read groups with Picard.
5. **Sort and index** with samtools.
6. **Call variants** with GATK `Mutect2`.
7. **Select SNPs** (`SelectVariants`) and **hard-filter** them (`VariantFiltration`: QD < 2, FS > 60, MQ < 40, SOR > 4, MQRankSum < -12.5, ReadPosRankSum < -8).
8. **Annotate** with `bcftools annotate` against the gene annotation BED.
9. **Extract cSNPs** with `find_cSNP.py` (see "Known gaps").

Single-end and paired-end data are both supported, and the script can be run under SLURM (`srun`) or on a standalone machine.

## Results

This repo contains the pipeline code only. No result files or figures are included. Methods and results are described in the author's undergraduate thesis (University of Central Florida; link to be added once published on UCF STARS).

## Reproducing

### Requirements

- Linux, Java 17 (Picard 3.0.0) and Java 8 (GATK 4.0.2.1 / 4.3.0.0 as invoked in the script)
- hg38 FASTA plus HISAT2 GRCh38 index ([download](https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz))
- GATK [4.0.2.1](https://github.com/broadinstitute/gatk/releases/tag/4.0.2.1) and 4.3.0.0

Most other tools can be installed with conda:

```bash
conda env create -f environment.yml
conda activate forensic-csnps
```

`picard.jar` is no longer committed; the environment provides Picard 3.0.0 (or download it from the [release page](https://github.com/broadinstitute/picard/releases/tag/3.0.0)).

### Running

```bash
bash preprocessing.sh <input_dir>/ <sample_name> <output_dir>/ <1|2>
```

| Argument | Meaning |
|---|---|
| `input_dir` | Directory containing the reference folder (`hg38/`) and `extractVariants/hdr.txt`. Keep the trailing `/`. |
| `sample_name` | Sample ID without `_R1`/`_R2` or extension, e.g. `SRRXXXXXX`. |
| `output_dir` | Where `<sample>_processing/` is created. Keep the trailing `/`. |
| `1` or `2` | `1` = single-end, `2` = paired-end. |

On SLURM: `srun bash preprocessing.sh /u/home/user/ SRRXXXXXX /u/home/user/output/ 1`

## Known gaps (not changed in this cleanup)

The cleanup did not alter pipeline logic, so these remain:

- `find_cSNP.py` is called by the last step but is not in this repo.
- `preprocessing.sh` has hard-coded cluster paths (`/home/al018041/utils/...`, `/share/apps/gatk-4.0.2.1`) and `module load` lines specific to the UCF cluster; the FASTQ location is hard-coded rather than taken from `input_dir`. Edit these for your system.
- `extractVariants/hdr.txt` and the annotation BED are referenced but not included.
- Tool versions for fastp, HISAT2, samtools and bcftools were not recorded.

## Repository structure

```
.
├── README.md
├── preprocessing.sh   # end-to-end FASTQ -> annotated SNP VCF -> cSNP list
├── environment.yml    # conda environment (partially pinned, see comments)
└── .gitignore         # blocks data, BAM/VCF/FASTQ, jars, secrets
```
