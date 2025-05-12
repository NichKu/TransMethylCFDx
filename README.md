# TransMethylCFDx - A pipeline for targeted methylation and genotype-based tissue-of-origin determination of cfDNA from transplant patients.

## Overview

* [Overview](#overview)
* [Introduction](#introduction)
* [Pipeline Summary](#pipeline-summary)
* [Requirements](#requirements)
* [Download](#download)
* [Usage Guide](#usage-guide)
* [Output Folder Structure](#output-folder-structure)

## Introduction

**TransMethylCFDx** is a pipeline for targeted methylation sequencing data analysis that deconvolutes cfDNA tissue-of-origin and simultaneously quantifies donor-derived cfDNA.

## Pipeline Summary

1. Inital fastq QC ([`FASTQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adpater/Quality trimming ([`Trim Galore`](https://github.com/FelixKrueger/TrimGalore))
3. Post-trimming QC ([`FASTQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Alignment ([`bwa-meth`](https://github.com/brentp/bwa-meth))
5. Deduplication ([`UMI-tools`](https://umi-tools.readthedocs.io/en/latest/index.html))/ Consensus Sequence Creation ([`fgbio`](http://fulcrumgenomics.github.io/fgbio/))
6. Alignment and Methlyation QC (various)
7. Pileup and dd-cfDNA calculation ([`samtools`](https://www.htslib.org))
8. Report Generation ([`MultiQC`](https://multiqc.info))


## Requirements

The these following dependencies need to be available:
- Conda >= 4.10
- Snakemake >= 7.17

> **Note**
> For the pipeline to run, a conda environment will be created with all necessary dependencies and will be done automatically.
These above requirements are only needed to launch the pipeline.

## Download

Clone the Github Repo:

```
https://github.com/NichKu/TransMethylCFDx.git
cd TransMethylCFDx/workflow
```

## Quick Start

To launch the analysis, run the following command:

```
snakemake --cores [n] -s [path/to/SnakemakeFile] --use-conda --conda-frontend conda
```

For detailed explanation of the Snakemake arguments please refer to the [User Guide](https://snakemake.readthedocs.io/en/stable/#).

> **Note**
> To test whether the pipeline would run, the argument --dry-run can be added preventing from acutally running the pipeline.

## Usage Guide

To set up and launch the pipline, the steps listed below should be followed:

1. Clone the github repo
2. Install the dependencies outlined in [Requirements](#requirements)
3. Go to /Config and configure the config.yaml file
4. cd 

#### Resources
The number of cores being allocated to each job is configured in the Snakemake file under # Threads. The pipline is recommended to be launched with at least 10 CPUs. The pipline has been tested using 4GB RAM per CPU.




## Output Folder Structure

Below is the structure of the output folder:

``` 
path/to/resultsfolder/name/
│  
├── results/
│   ├── 1_fastqc_init
│   ├── 2_trimmed
│   ├── 3_fastqc_trimmed
│   ├── 4_alignment
│   │   ├── bam
│   │   └── bam_filt
│   ├── 5_dedup_consensus
│   ├── 6_aligment_QC
│   │    ├── all_metrics
│   │    ├── flagstat
│   │    ├── gcbias_metrics
│   │    ├── idxstats
│   │    ├── meth_metrics
│   │    ├── meth_per_target_metrics
│   │    ├── preseq
│   │    ├── qualimap
│   │    ├── snp_metrics
│   │    └── snp_per_target_metrics
│   ├── 7_methylQC
│   ├── 8_dd-cfDNA
│   ├── 9_methyl_ctrl
│   │   └── CEREBIS
│   └── 10_multiqc
│
└── logs/  
    ├── 0_prepare_files
    ├── 1_fastqc_init
    ├── 2_trimmed
    ├── 3_fastqc_trimmed
    ├── 4_alignment
    ├── 5_dedup_consensus
    ├── 6_aligment_QC
    ├── 7_methylQC
    ├── 8_dd-cfDNA
    ├── 9_methyl_ctrl
    └── 10_multiqc
```
