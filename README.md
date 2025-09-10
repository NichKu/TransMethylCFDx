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

**TransMethylCFDx** is a pipeline for targeted methylation sequencing data analysis and subsequent cfDNA tissue-of-origin determination and donor-derived cfDNA quantification.

## Pipeline Summary

1. Initial fastq QC ([`FASTQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter/Quality trimming ([`Trim Galore`](https://github.com/FelixKrueger/TrimGalore))
3. Post-trimming QC ([`FASTQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Alignment ([`bwa-meth`](https://github.com/brentp/bwa-meth))
5. Deduplication ([`UMI-tools`](https://umi-tools.readthedocs.io/en/latest/index.html))/ Consensus Sequence Creation ([`fgbio`](http://fulcrumgenomics.github.io/fgbio/))
6. Alignment and Methylation QC (various)
7. Pileup Generation ([`samtools`](https://www.htslib.org))
8. Report Generation ([`MultiQC`](https://multiqc.info))


## Requirements

The following dependencies need to be available:
- Conda >= 4.10
- Snakemake >= 7.17

> **Note**
> A conda environment with all necessary dependencies will be created automatically.
These above requirements are only needed to launch the pipeline.

## Download

Clone the Github Repo:

```
git clone https://github.com/NichKu/TransMethylCFDx.git
```

## Quick Start

To launch the analysis, run the following command:

```
cd TransMethylCFDx/workflow
snakemake --cores <n> --use-conda -s ./Snakemake --conda-frontend conda
```

For a detailed explanation of the Snakemake arguments, please refer to the [User Guide](https://snakemake.readthedocs.io/en/stable/#).

> **Note**
> To test whether the pipeline will run, add the --dry-run argument to prevent it from actually running.

## Usage Guide

To set up and launch the pipeline, follow the steps below:

First, clone the github repository
Second, install the dependencies outlined in [Requirements](#requirements)
Third, go to config/ and configure the config.yaml file
```
cd ../workflow
snakemake --cores <n> --use-conda -s ./Snakemake --conda-frontend conda 
```

#### Resources
The number of cores allocated to each job is configured in the Snakemake file under the Threads setting.\
We recommend launching the pipeline with at least 10 CPUs. The pipeline has been tested with 4 GB of RAM per CPU.

## Downstream Analysis
### Tissue-of-Origin Deconvolution
To determine the tissue-of-origin, please use BAM files from 5_dedup_consensus and create PAT (.pat.gz) files using [`wgbstools`](https://github.com/nloyfer/wgbs_tools). The PAT files are then used with the reference atlas to perform the deconvolution using the [`UXM`](https://github.com/nloyfer/UXM_deconv) tool. The UXM tool also requires wgbstools.

The twist_TOO environment that was created to run the Snakemake workflow should be used for the wgbstools and UXM tool.

To perform the deconvolution starting from the BAM files:

First, activate the conda environment **twist_TOO**
```
export PATH=${PATH}:</path/to>/UXM_deconv:</path/to>/wgbs_tools
wgbstools init_genome hg19
# convert BAM to PAT files
wgbstools bam2pat -@ <n> --lbeta -o </path/to/resultsfolder>/results/pat_files -T <tmp_directory> </path/to/resultsfolder>/results/5_dedup_consensus/*cons.sort.bam
# perform deconvolution with the atlas and generate PAT files
uxm deconv --threads <n> --rlen 4 --atlas TransMethylCFDx/resources/Atlas.Submission3_1_final_full.l4.hg19.tsv -o <OUTPUT_NAME.csv> </path/to/resultsfolder>/results/pat_files/*pat.gz
```

> **Note**
The --rlen flag sets the number of CpGs that need to be on a single read for the read to be used for the deconvolution. --rlen 4 show the best LOD during validation.

### dd-cfDNA Calculation
The dd-cfDNA_calculation.py script in /scripts calculates dd-cfDNA from the pileup files saved to the 8_dd-cfDNA directory.

#### Dependencies
- Python 3.10 or higher
- pandas
- matplotlib
- scikit-learn

The dependencies are all available in the twist_TOO environment that was created to run the Snakemake pipeline.

#### Usage
Run the script with the required arguments as follows:
```
python dd-cfDNA_calculation.py -dir <input_directory> -o <output_file>
```
Use -h or --help for more information on more options.

## Output Folder Structure

Below is the structure of the output folder:

``` 
path/to/resultsfolder/
│  
├── results/
│   ├── 1_fastqc_init
│   ├── 2_trimmed
│   ├── 3_fastqc_trimmed
│   ├── 4_alignment
│   │   ├── bam
│   │   └── bam_filt
│   ├── 5_dedup_consensus
│   ├── 6_alignment_QC
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

## Cite
If you use this pipeline, please cite: \
***Kueng, N., Sandberg, F., Sidler, D., Banz, V., Berzigotti, A., Ng, C. K. Y., Largiader, C. R., & Amstutz, U. (2025). Integrated targeted deep sequencing reveals unique tissue-of-origin and donor-derived cell-free DNA signatures in organ transplant recipients. MedRxiv, 2025.04.29.25326125. doi:10.1101/2025.04.29.25326125.***

