![KIR*bloom logo](logo2.jpeg)

# KIR-BLOOM: KIR Allele Genotyping from Short-Read Sequencing

A pipeline for accurate KIR allele genotyping from short-read sequencing data.

---

## Getting Started: Clone the Repository

First, clone the KIR-BLOOM repository:

```bash
git clone https://github.com/YomnaGohar/KIR-BLOOM.git
cd KIR-BLOOM
```

This will download all the necessary pipeline files, including the Snakemake workflow and configuration templates.

---

## For First-Time Users

### Step 0: Download Reference Data

The KIR-BLOOM pipeline requires some reference files that must be placed in the `resources/` directory.
note that the reference files and indexes require approximately 60 GB of disk space
**Recommended: Use the automatic download script**

```bash
bash workflow/scripts/download_references.sh
```

This script will:
- Download all required FASTA files from Zenodo
- Extract and organize them in `resources/`
- Verify the download was successful

**Data Source**: https://zenodo.org/records/20532512

If the automatic script doesn't work, download manually from the Zenodo link above, then extract to the resources directory:

```bash
gunzip -c GRCh38_full_analysis_set_plus_decoy_hla.fa.gz > resources/GRCh38_full_analysis_set_plus_decoy_hla.fa
gunzip -c ref_with_utr_extended.fa.gz > resources/ref_with_utr_extended.fa
gunzip -c kir_gen_new_mod_with_utr_extended.fasta.gz > resources/kir_gen_new_mod_with_utr_extended.fasta
tar -xf msa_after_utr_extension.tar -C resources/
```

**Verify files are in place**:

```bash
ls -lh resources/
# Should show:
# GRCh38_full_analysis_set_plus_decoy_hla.fa
# ref_with_utr_extended.fa
# kir_gen_new_mod_with_utr_extended.fasta
# kir_allele_names_new_kir_only.bed
# annotations_mod_with_utr_extended_exon.bed
# annotations_mod_with_utr_extended_intron.bed
# annotations_mod_with_utr_extended_utr.bed
# Allelelist.txt
# kir_regions.bed
# allele_representatives_new_kir.pkl
# msa_after_utr_extension/
```

---

## Prerequisites Installation

### 1. Install Snakemake

```bash
# Using conda (recommended)
conda install -c bioconda snakemake=7.32

# Or using pip
pip install snakemake==7.32
```

### 2. Install BWA-MEM2

```bash
# Download pre-compiled binary
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xjf bwa-mem2-2.2.1_x64-linux.tar.bz2
# Add to PATH or note the path for config
```
### 3. Install Minimap2
conda install -c bioconda minimap2=2.28

### 4. Install Samtools
conda install -c bioconda samtools

### other Requirements

- **Python 3.9+**

---

## Quick Start Guide


### Step 1: Prepare Your Data

**FASTQ Requirements**:
- Paired-end reads
- Read 1 must end with `/1` and Read 2 with `/2`

### Step 2: Create Configuration File

**Create your own config file by copying from a template:**

```bash
cp config/example_config.yaml config/myconfig.yaml
```

Then edit `config/myconfig.yaml` with your paths.

### Step 3: Index Reference Files

Index the reference FASTA files and create BWA indices

```bash
# Inside KIR-BlOOM directory run 
snakemake --cores 4 all_index --configfile config/myconfig.yaml
```

**Verify indices were created**:
```bash
ls -lh resources/*.fai resources/*.amb
```


### Step 4: Preprocess Reads and Generate Pairwise Alignments

This step extracts reads mapping to the KIR region and the configured background region, retains their read pairs, and generates all possible pairwise alignments against the KIR allele reference set.


```bash
snakemake --cores 10  extract config/myconfig.yaml 
```
note that `--cores` specifies the maximum number of CPU cores available to the workflow, while `threads` in the configuration file specifies the maximum number of threads that can be used by an individual rule. If a rule requests more threads than are available through `--cores`, Snakemake will automatically scale the thread count down to the maximum number of available cores.
### Output Files

| File | Description |
|------|-------------|
| `mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam` | Contains reads mapping to background region. Used to estimate sequencing depth and error rates from the background region. |
| `paired_new_kir_sort_all4.bam` | Contains the corresponding read pairs on KIR allele reference sequences. |

## Step 5: Filter Candidate Alleles

This step filters for candidate alleles to reduce the solution space. These candidate alleles are carried forward to the copy number estimation and allele selection stages.


```bash
snakemake --cores 10 filter --configfile config/myconfig.yaml
```

### Output Files

| File | Description |
|------|-------------|
| `selected.txt` | Contains the candidate alleles selected after filtering. |
| `alleles_scored_and_grouped_by_genes.pkl` | scoring of alleles from EM algorithm which is used during CN estimation. |

## Step 6: Copy Number Inference

This step estimates the copy number for each KIR gene using the filtered candidate alleles and the background-based depth information.

To run the copy number inference step:

```bash
snakemake --cores 10 cn --configfile config/myconfig.yaml
```

### Output Files

| File | Description |
|------|-------------|
| `cn.tsv` | Contains the genes and their estimated copy numbers. |
| `coverage_plot_after_cn_inference.pdf` | coverage plot of the alleles selected during copy number inference step. Note that those alleles might changes after allele inference step. |

## Step 7: Five digit Allele Inference and variant calling

This step uses the previously inferred copy number as a constraint and selects the combination of alleles that best explains the observed sequencing data.

> **Note**
>
> KIR*BLOOM reports full-length allele sequences. However, the most reliable part of the prediction corresponds to the coding sequence (CDS). Intronic and other non-coding regions are not explicitly considered during allele selection and should therefore be interpreted with caution.


```bash
snakemake --cores 10 infer --configfile config/myconfig.yaml
```

### Output Files

| File | Description |
|------|-------------|
| `five_digit_allele_inference.tsv` | Contains the inferred alleles and their copy numbers prior to variant calling. |
| `five_digit_allele_inference.pdf` | Coverage plots for the alleles selected during the 5-digit allele inference step. |
| `kir_variants_in_exons.vcf` | VCF file containing the predicted variants in exonic regions. |
| `kir_mod.fa` | Modified allele sequences after incorporating the predicted variants. |
| `kir_mod_exon.bed` | BED file containing exon coordinates after incorporating the predicted variants. |



