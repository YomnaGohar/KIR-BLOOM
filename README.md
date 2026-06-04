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

Index the reference FASTA files and create BWA indices (replace `<your_config>` with your config file):

```bash
# Inside KIR-BlOOM directory run 
snakemake --cores 4 all_index --configfile config/myconfig.yaml
```

**Verify indices were created**:
```bash
ls -lh resources/*.fai resources/*.amb
```


### Step 4: Run the Pipeline

#### Dry Run

```bash
snakemake --np extract --configfile config/myconfig.yaml 

```

#### Full Run

```bash
snakemake --configfile config/<your_config>.yaml -j 4
snakemake --configfile config/<your_config>.yaml -j 10 --resources mem_mb=32000
```
