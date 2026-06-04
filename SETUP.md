# KIR-BLOOM Setup Guide for New Users

This guide walks you through setting up and running KIR-BLOOM from scratch.

## Prerequisites Installation

### 1. Install Snakemake

```bash
# Using conda (recommended)
conda install -c bioconda snakemake

# Or using pip
pip install snakemake
```

### 2. Install BWA-MEM2

```bash
# Download pre-compiled binary
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xjf bwa-mem2-2.2.1_x64-linux.tar.bz2
# Add to PATH or note the path for config

# Or compile from source
git clone https://github.com/bwa-mem2/bwa-mem2.git
cd bwa-mem2 && make
```

### 3. Get Immuannot

```bash
# Contact the Immuannot authors or download from their repository
# Note the path to immuannot.sh for your config
```

## Directory Setup

Create a clean project directory structure:

```bash
mkdir -p kir_project
cd kir_project

# Create subdirectories
mkdir -p resources data output config

# Your structure should look like:
# kir_project/
# ├── config/           (config files)
# ├── data/             (your FASTQ/BAM/CRAM files)
# ├── output/           (results will go here)
# ├── resources/        (reference files)
# └── workflow/         (copy from KIR-BLOOM)
```

## Step 1: Get Reference Files

All reference files must be placed in the `resources/` directory.

### Obtain Required Files

Contact the KIR-BLOOM developers or download from these sources:

```bash
cd resources

# 1. GRCh38 Reference with HLA and decoy
# Download from NCBI or your institution
# File: GRCh38_full_analysis_set_plus_decoy_hla.fa

# 2. KIR Reference Files (usually provided by KIR-BLOOM developers)
# These typically include:
mkdir -p kir_reference_2025
cd kir_reference_2025
# Download:
# - ref_with_utr_extended.fa
# - kir_gen_new_mod_with_utr_extended.fasta
# - kir_allele_names_new_kir_only.bed
# - annotations_mod_with_utr_extended_exon.bed
# - annotations_mod_with_utr_extended_intron.bed
# - annotations_mod_with_utr_extended_utr.bed
# - Allelelist.txt
# - allele_representatives_new_kir.pkl
# - msa_after_utr_extension/ (directory with MSA files)

cd ../..
```

### Verify Files Exist

```bash
ls -l resources/
# Should show:
# GRCh38_full_analysis_set_plus_decoy_hla.fa
# kir_reference_2025/

ls -l resources/kir_reference_2025/
# Should show all files mentioned above
```

## Step 2: Prepare Your Sequencing Data

Place FASTQ files in the `data/` directory:

```bash
# Copy paired-end FASTQ files
cp /your/data/sample_1_R1.fastq.gz data/
cp /your/data/sample_1_R2.fastq.gz data/
cp /your/data/sample_2_R1.fastq.gz data/
cp /your/data/sample_2_R2.fastq.gz data/

# Verify files
ls -lh data/
```

### FASTQ File Requirements

- **Paired-end**: Two files per sample (forward and reverse reads)
- **Naming**: Should distinguish between pairs. Examples:
  - `sample_R1.fastq.gz` and `sample_R2.fastq.gz`
  - `sample_1.fastq.gz` and `sample_2.fastq.gz`
  - (Snakemake will handle these patterns)
- **Format**: Can be gzipped (.gz) or uncompressed

## Step 3: Create Configuration

Copy the example config and customize it:

```bash
cp /path/to/KIR-BLOOM/config/example_config.yaml config/myconfig.yaml
```

Edit `config/myconfig.yaml`:

```yaml
# Edit these paths to match your setup:

Reference:
  fasta: "./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
  fasta_2: "./resources/kir_reference_2025/ref_with_utr_extended.fa"
  # ... (other fields already match if you followed directory structure)

Samples:
  samples_dir: "./output"
  sample_fastqs:
    "sample_1": ["./data/sample_1_R1.fastq.gz", "./data/sample_1_R2.fastq.gz"]
    "sample_2": ["./data/sample_2_R1.fastq.gz", "./data/sample_2_R2.fastq.gz"]

# If BWA or Immuannot are not in PATH:
bwa_path: "/path/to/bwa-mem2/bwa-mem2"
immuannot: "/path/to/Immuannot/scripts.pub.v3/immuannot.sh"
```

## Step 4: Validate Configuration

Before running the full pipeline, validate your config:

```bash
python /path/to/KIR-BLOOM/workflow/scripts/validate_config.py config/myconfig.yaml
```

This will check:
- All required files exist and are readable
- All paths are correct
- Configuration syntax is valid

**If validation passes**, you're ready to run!

## Step 5: Run the Pipeline

### Dry Run (Recommended First)

See what the pipeline will do without actually running:

```bash
snakemake --configfile config/myconfig.yaml --dryrun
```

This shows:
- All rules that will execute
- Input/output files
- Job scheduling

### Full Run

```bash
# Run with 4 threads
snakemake --configfile config/myconfig.yaml -j 4

# Run with more resources (e.g., 16 threads)
snakemake --configfile config/myconfig.yaml -j 16

# Run with detailed output (verbose)
snakemake --configfile config/myconfig.yaml -j 4 -v
```

### Run on HPC Clusters (Optional)

For SLURM-based clusters:

```bash
snakemake --configfile config/myconfig.yaml \
  --profile slurm \
  -j 100 \
  --cluster "sbatch -t {cluster.time} -c {cluster.threads}"
```

## Step 6: Check Results

After successful completion, results are in `output/`:

```bash
ls -l output/
# Should contain:
# sample_1/
# sample_2/
# ...

ls -l output/sample_1/
# Should contain:
# genotypes.txt (main results)
# alignments/ (BAM files)
# variants/ (VCF files)
```

## Troubleshooting

### Problem: Reference file not found

```
FileNotFoundError: ./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

**Solution**:
- Verify file exists: `ls -l resources/`
- Update path in config if file is in different location
- Check spelling and case sensitivity

### Problem: Validation fails

```
Error: validate_config.py - Missing required file
```

**Solution**:
1. Run validation to see what's missing:
   ```bash
   python workflow/scripts/validate_config.py config/myconfig.yaml
   ```
2. Download missing files to `resources/`
3. Update config with correct paths

### Problem: BWA-MEM not found

```
Error: [E::main] failed to execute bwa-mem2
```

**Solution**:
- Add to config: `bwa_path: "/full/path/to/bwa-mem2"`
- Or add to PATH: `export PATH="/path/to/bwa-mem2:$PATH"`

### Problem: Snakemake not found

```
command not found: snakemake
```

**Solution**:
```bash
# Activate conda environment if using conda
conda activate snakemake  # (or your environment name)

# Or install globally
pip install snakemake
```

### Problem: Out of memory

```
MemoryError or process killed
```

**Solution**:
- Increase memory allocation:
  ```bash
  snakemake --configfile config/myconfig.yaml -j 4 --resources mem_mb=64000
  ```
- Or update config:
  ```yaml
  parameters:
    memory_gb: 32
  ```

## Next Steps

- Read the main `README.md` for advanced options
- Check `config/config_template.yaml` for all available parameters
- Look at `workflow/rules/` to understand pipeline stages

## Getting Help

If you encounter issues:
1. Check validation output: `python workflow/scripts/validate_config.py config/myconfig.yaml`
2. Review Snakemake logs: `snakemake --configfile config/myconfig.yaml --debug`
3. Contact the KIR-BLOOM developers with:
   - Your config file (without sensitive paths)
   - Error messages
   - Output of validation script

