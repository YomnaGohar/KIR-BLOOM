#!/bin/bash
# KIR-BLOOM Reference Data Download Script
# This script downloads large reference files needed for KIR-BLOOM
#
# Usage: bash workflow/scripts/download_references.sh
#
# NOTE: Update the URLs below with your actual hosting location
# (Zenodo, Figshare, institutional server, Google Drive, etc.)

set -e  # Exit on error

RESOURCES_DIR="resources"
# Create directories if they don't exist
mkdir -p "$RESOURCES_DIR"

echo "=================================="
echo "KIR-BLOOM Reference Download Tool"
echo "=================================="
echo ""

# Using Zenodo dataset at https://zenodo.org/records/20532512
# These direct file URLs point to the files uploaded to that record.
FASTA_2_URL="https://zenodo.org/record/20532512/files/ref_with_utr_extended.fa.gz"
KIR_ALLELES_URL="https://zenodo.org/record/20532512/files/kir_gen_new_mod_with_utr_extended.fasta.gz"
GRCH38_URL="https://zenodo.org/record/20532512/files/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"
MSA_URL="https://zenodo.org/record/20532512/files/msa_after_utr_extension.tar"

# Function to download and verify file
download_file() {
    local url=$1
    local output_path=$2
    local description=$3

    if [ -f "$output_path" ] || [ -f "${output_path}.gz" ]; then
        echo "✓ $description already exists: $output_path"
        return 0
    fi

    echo "⬇ Downloading $description..."
    echo "   URL: $url"

    if command -v wget &> /dev/null; then
        wget -q --show-progress "$url" -O "$output_path.tmp"
    elif command -v curl &> /dev/null; then
        curl -# -L "$url" -o "$output_path.tmp"
    else
        echo "❌ Error: Neither wget nor curl is available"
        exit 1
    fi

    if [ $? -eq 0 ]; then
        mv "$output_path.tmp" "$output_path"
        echo "✓ Downloaded: $output_path"
    else
        echo "❌ Failed to download: $description"
        rm -f "$output_path.tmp"
        exit 1
    fi
}

# Function to decompress gzipped files
decompress_if_needed() {
    local file=$1
    if [[ "$file" == *.gz ]]; then
        echo "  Decompressing $file..."
        gunzip -f "$file"
        echo "  ✓ Done"
    fi
}

# Download each file
echo ""
echo "Starting downloads (this may take a while for large files)..."
echo ""

# 1. KIR reference FASTA
download_file "$FASTA_2_URL" "$RESOURCES_DIR/ref_with_utr_extended.fa.gz" "KIR reference sequence"
decompress_if_needed "$RESOURCES_DIR/ref_with_utr_extended.fa.gz"

# 2. KIR alleles FASTA
download_file "$KIR_ALLELES_URL" "$RESOURCES_DIR/kir_gen_new_mod_with_utr_extended.fasta.gz" "KIR alleles database"
decompress_if_needed "$RESOURCES_DIR/kir_gen_new_mod_with_utr_extended.fasta.gz"

# 3. GRCh38 reference
download_file "$GRCH38_URL" "$RESOURCES_DIR/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz" "GRCh38 reference genome"
decompress_if_needed "$RESOURCES_DIR/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"

# 4. MSA archive
if [ ! -d "$RESOURCES_DIR/msa_after_utr_extension" ]; then
    echo ""
    echo "⬇ Downloading MSA data..."
    download_file "$MSA_URL" "$RESOURCES_DIR/msa_after_utr_extension.tar" "MSA archive"

    echo "  Extracting MSA archive..."
    tar -xf "$RESOURCES_DIR/msa_after_utr_extension.tar" -C "$RESOURCES_DIR"
    rm "$RESOURCES_DIR/msa_after_utr_extension.tar"
    echo "  ✓ MSA extracted"
else
    echo "✓ MSA directory already exists"
fi

echo ""
echo "=================================="
echo "✓ Download complete!"
echo "=================================="
