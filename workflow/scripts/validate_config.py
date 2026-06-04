#!/usr/bin/env python3
"""
Configuration validator and processor for KIR-BLOOM.
Handles user-specified references, auto-generates KIR regions, and validates required files.
"""

import yaml
import sys
import os
import subprocess
from pathlib import Path

def load_config(config_path):
    """Load YAML configuration."""
    try:
        with open(config_path) as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"ERROR: Failed to load config: {e}")
        sys.exit(1)


def validate_file_exists(path, label, required=True):
    """Validate that a file exists."""
    if path is None:
        if required:
            print(f"ERROR: {label} is required but not specified")
            return False
        return True

    if not os.path.exists(path):
        print(f"ERROR: {label} not found: {path}")
        return False

    return True


def generate_kir_regions_if_needed(config):
    """Generate KIR regions BED file if not specified."""
    ref = config.get('Reference', {})
    fasta = ref.get('fasta')
    kir_bed = ref.get('KIR_regions_bed')

    if not fasta:
        print("ERROR: Reference fasta is required")
        return False

    if kir_bed and kir_bed != "auto":
        # User specified custom BED file
        if not validate_file_exists(kir_bed, "KIR_regions_bed"):
            return False
        print(f"Using user-specified KIR regions: {kir_bed}")
        return True

    # Auto-generate KIR regions
    print("Auto-generating KIR regions from reference...")

    script_path = os.path.join(os.path.dirname(__file__), "generate_kir_regions.py")
    output_bed = os.path.join(os.path.dirname(fasta), "kir_regions_auto.bed")

    try:
        subprocess.run(
            [sys.executable, script_path, fasta, output_bed],
            check=True
        )
        config['Reference']['KIR_regions_bed'] = output_bed
        print(f"Auto-generated KIR regions: {output_bed}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to generate KIR regions: {e}")
        return False
    except FileNotFoundError:
        print("ERROR: generate_kir_regions.py not found")
        return False


def validate_config(config_path):
    """
    Load and validate configuration.
    Auto-generates KIR regions if needed.
    Returns processed config or None if validation fails.
    """
    config = load_config(config_path)

    print("Validating configuration...")

    # Validate Reference section
    ref = config.get('Reference', {})

    required_refs = [
        ('fasta', 'Reference genome FASTA'),
        ('fasta_2', 'KIR reference FASTA'),
        ('KIR_alleles', 'KIR alleles database'),
        ('KIR_alleles_bed', 'KIR alleles BED'),
        ('exon_bed', 'Exon annotations'),
        ('intron_bed', 'Intron annotations'),
        ('utr_bed', 'UTR annotations'),
        ('msa_path', 'MSA path'),
        ('msa_mappings_path', 'MSA mappings path'),
        ('Allele_list', 'Allele list'),
        ('rep', 'Allele representatives'),
    ]

    valid = True
    for key, label in required_refs:
        if not validate_file_exists(ref.get(key), label, required=True):
            valid = False

    # Handle KIR regions (auto-generate if needed)
    if not generate_kir_regions_if_needed(config):
        valid = False

    # Validate Samples section
    samples = config.get('Samples', {})
    samples_dir = samples.get('samples_dir')

    if not samples_dir:
        print("ERROR: Samples.samples_dir is required")
        valid = False
    else:
        # Create samples directory if it doesn't exist
        os.makedirs(samples_dir, exist_ok=True)
        print(f"Output directory: {samples_dir}")

    # Check sample specification method
    sample_fastqs = samples.get('sample_fastqs')
    sample_csv = samples.get('sample_csv')

    if not sample_fastqs and not sample_csv:
        print("ERROR: Must specify either sample_fastqs or sample_csv")
        valid = False

    if sample_csv and not validate_file_exists(sample_csv, "Sample CSV"):
        valid = False

    # Validate tools (optional if in PATH)
    bwa = config.get('bwa_path')
    immuannot = config.get('immuannot')

    if bwa and not validate_file_exists(bwa, "BWA executable"):
        print("WARNING: BWA path specified but not found")

    if immuannot and not validate_file_exists(immuannot, "Immuannot script"):
        print("WARNING: Immuannot path specified but not found")

    if valid:
        print("✓ Configuration validation passed")
        return config
    else:
        print("✗ Configuration validation failed")
        return None


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: validate_config.py <config.yaml>")
        sys.exit(1)

    result = validate_config(sys.argv[1])
    sys.exit(0 if result else 1)

