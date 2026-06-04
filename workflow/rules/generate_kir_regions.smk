"""
Prepare KIR regions BED file.
NOTE: KIR regions BED is generated once and placed in samples_dir root
since it's shared across all samples. Each sample gets its own subdirectory
created by other rules where sample-specific outputs are placed.

Directory structure:
  samples_dir/
    ├── kir_regions.bed (shared across all samples)
    ├── sample1/ (created by other rules)
    ├── sample2/ (created by other rules)
    └── ...

If user specifies KIR_regions_bed, copy it to samples_dir.
Otherwise, generate it from user reference by matching against GRCh38.
This rule should be included first in the Snakefile.
"""

import os

# Get paths from config
REFERENCE_FASTA = config["Reference"]["fasta"]
SAMPLES_DIR = config["Samples"]["samples_dir"]
KIR_REGIONS_BED_OUTPUT = os.path.join(SAMPLES_DIR, "kir_regions.bed")
KIR_REGIONS_BED_INPUT = config["Reference"].get("KIR_regions_bed", None)

# Optional: path to GRCh38.fa.fai (if not in default location)
GRCHR38_FAI = config["Reference"].get("grchr38_fai", None)

if KIR_REGIONS_BED_INPUT:
    # User specified a BED file - copy it
    rule prepare_kir_regions_copy:
        """
        Copy user-specified KIR regions BED file to samples_dir root.
        Shared across all samples.
        """
        input:
            bed = KIR_REGIONS_BED_INPUT,
        output:
            bed = KIR_REGIONS_BED_OUTPUT,
        log:
            "logs/prepare_kir_regions.log"
        shell:
            """
            mkdir -p $(dirname {output.bed})
            cp {input.bed} {output.bed}
            echo "Copied KIR regions BED from {input.bed}" > {log}
            """
else:
    # User did not specify - generate from reference
    rule prepare_kir_regions_generate:
        """
        Generate KIR regions BED file by matching user reference sequences
        against GRCh38 KIR region sequences by length.
        Placed in samples_dir root, shared across all samples.
        """
        input:
            reference = REFERENCE_FASTA,
        output:
            bed = KIR_REGIONS_BED_OUTPUT,
        params:
            script = workflow.source_path("../scripts/generate_kir_regions.py"),
            grchr38_fai = GRCHR38_FAI,
        log:
            "logs/prepare_kir_regions.log"
        shell:
            """
            mkdir -p $(dirname {output.bed})
            if [ -z "{params.grchr38_fai}" ] || [ "{params.grchr38_fai}" = "None" ]; then
                # Use default location (resources/GRCh38.fa.fai)
                python {params.script} {input.reference} {output.bed} > {log} 2>&1
            else
                # Use specified GRCh38 FAI path
                python {params.script} {input.reference} {output.bed} {params.grchr38_fai} > {log} 2>&1
            fi
            """

# Update config with output BED path so downstream rules can use it
config["Reference"]["KIR_regions_bed"] = KIR_REGIONS_BED_OUTPUT

