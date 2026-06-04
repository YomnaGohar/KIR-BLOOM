"""
Generate KIR regions BED file from user-specified reference.
This rule should be included first in the Snakefile.
"""

# Get paths from config
REFERENCE_FASTA = config["Reference"]["fasta"]
SAMPLES_DIR = config["Samples"]["samples_dir"]
KIR_REGIONS_BED = os.path.join(SAMPLES_DIR, "kir_regions.bed")

# Optional: path to GRCh38.fa.fai (if not in default location)
GRCHR38_FAI = config["Reference"].get("grchr38_fai", None)

rule generate_kir_regions:
    """
    Generate KIR regions BED file by matching user reference sequences
    against GRCh38 KIR region sequences by length.
    """
    input:
        reference = REFERENCE_FASTA,
    output:
        bed = KIR_REGIONS_BED,
    params:
        script = workflow.source_path("../scripts/generate_kir_regions.py"),
        grchr38_fai = GRCHR38_FAI,
    log:
        "logs/generate_kir_regions.log"
    shell:
        """
        if [ -z "{params.grchr38_fai}" ] || [ "{params.grchr38_fai}" = "None" ]; then
            # Use default location (resources/GRCh38.fa.fai)
            python {params.script} {input.reference} {output.bed} > {log} 2>&1
        else
            # Use specified GRCh38 FAI path
            python {params.script} {input.reference} {output.bed} {params.grchr38_fai} > {log} 2>&1
        fi
        """

# Update config with generated BED path so downstream rules can use it
config["Reference"]["KIR_regions_bed"] = KIR_REGIONS_BED

