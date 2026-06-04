# KIR-BLOOM Reference Indexing Rules
# This rule file creates FASTA indices and BWA indices for reference sequences

import os

# Define paths to reference FASTA files
RESOURCES_DIR = "resources"
KIR_REF = os.path.join(RESOURCES_DIR, "ref_with_utr_extended.fa")
GRCH38_REF = os.path.join(RESOURCES_DIR, "GRCh38_full_analysis_set_plus_decoy_hla.fa")

# Get BWA path from config, default to "bwa-mem2" if not specified
BWA_PATH = config.get("bwa_path", "bwa-mem2")

# Target files
FASTA_INDEX_FILES = [
    f"{KIR_REF}.fai",
    f"{GRCH38_REF}.fai",
]

BWA_INDEX_FILES = [
    f"{KIR_REF}.amb",
    f"{KIR_REF}.ann",
    f"{KIR_REF}.bwt.2bit.64",
    f"{KIR_REF}.pac",
    f"{KIR_REF}.0123",
    f"{GRCH38_REF}.amb",
    f"{GRCH38_REF}.ann",
    f"{GRCH38_REF}.bwt.2bit.64",
    f"{GRCH38_REF}.pac",
    f"{GRCH38_REF}.0123",
]


rule all_index:
    input:
        FASTA_INDEX_FILES + BWA_INDEX_FILES


rule create_fasta_index:
    """
    Create FASTA index files using samtools for quick sequence lookup
    """
    input:
        fasta="{fasta}"
    output:
        index="{fasta}.fai"
    message:
        "Creating FASTA index for {input.fasta}..."
    shell:
        """
        samtools faidx {input.fasta}
        """


rule create_bwa_index:
    """
    Create BWA-MEM2 indices for sequence alignment

    This creates the following index files:
    - .amb: Ambiguous positions in reference
    - .ann: Index annotation file
    - .bwt.2bit.64: BWT index (2-bit encoded for 64-bit systems)
    - .pac: Packed sequence file
    - .0123: Alternative alphabet mapping
    """
    input:
        fasta="{fasta}"
    output:
        amb="{fasta}.amb",
        ann="{fasta}.ann",
        bwt="{fasta}.bwt.2bit.64",
        pac="{fasta}.pac",
        sa="{fasta}.0123"
    message:
        "Creating BWA-MEM2 index for {input.fasta}..."
    shell:
        """
        {BWA_PATH} index {input.fasta}
        """

