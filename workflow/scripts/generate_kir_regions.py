#!/usr/bin/env python3
"""
Generate KIR regions BED file from user-specified reference.
Indexes reference and auto-detects KIR regions based on sequence lengths.
"""

import subprocess
import sys
import os

# Sequences to remove (alternate chr19 representations)
SEQUENCES_TO_REMOVE = {
    "GL000209.2", "NT_113949.2",
    "GL949746.1", "NW_003571054.1",
    "GL949747.2", "NW_003571055.2",
    "GL949748.2", "NW_003571056.2",
    "GL949749.2", "NW_003571057.2",
    "GL949750.2", "NW_003571058.2",
    "GL949751.2", "NW_003571059.2",
    "GL949752.1", "NW_003571060.1",
    "GL949753.2", "NW_003571061.2",
    "KI270882.1", "NT_187636.1",
    "KI270883.1", "NT_187637.1",
    "KI270884.1", "NT_187638.1",
    "KI270885.1", "NT_187639.1",
    "KI270886.1", "NT_187640.1",
    "KI270887.1", "NT_187641.1",
    "KI270888.1", "NT_187642.1",
    "KI270889.1", "NT_187643.1",
    "KI270930.1", "NT_187684.1",
    "KI270931.1", "NT_187685.1",
    "KI270932.1", "NT_187686.1",
    "KI270933.1", "NT_187687.1",
    "KI270938.1", "NT_187693.1",
    "KV575246.1", "NW_016107300.1",
    "KV575247.1", "NW_016107301.1",
    "KV575248.1", "NW_016107302.1",
    "KV575249.1", "NW_016107303.1",
    "KV575250.1", "NW_016107304.1",
    "KV575251.1", "NW_016107305.1",
    "KV575252.1", "NW_016107306.1",
    "KV575253.1", "NW_016107307.1",
    "KV575254.1", "NW_016107308.1",
    "KV575255.1", "NW_016107309.1",
    "KV575256.1", "NW_016107310.1",
    "KV575257.1", "NW_016107311.1",
    "KV575258.1", "NW_016107312.1",
    "KV575259.1", "NW_016107313.1",
    "KV575260.1", "NW_016107314.1",
    "KI270890.1", "NT_187644.1",
    "KI270891.1", "NT_187645.1",
    "KI270914.1", "NT_187668.1",
    "KI270915.1", "NT_187669.1",
    "KI270916.1", "NT_187670.1",
    "KI270917.1", "NT_187671.1",
    "KI270918.1", "NT_187672.1",
    "KI270919.1", "NT_187673.1",
    "KI270920.1", "NT_187674.1",
    "KI270921.1", "NT_187675.1",
    "KI270922.1", "NT_187676.1",
    "KI270923.1", "NT_187677.1",
    "KI270929.1", "NT_187683.1"
}


def index_reference(fasta_path):
    """Create .fai index for reference FASTA if it doesn't exist."""
    fai_path = f"{fasta_path}.fai"

    if os.path.exists(fai_path):
        print(f"Index already exists: {fai_path}")
        return fai_path

    print(f"Creating index for {fasta_path}...")
    try:
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
        print(f"Index created: {fai_path}")
        return fai_path
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to index reference: {e}")
        sys.exit(1)
    except FileNotFoundError:
        print("ERROR: samtools not found. Please install samtools.")
        sys.exit(1)


def parse_fai(fai_path):
    """Parse FAI file and return dict of {sequence_name: length}."""
    sequences = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                name, length = parts[0], int(parts[1])
                sequences[name] = length
    return sequences


def build_reference_length_mapping(grchr38_fai_path):
    """
    Parse GRCh38.fa.fai and extract lengths for sequences in SEQUENCES_TO_REMOVE.
    Returns dict: {length: set of sequence names}
    This allows matching by length across different reference naming schemes.
    """
    reference_lengths = {}

    if not os.path.exists(grchr38_fai_path):
        print(f"ERROR: Reference FAI file not found: {grchr38_fai_path}")
        sys.exit(1)

    print(f"Parsing reference sequences from: {grchr38_fai_path}")

    with open(grchr38_fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                name = parts[0]
                length = int(parts[1])

                # Only include sequences that are in SEQUENCES_TO_REMOVE
                if name in SEQUENCES_TO_REMOVE:
                    if length not in reference_lengths:
                        reference_lengths[length] = set()
                    reference_lengths[length].add(name)

    print(f"Found {len(reference_lengths)} unique lengths for KIR regions")
    return reference_lengths


def generate_kir_regions_bed(fasta_path, output_bed_path, grchr38_fai_path=None):
    """
    Auto-generate KIR regions BED file from reference.
    Matches user sequences by length against GRCh38 reference sequences.

    Args:
        fasta_path: Path to user's reference FASTA
        output_bed_path: Output BED file path
        grchr38_fai_path: Path to GRCh38.fa.fai (if None, uses default location)
    """

    # Determine GRCh38 FAI path
    if grchr38_fai_path is None:
        # Look for it in resources directory relative to script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(os.path.dirname(script_dir))
        grchr38_fai_path = os.path.join(project_root, "resources", "GRCh38.fa.fai")

    # Index user reference
    fai_path = index_reference(fasta_path)
    user_sequences = parse_fai(fai_path)  # {name: length}

    # Parse GRCh38 reference lengths
    ref_lengths_map = build_reference_length_mapping(grchr38_fai_path)  # {length: set of names}

    # Create inverse mapping: {length: name} for user reference
    user_lengths = {}
    for user_name, length in user_sequences.items():
        user_lengths[length] = user_name

    print(f"Found {len(user_sequences)} sequences in user reference")
    print(f"Matching against {sum(len(names) for names in ref_lengths_map.values())} known KIR sequences")

    # Write BED file
    matched_count = 0
    matched_sequences = set()

    with open(output_bed_path, "w") as out:
        # Add chr19 region (standard KIR region)
        out.write("chr19\t52025634\t57084318\n")
        matched_count += 1

        # Match user sequences by length
        for length, ref_names in ref_lengths_map.items():
            # Check if this length exists in user reference
            if length in user_lengths:
                user_name = user_lengths[length]
                # Only add if not already added
                if user_name not in matched_sequences:
                    out.write(f"{user_name}\t0\t{length}\n")
                    matched_count += 1
                    matched_sequences.add(user_name)

    print(f"Generated KIR regions BED: {output_bed_path}")
    print(f"Included {matched_count} sequences (1 chr19 + {matched_count - 1} alts)")
    print(f"Matched {len(matched_sequences)} alternate chr19 representations by length")
    return output_bed_path


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: generate_kir_regions.py <reference.fa> [output.bed] [GRCh38.fa.fai]")
        sys.exit(1)

    fasta_path = sys.argv[1]

    if len(sys.argv) > 2:
        output_bed = sys.argv[2]
    else:
        # Remove extension and add _kir_regions.bed
        base_name = fasta_path.rsplit('.', 1)[0]
        output_bed = f"{base_name}_kir_regions.bed"

    grchr38_fai = sys.argv[3] if len(sys.argv) > 3 else None

    if not os.path.exists(fasta_path):
        print(f"ERROR: Reference file not found: {fasta_path}")
        sys.exit(1)

    generate_kir_regions_bed(fasta_path, output_bed, grchr38_fai)

