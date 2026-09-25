"""
Single pass over the background (Chr17) BAM:
  1. Assigns a unique ZP tag to every proper-pair read name
  2. Writes the tagged BAM to output.tagged_bam
  3. Computes insert-size statistics and saves them to output.stats (JSON)
"""
import json
import numpy as np
import pysam

n_threads    = snakemake.threads
bam_path     = snakemake.input.bam
chr17_prefix = snakemake.params.chr17

# ── Pass 1: collect insert sizes + build qname→tag map ───────────────────────
qname_to_tag = {}
tag_counter  = 1
insert_sizes = []

with pysam.AlignmentFile(bam_path, "rb", threads=n_threads) as infile:
    for read in infile.fetch(until_eof=True):
        if (not read.is_unmapped and not read.mate_is_unmapped
                and read.is_proper_pair
                and read.reference_name.startswith(chr17_prefix)):
            size = abs(read.template_length)
            if size > 0:
                insert_sizes.append(size)
        if not (read.is_unmapped or read.mate_is_unmapped or not read.is_proper_pair):
            qname = read.query_name
            if qname not in qname_to_tag:
                qname_to_tag[qname] = tag_counter
                tag_counter += 1

# ── Compute and save insert stats ────────────────────────────────────────────
data          = np.array(insert_sizes)
lower         = np.percentile(data, 5)
upper         = np.percentile(data, 95)
filtered      = data[(data >= lower) & (data <= upper)]
mean_insert   = float(np.mean(filtered))
std_insert    = float(np.std(filtered, ddof=1))

with open(snakemake.output.stats, "w") as f:
    json.dump({"mean_insert": mean_insert, "std_insert": std_insert}, f)

# ── Pass 2: write tagged BAM ──────────────────────────────────────────────────
with pysam.AlignmentFile(bam_path, "rb", threads=n_threads) as infile, \
     pysam.AlignmentFile(snakemake.output.tagged_bam, "wb", template=infile) as outfile:
    for read in infile.fetch(until_eof=True):
        if read.query_name in qname_to_tag:
            read.set_tag("ZP", qname_to_tag[read.query_name], value_type='i')
        outfile.write(read)
