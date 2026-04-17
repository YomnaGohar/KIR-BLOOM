import pysam
import matplotlib.pyplot as plt
import  numpy as np
from collections import defaultdict
import networkx as nx
from itertools import combinations
from collections import Counter
from Bio import SeqIO
import pandas as pd
import sys
def is_clipping_acceptable(read, max_clip_fraction=0.10):
    if read.cigartuples is None:
        return False
    clip_bases = 0
    for op, length in read.cigartuples:
        if op in {4, 5}:
            clip_bases += length
    read_length = 150
    return (clip_bases / read_length) <= max_clip_fraction
bam_path=snakemake.input.bam#sys.argv[4]#
def estimate_insert_size_from_Chr17q25(bam_path):
    insert_sizes = []
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for read in bamfile.fetch(until_eof=True):
            if (not read.is_unmapped and not read.mate_is_unmapped and read.is_proper_pair and read.reference_name.startswith(snakemake.params.chr17)):
                insert_size = abs(read.template_length)
                if insert_size > 0:
                    insert_sizes.append(insert_size)
    return insert_sizes
def tag_read_pairs_with_unique_id(bam_path, output_bam_path, tag_name="ZP"):
    tag_counter = 1
    read_buffer = defaultdict(dict)
    with pysam.AlignmentFile(bam_path, "rb") as infile, \
         pysam.AlignmentFile(output_bam_path, "wb", template=infile) as outfile:
        for read in infile.fetch(until_eof=True):
            if read.is_unmapped or read.mate_is_unmapped or not read.is_proper_pair:
                outfile.write(read)
                continue
            qname = read.query_name
            if read.is_read1:
                read_buffer[qname]['read1'] = read
            else:
                read_buffer[qname]['read2'] = read
            if 'read1' in read_buffer[qname] and 'read2' in read_buffer[qname]:
                read1 = read_buffer[qname]['read1']
                read2 = read_buffer[qname]['read2']
                read1.set_tag(tag_name, tag_counter, value_type='i')
                read2.set_tag(tag_name, tag_counter, value_type='i')
                outfile.write(read1)
                outfile.write(read2)
                del read_buffer[qname]
                tag_counter += 1
        for pair in read_buffer.values():
            for r in pair.values():
                outfile.write(r)
tag_read_pairs_with_unique_id(bam_path,snakemake.output.bam2, tag_name="ZP")    
insert_sizes=estimate_insert_size_from_Chr17q25(bam_path)
data = np.array(insert_sizes)
lower = np.percentile(data, 5)
upper = np.percentile(data, 95)
filtered_data = data[(data >= lower) & (data <= upper)]
mean_insert = np.mean(filtered_data)
std_insert = np.std(filtered_data, ddof=1) 
clippe=0
bam1_path=snakemake.input.bam1#sys.argv[1]#
bam2_path=snakemake.input.bam2#sys.argv[2]#
def collect_read_info(bam_path):
    read_info_dict = defaultdict(list)
    clipped = 0  # Counter for excluded reads
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not read.is_unmapped and (read.reference_name.startswith("KIR") or read.reference_name.startswith("IAG")): # 
                    contig = read.reference_name
                    start = read.reference_start
                    stop = read.reference_end
                    strand = '-' if read.is_reverse else '+'
                    mismatches = read.get_tag('NM') if read.has_tag('NM') else 0
                    cigar = read.cigarstring
                    seq = read.query_sequence
                    qual = read.query_qualities
                    name = read.query_name
                    name = name.removesuffix("/1").removesuffix("/2")
                    read_info_dict[name].append({
                        "contig": contig,
                        "start": start,
                        "stop": stop,
                        "strand": strand,
                        "NM": mismatches,
                        "cigar": cigar,
                        "seq": seq,
                        "qual": qual,
                        "clipping":is_clipping_acceptable(read, max_clip_fraction=0.10)
                    })


    return read_info_dict
 
    
original_seqs = {rec.id: str(rec.seq).replace("-", "") for rec in SeqIO.parse(snakemake.input.fasta, "fasta")}

def start_or_end(read_info, original_seqs):
    start = read_info["start"]
    contig_len = len(original_seqs[read_info["contig"]])
    end = contig_len - start
    return "beginning" if end > start else "end"
    
    
read1_info = collect_read_info(bam1_path)
read2_info = collect_read_info(bam2_path)

paired_reads = []
for qname in read1_info:
    r1_entries = read1_info[qname]
    if qname in read2_info:
        r2_entries = read2_info[qname]
        for r1 in r1_entries:
            for r2 in r2_entries:
                if r1["contig"] != r2["contig"]:
                    continue
                if r1["strand"] == r2["strand"]:
                    continue
                r1_5p = r1["start"] if r1["strand"] == '+' else r1["stop"] - 1
                r2_5p = r2["start"] if r2["strand"] == '+' else r2["stop"] - 1
                inward = (
                    (r1["strand"] == '+' and r2["strand"] == '-' and r1_5p < r2_5p) or
                    (r1["strand"] == '-' and r2["strand"] == '+' and r2_5p < r1_5p))
                if not inward:
                    continue
                insert_size = max(r1["stop"], r2["stop"]) - min(r1["start"], r2["start"])
                if not (mean_insert - 4 * std_insert <= insert_size <= mean_insert + 4 * std_insert):
                    continue
                paired_reads.append((qname, r1, r2))

with pysam.AlignmentFile(bam1_path, "rb") as bam_in:
    header = bam_in.header
output_bam=snakemake.output.bam1#sys.argv[3]#
with pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:
    tag_counter = 1 
    for qname, r1, r2 in paired_reads:
        if r1 and r2:
            a1 = pysam.AlignedSegment()
            a1.query_name = qname
            a1.reference_id = bam_out.get_tid(r1["contig"])
            a1.reference_start = r1["start"]
            a1.cigarstring = r1["cigar"]
            a1.query_sequence = r1["seq"]
            a1.query_qualities = r1["qual"]
            a1.flag = 99 if r1["strand"] == '+' else 83  # Adjust if needed
            a1.set_tag("NM", r1["NM"])
            tag = f"TAG{tag_counter}"
            a1.set_tag("TP", tag) 
            
            a2 = pysam.AlignedSegment()
            a2.query_name = qname
            a2.reference_id = bam_out.get_tid(r2["contig"])
            a2.reference_start = r2["start"]
            a2.cigarstring = r2["cigar"]
            a2.query_sequence = r2["seq"]
            a2.query_qualities = r2["qual"]
            a2.flag = 147 if r2["strand"] == '-' else 163  # Adjust if needed
            a2.set_tag("NM", r2["NM"])
            a1.is_read1          = True
            a1.is_read2          = False
            a1.is_paired = True
            a2.is_paired = True
            a1.is_proper_pair = True
            a2.is_proper_pair = True
            a1.mate_is_unmapped = False
            a2.mate_is_unmapped = False
            a1.next_reference_id = a2.reference_id
            a2.next_reference_id = a1.reference_id
            a1.next_reference_start = a2.reference_start
            a2.next_reference_start = a1.reference_start
            a2.is_read1          = False
            a2.is_read2          = True    
            tlen = max(r1["stop"], r2["stop"]) - min(r1["start"], r2["start"])
            a1.template_length = tlen if a1.reference_start < a2.reference_start else -tlen
            a2.template_length = -a1.template_length
            a2.set_tag("TP", tag)
            bam_out.write(a1)
            bam_out.write(a2)
        elif r1:    
            a1 = pysam.AlignedSegment()
            a1.query_name = qname
            a1.reference_id = bam_out.get_tid(r1["contig"])
            a1.reference_start = r1["start"]
            a1.cigarstring = r1["cigar"]
            a1.query_sequence = r1["seq"]
            a1.query_qualities = r1["qual"]
            a1.flag = 0 if r1["strand"] == '+' else 16  # unpaired strand flags
            a1.set_tag("NM", r1["NM"])
            a1.is_read1          = True
            a1.is_read2          = False 
            a1.is_paired = True
            a1.is_proper_pair = False
            a1.mate_is_unmapped = True
            a1.template_length = 0  # No mate, so no TLEN
            tag = f"TAG{tag_counter}"
            a1.set_tag("TP", tag)
            bam_out.write(a1)
        elif r2:    
            a1 = pysam.AlignedSegment()
            a1.query_name = qname
            a1.reference_id = bam_out.get_tid(r2["contig"])
            a1.reference_start = r2["start"]
            a1.cigarstring = r2["cigar"]
            a1.query_sequence = r2["seq"]
            a1.query_qualities = r2["qual"]
            a1.flag = 0 if r2["strand"] == '+' else 16  # unpaired strand flags
            a1.set_tag("NM", r2["NM"])
            a1.is_read1          = False
            a1.is_read2          = True 
            a1.is_paired = True
            a1.is_proper_pair = False
            a1.mate_is_unmapped = True
            a1.template_length = 0  # No mate, so no TLEN
            tag = f"TAG{tag_counter}"
            a1.set_tag("TP", tag)
            bam_out.write(a1)    
        tag_counter += 1
