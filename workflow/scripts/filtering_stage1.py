#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 13:53:57 2025

@author: yomna
"""
import pysam
from collections import Counter
import re
import pickle
from collections import defaultdict
import pandas as pd
import numpy as np
from Bio import SeqIO
def get_paired_by_mate_info_fast(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    read1_cands = defaultdict(list)
    read2_cands = defaultdict(list)
    initial_unpaired = defaultdict(list)
    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped:
            continue
        if not (aln.reference_name.startswith("KIR") or aln.reference_name.startswith("IAG")):
            continue
        qn = aln.query_name
        if aln.next_reference_id == -1 or aln.mate_is_unmapped:
            initial_unpaired[qn].append((aln, None))
            continue 
        if aln.is_read1:
            if aln.next_reference_id == aln.reference_id:
                read1_cands[qn].append((aln, aln.next_reference_start))
            else:
                initial_unpaired[qn].append((aln, None))
        else:
            if aln.next_reference_id == aln.reference_id:
                read2_cands[qn].append(aln)
            else:
                initial_unpaired[qn].append((aln,None))
    bam.close()
    paired_dict = defaultdict(list)
    all_qns = set(read1_cands) | set(read2_cands)
    for qn in all_qns:
        r1_list = read1_cands.get(qn, [])
        r2_list = read2_cands.get(qn, [])
        for r1, mate_start in r1_list:
            for r2 in r2_list:
                if r2.reference_id == r1.reference_id and r2.reference_start == mate_start:
                    paired_dict[qn].append((r1, r2))
                    break
    for entry in initial_unpaired:
        paired_dict[entry].extend(initial_unpaired[entry])
    return paired_dict
bam_path=snakemake.input.sort
paired_dict= get_paired_by_mate_info_fast(bam_path)   
partition = set()
for qname, alignments in paired_dict.items():
    for r1, r2 in alignments:
        partition.add(r1.reference_name)
fasta_file = snakemake.input.fasta
allele_info = {}
utr_bed=snakemake.input.utr
exon_bed=snakemake.input.exon
intron_bed=snakemake.input.intron
utr_df = pd.read_csv(utr_bed, sep="\t", header=None, names=["allele", "start", "end", "feature"])
exon_df = pd.read_csv(exon_bed, sep="\t", header=None, names=["allele", "start", "end", "feature"])
intron_df = pd.read_csv(intron_bed, sep="\t", header=None, names=["allele", "start", "end", "feature"])
def build_ranges(df):
    ranges = {}
    for allele, group in df.groupby("allele"):
        ranges[allele] = list(zip(group["start"].astype(int), group["end"].astype(int)))
    return ranges
utr_ranges = build_ranges(utr_df)
exon_ranges = build_ranges(exon_df)
intron_ranges = build_ranges(intron_df)
original_seqs = {rec.id: str(rec.seq).replace("-", "") for rec in SeqIO.parse(snakemake.input.fasta, "fasta")}
allele_lengths = {record: len(original_seqs[record])for record in original_seqs}
for record in SeqIO.parse(fasta_file, "fasta"):
    if record.id in partition:
        allele_name = record.id
        allele_len = allele_lengths[allele_name]
        left_trim = int(allele_len * 0.15)
        right_trim = allele_len - int(allele_len * 0.15)
        utr_list = utr_ranges.get(allele_name, [])
        utr_start_ext = left_trim 
        utr_end_ext = right_trim 
        if utr_list:
            utr_start_end = max(e for s, e in utr_list if s == 0) if any(s == 0 for s, e in utr_list) else left_trim
            utr_end_start = min(s for s, e in utr_list if e == allele_len) if any(e == allele_len for s, e in utr_list) else right_trim
            utr_start_ext = max(left_trim, utr_start_end)
            utr_end_ext = min(right_trim, utr_end_start)
        allele_info[allele_name] = {
            "length": allele_len,
            "valid_start": utr_start_ext,
            "valid_end": utr_end_ext,
            "read_count": 0,
            "total_mismatches": 0,
            "covered_positions": set(),
            "total_mismatches/read_count/length": 1,
            "cluster": np.nan }
for qname, pairs in paired_dict.items():
    for r1, r2 in pairs:
        if r1 is None:
            continue
        allele = r1.reference_name
        if allele not in allele_info:
            continue
        info = allele_info[allele]
        start_cut, end_cut = info["valid_start"], info["valid_end"]
        info["read_count"] += 1
        if r1.reference_start is not None and r1.reference_end is not None:
            start = max(r1.reference_start, start_cut)
            end = min(r1.reference_end, end_cut)
            if start < end:
                info["covered_positions"].update(range(start, end))
        if r2 is not None and r2.reference_start is not None and r2.reference_end is not None:
            start = max(r2.reference_start, start_cut)
            end = min(r2.reference_end, end_cut)
            if start < end:
                info["covered_positions"].update(range(start, end))
for allele, info in allele_info.items():
    covered = info["covered_positions"]
    total_len = info["length"]
    start_cut, end_cut = info["valid_start"], info["valid_end"]
    def region_coverage(ranges):
        valid_positions = set()
        for s, e in ranges.get(allele, []):
            s = max(s, start_cut)
            e = min(e, end_cut)
            if s < e:
                valid_positions.update(range(s, e))
        if not valid_positions:
            return np.nan
        covered_here = valid_positions & covered
        return len(covered_here) / len(valid_positions) * 100
    exon_cov = region_coverage(exon_ranges)
    intron_cov = region_coverage(intron_ranges)
    info["exon_coverage_percent"] = round(exon_cov, 3) if not np.isnan(exon_cov) else np.nan
    info["intron_coverage_percent"] = round(intron_cov, 3) if not np.isnan(intron_cov) else np.nan
    info["exon_fully_covered"] = exon_cov >= 100 if not np.isnan(exon_cov) else False
    info["intron_95percent_covered"] = intron_cov >= 100 if not np.isnan(intron_cov) else False
    del info["covered_positions"]
allele_fully_covered = {
    a: "0" for a, info in allele_info.items()
    if info["exon_fully_covered"] and info["intron_95percent_covered"] }
with open(snakemake.output.dict, "wb") as f:
    pickle.dump(allele_fully_covered, f)     
