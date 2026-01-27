#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 15:49:49 2025

@author: yomna
"""
from itertools import combinations
import pysam
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt
import re
import numpy as np
import pandas as pd
from scipy.stats import nbinom, poisson
import pickle
from scipy.stats import binom
from Bio import SeqIO, AlignIO
import random
import os
from matplotlib.patches import Patch
import math
import scipy.sparse as sp
from scipy.sparse import coo_matrix

l=snakemake.input.list
with open(l) as f:
    allele_dict = {
        "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] : line.strip().split(",")[1]
        for line in f if not line.startswith("#") and line.strip()}
with open(l) as f:
    allele_dict_2 = {
     line.strip().split(",")[1] :  "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] 
        for line in f if not line.startswith("#") and line.strip()}  
Intron_regions=defaultdict(list)
with open(snakemake.input.intron) as bed:
    for line in bed:
        if line.startswith('#') or not line.strip():
            continue
        chrom, start, end = line.strip().split()[:3]
        Intron_regions[chrom].append((int(start), int(end)))     
UTR_regions=defaultdict(list)
with open(snakemake.input.utr) as bed: 
    for line in bed:
        if line.startswith('#') or not line.strip():
            continue
        chrom, start, end = line.strip().split()[:3]
        UTR_regions[chrom].append((int(start), int(end)))
original_seqs = {rec.id: str(rec.seq).replace("-", "") for rec in SeqIO.parse(snakemake.input.fa, "fasta")}
allele_lengths = {record: len(original_seqs[record])for record in original_seqs}
msa_files = snakemake.input.msa
seq_dict = {}
for msa_file in msa_files:
    alignment = AlignIO.read(msa_file, "fasta")
    for rec in alignment:
        seq_dict[rec.id.replace("KIR","")] = str(rec.seq)
_md_re = re.compile(r'(\d+)|([A-Z]|\^[A-Z]+)')   
fasta_path=snakemake.input.ref
def compute_error_Pc_chr17(bam_path, chrom="17"):
    md_re = re.compile(r'(\d+)|([A-Z]|\^[A-Z]+)')
    counts = Counter()
    c=0
    bam = pysam.AlignmentFile(bam_path, "rb")
    for aln in bam.fetch(chrom):
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        if aln.reference_name != chrom:
            continue
        c +=1
        for op, length in aln.cigartuples or []:
            if op == 0:   # M = alignment match (can be mismatch or match)
                counts['M'] += length
            elif op == 1:   # I
                counts['I'] += length
            elif op == 2: # D
                counts['D'] += length
            elif op == 4 or op == 5 : # S
                counts['S'] += length
        try:
            md = aln.get_tag("MD")
            
        except KeyError:
            continue
        mismatches = 0
        matches = 0
        for tok in md_re.finditer(md):
            num, var = tok.groups()
            if num:
                matches += int(num)
            if var and not var.startswith('^'):
                # a single-base substitution
                mismatches += 1
        counts['X'] += mismatches
        counts['M'] -= mismatches
    bam.close()
    cats = ['M', 'X', 'I', 'D', 'S']
    total = sum(counts[c] for c in cats)
    if total == 0:
        raise RuntimeError(f"No error operations found on chr{chrom}")
    P_c = {c: counts[c]/total for c in cats}
    return P_c     
def count_mismatches(aln):
    try:
        md = aln.get_tag("MD")
    except KeyError:
        mismatches = 0
        for op, length in aln.cigartuples or []:
            if op == 8:  # CIGAR 'X' = mismatch
                mismatches += length
        return mismatches
    mismatches = 0
    for num, var in _md_re.findall(md):
        if var and not var.startswith('^'):
            mismatches += 1
    return mismatches
def process_cigar(aln, P_c, allele):
    clip = 0
    Del = 0
    ins = 0
    m = 0
    mm = 0
    read_ll = 0
    if aln.cigartuples is not None:
        for i, (op, length) in enumerate(aln.cigartuples):
            if op in (4, 5):  # soft or hard clip
                if aln.reference_start == 0 and i == 0:
                    continue
                if aln.reference_end == allele_lengths[allele] and i == len(aln.cigartuples) - 1:
                    continue
                clip += length
            elif op == 1:  # Insertion
                ins += length
            elif op == 2:  # Deletion
                Del += length
            elif op == 7:  # Match
                m += length
    mm += count_mismatches(aln)
    read_ll += mm * np.log(P_c["X"]) + ins * np.log(P_c["I"]) + Del * np.log(P_c["D"])  + clip * np.log(P_c["S"]) + m * np.log(P_c["M"])
    return read_ll

  
def get_paired_by_mate_info_fast_kir_with_tag(kir_bam,selected_alleles):
    bam = pysam.AlignmentFile(kir_bam, "rb")
    paired_dict = {}
    read_to_tag_name=defaultdict(set)
    tag_to_read_name={}
    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped:
            continue
        if not aln.reference_name in selected_alleles :
            continue
        if aln.has_tag("TP"):
           tp_tag = aln.get_tag("TP")
           qn = aln.query_name
           read_to_tag_name[qn].add(tp_tag)
           tag_to_read_name[tp_tag]=qn
           allele= aln.reference_name
           if allele not in paired_dict:
              paired_dict[allele] = {}
           if tp_tag not in  paired_dict[allele]:
              paired_dict[allele][tp_tag] =0
              paired_dict[allele][tp_tag] +=  np.log(1 / allele_lengths[allele]) 
              paired_dict[allele][tp_tag] += process_cigar(aln,P_c,allele)
           else:   
              paired_dict[allele][tp_tag] += process_cigar(aln,P_c,allele)
    bam.close()
    return paired_dict,read_to_tag_name,tag_to_read_name
bed_file_exon=snakemake.input.exon
def parse_exons_bed(bed_file):
    exons = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            chrom, start, end = line.strip().split()[:3]
            exons[chrom].append((int(start), int(end)))
    lumped_exons = {}
    for chrom, intervals in exons.items():
        if not intervals:
            lumped_exons[chrom] = []
            continue
        sorted_intervals = sorted(intervals)
        merged = [list(sorted_intervals[0])]
        for s, e in sorted_intervals[1:]:
            last = merged[-1]
            if s - last[1] <= 1:  
                last[1] = max(last[1], e)
            else:
                merged.append([s, e])
        lumped_exons[chrom] = [tuple(m) for m in merged]            
    return lumped_exons
exon_ranges =parse_exons_bed(bed_file_exon)

def compute_tag_posterior_given_read(
    paired_dict_with_tag,        # {allele: {tag: log-likelihood}}
    read_to_tag_name,            # {read: [tags]}
    selected_alleles_with_cn,     # {allele: copy_number}
    min_cn=10**-30):
    tag_to_allele = {
        tag: allele
        for allele, tag_dict in paired_dict_with_tag.items()
        for tag in tag_dict
    }
    tag_posteriors = {}  # Final result: {tag: posterior probability given read}
    for read, tags in read_to_tag_name.items():
        log_numerators = {}
        for tag in tags:
            allele = tag_to_allele.get(tag)
            if allele is None:
                continue
            cn = selected_alleles_with_cn.get(allele, 0)
            if cn == 0:
                continue
            log_likelihood = paired_dict_with_tag[allele][tag]
            cn = min_cn  if cn == 0.01 else cn # penalize very small cn to force the reads to map on cn >= 1
            log_prob = np.log(cn) + log_likelihood
            log_numerators[tag] = log_prob
        if not log_numerators:
            continue
        max_logp = max(log_numerators.values())
        exp_shifted = {tag: np.exp(logp - max_logp) for tag, logp in log_numerators.items()}
        total = sum(exp_shifted.values())
        for tag, val in exp_shifted.items():
            tag_posteriors[tag] = val / total
    return tag_posteriors  # {tag: P(tag | read)}
def compute_weighted_stats_from_tag_posteriors(tag_posteriors, read_stats_by_allele,candidate,tag_to_read_name):
    stats_by_allele = {}
    for allele, tags in read_stats_by_allele.items():
        if candidate[allele] == 0:
           continue
        real_depth = defaultdict(float)
        depth = defaultdict(float)
        matches = defaultdict(float)
        mismatches = defaultdict(float)
        insertions = defaultdict(float)
        deletions = defaultdict(float)
        for tag, pos_dict in tags.items():
            prob = tag_posteriors.get(tag, 0)
            for pos, (m, mm, ins, dele) in pos_dict.items():
                total = m + mm + ins + dele
                if total == 0:
                    continue
                read_name=tag_to_read_name[tag] 
                if allele in read_name:
                    real_depth[pos] += total
                matches[pos] += prob * m
                mismatches[pos] += prob * mm
                insertions[pos] += prob * ins
                deletions[pos] += prob * dele
                depth[pos] += prob * total
        stats_by_allele[allele] = {
            'depth': dict(depth),
            'matches': dict(matches),
            'mismatches': dict(mismatches),
            'insertions': dict(insertions),
            'deletions': dict(deletions),
            "real_depth":  dict(real_depth)}
    return stats_by_allele


def allele_to_gene(allele):
    gene=allele.split('*', 1)[0]
    if gene == "2DL5A" or gene== "2DL5B":
       return "2DL5"
    else:
       return gene
def _normalize_segments(segments, L):
    out = []
    for s,e in segments:
        s = max(0, int(s)); e = min(L-1, int(e))
        if e >= s:
            out.append((s,e))
    if not out:
        return []
    out.sort()
    merged = [out[0]]
    for s,e in out[1:]:
        ps,pe = merged[-1]
        if s <= pe + 1:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s,e))
    return merged
def _complement_segments(exons, L):
    if not exons:
        return [(0, L-1)] if L > 0 else []
    introns, cur = [], 0
    for s,e in exons:
        if s > cur:
            introns.append((cur, s-1))
        cur = e + 1
    if cur <= L-1:
        introns.append((cur, L-1))
    return introns
def _total_len(segments):
    return sum(e - s + 1 for s,e in segments)
def _lump_segments_into_windows(segments, window_size, min_block_size):
    if not segments:
        return []
    windows, cur_block, cur_size = [], [], 0
    seg_idx, pos = 0, segments[0][0]
    while seg_idx < len(segments):
        s, e = segments[seg_idx]
        pos = max(pos, s)
        if pos > e:
            seg_idx += 1
            if seg_idx < len(segments):
                pos = segments[seg_idx][0]
            continue
        space_left = window_size - cur_size
        take = min(space_left, e - pos + 1)
        cur_block.append((pos, pos + take - 1))
        cur_size += take
        pos += take
        if cur_size == window_size:
            windows.append(cur_block); cur_block, cur_size = [], 0
        if pos > e:
            seg_idx += 1
            if seg_idx < len(segments):
                pos = segments[seg_idx][0]
    if cur_block:
        windows.append(cur_block)
    if len(windows) >= 2:
        last_size = sum(e - s + 1 for s,e in windows[-1])
        if last_size < min_block_size:
            windows[-2] = windows[-2] + windows[-1]
            windows.pop()
    return windows
def _slice_segments_into_K_windows(segments, K):
    """Split segments into exactly K consecutive windows with near-equal total sizes (quantile-style)."""
    Ltot = _total_len(segments)
    if Ltot == 0 or K <= 0:
        return []
    base = Ltot // K
    r = Ltot % K
    targets = [base + (1 if i < r else 0) for i in range(K)]  # sums to Ltot, min = floor(Ltot/K)
    windows, seg_i = [], 0
    pos = segments[0][0] if segments else 0
    for t in targets:
        remain, block = t, []
        while remain > 0 and seg_i < len(segments):
            s, e = segments[seg_i]
            pos = max(pos, s)
            if pos > e:
                seg_i += 1
                if seg_i < len(segments):
                    pos = segments[seg_i][0]
                continue
            take = min(remain, e - pos + 1)
            block.append((pos, pos + take - 1))
            pos += take
            remain -= take
            if pos > e:
                seg_i += 1
                if seg_i < len(segments):
                    pos = segments[seg_i][0]
        windows.append(block)
    return windows

def lump_exon_intron_windows_equalized_exons(
    selected_alleles,
    allele_lengths,                 # {allele: L}
    exon_ranges,                    # {allele: [(s,e_incl), ...]}
    window_size_hint=150,
    min_exon_window_size=100
):
    exons_norm, introns_norm = {}, {}
    for a in selected_alleles:
        L = int(allele_lengths[a])
        ex = _normalize_segments(exon_ranges.get(a, []), L)
        intr = _complement_segments(ex, L)
        exons_norm[a] = ex
        introns_norm[a] = intr

    gene_to_alleles = defaultdict(list)
    for a in selected_alleles:
        gene_to_alleles[allele_to_gene(allele_dict[a])].append(a)

    # Exon window count per gene (as before: based on shortest exon length)
    K_exon_by_gene = {}
    for g, alleles in gene_to_alleles.items():
        a_short = min(alleles, key=lambda a: _total_len(exons_norm[a]))
        wins_short = _lump_segments_into_windows(exons_norm[a_short], window_size_hint, min_exon_window_size)
        K_exon_by_gene[g] = len(wins_short)

    # NEW: Intron window count per gene = MIN number of intron windows across alleles (independent of exons)
    K_intron_by_gene = {}
    for g, alleles in gene_to_alleles.items():
        # build intron windows for each allele with the usual lumping (min size = 100)
        counts = []
        for a in alleles:
            wins = _lump_segments_into_windows(introns_norm[a], window_size_hint, min_exon_window_size)
            counts.append(len(wins))
        K_intron_by_gene[g] = min(counts) if counts else 0

    result = defaultdict(dict)
    window_type = defaultdict(dict)
    pos_type = defaultdict(lambda: np.array([], dtype="<U10"))
    lumped_sizes = set()

    for a in selected_alleles:
        g = allele_to_gene(allele_dict[a])
        idx = 0
        L = int(allele_lengths[a])

        # Per-base labels
        base_labels = np.full(L, "intergenic", dtype="<U10")
        for s, e in exons_norm[a]:
            base_labels[s:e+1] = "exon"
        for s, e in introns_norm[a]:
            base_labels[s:e+1] = "intron"
        pos_type[a] = base_labels

        # ----- EXON windows (unchanged: equalized per gene) -----
        Kx = K_exon_by_gene.get(g, 0)
        if Kx > 0 and _total_len(exons_norm[a]) > 0:
            exon_wins = _slice_segments_into_K_windows(exons_norm[a], Kx)
            for w in exon_wins:
                result[a][idx] = w
                window_type[a][idx] = "exon"
                lumped_sizes.add(sum(e - s + 1 for s, e in w))
                idx += 1

        # ----- INTRON windows (NEW: same count per gene = min across alleles) -----
        Ki = K_intron_by_gene.get(g, 0)
        if Ki > 0 and _total_len(introns_norm[a]) > 0:
            intron_wins = _slice_segments_into_K_windows(introns_norm[a], Ki)
        else:
            intron_wins = []  # no intron windows requested for this gene

        for w in intron_wins:
            result[a][idx] = w
            window_type[a][idx] = "intron"
            lumped_sizes.add(sum(e - s + 1 for s, e in w))
            idx += 1
    return result, window_type, lumped_sizes, K_exon_by_gene, pos_type


def tp_to_read(kir_bam,exon_windows):
    tp_to_reads = defaultdict(lambda: defaultdict(dict))
    with pysam.AlignmentFile(kir_bam, "rb") as bam:
        for allele, window_indices in exon_windows.items():
            allele_start = 0
            allele_end = allele_lengths[allele]
            for read in bam.fetch(allele, allele_start, allele_end + 1):
                if read.has_tag("TP") and not read.is_unmapped:
                    tp_tag = read.get_tag("TP")
                    if read.is_read1:
                        tp_to_reads[allele][tp_tag]['read1'] = read
                    elif read.is_read2:
                        tp_to_reads[allele][tp_tag]['read2'] = read  
    return tp_to_reads                
def bam_to_read_stats_by_allele_tag(
    kir_bam,
    selected_alleles,
    original_seqs,
    result, 
    window_type,
    min_base_quality=0):
    read_stats_by_allele_exon = {}    
    read_stats_by_allele_intron = {} 
    with pysam.AlignmentFile(kir_bam, "rb") as bamfile:
        for chrom in bamfile.references:
            if chrom not in selected_alleles:
                continue
            ref_seq = original_seqs[chrom].upper()
            exon_positions = set()
            intron_positions = set()
            for i in  result[chrom]:
                if window_type[chrom][i] == "exon":
                    for  start, end  in result[chrom][i]:
                          exon_positions.update(range(start, end+1))
                else:
                     for  start, end  in result[chrom][i]:
                           intron_positions.update(range(start, end+1))                    
            allele_data = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0]))  # match, mismatch, ins, del            
            for pileupcolumn in bamfile.pileup(
                chrom,
                0,
                len(ref_seq),
                truncate=True,
                stepper='all',
                min_base_quality=min_base_quality):
                pos = pileupcolumn.reference_pos
                if pos  in exon_positions or pos in intron_positions:
                    ref_base = ref_seq[pos]
                    for pileupread in pileupcolumn.pileups:
                        alignment = pileupread.alignment  # Store the alignment object in a local variable
                        #read_id = alignment.query_name
                        tp_value = alignment.get_tag("TP")
                        if pileupread.is_refskip:
                            continue
                        if pileupread.is_del:
                            allele_data[tp_value][pos][3] += 1  # deletion
                            continue
                        if pileupread.indel > 0:
                            allele_data[tp_value][pos][2] += 1  # insertion
                            continue
                        query_pos = pileupread.query_position
                        if query_pos is None:
                            continue                    
                        read_base = alignment.query_sequence[query_pos].upper()
                        if read_base == ref_base:
                            allele_data[tp_value][pos][0] += 1  # match
                        else:
                            allele_data[tp_value][pos][1] += 1  # mismatch                    
            read_stats_by_allele_exon[chrom] = {
                    tp_value: {pos: tuple(counts) for pos, counts in pos_dict.items()}
                    for tp_value, pos_dict in allele_data.items()}
    return read_stats_by_allele_exon

def build_sparse_counts(tp_to_reads, allele_to_idx, tag_to_idx, allele_lengths):
    rows, cols, data = [], [], []
    n_tags = len(tag_to_idx)
    for allele, tag_dict in tp_to_reads.items():
        if allele not in allele_to_idx:
            continue
        ai = allele_to_idx[allele]
        L = allele_lengths[allele]
        for tag, reads in tag_dict.items():
            if tag not in tag_to_idx:
                continue
            ti = tag_to_idx[tag]
            for rname in ("read1", "read2"):
                aln = reads.get(rname)
                if aln is None or aln.cigartuples is None:
                    continue
                refpos = aln.reference_start
                for op, ln in aln.cigartuples:
                    if op in (0, 7):  # match
                        feat = 0
                        for p in range(refpos, refpos + ln):
                            if p < L:
                                rows.append((ai*n_tags*4) + (ti*4) + feat)
                                cols.append(p)
                                data.append(1)
                        refpos += ln
                    elif op == 8:  # mismatch
                        feat = 1
                        for p in range(refpos, refpos + ln):
                            if p < L:
                                rows.append((ai*n_tags*4) + (ti*4) + feat)
                                cols.append(p)
                                data.append(1)
                        refpos += ln
                    elif op == 1:  # insertion I done we nned insertion to be included in the depth
                         feat = 2
                         if refpos < L:
                             rows.append((ai*n_tags*4) + (ti*4) + feat)
                             cols.append(refpos)
                             data.append(1) # whatever the size of the insertion they are considered one
                    elif op == 2:  # deletion
                        feat = 3
                        for p in range(refpos, refpos + ln):
                            if p < L:
                                rows.append((ai*n_tags*4) + (ti*4) + feat)
                                cols.append(p)
                                data.append(1)
                        refpos += ln
                    else:
                        pass
    counts_sparse = sp.csr_matrix((data, (rows, cols)))
    return counts_sparse
def expected_counts_with_features_variable(
    counts_sparse, tag_posteriors, tag_to_idx, allele_lengths, n_tags, n_features, allele_to_idx,selected_alleles
):
    idx_to_allele = {idx: allele for allele, idx in allele_to_idx.items()}
    tag_probs = np.zeros(n_tags, dtype=float)
    for tag, prob in tag_posteriors.items():
        tag_probs[tag_to_idx[tag]] = prob
    counts_coo = counts_sparse.tocoo()
    rows, cols, vals = counts_coo.row, counts_coo.col, counts_coo.data
    alleles  = rows // (n_tags * n_features)
    tags     = (rows // n_features) % n_tags
    features = rows % n_features
    weights  = tag_probs[tags]
    data = vals * weights
    max_pos = max(allele_lengths.values())
    row_idx = alleles * n_features * max_pos + features * max_pos + cols
    agg = coo_matrix((data, (row_idx, np.zeros_like(row_idx))),
                     shape=(len(allele_to_idx)*n_features*max_pos, 1)).toarray().ravel()

    expected = {}
    for a, allele in idx_to_allele.items():
        L = allele_lengths[allele]
        mat = np.zeros((n_features, L), dtype=float)
        for f in range(n_features):
            start = a*n_features*max_pos + f*max_pos
            mat[f, :] = agg[start:start+L]
        expected[allele] = mat
    return expected
def compute_avg_depth_per_window(expected, result, window_type):
    avg_depth = {}
    avg_depth_exon = {}
    for allele, feats in expected.items():
        depth = feats[[0, 1, 3], :].sum(axis=0) # total depth per position exclusing insertion
        avg_depth[allele] = {}
        avg_depth_exon[allele] = {}
        for win_idx, ranges in result[allele].items():
            total = 0.0
            length = 0
            for s, e in ranges:
                total += depth[s:e+1].sum()
                length += (e - s + 1)
            avg = total / length if length > 0 else 0.0
            avg_depth[allele][win_idx] = avg
            if window_type[allele][win_idx] == "exon":
                avg_depth_exon[allele][win_idx] = avg
    return avg_depth, avg_depth_exon
with open(snakemake.input.allele_rep, "rb") as f:
    allele_representatives = pickle.load(f)

def ll_column_wise_mapped_divergent(candidate, stats_by_allele, f, seq_dict, pos_type):
    msa_len = len(next(iter(f.values())))
    # depth stays numeric, matches will be per-allele values
    position_sums = {
        'depth': [0.0] * msa_len,
        'matches': [[] for _ in range(msa_len)],
        'positions': [[] for _ in range(msa_len)]  # new key for original positions
    }
    for allele in candidate:
        depth_d = stats_by_allele[allele]["depth"]
        m_d = stats_by_allele[allele]["matches"]
        for msa_pos, orig_pos in enumerate(f[allele]):
            if orig_pos is None:
                continue
            position_sums['depth'][msa_pos] += depth_d.get(orig_pos, 0)
            position_sums['matches'][msa_pos].append(m_d.get(orig_pos, 0))
            position_sums['positions'][msa_pos].append(orig_pos)
    divergent_positions = []
    for msa_pos in range(msa_len):
        bases = set()
        for allele in candidate:
            seq = seq_dict[allele_dict[allele].upper()]
            if msa_pos < len(seq):
                base = seq[msa_pos]
                bases.add(base)
        if len(bases) > 1 and "-" not in bases:
            divergent_positions.append(msa_pos)
    exons = {k: [] for k in position_sums}
    others = {k: [] for k in position_sums}

    for msa_pos in divergent_positions:
        orig_pos = None
        allele_for_type = None
        for allele in candidate:
            if f[allele][msa_pos] is not None:
                orig_pos = f[allele][msa_pos]
                allele_for_type = allele
                break
        if orig_pos is None or allele_for_type is None:
            region_type = "other"
        else:
            region_type = pos_type[allele_for_type][orig_pos]

        if region_type == "exon":
            for key in position_sums:
                exons[key].append(position_sums[key][msa_pos])
        else:
            for key in position_sums:
                others[key].append(position_sums[key][msa_pos])
    exons['matches'] = [tuple(x) for x in exons['matches']]
    others['matches'] = [tuple(x) for x in others['matches']]
    exons['positions'] = [tuple(x) for x in exons['positions']]
    others['positions'] = [tuple(x) for x in others['positions']]
    return exons, others
def map_alignment_to_original_positions(seq_dict, alleles):
    allele_to_positions = {allele: [] for allele in alleles}
    for allele in alleles:
        aligned_seq = seq_dict[allele_dict[allele].upper()]
        orig_pos = -1
        for base in aligned_seq:
            if base != '-':
                orig_pos += 1
                allele_to_positions[allele].append(orig_pos)
            else:
                allele_to_positions[allele].append(None)
    return allele_to_positions
def divergent_exons_to_vcf(alleles, divergent_sum_e, original_seqs, f, min_frac=0.75):
    vcf_lines = []
    for depth,matches, positions in zip(divergent_sum_e["depth"],divergent_sum_e["matches"], divergent_sum_e["positions"]):
        total =depth
        frac = [m / total for m in matches]
        if max(frac) < min_frac:
           continue
        donor_idx = frac.index(max(frac))     # allele with strong support
        donor_allele = alleles[donor_idx]
        donor_pos = positions[donor_idx]
        for idx, fval in enumerate(frac):
            if idx == donor_idx:
                continue
            if fval <= (1 - min_frac):  # e.g., donor ≥75%, target ≤25%
                target_allele =alleles[idx]
                target_pos = positions[idx]
                ref_base = original_seqs[target_allele][target_pos]
                alt_base = original_seqs[donor_allele][donor_pos]
                chrom = target_allele
                pos = target_pos + 1  # VCF is 1-based
                vid = "."
                qual = "."
                flt = "PASS"
                info = f"Donor={allele_dict[donor_allele]};Donor_matches_freq={max(frac)}"
                fmt = "GT"
                sample = "1/1"
                vcf_line = f"{ allele_dict[chrom]}\t{pos}\t{vid}\t{ref_base}\t{alt_base}\t{qual}\t{flt}\t{info}\t{fmt}\t{sample}"
                vcf_lines.append(vcf_line)
    return vcf_lines
def check_if_there_is_more_than_x_perc_errors(
    g_nonzero, expected, pos_type,
    kir_bam, ref_fasta,
    which="exon", thresh=0.75,
): 
    vcf_lines = []
    bamfile = pysam.AlignmentFile(kir_bam, "rb", reference_filename=ref_fasta)

    for allele in g_nonzero:
        mat = expected[allele]      # shape (4, L)
        L = mat.shape[1]

        region_mask = (pos_type[allele][:L] == which)
        if not np.any(region_mask):
            continue

        depth     = mat.sum(axis=0)
        mm_vals   = mat[1, :]
        ins_vals  = mat[2, :]
        del_vals  = mat[3, :]

        err_tot   = mm_vals + ins_vals + del_vals
        with np.errstate(divide="ignore", invalid="ignore"):
            err_frac = np.where(depth > 0, err_tot / depth, 0.0)

        mask = (err_frac > thresh) & region_mask & (depth > 0)
        if not np.any(mask):
            continue

        positions = np.where(mask)[0]
        depths    = depth[mask].astype(int)

        for pos, d in zip(positions, depths):
            if d == 0:
                continue

            chrom = allele

            # reference base from the original sequence
            ref_base = original_seqs[allele][pos]

            # collect alternative alleles
            alt_candidates = Counter()

            for pileupcolumn in bamfile.pileup(chrom, pos, pos+1, truncate=True):
                if pileupcolumn.pos != pos:
                    continue

                for read in pileupcolumn.pileups:

                    # INSERTION ---------------------
                    if read.indel > 0:
                        ins = read.alignment.query_sequence[
                            read.query_position + 1 :
                            read.query_position + 1 + read.indel
                        ]
                        ALT = ref_base + ins   # VCF ALT
                        REF = ref_base         # VCF REF
                        alt_candidates[(REF, ALT,"INS")] += 1

                    # DELETION ----------------------
                    elif read.is_del:
                        # Count how many reads support deletion at this position
                        del_len = 1
                        deleted = original_seqs[allele][pos-1 : pos + del_len]
                        REF = deleted
                        ALT = deleted[0]    # VCF rule: ALT = first base of REF
                        alt_candidates[(REF, ALT,"DEL")] += 1

                    # MISMATCH ----------------------
                    elif not read.is_del and not read.is_refskip:
                        base = read.alignment.query_sequence[read.query_position]
                        if base.upper() != ref_base.upper():
                            REF = ref_base
                            ALT = base
                            alt_candidates[(REF, ALT,"SNP")] += 1

            if not alt_candidates:
                continue

            # pick the most common (REF, ALT) pair
            (REF, ALT,T), _ = alt_candidates.most_common(1)[0]

            # construct VCF line
            if T == "DEL":
                vcf_line = (
                    f"{allele_dict[chrom]}\t{pos}\t.\t{REF}\t{ALT}\t.\t"
                    f"FAIL_ERRFRAC\tErrorFrac={err_frac[pos]:.3f};Depth={d}\tGT\t./."
                )
            else:
                vcf_line = (
                    f"{allele_dict[chrom]}\t{pos+1}\t.\t{REF}\t{ALT}\t.\t"
                    f"FAIL_ERRFRAC\tErrorFrac={err_frac[pos]:.3f};Depth={d}\tGT\t./."
                )                
            vcf_lines.append(vcf_line)

    bamfile.close()
    return vcf_lines


def write_vcf(vcf_path, all_vcf_lines, sample_name="SAMPLE"):
    """
    Write VCF lines into a valid VCF file.
    """
    contigs = sorted({line.split("\t")[0] for line in all_vcf_lines})
    with open(vcf_path, "w") as out:
        # Mandatory metadata
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=KIR_variant_checker\n")
        for c in contigs:
            out.write(f"##contig=<ID={c}>\n")
        out.write("##INFO=<ID=ErrorFrac,Number=1,Type=Float,Description=\"Observed error fraction\">\n")
        out.write("##INFO=<ID=Depth,Number=1,Type=Integer,Description=\"Read depth\">\n")
        out.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        out.write("##FILTER=<ID=FAIL_ERRFRAC,Description=\"Error fraction too high\">\n")
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name + "\n")
        for line in all_vcf_lines:
            out.write(line + "\n")


i=snakemake.params.sample
with open(snakemake.input.a, "rb") as f:
        gene_copy_number_dict=pickle.load(f)
with open(snakemake.input.b, "rb") as f:
            gene_allele_level=pickle.load(f) 
selected_alleles=[a for a in gene_allele_level if gene_allele_level[a] > 0.01 ]  
allele_to_gene_map={}   
ref_fasta=snakemake.input.fa
for a in selected_alleles:
    gene=allele_dict[a].split("*") [0]
    allele_to_gene_map[a] = "2DL5" if gene == "2DL5A" or gene == "2DL5B" else gene             
result, window_type, lumped_sizes, K_exon_by_gene,pos_type = lump_exon_intron_windows_equalized_exons(
    selected_alleles,
    allele_lengths,                
    exon_ranges, 
    window_size_hint=150,
    min_exon_window_size=100)            
kir_bam = snakemake.input.paired
_md_re = re.compile(r'(\d+)|([A-Z]|\^[A-Z]+)')
output=    snakemake.output.vcf
bam_path= snakemake.input.bam
P_c = compute_error_Pc_chr17(bam_path, chrom=snakemake.params.chr17)  
tp_to_reads=tp_to_read(kir_bam,result)
all_tags=[t for a in tp_to_reads for t in tp_to_reads[a] ]
allele_to_idx = {allele: i for i, allele in enumerate(selected_alleles)}
tag_to_idx    = {tag: i for i, tag in enumerate(all_tags)}
counts_sparse=build_sparse_counts(tp_to_reads, allele_to_idx, tag_to_idx, allele_lengths)
paired_dict_with_tag,read_to_tag_name,tag_to_read_name=get_paired_by_mate_info_fast_kir_with_tag(kir_bam,selected_alleles)
selected_alleles_with_cn = { allele: gene_allele_level.get(allele, 0.01) for allele in selected_alleles}                   
tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, selected_alleles_with_cn,min_cn=10**-300) 
expected=expected_counts_with_features_variable( counts_sparse, tag_posteriors, tag_to_idx, allele_lengths, len(all_tags), 4, allele_to_idx,selected_alleles)
avg_depth,avg_depth_exon=compute_avg_depth_per_window(expected, result, window_type) 
read_stats_by_allele_tag_total=bam_to_read_stats_by_allele_tag(kir_bam,selected_alleles, original_seqs,result, window_type,min_base_quality=0)
stats_by_allele_c=compute_weighted_stats_from_tag_posteriors(tag_posteriors, read_stats_by_allele_tag_total,selected_alleles_with_cn,tag_to_read_name)
all_vcf_lines = []
for gene, CN_g in gene_copy_number_dict.items():
    g_alleles_t = [a for a in selected_alleles if allele_to_gene_map[a] == gene]
    gg=defaultdict(list)
    for allele in g_alleles_t:
        gg[allele_dict[allele].split("*")[0]].append(allele)
    for gene in gg:
        g_alleles = gg[gene]
        f = map_alignment_to_original_positions(seq_dict, g_alleles)
        divergent_sum_e, divergent_sum_o = ll_column_wise_mapped_divergent(
            g_alleles, stats_by_allele_c, f, seq_dict, pos_type)
        err_lines = check_if_there_is_more_than_x_perc_errors(
            g_alleles, expected, pos_type, kir_bam, ref_fasta,
            which="exon", thresh=0.75)
        all_vcf_lines.extend(err_lines)
        div_lines = divergent_exons_to_vcf(g_alleles, divergent_sum_e, original_seqs, f, min_frac=0.85) #0.8
        all_vcf_lines.extend(div_lines)    
write_vcf(output, all_vcf_lines, sample_name=i)
