
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 15:49:49 2025

@author: yomna
"""
from itertools import combinations
import pysam
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from collections import Counter
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
import scipy.sparse as sp
from scipy.sparse import coo_matrix
import math

l=snakemake.input.list
with open(l) as f:
    allele_dict = {
        "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] : line.strip().split(",")[1]
        for line in f if not line.startswith("#") and line.strip()}
with open(l) as f:
    allele_dict_2 = {
     line.strip().split(",")[1] :  "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] 
        for line in f if not line.startswith("#") and line.strip()}  
original_seqs = {rec.id: str(rec.seq).replace("-", "") for rec in SeqIO.parse(snakemake.input.fa, "fasta")}
allele_lengths = {record: len(original_seqs[record])for record in original_seqs}
_md_re = re.compile(r'(\d+)|([A-Z]|\^[A-Z]+)')   
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
def plot_depth_full_length(
    stats_by_allele_c_t,
    stats_by_allele_w_t,
    allele_list,
    output,
    allele_lengths,
    windows_by_allele=None,          # {allele: {win_idx: [(s,e), ...], ...}} in genomic coords
    copy_number_dict=None,
    window_assignments=None,          # {cn: {size: (mean_val, p, ttype)}}
    allele_dict=None,                 # pretty-name map
    exon_ranges=None,                 # {allele: [(start,end), ...]} in genomic coords
    show_exon_labels=True
):
    num_alleles = len(allele_list)
    fig, axes = plt.subplots(num_alleles, 1, figsize=(12, 4 * num_alleles), sharex=False)
    if num_alleles == 1:
        axes = [axes]

    for ax, allele in zip(axes, allele_list):
        L = allele_lengths.get(allele)
        if allele not in stats_by_allele_c_t or L is None:
            ax.set_title(f"{allele} (no data)")
            ax.set_ylabel("Count")
            ax.grid(True)
            continue

        stats = stats_by_allele_c_t[allele]
        # If your per-position dicts are 1-based, switch to range(1, L+1)
        flat_positions, flat_depth, flat_matches, flat_mismatches, flat_ins, flat_del, flat_true = [], [], [], [], [], [],[]
        for pos in range(L):
            if pos in stats['depth']:
                flat_positions.append(pos)
                flat_depth.append(stats['depth'][pos])
                flat_matches.append(stats['matches'].get(pos, 0))
                flat_mismatches.append(stats['mismatches'].get(pos, 0))
                flat_ins.append(stats['insertions'].get(pos, 0))
                flat_del.append(stats['deletions'].get(pos, 0))
                flat_true.append(stats['real_depth'].get(pos, 0))
        if flat_positions:
            ax.plot(flat_positions, flat_depth, label='Depth', color='blue')
            ax.plot(flat_positions, flat_matches, label='Matches', color='green')
            ax.plot(flat_positions, flat_mismatches, label='Mismatches', color='red')
            ax.plot(flat_positions, flat_ins, label='Insertions', color='purple')
            ax.plot(flat_positions, flat_del, label='Deletions', color='orange')
            ax.plot(flat_positions, flat_true, label='True_depth', color='brown')
        if exon_ranges and allele in exon_ranges:
            ymax = ax.get_ylim()[1] if flat_positions else 1.0
            for i, (s, e) in enumerate(sorted(exon_ranges[allele], key=lambda x: x[0]), start=1):
                # Clamp to allele length just in case
                s = max(0, min(s, L - 1))
                e = max(0, min(e, L - 1))
                if e < s:
                    continue
                ax.axvspan(s, e, alpha=0.3, color='lightgreen')
                if show_exon_labels:
                    ax.text((s + e) / 2, 0.95 * ymax, f"Exon {i}",
                            ha='center', va='top', fontsize=8, color='darkgreen')
        if windows_by_allele and allele in windows_by_allele and allele in stats_by_allele_w_t:
            window_counts = stats_by_allele_w_t[allele]      # {win_idx: count}
            first_line = True
            ymax = ax.get_ylim()[1] if flat_positions else 1.0

            for window_index, count in window_counts.items():
                segs = windows_by_allele[allele].get(window_index)
                if not segs:
                    continue

                size = 0
                last_s = last_e = None
                for s, e in segs:
                    s = max(0, min(s, L - 1))
                    e = max(0, min(e, L - 1))
                    if e < s:
                        continue
                    size += (e - s + 1)
                    ax.hlines(count, s, e, color="black",linestyles='-', linewidth=2,
                              label='Window Depth' if first_line else "")
                    last_s, last_e = s, e
                    first_line = False
                if copy_number_dict and allele in copy_number_dict and window_assignments:
                    cn = copy_number_dict[allele]
                    if cn in window_assignments and size in window_assignments[cn]:
                        mean_val, p, ttype = window_assignments[cn][size]
                        if ttype == "NegativeBinomial":
                            mean_val = mean_val * (1 - p) / p
                        ax.hlines(mean_val, last_s, last_e, linestyles='dashed', linewidth=1.5,
                                  label='Expected Mean' if window_index == 0 else "")

        cn_str = f" CN={copy_number_dict[allele]}" if (copy_number_dict and allele in copy_number_dict) else ""
        pretty = allele_dict.get(allele, allele) if allele_dict else allele
        ax.set_title(f"Allele: {pretty}{cn_str}")
        ax.set_ylabel("Count")
        ax.grid(True)
        ax.legend(loc='upper right', fontsize=8)

    axes[-1].set_xlabel("Allele position (bp)")
    plt.tight_layout()
    plt.savefig(output)


def ll_column_wise_and_windows(candidate,stats_by_allele_w, P_c, result, window_assignments,window_type):#stats_by_allele_c
    prob_d = 0.0
    prob_mm = 0.0
    prob_ins = 0.0
    prob_del = 0.0
    total = 0.0
    w_exon = 1.0
    w_intron = 1.0
    for allele, cn in candidate.items():
        if cn == 0:
            continue 
        params_for_w= window_assignments[cn] #window_assignments[1]  if cn >= 1 else window_assignments[cn] 
        for win_idx in result[allele]:
            y_hat = float(stats_by_allele_w[allele][win_idx])
            win_len = sum((end - start + 1) for start, end in result[allele][win_idx])
            depth_param = params_for_w[win_len]
            kind = depth_param[-1]    
            if window_type[allele][win_idx] == 'exon':
                w_win = w_exon
            elif window_type[allele][win_idx] == 'intron':
                w_win = w_intron
            if kind == "Poisson":
                prob_d +=  w_win* poisson.logpmf(round(y_hat), depth_param[0])# if cn >=1 else   w_win* poisson.logpmf(round(y_hat), depth_param[0]) 
    total +=   prob_d + prob_mm+ prob_ins +  prob_del     
    return total
def forward_backward_selection(paired_dict_with_tag,read_to_tag_name,allele_lengths,exon_windows,P_c ,tp_to_reads,window_assignments,output ,window_type,counts_sparse):
    selected_alleles_with_cn = {allele: 0.01 for allele in selected_alleles}
    current_ll = -np.inf
    iteration = 1
    c=0
    for cn_cap in [4]:
        with open(snakemake.output.log, "a") as f:
            print(f"\n=== Starting forward/backward pass with CN cap = {cn_cap} ===",file=f)
        forward_improvement = True
        while forward_improvement:
            forward_improvement = False
            best_allele = None
            best_ll = current_ll
            # Forward step: try incrementing each allele
            for allele in selected_alleles:
                candidate = selected_alleles_with_cn.copy()
                candidate[allele] += 1
                candidate[allele]  = int(candidate[allele])
                if candidate[allele] >=cn_cap:
                    continue
                tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, candidate, min_cn=10**-4)
                expected=expected_counts_with_features_variable( counts_sparse, tag_posteriors, tag_to_idx, allele_lengths, len(all_tags), 4, allele_to_idx,selected_alleles)
                avg_depth,avg_depth_exon=compute_avg_depth_per_window(expected, result, window_type)
                total_ll = ll_column_wise_and_windows(candidate,avg_depth, P_c, result, window_assignments,window_type) 
                if total_ll > best_ll:
                    best_ll = total_ll
                    best_allele = allele
                    with open(snakemake.output.log, "a") as f:
                        print("Forward candidate:", allele_dict[best_allele], best_ll, file=f)
            if best_ll > current_ll:
                selected_alleles_with_cn[best_allele] += 1
                selected_alleles_with_cn[best_allele]= int(selected_alleles_with_cn[best_allele])
                current_ll = best_ll
                forward_improvement = True
                with open(snakemake.output.log, "a") as f:
                    print(f"[Iteration {iteration}] Forward: incremented {allele_dict[best_allele]} to CN={selected_alleles_with_cn[best_allele]} | Log-likelihood: {current_ll:.2f}",file=f)
                iteration += 1
                continue     
            gene_copy_number_dict=defaultdict(int)
            for a in selected_alleles_with_cn:
                if selected_alleles_with_cn[a] > 0.01:
                    g= allele_dict[a].split("*")[0]
                    if g == "2DL5A" or g == "2DL5B":
                       gene_copy_number_dict["2DL5"] += selected_alleles_with_cn[a] 
                    else:
                        gene_copy_number_dict[g] += selected_alleles_with_cn[a]   
    gene_allele_level=      selected_alleles_with_cn
    return {allele: cn for allele, cn in selected_alleles_with_cn.items() if cn > 0},gene_allele_level,gene_copy_number_dict

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
    if not exons or len(exons) < 2:
        return []  # no introns between exons
    introns = []
    exons = sorted(exons)
    for (s1, e1), (s2, e2) in zip(exons[:-1], exons[1:]):
        if s2 > e1 + 1:
            introns.append((e1 + 1, s2 - 1))
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
    K_exon_by_gene = {}
    for g, alleles in gene_to_alleles.items():
        # build intron windows for each allele with the usual lumping (min size = 100)
        counts = []
        for a in alleles:
            wins = _lump_segments_into_windows(exons_norm[a], window_size_hint, min_exon_window_size)
            counts.append(len(wins))
        K_exon_by_gene[g] = min(counts) if counts else 0
    # NEW: Intron window count per gene = MIN number of intron windows across alleles (independent of exons)
    K_intron_by_gene = {}
    for g, alleles in gene_to_alleles.items():
        # build intron windows for each allele with the usual lumping (min size = 100)
        counts = []
        for a in alleles:
            wins = _lump_segments_into_windows(introns_norm[a], window_size_hint, min_exon_window_size)
            counts.append(len(wins))
        K_intron_by_gene[g] = min(counts) if counts else 0
        positive_counts = [c for c in counts if c > 0]
        if positive_counts:
           K_intron_by_gene[g] = min(positive_counts)
        else:
           K_intron_by_gene[g] = 0
    result = defaultdict(dict)
    window_type = defaultdict(dict)
    pos_type = defaultdict(lambda: np.array([], dtype="<U10"))
    lumped_sizes = set()
    for a in selected_alleles:
        g = allele_to_gene(allele_dict[a])
        idx = 0
        L = int(allele_lengths[a])

        # Per-base labels
        base_labels = np.full(L, "utr", dtype="<U10")
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



def make_windows(region_start, region_end, window_size):
    """Return non-overlapping windows [(start, end), ...] covering [start, end)."""
    windows = []
    n_bp = region_end - region_start
    num_windows = math.ceil(n_bp / window_size)
    for i in range(num_windows):
        w_start = region_start + i * window_size
        w_end   = min(region_end - 1, w_start + window_size - 1)
        windows.append((w_start, w_end))
    return windows

def background_avg_depth_on_region(
    bam_path,
    contig="NC_000017.11",
    region_start=74_000_000,
    region_end=76_000_000,
    min_base_quality=0
):
    total_depth = 0
    region_length = region_end - region_start +1
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for col in bam.pileup(
            contig,
            region_start,
            region_end,
            truncate=True,
            stepper="all",
            min_base_quality=min_base_quality):
            pos = col.reference_pos
            if pos < region_start or pos >= region_end:
                continue
            depth = 0
            for pr in col.pileups:
                depth += 1
            total_depth += depth
    avg_depth = total_depth / region_length 
    return avg_depth


def tp_to_read(kir_bam,exon_windows):
    tp_to_reads = defaultdict(lambda: defaultdict(dict))
    with pysam.AlignmentFile(kir_bam, "rb") as bam:
        for allele, _ in exon_windows.items():
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
_REF_CONSUME = {0, 2, 7, 8}  # M, D, =, X
def _add_spans_to_diff(diff_arr, read, w, include_refskip=False):
    """Accumulate posterior weight w into a difference array via CIGAR spans.
       Counts deletions (D). Optionally counts refskip (N) if include_refskip=True."""
    ct = read.cigartuples
    if not ct:
        return
    L = diff_arr.shape[0]
    r = read.reference_start
    for op, ln in ct:
        if op in _REF_CONSUME or (include_refskip and op == 3):  # consume reference
            s = r
            e = r + ln  # half-open
            # clamp
            if e > 0 and s < L:
                s0 = 0 if s < 0 else s
                e0 = L if e > L else e
                if e0 > s0:
                    diff_arr[s0] += w
                    if e0 < L:
                        diff_arr[e0] -= w
            r += ln
        elif op == 3:  # N, consume ref but not counted unless include_refskip
            r += ln
        else:
            pass

def get_window_depth_counts_avg_with_deletions_fast(
    result,             # {allele: {win_idx: [(s,e_inclusive), ...]}}
    tag_posteriors,     # {tp_tag: weight}
    tp_to_reads,        # {allele: {tp_tag: {'read1': Aln, 'read2': Aln}}}
    allele_lengths,     # {allele: length}
    candidate,          # {allele: CN}
    include_refskip=False
):
    diff = {}
    for allele, L in allele_lengths.items():
        if candidate.get(allele, 0) > 0:
            diff[allele] = np.zeros(L, dtype=np.float32)
    for allele, tp_tags in tp_to_reads.items():
        if allele not in diff:
            continue
        darr = diff[allele]
        for tp_tag, pair in tp_tags.items():
            w = tag_posteriors.get(tp_tag)
            if not w:
                continue
            w = float(w)
            r1 = pair.get('read1'); r2 = pair.get('read2')
            if r1 is not None:
                _add_spans_to_diff(darr, r1, w, include_refskip=include_refskip)
            if r2 is not None:
                _add_spans_to_diff(darr, r2, w, include_refskip=include_refskip)
    depth_by_pos = {}
    for allele, darr in diff.items():
        depth_by_pos[allele] = np.cumsum(darr, dtype=np.float32)
    avg_depth_by_window = defaultdict(dict)
    for allele, windows in result.items():
        arr = depth_by_pos.get(allele)
        if arr is None:
            for i, ranges in windows.items():
                total_len = sum(max(0, min(allele_lengths[allele]-1, e) - max(0, s) + 1) for s, e in ranges)
                avg_depth_by_window[allele][i] = 0.0 if total_len > 0 else 0.0
            continue
        csum = np.empty(arr.size + 1, dtype=np.float64)
        csum[0] = 0.0
        np.cumsum(arr, out=csum[1:], dtype=np.float64)
        for i, ranges in windows.items():
            total_len = 0
            total_sum = 0.0
            L = arr.size
            for s, e_incl in ranges:
                s0 = 0 if s < 0 else s
                e0 = L-1 if e_incl >= L else e_incl
                if e0 < s0:
                    continue
                # inclusive -> half-open
                e_ex = e0 + 1
                total_len += (e_ex - s0)
                total_sum += (csum[e_ex] - csum[s0])
            avg_depth_by_window[allele][i] = (total_sum / total_len) if total_len > 0 else 0.0
    return avg_depth_by_window, depth_by_pos
def simulate_cn3_nb(window_assignments, window_sizes, num_samples=10000):
    cn3_assignments = {}
    n1 = window_assignments[1][list(window_sizes)[0]][0] 
    n2= window_assignments[2][list(window_sizes)[0]][0]
    samples1 = poisson.rvs(n1, size=num_samples)
    samples2 = poisson.rvs(n2, size=num_samples)            
    combined_samples = samples1 + samples2
    mu = np.mean(combined_samples)
    for w in window_sizes:
        cn3_assignments[w] = (mu,"","Poisson") #(n3, p3,t)
    return cn3_assignments
def simulate_cn4_nb(window_assignments, window_sizes, num_samples=10000):
    cn3_assignments = {}
    n1 = window_assignments[2][list(window_sizes)[0]][0] 
    n2= window_assignments[2][list(window_sizes)[0]][0] 
    samples1 = poisson.rvs(n1, size=num_samples)
    samples2 = poisson.rvs(n2, size=num_samples)            
    combined_samples = samples1 + samples2
    mu = np.mean(combined_samples)
    for w in window_sizes:
        cn3_assignments[w] = (mu,"","Poisson")
    return cn3_assignments                 
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
        depth = feats.sum(axis=0)  # total depth per position
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
def filter_by_threshold(allele_scores, frac=0.1):
    """
    Filter alleles per gene with 10% threshold relative to top.
    Returns dict: {gene: {allele: True/False}}.
    """
    filtered = defaultdict(list)
    for gene, scores in allele_scores.items():
        # linear scale
        scores_lin = [(a, np.exp(s)) for a, s in scores]
        if not scores_lin:
            continue
        scores_lin.sort(key=lambda x: x[1], reverse=True)

        top_val = scores_lin[0][1]
        thresh = ( frac) * top_val

        for a, v in scores_lin:
            if  (v >= thresh):
                filtered[gene.replace("KIR","")].append(a.replace("KIR","") )
    return filtered
  
for i in [ snakemake.params.sample ]:
    with open(snakemake.output.log, "w") as f:
          print(i, file=f)
    gene_copy_number_dict={}
    gene_allele_level={}
    bam_path=snakemake.input.bam2
    P_c = compute_error_Pc_chr17(bam_path, chrom=snakemake.params.chr17)
    sele = snakemake.input.allele_representatives_all_comp
    with open(sele, "rb") as f:
        Top_4_Predicted_Alleles_all = pickle.load(f)
    Top_4_Predicted_Alleles_all_mod=defaultdict(list)
    for gene, alleles in Top_4_Predicted_Alleles_all.items():
        if gene== "2DL5A" or gene=="2DL5B":
            for allele in alleles:
                Top_4_Predicted_Alleles_all_mod["2DL5"].append(allele)
        else:
            for allele in alleles:
                Top_4_Predicted_Alleles_all_mod[gene].append(allele)    
    filtered_all   = filter_by_threshold(Top_4_Predicted_Alleles_all_mod, frac=0.1)
    selected_alleles = []            
    for gene, allele in filtered_all.items(): 
        for a in allele:
            selected_alleles.append(allele_dict_2[a])      
    allele_to_gene_map={}        
    for a in selected_alleles:
        gene=allele_dict[a].split("*") [0]
        allele_to_gene_map[a] = "2DL5" if gene == "2DL5A" or gene == "2DL5B" else gene               
    result, window_type, lumped_sizes, K_exon_by_gene,pos_type = lump_exon_intron_windows_equalized_exons(
        selected_alleles,
        allele_lengths,
        exon_ranges,  
        window_size_hint=150,
        min_exon_window_size=100) 
    cn=[0.01,1,2,3,4]
    avg_depth_by_window_= background_avg_depth_on_region(
            snakemake.input.bam2,
            contig=snakemake.params.chr17, #"NC_000017.11"
            region_start=74_000_000,
            region_end=76_000_000,
            min_base_quality=0)
    window_assignments=defaultdict(dict)
    for c in cn:
        avg_depth_by_window=avg_depth_by_window_* c/2
        for w in lumped_sizes:
            window_assignments[c][w]= (avg_depth_by_window,"","Poisson")#get_nb_or_binomial_params(mu, var,c)
        with open(snakemake.output.log, "a") as f:
             print(c,avg_depth_by_window,file=f)
    kir_bam = snakemake.input.paired
    _md_re = re.compile(r'(\d+)|([A-Z]|\^[A-Z]+)')
    output=    snakemake.params.out# f"{path}/kir/real_data/{i}/cn_optimization/29"     
    tp_to_reads=tp_to_read(kir_bam,result)
    all_tags=[t for a in tp_to_reads for t in tp_to_reads[a] ]
    allele_to_idx = {allele: i for i, allele in enumerate(selected_alleles)}
    tag_to_idx    = {tag: i for i, tag in enumerate(all_tags)}
    counts_sparse=build_sparse_counts(tp_to_reads, allele_to_idx, tag_to_idx, allele_lengths)
    paired_dict_with_tag,read_to_tag_name,tag_to_read_name=get_paired_by_mate_info_fast_kir_with_tag(kir_bam,selected_alleles)
    infer,gene_allele_level,gene_copy_number_dict= forward_backward_selection(paired_dict_with_tag,read_to_tag_name,allele_lengths,result,P_c ,tp_to_reads,window_assignments,output,window_type ,counts_sparse)
    tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, infer, min_cn=10**-4)
    avg_depth_by_window, depth_by_pos=get_window_depth_counts_avg_with_deletions_fast(
            result,             # {allele: {win_idx: [(s,e_inclusive), ...]}}
            tag_posteriors,     # {tp_tag: weight}
            tp_to_reads,        # {allele: {tp_tag: {'read1': Aln, 'read2': Aln}}}
            allele_lengths,     # {allele: length}
            infer )
    read_stats_by_allele_tag_total=bam_to_read_stats_by_allele_tag(kir_bam,selected_alleles, original_seqs,result, window_type, min_base_quality=0)
    stats_by_allele_c=compute_weighted_stats_from_tag_posteriors(tag_posteriors, read_stats_by_allele_tag_total,infer,tag_to_read_name)
    plot_depth_full_length(
            stats_by_allele_c,
            avg_depth_by_window,
            selected_alleles,
            output=output+f"/coverage_plot_after_cn_inference.pdf",
            allele_lengths= allele_lengths,
            windows_by_allele=result, 
            copy_number_dict=infer,
            window_assignments=window_assignments, 
            allele_dict=allele_dict,
            exon_ranges=exon_ranges)  
    with open(snakemake.output.cn, "wb") as f:
            pickle.dump(gene_copy_number_dict, f)   
    with open(snakemake.output.allele, "wb") as f:
            pickle.dump(gene_allele_level, f)  
