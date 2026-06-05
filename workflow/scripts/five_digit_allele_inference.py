#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 15:49:49 2025

@author: yomna
"""
import scipy.sparse as sp
from itertools import combinations
import pysam
import numpy as np
from collections import defaultdict
from collections import Counter
import matplotlib.pyplot as plt
import re
import pandas as pd
from scipy.stats import nbinom, poisson
import pickle
from scipy.stats import binom
from Bio import SeqIO, AlignIO
import random
import os
from matplotlib.patches import Patch
import math
from scipy.sparse import coo_matrix
from scipy.stats import multinomial
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

def ll_column_wise_and_windows_backward(candidate,avg_depth, P_c, result, window_assignments,window_type,expected,pos_type,which="exon", w_exon = 1.0,w_intron = 0.0,threshould=0.2):
    total = 0.0    
    for allele, cn in candidate.items():
        prob_d = 0.0
        prob_mm = 0.0
        prob_ins = 0.0
        prob_del = 0.0
        prob_total_err= 0.0
        if cn > 0.01:
            mat = expected[allele]      # shape (4, L)
            L = mat.shape[1]
            exon_mask = (pos_type[allele][:L] == which )   # only exon bases
            if not np.any(exon_mask):
                continue
            depth = mat.sum(axis=0)        # per position depth
            mm_vals  = mat[1, :]
            ins_vals = mat[2, :]
            del_vals = mat[3, :]
            err_tot = mm_vals + ins_vals + del_vals
            with np.errstate(divide="ignore", invalid="ignore"):
                err_frac = np.where(depth > 0, err_tot / depth, 0.0)
            mask = (err_frac > threshould) & exon_mask
            n = np.rint(depth[mask]).astype(int)
            k_mm  = np.rint(mm_vals[mask]).astype(int)
            k_ins = np.rint(ins_vals[mask]).astype(int)
            k_del = np.rint(del_vals[mask]).astype(int)
            k_mm  = np.clip(k_mm,  0, n)
            k_ins = np.clip(k_ins, 0, n)
            k_del = np.clip(k_del, 0, n)
            total_err= k_mm + k_ins + k_del
            prob_mm  += np.sum(binom.logpmf(k_mm,  n, P_c["X"]))
            pos_idx = np.where(mask)[0] 
            prob_ins += np.sum(binom.logpmf(k_ins, n, P_c["I"]))
            prob_del += np.sum(binom.logpmf(k_del, n, P_c["D"]))
            prob_total_err += np.sum(binom.logpmf(total_err, n, P_c["D"]+ P_c["I"]+ P_c["X"]))
        params_for_w=  window_assignments[cn]
        for win_idx in result[allele]:
            y_hat = float(avg_depth[allele][win_idx])
            win_len = sum((end - start + 1) for start, end in result[allele][win_idx])
            depth_param = params_for_w[win_len]
            kind = depth_param[-1]    
            if window_type[allele][win_idx] == 'exon':
                w_win = w_exon
            else:
                w_win = w_intron
            if kind == "Poisson":
                prob_d += w_win*poisson.logpmf(round(y_hat), depth_param[0])
        total +=   prob_d  +prob_total_err#(prob_mm+ prob_ins +  prob_del ) #logL_error
    return total


def check_if_there_is_more_than_20_perc_errors(
    g_nonzero, avg_depth, P_c, result, window_assignments,
    window_type, expected, pos_type, selected_alleles_with_cn, which="exon",threshold=0.2):
    more_than_20_errors = False
    significant_exon_depth = False
    for allele in g_nonzero:
        mat = expected[allele] 
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
        mask = (err_frac >= threshold) & region_mask
        if np.any(mask):
            more_than_20_errors = True
            positions = np.where(mask)[0]
            depths    = depth[mask].astype(int)
            with open(snakemake.output.log, "a") as f:
               print(f"[{allele_dict[allele]}] >{threshold} errors at positions: {positions}, depths: {depths}",file=f)
        cn = selected_alleles_with_cn[allele]
        for win_idx, ranges in result[allele].items():
            if window_type[allele][win_idx] != "exon":
                continue
            obs = sum(depth[s:e+1].sum() for s, e in ranges)
            win_len = sum((e - s + 1) for s, e in ranges)
            obs= obs/win_len
            params_for_w = window_assignments[cn]
            depth_param = params_for_w[win_len]
            lam = depth_param[0]  # Poisson λ
            if obs <= lam:
                pval = poisson.cdf(obs, lam)
            else:
                pval = poisson.sf(obs - 1, lam)
            pval_two_sided = 2 * min(pval, 1 - pval)
            if pval_two_sided < 0.01:
                significant_exon_depth = True
                with open(snakemake.output.log, "a") as f:
                    print(f"[{allele_dict[allele]}] Window {win_idx} depth {obs} significantly different from expectation "
                    f"(λ={lam:.2f}, p={pval_two_sided:.3g})",file=f)
    return more_than_20_errors, significant_exon_depth
def check_the_same_exon(allele1,allele2):
    # code1=allele_dict[allele1].split("*")[-1]
    # code2=allele_dict[allele2].split("*")[-1]
    # same_exon=False
    # if len(code1) == len(code2) and len(code1) == 7 : This method is a bug
    #     if code1[:-2] == code2[:-2]:
    #         same_exon=True
    code1=allele_representatives["KIR"+allele_dict[allele1]]
    code2=allele_representatives["KIR"+allele_dict[allele2]]
    same_exon=False
    if code1 == code2:
        same_exon=True            
    return same_exon
def forward_backward_selection(paired_dict_with_tag,read_to_tag_name,allele_lengths,exon_windows,P_c ,tp_to_reads,window_assignments,output ,window_type,counts_sparse,gene_copy_number_dict,gene_allele_level,threshold):
    selected_alleles_with_cn = {
        allele: gene_allele_level.get(allele, 0.01) 
        for allele in selected_alleles}
    for g in gene_copy_number_dict:
        if gene_copy_number_dict[g]==1:
            a=sorted(Top_4_Predicted_Alleles_all_mod[g], key=lambda x: x[1], reverse=True)[0][0]
            g_alleles = [a for a in selected_alleles if allele_to_gene_map[a] == g]
            selected_alleles_with_cn[allele_dict_2[a]]=1
            for other in g_alleles:
                if other != allele_dict_2[a]:
                    selected_alleles_with_cn[other]=0.01
    #3    
    tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, selected_alleles_with_cn,min_cn=10**-300) 
    expected=expected_counts_with_features_variable( counts_sparse, tag_posteriors, tag_to_idx, allele_lengths, len(all_tags), 4, allele_to_idx,selected_alleles)
    avg_depth,avg_depth_exon=compute_avg_depth_per_window(expected, result, window_type) 
    current_ll = ll_column_wise_and_windows_backward(selected_alleles_with_cn,avg_depth, P_c, result, window_assignments,window_type,expected,pos_type,which="exon", w_exon = 1.0,w_intron = 0.0,threshould=threshold)
    backward_improvement = True
    genes_to_process=[]
    for gene, CN_g in gene_copy_number_dict.items():
        g_alleles = [a for a in selected_alleles if allele_to_gene_map[a] == gene]
        if not g_alleles:
            continue
        g_nonzero = [a for a in g_alleles if int(selected_alleles_with_cn[a]) > 0]  # candidates to decrease
        uu,ff=check_if_there_is_more_than_20_perc_errors( g_nonzero, avg_depth, P_c, result, window_assignments, window_type, expected, pos_type,selected_alleles_with_cn, which="exon",threshold=threshold)
        if uu or ff:
            genes_to_process.append(gene) 
    with open(snakemake.output.log, "a") as f:    
         print("gene to process in one-allele-at-a-time swaping movement",genes_to_process,file=f)        
    while backward_improvement:
        backward_improvement = False
        back_best_ll = current_ll
        back_best_move = None  # (gene, dec_tuple, inc_tuple)
        for gene in genes_to_process:
            checked_comp=[]
            g_alleles = [a for a in selected_alleles if allele_to_gene_map[a] == gene]
            if not g_alleles:
                continue
            g_nonzero = [a for a in g_alleles if int(selected_alleles_with_cn[a]) > 0]  # candidates to decrease
            g_zero    = [a for a in g_alleles if int(selected_alleles_with_cn[a]) >= 0] # candidates to increase
            for dec_tuple in combinations(g_nonzero, 1):
                for inc_tuple in combinations(g_zero, 1):
                    if inc_tuple[0] == dec_tuple[0]:
                        continue
                    kkk= (allele_representatives["KIR"+allele_dict[inc_tuple[0]]], allele_representatives["KIR"+allele_dict[dec_tuple[0]]])
                    if kkk in checked_comp:
                        continue
                    else:
                        checked_comp.append(kkk)
                    candidate = selected_alleles_with_cn.copy()
                    if check_the_same_exon(dec_tuple[0],inc_tuple[0]):
                        continue
                    for a in dec_tuple:
                        candidate[a] = int(candidate[a]) - 1
                        if candidate[a] ==0:
                            candidate[a]= 0.01
                    for a in inc_tuple:
                        candidate[a] = int(candidate[a]) + 1
                    tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, candidate, min_cn=10**-300)
                    expected=expected_counts_with_features_variable( counts_sparse, tag_posteriors, tag_to_idx, allele_lengths, len(all_tags), 4, allele_to_idx,selected_alleles)
                    avg_depth,avg_depth_exon=compute_avg_depth_per_window(expected, result, window_type)
                    total_ll = ll_column_wise_and_windows_backward(candidate,avg_depth, P_c, result, window_assignments,window_type,expected,pos_type,which="exon", w_exon = 1.0,w_intron = 0.0,threshould=threshold)
                    if total_ll > back_best_ll:
                        back_best_ll = total_ll
                        back_best_move = (gene, dec_tuple, inc_tuple)
                        with open(snakemake.output.log, "a") as f:
                            print(f"Backward candidate (gene {gene}): "
                            f"dec {[allele_dict[d] for d in dec_tuple]} "
                            f"inc {[allele_dict[i] for i in inc_tuple]} | LL: {back_best_ll:.3f}",file=f)      
        if back_best_move:
            gene, dec_tuple, inc_tuple = back_best_move
            for a in dec_tuple:
                selected_alleles_with_cn[a] -= 1
                if selected_alleles_with_cn[a] ==0:
                    selected_alleles_with_cn[a] = 0.01
            for a in inc_tuple:
                selected_alleles_with_cn[a] += 1
                selected_alleles_with_cn[a] = int(selected_alleles_with_cn[a])
            current_ll = back_best_ll
            with open(snakemake.output.log, "a") as f:
                print(f"Performed Backward k-swap (gene {gene}): "
                f"dec {[allele_dict[d] for d in dec_tuple]}, "
                f"inc {[allele_dict[i] for i in inc_tuple]} | LL: {current_ll:.3f}",file=f)
            backward_improvement = True
        else:
            with open(snakemake.output.log, "a") as f:
                print("No further backward improvement. Entering two-allele-at-a-time swaping movement.",file=f)
            break  
             #4   
    backward_improvement = True
    tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, selected_alleles_with_cn, min_cn=10**-300)
    expected=expected_counts_with_features_variable( counts_sparse, tag_posteriors, tag_to_idx, allele_lengths, len(all_tags), 4, allele_to_idx,selected_alleles)
    avg_depth,avg_depth_exon=compute_avg_depth_per_window(expected, result, window_type)
    current_ll = ll_column_wise_and_windows_backward(selected_alleles_with_cn,avg_depth, P_c, result, window_assignments,window_type,expected,pos_type,which="exon", w_exon = 1.0,w_intron = 0.0,threshould=threshold)            
    genes_to_process=[]
    for gene, CN_g in gene_copy_number_dict.items():
        if CN_g <=1:
            continue
        g_alleles = [a for a in selected_alleles if allele_to_gene_map[a] == gene]
        if not g_alleles:
            continue
        g_nonzero = [a for a in g_alleles if int(selected_alleles_with_cn[a]) > 0]  # candidates to decrease
        uu,ff=check_if_there_is_more_than_20_perc_errors( g_nonzero, avg_depth, P_c, result, window_assignments, window_type, expected, pos_type,selected_alleles_with_cn, which="exon",threshold=threshold)
        if uu or ff:
            genes_to_process.append(gene) 
    with open(snakemake.output.log, "a") as f:    
        print("gene to process in two-alleles-at-a-time swaping movement",genes_to_process,file=f)                    
    while backward_improvement:
        backward_improvement = False
        back_best_ll = current_ll
        back_best_move = None  # (gene, dec_tuple, inc_tuple)
        for gene in genes_to_process:
            g_alleles = [a for a in selected_alleles if allele_to_gene_map[a] == gene]
            if not g_alleles:
                continue
            g_nonzero = [a for a in g_alleles if int(selected_alleles_with_cn[a]) > 0]  # candidates to decrease
            g_zero    = [a for a in g_alleles if int(selected_alleles_with_cn[a]) == 0] # candidates to increase
            checked_comp=[]
            trial= min(2, len(g_nonzero),len(g_zero) )
            for dec_tuple in combinations(g_nonzero, trial):
                ttt= sorted([allele_representatives["KIR"+allele_dict[a]] for a in dec_tuple])
                for inc_tuple in combinations(g_zero, trial):
                    zzz= sorted([allele_representatives["KIR"+allele_dict[a]] for a in inc_tuple])
                    aaa= any ([True for a in ttt if a in zzz])
                    if (ttt,zzz) in checked_comp or aaa:
                        continue
                    else:
                        checked_comp.append((ttt,zzz))
                    if sorted(dec_tuple) == sorted(inc_tuple) or ttt == zzz:
                        continue 
                    candidate = selected_alleles_with_cn.copy()
                    for a in dec_tuple:
                        candidate[a] = int(candidate[a]) - 1
                        if candidate[a] ==0:
                            candidate[a]= 0.01    
                    for a in inc_tuple:
                        candidate[a] = int(candidate[a]) + 1
                    tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, candidate, min_cn=10**-300)
                    expected=expected_counts_with_features_variable( counts_sparse, tag_posteriors, tag_to_idx, allele_lengths, len(all_tags), 4, allele_to_idx,selected_alleles)
                    avg_depth,avg_depth_exon=compute_avg_depth_per_window(expected, result, window_type)
                    total_ll = ll_column_wise_and_windows_backward(candidate,avg_depth, P_c, result, window_assignments,window_type,expected,pos_type,which="exon", w_exon = 1.0,w_intron = 0.0,threshould=threshold)
                    if total_ll > back_best_ll:
                        back_best_ll = total_ll
                        back_best_move = (gene, dec_tuple, inc_tuple)
                        with open(snakemake.output.log, "a") as f:
                            print(f"Backward candidate (gene {gene}): "
                            f"dec {[allele_dict[d] for d in dec_tuple]} "
                            f"inc {[allele_dict[i] for i in inc_tuple]} | LL: {back_best_ll:.3f}",file=f)
  
        if back_best_move:
            gene, dec_tuple, inc_tuple = back_best_move
            for a in dec_tuple:
                selected_alleles_with_cn[a] -= 1
                if  selected_alleles_with_cn[a] == 0:
                    selected_alleles_with_cn[a] =0.01
            for a in inc_tuple:
                selected_alleles_with_cn[a] += 1
                selected_alleles_with_cn[a] = int(selected_alleles_with_cn[a] )
            current_ll = back_best_ll
            with open(snakemake.output.log, "a") as f:
                print(f"Performed Backward k-swap (gene {gene}): "
                f"dec {[allele_dict[d] for d in dec_tuple]}, "
                f"inc {[allele_dict[i] for i in inc_tuple]} | LL: {current_ll:.3f}",file=f)
            backward_improvement = True
            backward_improvement = True
            tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, selected_alleles_with_cn, min_cn=10**-300)
            expected=expected_counts_with_features_variable( counts_sparse, tag_posteriors, tag_to_idx, allele_lengths, len(all_tags), 4, allele_to_idx,selected_alleles)
            avg_depth,avg_depth_exon=compute_avg_depth_per_window(expected, result, window_type)
            current_ll = ll_column_wise_and_windows_backward(selected_alleles_with_cn,avg_depth, P_c, result, window_assignments,window_type,expected,pos_type,which="exon", w_exon = 1.0,w_intron = 0.0,threshould=threshold)  
            for gene in genes_to_process:
                g_alleles = [a for a in selected_alleles if allele_to_gene_map[a] == gene]
                if not g_alleles:
                    continue
                g_nonzero = [a for a in g_alleles if int(selected_alleles_with_cn[a]) > 0]  # candidates to decrease
                if not check_if_there_is_more_than_20_perc_errors( g_nonzero, avg_depth, P_c, result, window_assignments, window_type, expected, pos_type,selected_alleles_with_cn, which="exon",threshold=threshold):
                    genes_to_process.remove(gene)   
            if  not genes_to_process:
                break
                
        else:
            with open(snakemake.output.log, "a") as f:
                print("No further improvement. Selection complete.",file=f)
            break                 
    gene_allele_level=      selected_alleles_with_cn
    return {allele: cn for allele, cn in selected_alleles_with_cn.items() if cn > 0},gene_allele_level

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
def write_selected_alignments_to_bam_all(tag_posteriors, read_to_tag_name, tp_to_reads, kir_bam,output_bam_2):
    flat_tp_to_reads = {}
    for allele_dict in tp_to_reads.values():
        for tag, read_pair in allele_dict.items():
            flat_tp_to_reads[tag] = read_pair
    with pysam.AlignmentFile(kir_bam, "rb") as bam_in:
        header = bam_in.header
    with pysam.AlignmentFile(output_bam_2, "wb", header=header) as bam_out:
        for read_name, tag_list in read_to_tag_name.items():
            valid_tags = [tag for tag in tag_list if( tag in tag_posteriors)]
            if not valid_tags:
                continue
            probs = np.array([tag_posteriors[tag] for tag in valid_tags])
            best_tag = np.random.choice(valid_tags, p=probs)
            if best_tag not in flat_tp_to_reads:
                continue
            read_pair = flat_tp_to_reads[best_tag]
            for rtype in ['read1', 'read2']:
                if rtype in read_pair:
                    read = read_pair[rtype]
                    bam_out.write(read)                    
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

gene_copy_number_dict={}
gene_allele_level={}
with open(snakemake.output.log, "w") as f:
     print("",file=f)
i=snakemake.params.sample
with open(snakemake.input.cn, "rb") as f:
        gene_copy_number_dict=pickle.load(f)
with open(snakemake.input.allele, "rb") as f:
            gene_allele_level=pickle.load(f)        
bam_path=snakemake.input.bam2
P_c = compute_error_Pc_chr17(bam_path, chrom=snakemake.params.chr17)
with open(snakemake.input.allele_representatives_all_comp, "rb") as f:
    comp = pickle.load(f)
Top_4_Predicted_Alleles_all_mod=defaultdict(list)
for gene, alleles in comp.items():
    if gene== "2DL5A" or gene=="2DL5B":
        for allele in alleles:
            Top_4_Predicted_Alleles_all_mod["2DL5"].append(allele)
    else:
        for allele in alleles:
            Top_4_Predicted_Alleles_all_mod[gene].append(allele)    
selected_alleles=[allele_dict_2[g[0]] for a in Top_4_Predicted_Alleles_all_mod if a in gene_copy_number_dict for g in Top_4_Predicted_Alleles_all_mod[a]  ]  
allele_to_gene_map={}        
for a in selected_alleles:
    gene=allele_dict[a].split("*") [0]
    allele_to_gene_map[a] = "2DL5" if gene == "2DL5A" or gene == "2DL5B" else gene             
result, window_type, lumped_sizes, K_exon_by_gene,pos_type = lump_exon_intron_windows_equalized_exons(
    selected_alleles,
    allele_lengths,                 # {allele: L}
    exon_ranges,          # {allele: [(s,e_incl), ...]}
    window_size_hint=150,
    min_exon_window_size=100)
cn=[0.01,1,2,3,4]
avg_depth_by_window_= background_avg_depth_on_region(
            snakemake.input.bam2,
            contig=snakemake.params.chr17, #"NC_000017.11"
            region_start=snakemake.params.start,
            region_end=snakemake.params.end,
            min_base_quality=0)
window_assignments=defaultdict(dict)
for c in cn:
    avg_depth_by_window=avg_depth_by_window_* c/2
    for w in lumped_sizes:
        window_assignments[c][w]= (avg_depth_by_window,"","Poisson") #get_nb_or_binomial_params(mu, var,c)
    with open(snakemake.output.log, "a") as f:
        print(c, avg_depth_by_window,file=f)
kir_bam = snakemake.input.paired
_md_re = re.compile(r'(\d+)|([A-Z]|\^[A-Z]+)')
output=    snakemake.params.out 
tp_to_reads=tp_to_read(kir_bam,result)
all_tags=[t for a in tp_to_reads for t in tp_to_reads[a] ]
allele_to_idx = {allele: i for i, allele in enumerate(selected_alleles)}
tag_to_idx    = {tag: i for i, tag in enumerate(all_tags)}
counts_sparse=build_sparse_counts(tp_to_reads, allele_to_idx, tag_to_idx, allele_lengths)
paired_dict_with_tag,read_to_tag_name,tag_to_read_name=get_paired_by_mate_info_fast_kir_with_tag(kir_bam,selected_alleles)
threshold= 0.2
infer,gene_allele_level_new= forward_backward_selection(paired_dict_with_tag,read_to_tag_name,allele_lengths,result,P_c ,tp_to_reads,window_assignments,output,window_type ,counts_sparse,gene_copy_number_dict,gene_allele_level,threshold)
infer={a: infer[a] for a in infer if infer[a] > 0.01}
infer={a: infer[a] for a in infer if infer[a] > 0.01}
tag_posteriors=compute_tag_posterior_given_read(  paired_dict_with_tag,   read_to_tag_name, infer, min_cn=10**-300)
avg_depth_by_window, depth_by_pos=get_window_depth_counts_avg_with_deletions_fast(
        result,             # {allele: {win_idx: [(s,e_inclusive), ...]}}
        tag_posteriors,     # {tp_tag: weight}
        tp_to_reads,        # {allele: {tp_tag: {'read1': Aln, 'read2': Aln}}}
        allele_lengths,     # {allele: length}
        infer )
read_stats_by_allele_tag_total=bam_to_read_stats_by_allele_tag(kir_bam,list(infer.keys()), original_seqs,result, window_type,min_base_quality=0)
stats_by_allele_c=compute_weighted_stats_from_tag_posteriors(tag_posteriors, read_stats_by_allele_tag_total,infer,tag_to_read_name)
output_bam_2=snakemake.output.output_bam_2#f"{path}/kir/real_data/{i}/cn_optimization/29/assignment_all_with_intron.bam"
write_selected_alignments_to_bam_all(tag_posteriors, read_to_tag_name, tp_to_reads, kir_bam,output_bam_2)
plot_depth_full_length(
        stats_by_allele_c,
        avg_depth_by_window,
        list(infer.keys()),
        output=snakemake.output.pdf,
        allele_lengths= allele_lengths,
        windows_by_allele=result,          # {allele: {win_idx: [(s,e), ...], ...}} in genomic coords
        copy_number_dict=infer,
        window_assignments=window_assignments,          # {cn: {size: (mean_val, p, ttype)}}
        allele_dict=allele_dict,                 # pretty-name map
        exon_ranges=exon_ranges)                 # {allele: [(start,end), ...]} in genomic coord)
with open(snakemake.output.allele_local, "wb") as f: #f"{path}/kir/real_data/{i}/cn_optimization/29/29_allele_local.pkl"
        pickle.dump(gene_allele_level_new, f)  
