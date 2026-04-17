#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 13:53:57 2025

@author: yomna
"""
import pysam
from collections import Counter
import re
from collections import defaultdict
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import pickle
with open(snakemake.input.pkl, "rb") as f:
    allele_representatives = pickle.load(f)
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
bam_path1=snakemake.input.bam_path1
P_c=compute_error_Pc_chr17(bam_path1, chrom=snakemake.params.chr17)
bam_path2=snakemake.input.bam_path2
l=snakemake.input.list
paired_dict= get_paired_by_mate_info_fast(bam_path2)
bed_file_exon= snakemake.input.exon
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
with open(l) as f:
    allele_dict = {
        "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] : line.strip().split(",")[1]
        for line in f if not line.startswith("#") and line.strip()}
with open(l) as f:
    allele_dict_2 = {
     line.strip().split(",")[1] :  "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] 
        for line in f if not line.startswith("#") and line.strip()}    
original_seqs = {rec.id: str(rec.seq).replace("-", "") for rec in SeqIO.parse(snakemake.input.fasta, "fasta")}
allele_lengths = {allele_dict[record]: len(original_seqs[record]) for record in original_seqs}        
sele=snakemake.input.number_of_other
top_n_per_cluster_new=defaultdict(list)
with open(sele) as f: 
    for line in f:
         s = line.strip()
         gene_name = s
         if gene_name in ["2DL5A", "2DL5B"]:
            gene_name = "2DL5"
         top_n_per_cluster_new[gene_name].append(s)
def group_and_filter_pairs_by_cluster(d, allele_dict, top_n_per_cluster_new):
    cluster_pairs = defaultdict(lambda: defaultdict(list))
    filtered_count = 0
    allele_to_gene = {allele: ("2DL5" if gene in ["2DL5A", "2DL5B"] else gene)for gene, alleles in top_n_per_cluster_new.items() for allele in alleles }
    for read_name, pairs in d.items():
        for aln1, aln2 in pairs:
            allele1 = allele_dict.get(aln1.reference_name) if aln1 else  allele_dict.get(aln2.reference_name)
            gene1 = allele_to_gene.get(allele1)
            if gene1:
                cluster_pairs[gene1][read_name].append((aln1, aln2))
                filtered_count += 1
    return cluster_pairs, filtered_count
cluster_pairs,filtered_count = group_and_filter_pairs_by_cluster(paired_dict, allele_dict, top_n_per_cluster_new)    
_md_re = re.compile(r'(\d+)|([A-Z]|\^[A-Z]+)')
def count_mismatches(aln):
    try:
        md = aln.get_tag("MD")
    except KeyError:
        mismatches = 0
        matches = 0
        for op, length in aln.cigartuples or []:
            if op == 8:  # CIGAR X (mismatch)
                mismatches += length
            if op ==7:
                matches +=length 
        return (mismatches, matches)
    mismatches = 0
    matches = 0
    for num, var in _md_re.findall(md):
        if var and not var.startswith('^'):
            mismatches += 1
        if num:
           matches += int(num)    
    return (mismatches, matches)
def compute_read_cluster_likelihoods_no_norm(cluster_pairs, P_c, top_n_per_cluster_new, allele_dict, allele_lengths):
    read_likelihoods = defaultdict(dict)
    for gene, reads_dict in cluster_pairs.items():
        norm_gene = "2DL5" if gene in ["2DL5A", "2DL5B"] else gene
        for read_name, aln_list in reads_dict.items():
            best_p = -np.inf
            for aln1, aln2 in aln_list:
                aln = aln1 if aln1 is not None else aln2
                ref_name = aln.reference_name
                allele_length = allele_lengths.get(allele_dict[ref_name])
                counts = {'X': 0, 'I': 0, 'D': 0,"M":0, 'S': 0} #,
                for a in (aln1, aln2):
                    if a is None:
                        continue
                    cigars = a.cigartuples or []
                    ref_start = a.reference_start
                    ref_end = a.reference_end
                    ref_len = allele_length
                    for i, (op, length) in enumerate(cigars):
                        if op == 1:
                            counts['I'] += length
                        elif op == 2:
                            counts['D'] += length
                        elif op in (4, 5):  # soft/hard clip
                            is_first = (i == 0)
                            is_last = (i == len(cigars) - 1)
                            if is_first and ref_start == 0:
                                continue
                            if is_last and ref_end == ref_len:
                                continue
                            counts['S'] += length    
                xm1, m1 = count_mismatches(aln1) if aln1 else (0,0)
                xm2, m2 = count_mismatches(aln2) if aln2 else (0,0)
                counts['X'] = xm1 + xm2
                counts['M'] = m1 + m2
                p = 0.0
                #p = 1
                for c, nc in counts.items():
                    if nc:
                        p += nc*np.log(P_c[c])
                normalized_p = p - np.log( allele_length)
                best_p = max(best_p, normalized_p)
            read_likelihoods[read_name][norm_gene] = best_p
    return read_likelihoods     
read_likelihoods=    compute_read_cluster_likelihoods_no_norm(cluster_pairs, P_c, top_n_per_cluster_new,allele_dict, allele_lengths)
read_to_alleles_to_cigar={}
for read in read_likelihoods:
     read_to_alleles_to_cigar[read]={}
     for allele in read_likelihoods[read]:
         read_to_alleles_to_cigar[read][allele]=[]
         for al in cluster_pairs[allele][read]:
              read1_cigar = tuple(al[0].cigar) if al[0] else None
              read2_cigar = tuple(al[1].cigar) if al[1] else None    
              read_to_alleles_to_cigar[read][allele].append((read1_cigar,read2_cigar))
identical_cigar_reads = []
for read, allele_to_cigars in read_to_alleles_to_cigar.items():
    cigar_sets = set()
    for cigars in allele_to_cigars.values():
        cigar_tuple = tuple(end for end in cigars)
        cigar_sets.add(cigar_tuple)
    if len(cigar_sets) == 1:
        identical_cigar_reads.append(read)     
def  read_likelihoods_all_mod(read_likelihoods_all):
     read_likelihoods_all_mod={}
     for read, genes in read_likelihoods_all.items():
         read_likelihoods_all_mod[read]={}
         for gene, alleles in genes.items():
             for allele in alleles:
                 read_likelihoods_all_mod[read][allele_dict[allele]]=read_likelihoods_all[read][gene][allele][0]
     return read_likelihoods_all_mod  
def run_em(read_likelihoods, interesting_reads_post=None,tol=0.01, max_iters=1000):
    zeft=[]
    allele_set = set()
    for read, alleles in read_likelihoods.items():
        for allele in alleles:
            allele_set.add(allele)
    F = {allele: np.log(1 / len(allele_set)) for allele in allele_set}
    zeft.append(F.copy())
    reads = list(read_likelihoods.keys())
    R = len(reads)
    LL = [-np.inf]
    for iteration in range(max_iters):
        read_given_allele_times_abundance = defaultdict(dict)
        allele_given_read = defaultdict(dict)
        for read, alleles_likelihoods in read_likelihoods.items():
            scores = [F[a] + alleles_likelihoods[a] for a in alleles_likelihoods]
            max_score = max(scores)
            exp_scores = {a: np.exp(F[a] + alleles_likelihoods[a] - max_score) for a in alleles_likelihoods}
            denom = sum(exp_scores.values())
            for a in exp_scores:
                allele_given_read[read][a] = exp_scores[a] / denom    
        for a in F:
            total = np.log(sum(allele_given_read[read][a] if a in allele_given_read[read] else 0 for read in allele_given_read ) ) - np.log( R)
            F[a] = total
        zeft.append(F.copy())
        temp = {
            read: sum(allele_given_read[read][a] * np.exp(F[a]) for a in allele_given_read[read])
            for read in allele_given_read
        }
        L = sum(np.log(temp[read]) for read in temp)
        if L - LL[-1] > tol:
            LL.append(L)
        else:
            break
    return F, allele_given_read, LL,zeft
alleles = sorted({a for alleles in read_likelihoods.values() for a in alleles})
reads   = list(read_likelihoods.keys())
allele_to_idx = {a:i for i,a in enumerate(alleles)}
read_to_idx   = {r:i for i,r in enumerate(reads)}
LL_mat = np.full((len(reads), len(alleles)), -np.inf)
for r, alleles_ll in read_likelihoods.items():
    ri = read_to_idx[r]
    for a, ll in alleles_ll.items():
        ai = allele_to_idx[a]
        LL_mat[ri, ai] = ll
F_dict, posteriors, LL,zeft =run_em(read_likelihoods, tol=0.01, max_iters=100)
zeft_dict = defaultdict(list)
for f_iter in zeft:
    for a in f_iter:
        zeft_dict[a].append(f_iter[a])
gene_to_alleles = defaultdict(list)
for allele, freq in F_dict.items():
    gene = allele.split("*")[0]
    gene_to_alleles[gene].append((allele, freq))
for gene in gene_to_alleles:
    gene_to_alleles[gene] = sorted(gene_to_alleles[gene], key=lambda x: x[1], reverse=True)
with open(snakemake.output.allele_representatives_all_comp, "wb") as f:
    pickle.dump(gene_to_alleles, f)    
    
