#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 01:34:47 2025
"""
import numpy as np
from Bio import SeqIO, AlignIO
import pysam
from collections import defaultdict, Counter
import pickle
from scipy.stats import poisson
import string
with open(snakemake.input.dict, "rb") as f:
    allele_fully_covered = pickle.load(f)
alleles_to_keep = set(allele_fully_covered.keys())
msa_files = snakemake.input.msa
seq_dict = {}
for msa_file in msa_files:
    alignment = AlignIO.read(msa_file, "fasta")
    for rec in alignment:
        seq_dict[rec.id.replace("KIR","")] = str(rec.seq)

original_seqs = {rec.id: str(rec.seq).replace("-", "") for rec in SeqIO.parse(snakemake.input.fasta, "fasta")}

bam_path = snakemake.input.bam1

mapping_paths = list(snakemake.input.map)  # make sure it's a list, not a string
all_mappings = {}
for mapping_path in mapping_paths:
    with open(mapping_path, "rb") as f:
        msa_map = pickle.load(f)
    all_mappings.update(msa_map)

def get_bases_at_position_nodup(bam_path, ref, pos):
    bases = []
    seen_templates = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for col in bam.pileup(
            ref, pos, pos+1,
            truncate=True,
            stepper="all",
            min_base_quality=0,
            min_mapping_quality=0,
            flag_filter=0
        ):
            if col.reference_pos != pos:
                continue
            for pr in col.pileups:
                read = pr.alignment
                qname = read.query_name
                if qname in seen_templates:
                    continue
                seen_templates.add(qname)
                if pr.is_del:
                    bases.append('-')
                elif pr.is_refskip:
                    bases.append('<')
                elif pr.query_position is not None:
                    bases.append(read.query_sequence[pr.query_position])
    return bases


# Extract divergent positions and assess presence based on Poisson p-value
l=snakemake.input.list
with open(l) as f:
    allele_dict = {
        "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] : line.strip().split(",")[1]
        for line in f if not line.startswith("#") and line.strip()}
with open(l) as f:
    allele_dict_2 = {
     line.strip().split(",")[1] :  "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] 
        for line in f if not line.startswith("#") and line.strip()} 
    
allele_fully_covered_names = {allele_dict[a]: cluster for a, cluster in allele_fully_covered.items()}
exon_regions = defaultdict(list)
with open(snakemake.input.bed_exon) as bed:
    for line in bed:
        if line.startswith('#') or not line.strip():
            continue
        chrom, start, end = line.strip().split()[:3]
        exon_regions[chrom].append((int(start), int(end)))

def in_exon(chrom, pos):
    """Return True if 0-based pos is within any exon interval on chrom"""
    for s, e in exon_regions.get(chrom, []):
        if s <= pos < e:
            return True
    return False
Intron_regions=defaultdict(list)
with open(snakemake.input.bed_intron) as bed:
    for line in bed:
        if line.startswith('#') or not line.strip():
            continue
        chrom, start, end = line.strip().split()[:3]
        Intron_regions[chrom].append((int(start), int(end)))
def in_intron(chrom, pos):
    """Return True if 0-based pos is within any exon interval on chrom"""
    for s, e in Intron_regions.get(chrom, []):
        if s <= pos < e:
            return True
    return False        
UTR_regions=defaultdict(list)
with open(snakemake.input.bed_utr) as bed:
    for line in bed:
        if line.startswith('#') or not line.strip():
            continue
        chrom, start, end = line.strip().split()[:3]
        UTR_regions[chrom].append((int(start), int(end)))
def in_utr(chrom, pos):
    """Return True if 0-based pos is within any exon interval on chrom"""
    for s, e in UTR_regions.get(chrom, []):
        if s <= pos < e:
            return True
    return False         
def map_pos_between_alleles(allele1, allele2, target_pos1):
    aln1 = seq_dict[allele1.upper()]
    aln2 = seq_dict[allele2.upper()]
    pos1 = pos2 = 0
    for i in range(len(aln1)):
        # if we've reached the target in allele1's ungapped seq
        if pos1 == target_pos1:
            return pos2 if pos2 < len(orig2) else None
        if aln1[i] != '-':
            pos1 += 1
        if aln2[i] != '-':
            pos2 += 1
    return None
from collections import defaultdict
def parse_pileup(pileup_path):
    pileup_data = defaultdict(dict)
    with open(pileup_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            allele = allele_dict[parts[0]]
            pos = int(parts[1])-1
            ref_base = parts[2].upper()
            counts_str = parts[3]

            ref_count = 0
            for token in counts_str.split(";"):
                if "=" not in token:
                    continue

                base_raw, count_str = token.split("=")
                if not base_raw or not count_str:
                    continue

                base_clean = base_raw[0].upper()  # first character = actual base
                try:
                    count = int(count_str)
                except ValueError:
                    continue
                # add only counts matching reference base (case-insensitive)
                if base_clean == ref_base:
                    ref_count += count
            pileup_data[allele][pos] = ref_count
    return pileup_data


pile=parse_pileup(snakemake.input.pileup)
def regroup_clusters_by_gene(cluster_to_alleles, extract_gene_name):
    # Step 1: group alleles by gene across ALL clusters
    gene_to_alleles = defaultdict(list)
    for cluster, alleles in cluster_to_alleles.items():
        for allele in alleles:
            gene = extract_gene_name(allele)
            gene_to_alleles[gene].append(allele)

    # Step 2: assign new cluster IDs, one per gene
    new_allele_to_cluster = {}
    new_cluster_to_alleles = {}

    for i, (gene, alleles) in enumerate(gene_to_alleles.items(), start=1):
        cluster_id = str(i)  # could also keep gene name if you prefer
        new_cluster_to_alleles[cluster_id] = alleles
        for allele in alleles:
            new_allele_to_cluster[allele] = cluster_id

    return new_allele_to_cluster, new_cluster_to_alleles

def extract_gene_name(allele_name):
    return allele_name.split("*")[0].split(":")[-1][:5]
cluster_to_alleles = defaultdict(list)
for allele, cluster in allele_fully_covered_names.items():
    cluster_to_alleles[cluster].append(allele)
new_allele_to_cluster, new_cluster_to_allele=regroup_clusters_by_gene(cluster_to_alleles, extract_gene_name)
significance_threshold = 0.05
c=0
allele_position_presence_exons = defaultdict(dict)
allele_position_pvalues_exons = defaultdict(dict)
allele_position_presence_introns = defaultdict(dict)
allele_position_pvalues_intons = defaultdict(dict)
allele_position_presence_utrs = defaultdict(dict)
allele_position_pvalues_utrs = defaultdict(dict)
seen_positions = set()
from itertools import combinations
for cluster, alleles in new_cluster_to_allele.items():
    for allele1, allele2 in combinations(alleles, 2):
            pair_key = tuple(sorted([allele1, allele2]))
            if pair_key in seen_positions:
                continue
            seq1 = seq_dict[allele1.upper()]
            seq2 = seq_dict[allele2.upper()]
            orig1 = original_seqs.get(allele_dict_2[allele1], "")
            orig2 = original_seqs.get(allele_dict_2[allele2], "")
            beg10_1 = len(orig1) * 0.15
            end90_1 = len(orig1) - len(orig1) * 0.15
            beg10_2 = len(orig2) * 0.15
            end90_2 = len(orig2) - len(orig2) * 0.15
            s1 = np.frombuffer(seq1.encode("ascii"), dtype="S1")
            s2 = np.frombuffer(seq2.encode("ascii"), dtype="S1")
            valid_mask = (s1 != b'-') & (s2 != b'-')
            diff_mask = (s1 != s2) & valid_mask
            divergent_positions = np.where(diff_mask)[0]           
            for i in divergent_positions:
                pos1= all_mappings[allele1.upper()][i]
                pos2= all_mappings[allele2.upper()][i]
                nt1=orig1[pos1]
                nt2=orig2[pos2]
                if (pos1 < len(orig1)) and (pos2 < len(orig2)):
                        if in_exon(allele_dict_2[allele1], pos1) and in_exon(allele_dict_2[allele2], pos2):
                            tag1 = f"pos_{pos1}_{nt1}"
                            if tag1 not in allele_position_presence_exons[allele1] and orig1[pos1].isupper():
                                if pos1 in pile[allele1]:
                                   count = int(pile[allele1][pos1])
                                else:
                                    count= 0                          
                                p_val = 1 - poisson.cdf(count - 1, mu=1)
                                allele_position_presence_exons[allele1][tag1] = count
                                allele_position_pvalues_exons[allele1][tag1] = int(p_val < 0.05) if beg10_1 < pos1 < end90_1 else  int(p_val < 0.1) 
                            tag2 = f"pos_{pos2}_{nt2}"
                            if tag2 not in allele_position_presence_exons[allele2] and orig2[pos2].isupper():
                                if pos2 in pile[allele2]:
                                    count = int(pile[allele2][pos2])
                                else: 
                                    count = 0    
                                p_val = 1 - poisson.cdf(count - 1, mu=1)
                                allele_position_presence_exons[allele2][tag2] = count
                                allele_position_pvalues_exons[allele2][tag2] = int(p_val < 0.05) if beg10_2 < pos2 < end90_2 else  int(p_val < 0.1)
                        elif in_intron(allele_dict_2[allele1], pos1) and in_intron(allele_dict_2[allele2], pos2):
                            tag1 = f"pos_{pos1}_{nt1}"
                            if tag1 not in allele_position_presence_introns[allele1] and orig1[pos1].isupper():
                                if pos1 in pile[allele1]:
                                   count = int(pile[allele1][pos1])
                                else:
                                    count= 0   
                                p_val = 1 - poisson.cdf(count - 1, mu=1)
                                allele_position_presence_introns[allele1][tag1] = count
                                allele_position_pvalues_intons[allele1][tag1] = int(p_val < 0.05) if beg10_1 < pos1 < end90_1 else  int(p_val < 0.1)
                            tag2 = f"pos_{pos2}_{nt2}"
                            if tag2 not in allele_position_presence_introns[allele2] and orig2[pos2].isupper():
                                if pos2 in pile[allele2]:
                                    count = int(pile[allele2][pos2])
                                else: 
                                    count = 0 
                                p_val = 1 - poisson.cdf(count - 1, mu=1)
                                allele_position_presence_introns[allele2][tag2] = count
                                allele_position_pvalues_intons[allele2][tag2] = int(p_val < 0.05) if beg10_2 < pos2 < end90_2 else  int(p_val < 0.1)
                        elif in_utr(allele_dict_2[allele1], pos1) and in_utr(allele_dict_2[allele2], pos2):
                            tag1 = f"pos_{pos1}_{nt1}"
                            if tag1 not in allele_position_presence_utrs[allele1] and orig1[pos1].isupper():
                                if pos1 in pile[allele1]:
                                   count = int(pile[allele1][pos1])
                                else:
                                    count= 0                                  
                                p_val = 1 - poisson.cdf(count - 1, mu=1)
                                allele_position_presence_utrs[allele1][tag1] = count
                                allele_position_pvalues_utrs[allele1][tag1] =  int(p_val < 0.05) if beg10_1 < pos1 < end90_1 else  int(p_val < 0.1)
                            tag2 = f"pos_{pos2}_{nt2}"
                            if tag2 not in allele_position_presence_utrs[allele2] and orig2[pos2].isupper():
                                if pos2 in pile[allele2]:
                                    count = int(pile[allele2][pos2])
                                else: 
                                    count = 0 
                                p_val = 1 - poisson.cdf(count - 1, mu=1)
                                allele_position_presence_utrs[allele2][tag2] = count
                                allele_position_pvalues_utrs[allele2][tag2] = int(p_val < 0.05) if beg10_2 < pos2 < end90_2 else  int(p_val < 0.1)                                     
            seen_positions.add(pair_key)

loses = {}
for allele, pvals in allele_position_pvalues_exons.items():
    count_exons = sum(1 for pv in pvals.values() if pv == 0)
    count_intron= sum(1 for pv in allele_position_pvalues_intons[allele] if allele_position_pvalues_intons[allele][pv] == 0) if  allele_position_pvalues_intons[allele] else None
    count_utrs= sum(1 for pv in allele_position_pvalues_utrs[allele] if allele_position_pvalues_utrs[allele][pv] == 0) if allele_position_pvalues_utrs[allele] else None
    loses[allele] = (count_exons,  count_intron,count_utrs)



selected = []
gene_to_alleles = defaultdict(list)
for allele in loses:
    gene = extract_gene_name(allele)
    gene_to_alleles[gene].append(allele)
for gene, alleles in gene_to_alleles.items():
    min_exon = min(loses[a][0] for a in alleles)
    candidates = [a for a in alleles if loses[a][0] == min_exon]
    code_to_alleles = defaultdict(list)
    for a in candidates:
        code = a.split("*")[1][:5]
        code_to_alleles[code].append(a)
    for code, group in code_to_alleles.items():
        exon_only = [a for a in group if loses[a][1] is None]
        selected.extend(exon_only)
        intronal = [a for a in group if loses[a][1] is not None]
        if intronal:
            best_intron = min( loses[a][1] for a in intronal )
            a_ = [a for a in intronal if loses[a][1] == best_intron]
            selected.extend(a_) 

for a in new_cluster_to_allele:
    if len(new_cluster_to_allele[a]) ==1:
        selected.extend(new_cluster_to_allele[a])
        
out_path =snakemake.output.number_of_other
with open(out_path, "w") as out:
    for allele in selected:
        out.write(f"{allele}\n")
 
