#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
from Bio import SeqIO
import pickle
from collections import defaultdict

open(snakemake.output.bed, "w").close()
# 1) Load allele_dict_2
l = snakemake.input.list
with open(l) as f:
    allele_dict = {
        "KIR:" + line.split(",")[0] if "KIR" in  line.split(",")[0] else line.split(",")[0] : line.strip().split(",")[1]
        for line in f if not line.startswith("#") and line.strip()}
with open(l) as f:
    allele_dict_2 = {
        line.strip().split(",")[1]:
        ("KIR:" + line.split(",")[0] if "KIR" in line.split(",")[0] else line.split(",")[0])
        for line in f if not line.startswith("#") and line.strip()
    }
def update_exon_bed_coordinates(exons_bed, allele_name,variants, allele_dict_2,nn, output_bed):
    with open(output_bed, "a") as out:
            mapped_allele_name = allele_dict_2[allele_name]
            variants = sorted(variants, key=lambda r: r.pos)

            updated_intervals = []
            cumulative_shift = 0
            var_idx = 0

            original_intervals = sorted(exons_bed[mapped_allele_name])

            for start, end in original_intervals:
                # add all indel effects before exon start
                while var_idx < len(variants) and (variants[var_idx].pos - 1) < start:
                    rec = variants[var_idx]
                    ref = rec.ref
                    alt = rec.alts[0]
                    if len(ref) != len(alt):   # only indels shift coordinates
                        cumulative_shift += len(alt) - len(ref)
                    var_idx += 1

                new_start = start + cumulative_shift
                new_end = end + cumulative_shift
                temp_idx = var_idx
                internal_shift = 0
                while temp_idx < len(variants) and (variants[temp_idx].pos - 1) < end:
                    rec = variants[temp_idx]
                    ref = rec.ref
                    alt = rec.alts[0]
                    if len(ref) != len(alt):
                        internal_shift += len(alt) - len(ref)
                    temp_idx += 1
                new_end += internal_shift
                updated_intervals.append((new_start, new_end))
            for new_start, new_end in updated_intervals:
                out.write(f"{nn}\t{new_start}\t{new_end}\n")   
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
exons_bed =parse_exons_bed(bed_file_exon)                  
vcf_path = snakemake.input.vcf
msa_path =  snakemake.input.ref
ref_sequences = {rec.id: rec for rec in SeqIO.parse(msa_path, "fasta")}
vcf = pysam.VariantFile(vcf_path)

allele_variants = {}
for rec in vcf:
    allele_name = rec.chrom
    allele_variants.setdefault(allele_name, []).append(rec)

with open(snakemake.input.b,"rb") as f: #29_allele_local
        gene_allele_level_new=pickle.load( f)  
genotype={allele_dict[a]:gene_allele_level_new[a] for a in    gene_allele_level_new  if gene_allele_level_new[a] > 0.01 }    
output_fasta = snakemake.output.mod#f"{mapped_allele_name}_modified.fasta"
observed=[]
with open(output_fasta, "w") as out:
    for allele_name, variants in allele_variants.items():
        mapped_allele_name = allele_dict_2[allele_name]
        observed.append(allele_name)
        sequence = list(str(ref_sequences[mapped_allele_name].seq))
        for rec in variants:
            pos = rec.pos - 1
            ref = rec.ref
            alt = rec.alts[0]
            if len(ref) == 1 and len(alt) == 1:        # SNP
                sequence[pos] = alt

            elif len(ref) == 1 and len(alt) > 1:       # Insertion
                sequence[pos] = alt[0]
                sequence[pos+1:pos+1] = list(alt[1:])

            elif len(ref) > 1 and len(alt) == 1:       # Deletion
                sequence[pos] = alt
                del sequence[pos+1 : pos + len(ref)]
        for count in range(genotype[allele_name]):
            modified_seq = "".join(sequence)
            nn= f"{allele_name}_modified_{count+1}"
            out.write(f">{nn}\n")
            out.write(modified_seq + "\n")
            update_exon_bed_coordinates(exons_bed, allele_name,variants, allele_dict_2,nn, snakemake.output.bed)
    for a in   genotype:
        count= genotype[a]
        if a not in observed:
            sequence = str(ref_sequences[allele_dict_2[a]].seq)
            for c in range(count):
                nn=f"{a}_modified_{c+1}"
                out.write(f">{nn}\n")
                out.write(sequence + "\n")
                update_exon_bed_coordinates(exons_bed, a,[], allele_dict_2,nn, snakemake.output.bed)
        
