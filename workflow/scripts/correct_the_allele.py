#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
from Bio import SeqIO
import pickle
# 1) Load allele_dict_2
l = snakemake.input.list#""/home/yomna/hpc_project/kir/resources/Allelelist_new.txt"#
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
        
vcf_path = snakemake.input.vcf# "/home/yomna/hpc_project/kir/real_data/ERR3989060/cn_optimization/29_new_kir_all_no_utr4/kir_variants_in_exons.vcf"
msa_path =  snakemake.input.ref#"/home/yomna/hpc_project/kir/resources/kir_reference/kir_gen_new_mod_with_utr_extended.fasta"#
ref_sequences = {rec.id: rec for rec in SeqIO.parse(msa_path, "fasta")}
vcf = pysam.VariantFile(vcf_path)

allele_variants = {}
for rec in vcf:
    allele_name = rec.chrom
    allele_variants.setdefault(allele_name, []).append(rec)

# if len(allele_variants) == 0:
#     output_fasta = snakemake.output.mod
#     with open(output_fasta, "w") as out:
#         pass   # empty file
#     print("VCF empty → wrote empty FASTA:", output_fasta)
#     exit(0)
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
            out.write(f">{allele_name}_modified_{count+1}\n")
            out.write(modified_seq + "\n")
    for a in   genotype:
        count= genotype[a]
        if a not in observed:
            sequence = str(ref_sequences[allele_dict_2[a]].seq)
            for c in range(count):
                out.write(f">{a}_modified_{c+1}\n")
                out.write(sequence + "\n")