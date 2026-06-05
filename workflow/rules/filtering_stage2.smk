rule filter2:
    input:
        expand("{DATA_DIR}/{sample}/selected_alleles.txt", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

GENES = [
    "KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5A","KIR2DL5B",
    "KIR2DS1","KIR2DS2","KIR2DP1","KIR2DS3","KIR2DS4","KIR2DS5",
    "KIR3DL1","KIR3DL2","KIR3DL3","KIR3DP1","KIR3DS1"
]      
rule pileup4:
  input:
   bam="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
   ref=config["Reference"]["KIR_alleles"]
  output:
    "{DATA_DIR}/{sample}/mpileup4.txt"
  params:
    samtools_bin=config["samtools_path"]
  threads: min(config["threads"], 5)
  shell:  
   """
   perl workflow/scripts/extractPileUpFrequencies_v2.pl --samtools_bin {params.samtools_bin} --BAM {input.bam} --outputFile {output} --minMappingQuality 0 --minBaseQuality 0 --referenceGenome {input.ref}    
   """
rule filter_alleles_with_no_variant_position_support:
     input:
          dict= "{DATA_DIR}/{sample}/alleles_with_no_gaps.pkl",
          bam1="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
          index="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam.bai",
          fasta =config["Reference"]["KIR_alleles"],
          msa = expand(config["Reference"]["msa_path"] + "/{gene}_updated_alignment.fasta", gene=GENES),
          map = expand(config["Reference"]["msa_mappings_path"] + "/{gene}_msa_to_original.pkl", gene=GENES),
          list= config["Reference"]["Allele_list"],
          bed_exon=config["Reference"]["exon_bed"],
          pileup="{DATA_DIR}/{sample}/mpileup4.txt",
          bed_intron=config["Reference"]["intron_bed"],
          bed_utr=config["Reference"]["utr_bed"],
     output:
          number_of_other="{DATA_DIR}/{sample}/selected_alleles.txt",
          #intron="{DATA_DIR}/{sample}/intron_count_new_kir_all3.pkl",
          #exon="{DATA_DIR}/{sample}/exon_count_new_kir_all3.pkl",
     threads: min(config["threads"], 5)  
     script:
          "../scripts/filtering_stage2.py"    
          