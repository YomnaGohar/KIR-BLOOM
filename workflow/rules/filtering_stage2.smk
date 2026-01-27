rule filter2:
    input:
        expand("{data_dir}/{sample}/selected_alleles.txt", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

GENES = [
    "KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5A","KIR2DL5B",
    "KIR2DS1","KIR2DS2","KIR2DP1","KIR2DS3","KIR2DS4","KIR2DS5",
    "KIR3DL1","KIR3DL2","KIR3DL3","KIR3DP1","KIR3DS1"
]        
rule pileup4:
  input:
   bam="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
   ref='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta'
  output:
    "{DATA_DIR}/{sample}/mpileup4.txt"
  params:
    samtools_bin="/home/yogah100/miniforge3/envs/medaka2/bin/samtools"
  threads: 10 
  shell:  
   """
   perl workflow/scripts/extractPileUpFrequencies_v2.pl --samtools_bin {params.samtools_bin} --BAM {input.bam} --outputFile {output} --minMappingQuality 0 --minBaseQuality 0 --referenceGenome {input.ref}    
   """
rule filter_alleles_with_no_variant_position_support:
     input:
          dict= "{DATA_DIR}/{sample}/alleles_with_no_gaps.pkl",
          bam1="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
          index="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam.bai",
          fasta = ('/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta'),
          msa=expand('/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/msa_after_utr_extension/{gene}_updated_alignment.fasta', gene=GENES),
          map=expand('/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/msa_after_utr_extension/msa_mappings/{gene}_msa_to_original.pkl', gene=GENES),
          list="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/Allelelist.txt",
          bed_exon="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_exon.bed",
          pileup="{DATA_DIR}/{sample}/mpileup4.txt",
          bed_intron='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_intron.bed',
          bed_utr='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_utr.bed'
     output:
          number_of_other="{DATA_DIR}/{sample}/selected_alleles.txt"
     threads: 20     
     script:
          "../scripts/filtering_stage2.py"    
          