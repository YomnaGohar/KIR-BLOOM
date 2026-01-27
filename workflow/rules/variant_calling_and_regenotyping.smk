rule var:
     input:
          expand("{data_dir}/{sample}/kir_mod_immu.gtf.gz", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
          expand("{data_dir}/{sample}/kir_mod.fa.fai", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

rule detect_var_new_kir2_no_utr4_per_error:
     input:
        a="{DATA_DIR}/{sample}/cn.pkl",
        b="{DATA_DIR}/{sample}/five_digit_allele_inference.pkl",
        list="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/Allelelist.txt",
        fa='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta',
        msa=expand('/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/msa_after_utr_extension/{gene}_updated_alignment.fasta', gene=GENES),    
        ref='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/ref_with_utr_extended.fa',
        intron='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_intron.bed',
        utr='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_utr.bed',    
        paired="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
        bam="{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam",
        exon="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_exon.bed",
        allele_rep="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/allele_representatives_new_kir.pkl"                       
     output:
        vcf="{DATA_DIR}/{sample}/kir_variants_in_exons.vcf" 
     params:
        sample = "{sample}" ,
        chr17="chr17"    
     script:   
       "../scripts/variant_detection_new_kir.py"   

rule modify_allele4_per_error:
     input:
        b="{DATA_DIR}/{sample}/five_digit_allele_inference.pkl",
        list="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/Allelelist.txt",
        ref='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta',
        vcf="{DATA_DIR}/{sample}/kir_variants_in_exons.vcf" 
     output:
        mod="{DATA_DIR}/{sample}/kir_mod.fa" 
     params:
        sample = "{sample}"   
     script:   
       "../scripts/correct_the_allele.py"  
      
rule index_kir_mod_4_per_error:
     input:
        mod="{DATA_DIR}/{sample}/kir_mod.fa" 
     output:
        mod="{DATA_DIR}/{sample}/kir_mod.fa.fai" 
     shell:
         """
         samtools faidx {input.mod}
         """     
rule immuannot_4_per_error:
      input:
        mod="{DATA_DIR}/{sample}/kir_mod.fa" 
      output:
        mod="{DATA_DIR}/{sample}/kir_mod_immu.gtf.gz"        
      params:
        sample = "{DATA_DIR}/{sample}/kir_mod_immu"   
      shell:  
       """
      bash /gpfs/project/yogah100/Software/Immuannot/scripts.pub.v3/immuannot.sh -c {input.mod} -r /gpfs/project/yogah100/Software/Immuannot/scripts.pub.v3/prepare-reference/Data-2025Dez28/ -o {params.sample}    
       """  