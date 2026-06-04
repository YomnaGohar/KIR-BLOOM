rule var:
     input:
          expand("{DATA_DIR}/{sample}/kir_variants_in_exons.vcf", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
          #expand("{DATA_DIR}/{sample}/kir_mod_immu.gtf.gz", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
          expand("{DATA_DIR}/{sample}/kir_mod.fa.fai", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
          expand("{DATA_DIR}/{sample}/kir_mod_exon.bed", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
GENES = [
    "KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5A","KIR2DL5B",
    "KIR2DS1","KIR2DS2","KIR2DP1","KIR2DS3","KIR2DS4","KIR2DS5",
    "KIR3DL1","KIR3DL2","KIR3DL3","KIR3DP1","KIR3DS1"
]  
rule detect_var_new_kir2_no_utr4_per_error:
     input:
        a="{DATA_DIR}/{sample}/cn.pkl",#rules.CN_inference.output.cn,#"
        b="{DATA_DIR}/{sample}/five_digit_allele_inference.pkl",#rules.five_digit_inf.output.allele_local,#
        list=config["Reference"]["Allele_list"],
        fa=config["Reference"]["KIR_alleles"], 
        msa=expand(config["Reference"]["msa_path"] + "/{gene}_updated_alignment.fasta", gene=GENES),
        paired="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
        bam="{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam",
        exon=config["Reference"]["exon_bed"],
        allele_rep=config["Reference"]["rep"],                           
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
        list=config["Reference"]["Allele_list"],
        ref=config["Reference"]["KIR_alleles"], 
        vcf="{DATA_DIR}/{sample}/kir_variants_in_exons.vcf",
        exon=config["Reference"]["exon_bed"],
     output:
        mod="{DATA_DIR}/{sample}/kir_mod.fa" ,
        bed= "{DATA_DIR}/{sample}/kir_mod_exon.bed"
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
#rule immuannot_4_per_error:
#      input:
#        mod="{DATA_DIR}/{sample}/kir_mod.fa" 
#      output:
#        mod="{DATA_DIR}/{sample}/kir_mod_immu.gtf.gz"        
#      params:
#        sample = "{DATA_DIR}/{sample}/kir_mod_immu",
#        immu=   config["immuannot"] 
#      shell:  
#       """
#      bash {params.immu} -c {input.mod} -r /gpfs/project/yogah100/Software/Immuannot/scripts.pub.v3/prepare-reference/Data-2025Dez28/ -o {params.sample}    
#       """  