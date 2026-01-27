rule cn:
    input:
        expand("{data_dir}/{sample}/cn.pkl", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

rule CN_inference:
    input:
        intron='/gpfs/project/yogah100/kir/resources/kir_reference_2025/annotations_mod_with_utr_extended_intron.bed',
        utr='/gpfs/project/yogah100/kir/resources/kir_reference_2025/annotations_mod_with_utr_extended_utr.bed', 
        fa='/gpfs/project/yogah100/kir/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta',   
        ref='/gpfs/project/yogah100/kir/resources/kir_reference_2025/ref_with_utr_extended.fa',
        exon="/gpfs/project/yogah100/kir/resources/kir_reference_2025/annotations_mod_with_utr_extended_exon.bed",
        bam2 = "{data_dir}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam",
        allele_representatives_all_comp= "{data_dir}/{sample}/alleles_cored_and_grouped_by_genes.pkl",
        bam_0_sort = "{data_dir}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.01_sort_all4.bam",
        bam_half_sort = "{data_dir}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.5_sort_all4.bam",
        list="/gpfs/project/yogah100/kir/resources/kir_reference_2025/Allelelist.txt",
        paired="{data_dir}/{sample}/paired_new_kir_sort_all4.bam",
        allele_rep="/gpfs/project/yogah100/kir/resources/kir_reference_2025/allele_representatives_new_kir.pkl"
    output:
        cn="{data_dir}/{sample}/cn.pkl",
        allele="{data_dir}/{sample}/cn_and_allele.pkl",
        log="{data_dir}/{sample}/script.log"     
    params:
        sample = "{sample}",
        out= "{data_dir}/{sample}",
        chr17="chr17" 
    script:
        "../scripts/Infer_CN.py" 