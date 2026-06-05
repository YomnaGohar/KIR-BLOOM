rule infer0:
     input:
        expand("{DATA_DIR}/{sample}/five_digit_allele_inference.pdf" ,DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"])
rule  five_digit_inf:         
    input:
        cn="{DATA_DIR}/{sample}/cn.pkl",
        allele="{DATA_DIR}/{sample}/cn_and_allele.pkl",
        fa=config["Reference"]["KIR_alleles"], 
        exon=config["Reference"]["exon_bed"],
        allele_representatives_all_comp= "{DATA_DIR}/{sample}/alleles_scored_and_grouped_by_genes.pkl",
        bam2 = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam",
        list=config["Reference"]["Allele_list"],
        paired="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
        allele_rep=config["Reference"]["rep"],             
    output:
        allele_local=temp("{DATA_DIR}/{sample}/five_digit_allele_inference.pkl"),
        tab="{DATA_DIR}/{sample}/five_digit_allele_inference.tsv",
        output_bam_2="{DATA_DIR}/{sample}/read_assignment.bam",
        pdf="{DATA_DIR}/{sample}/five_digit_allele_inference.pdf",
        log="{DATA_DIR}/{sample}/five_digit_allele_inference.log"
    params:
        sample = "{sample}",
        out= "{DATA_DIR}/{sample}/" ,
        chr17=config["background_region"].split(":")[0],
        start=int(config["background_region"].split(":")[1].split("-")[0]),
        end=int(config["background_region"].split(":")[1].split("-")[1])
    threads:  min(config["threads"], 4)  
    script:
        "../scripts/five_digit_allele_inference.py" 