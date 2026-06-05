rule cn:
    input:
        expand("{DATA_DIR}/{sample}/cn.pkl", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

rule CN_inference:
    input:
        fa=config["Reference"]["KIR_alleles"],   
        exon=config["Reference"]["exon_bed"],
        bam2 = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam",
        allele_representatives_all_comp= "{DATA_DIR}/{sample}/alleles_scored_and_grouped_by_genes.pkl",
        list=config["Reference"]["Allele_list"],
        paired="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
        allele_rep=config["Reference"]["rep"],
    output:
        cn="{DATA_DIR}/{sample}/cn.pkl",
        tab="{DATA_DIR}/{sample}/cn.tsv",
        allele="{DATA_DIR}/{sample}/cn_and_allele.pkl",
        log="{DATA_DIR}/{sample}/cn.log"     
    params:
        sample = "{sample}",
        out= "{DATA_DIR}/{sample}",
        chr17=config["background_region"].split(":")[0],
        start=int(config["background_region"].split(":")[1].split("-")[0]),
        end=int(config["background_region"].split(":")[1].split("-")[1])
    script:
        "../scripts/Infer_CN.py" 