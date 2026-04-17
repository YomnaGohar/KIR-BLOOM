rule EM:
    input:
        expand("{DATA_DIR}/{sample}/alleles_cored_and_grouped_by_genes.pkl", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

rule EMU_with_new_KIR4:
    input:
        bam_path1='{DATA_DIR}/{sample}/mapped_filt_Chr17q25.bam',
        bam_path2="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
        list= config["Reference"]["Allele_list"],
        number_of_other="{DATA_DIR}/{sample}/selected_alleles.txt",
        exon=config["Reference"]["exon_bed"],
        fasta = config["Reference"]["KIR_alleles"],
        pkl= config["Reference"]["rep"]
    params:
        chr17="chr17"    
    output:
        allele_representatives_all_comp ="{DATA_DIR}/{sample}/alleles_scored_and_grouped_by_genes.pkl",
    script:
        "../scripts/EMU.py"  