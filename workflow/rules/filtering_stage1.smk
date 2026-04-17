rule filter1:
    input:
        expand("{DATA_DIR}/{sample}/alleles_with_no_gaps.pkl", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

 
rule filter_alleles_with_gaps:
    input:
        sort="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
        bam  = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25.bam",
        fasta = config["Reference"]["KIR_alleles"],
        intron=config["Reference"]["intron_bed"],
        utr=config["Reference"]["utr_bed"],
        exon= config["Reference"]["exon_bed"],
    output:
        dict= "{DATA_DIR}/{sample}/alleles_with_no_gaps.pkl"
    threads: 72
    script:
        "../scripts/filtering_stage1.py"  