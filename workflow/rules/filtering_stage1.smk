rule filter1:
    input:
        expand("{data_dir}/{sample}/alleles_with_no_gaps.pkl", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

 
rule filter_alleles_with_gaps:
    input:
        sort="{data_dir}/{sample}/paired_new_kir_sort_all4.bam",
        bam  = "{data_dir}/{sample}/mapped_filt_Chr17q25.bam",
        fasta = '/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta',
        intron='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_intron.bed',
        utr='/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_utr.bed',
        exon="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_exon.bed",
    output:
        dict= "{data_dir}/{sample}/alleles_with_no_gaps.pkl"
    threads: 72
    script:
        "../scripts/filtering_stage1.py"  