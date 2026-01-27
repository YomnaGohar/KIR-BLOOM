rule EM:
    input:
        expand("{data_dir}/{sample}/alleles_cored_and_grouped_by_genes.pkl", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),

rule EMU_with_new_KIR4:
    input:
        bam_path1='{DATA_DIR}/{sample}/mapped_filt_Chr17q25.bam',
        bam_path2="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
        list="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/Allelelist.txt",
        number_of_other="{DATA_DIR}/{sample}/selected_alleles.txt",
        bed="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended.bed",
        exon="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/annotations_mod_with_utr_extended_exon.bed",
        fasta = ('/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta'),
        pkl="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/allele_representatives_new_kir.pkl"
    params:
        d=directory("{DATA_DIR}/{sample}/cluster"),
        chr17="chr17"    
    output:
        #allele_representatives_all="{DATA_DIR}/{sample}/allele_representatives_all_new_kir_all4.pkl",
        allele_representatives_all_comp ="{DATA_DIR}/{sample}/alleles_cored_and_grouped_by_genes.pkl",
    script:
        "../scripts/EMU.py"  