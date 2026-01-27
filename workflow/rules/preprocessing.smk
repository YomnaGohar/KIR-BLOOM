rule extract:
    input:
        expand("{data_dir}/{sample}/mapped.bam", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
        expand("{data_dir}/{sample}/mapped_filt.read1_mod.fastq", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
        expand("{data_dir}/{sample}/mapped_filt.read1.fastq", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
        expand("{data_dir}/{sample}/paired_new_kir_sort_all4.bam.bai", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
        expand("{data_dir}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.0001_sort_all4.bam", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
        expand("{data_dir}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.01_sort_all4.bam", data_dir= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"])


rule extract_reads:
    input:
        reads=lambda wildcards: config["Samples"]["sample_fastqs"][wildcards.sample],
        ref=config["Reference"]["fasta"],
        bed=config["Reference"]["KIR_regions_bed"]
    output:
        bam="{data_dir}/{sample}/mapped.bam"
    threads: 30
    params:
        tmpdir="{data_dir}/{sample}/tmp"    
    run:
        import os

        reads = input.reads
        if reads[0].endswith(".cram"):
            shell("""
                mkdir -p {params.tmpdir}
                # 1. Reads overlapping KIR regions
                samtools view -@ {threads} -b -T {input.ref} -L {input.bed} {reads} \
                > {output.bam}.kir.tmp
                samtools index {output.bam}.kir.tmp

                # 2. Read names from KIR-overlapping reads
                samtools view {output.bam}.kir.tmp \
                | cut -f1 | sort -u > {output.bam}.names.tmp

                # 3. Extract all reads with those names
                samtools view -@ {threads} -b -T {input.ref} \
                -N {output.bam}.names.tmp {reads} > {output.bam}.kir.named.bam

                # 4. Extract unmapped reads
                samtools view -@ {threads} -b -T {input.ref} -f 4 {reads} \
                > {output.bam}.unmapped.bam

                # 5. Merge + sort
                samtools merge -@ {threads} -f {output.bam}.merged.bam \
                {output.bam}.kir.named.bam {output.bam}.unmapped.bam

                samtools sort -@ {threads} -T {params.tmpdir} -o {output.bam} {output.bam}.merged.bam

                rm -f {output.bam}.kir.tmp* {output.bam}.names.tmp \
                    {output.bam}.kir.named.bam {output.bam}.unmapped.bam \
                    {output.bam}.merged.bam
            """)
        elif reads[0].endswith(".fastq") or reads[0].endswith(".fq") or reads[0].endswith(".fastq.gz"):
            bwa = config["bwa_path"]
            shell("""
                set -euo pipefail
                mkdir -p {params.tmpdir}
                # 0. Align reads
                {bwa} mem -t {threads} {input.ref} {reads} > {output.bam}.align.sam
                samtools view -@ {threads} -b {output.bam}.align.sam \
                > {output.bam}.align.tmp

                # 1. Reads overlapping KIR regions
                samtools view -@ {threads} -b -L {input.bed} {output.bam}.align.tmp \
                > {output.bam}.kir.tmp
                samtools sort -T {params.tmpdir} {output.bam}.kir.tmp > {output.bam}.sort.kir.tmp
                samtools index {output.bam}.sort.kir.tmp

                # 2. Read names from KIR-overlapping reads
                samtools view {output.bam}.sort.kir.tmp | cut -f1 > {output.bam}.names.raw.tmp
                sort -u {output.bam}.names.raw.tmp > {output.bam}.names.tmp

                # 3. Extract all reads with those names
                samtools view -@ {threads} -hb \
                -N {output.bam}.names.tmp {output.bam}.sort.kir.tmp \
                > {output.bam}.kir.named.bam

                # 4. Extract unmapped reads
                samtools view -@ {threads} -hb -f 4 {output.bam}.align.tmp \
                > {output.bam}.unmapped.bam

                # 5. Merge + sort
                samtools merge -@ {threads} -f {output.bam}.merged.bam \
                {output.bam}.kir.named.bam {output.bam}.unmapped.bam

                samtools sort -@ {threads} -T {params.tmpdir} -o {output.bam} {output.bam}.merged.bam
                samtools index {output.bam}

                # 6. Cleanup
                rm -f {output.bam}.align.sam \
                    {output.bam}.align.tmp \
                    {output.bam}.kir.tmp* \
                    {output.bam}.names.raw.tmp \
                    {output.bam}.names.tmp \
                    {output.bam}.kir.named.bam \
                    {output.bam}.unmapped.bam \
                    {output.bam}.merged.bam \
                    {output.bam}.sort.kir.tmp \
                    {output.bam}.sort.kir.tmp
            """)
rule remove_dup_in_bam_files: 
    input:
       bam="{data_dir}/{sample}/mapped.bam"
    output:
       bam="{data_dir}/{sample}/dedup_mapped_noRG.bam",
    shell:
        """
        samtools view -h {input.bam}|  awk 'BEGIN{{OFS="\t"}} /^@/ {{print; next}} !seen[$1,$2,$3,$4,$5,$6,$7,$8]++'  | sed 's/\tRG:Z:[^\t]*//g' | sed 's/\PG:Z:[^\t]*//g' | samtools view -b -o {output.bam}
        """       
rule estimate_insert_size_new_kir:
     input:
          bam="{DATA_DIR}/{sample}/dedup_mapped_noRG.bam",
     output:
         bam="{DATA_DIR}/{sample}/mapped_filt_Chr17q25.bam",    
         kir="{DATA_DIR}/{sample}/mapped_filt_noChr17.bam",
     shell:
         """
         samtools index {input.bam}
         samtools view -b {input.bam} chr17:74000000-76000000 > {output.bam}
         samtools index {output.bam}
         samtools view -b -e 'rname!="chr17"' {input.bam} > {output.kir}
         samtools index {output.kir}
         """          
rule bam_to_fastq:
    input:
        bam="{DATA_DIR}/{sample}/mapped_filt_noChr17.bam"
    output:
        read1="{DATA_DIR}/{sample}/mapped_filt.read1.fastq",
        read2="{DATA_DIR}/{sample}/mapped_filt.read2.fastq",
        sing="{DATA_DIR}/{sample}/singletons.fastq"
    threads: 20
    shell:
        """
        picard SamToFastq I={input.bam} F={output.read1}  F2={output.read2} FU={output.sing}  VALIDATION_STRINGENCY=SILENT
        """
rule modify_fastq_naming:
    input:
        read1="{DATA_DIR}/{sample}/mapped_filt.read1.fastq",
        read2="{DATA_DIR}/{sample}/mapped_filt.read2.fastq",
    output:
        read1="{DATA_DIR}/{sample}/mapped_filt.read1_mod.fastq",
        read2="{DATA_DIR}/{sample}/mapped_filt.read2_mod.fastq",
    shell:
        """
        sed '1~4s/\\/1\\b//' {input.read1} > {output.read1}
        sed '1~4s/\\/2\\b//' {input.read2} > {output.read2}
        """
rule filter_for_the_second_time:
    input:
        fastq1="{data_dir}/{sample}/mapped_filt.read1_mod.fastq",
        fastq2="{data_dir}/{sample}/mapped_filt.read2_mod.fastq",  
        reference="/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/ref_with_utr_extended.fa",
    output:
        sam=temp("{data_dir}/{sample}/mapped_the_second_time.sam")
    threads: 72
    shell:
        """
        bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} > {output.sam}
        """  
rule sam_to_bam2:
    input:
         "{data_dir}/{sample}/mapped_the_second_time.sam"
    output:
        "{data_dir}/{sample}/mapped_the_second_time.bam"
    threads: 72
    shell:
        "samtools view -Sb {input} > {output}"

rule filter_unmapped_new_kir2_keep_utrs:
    input:
        bam = "{data_dir}/{sample}/mapped_the_second_time.bam",
        bed = "/gpfs/project/yogah100/kir/resources/kir_reference_2025/kir_allele_names_new_kir_only.bed",   # ROI
    output:
        bam = "{data_dir}/{sample}/mapped_filt_new_kir2.bam",
        roi = temp("{data_dir}/{sample}/roi_tmp2.bam"),
    threads: 72
    shell:
        """
        samtools view -b -L {input.bed} {input.bam} > {output.roi}
        samtools sort -@ {threads} -o {output.bam} {output.roi}
        samtools index -@ {threads} {output.bam}
        """        
rule extract_mapped_reads_with_seqkit_new_kir2:
    input:
        bam="{data_dir}/{sample}/mapped_filt_new_kir2.bam",
        fq1="{data_dir}/{sample}/mapped_filt.read1_mod.fastq",
        fq2="{data_dir}/{sample}/mapped_filt.read2_mod.fastq"
    output:
        read1="{data_dir}/{sample}/mapped_filt.read1_new_kir2.fastq",
        read2="{data_dir}/{sample}/mapped_filt.read2_new_kir2.fastq",
        names=temp("{data_dir}/{sample}/mapped_names_new_kir2.txt"),
    threads: 20
    shell:
        """
        samtools view {input.bam} | awk '{{print $1}}' | sort -u > {output.names}
        seqkit grep -f {output.names} --threads {threads} {input.fq1} > {output.read1}
        seqkit grep -f {output.names} --threads {threads} {input.fq2} > {output.read2}
        """                     
rule proper_mapping_with_new_KIR_4:
    input:
        read1="{DATA_DIR}/{sample}/mapped_filt.read1_new_kir2.fastq",
        read2="{DATA_DIR}/{sample}/mapped_filt.read2_new_kir2.fastq",
        mmi='/gpfs/project/yogah100/kir/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta',
    output:
        bam1="{DATA_DIR}/{sample}/remapped1_new_kir4.bam",
        bam2="{DATA_DIR}/{sample}/remapped2_new_kir4.bam"
    threads: 72
    shell:
        """
        minimap2 -ax sr --secondary-seq --MD --eqx --secondary=yes -N 2000 -t {threads} {input.mmi} {input.read1} | samtools view -b | samtools view -hF 4 | samtools sort > {output.bam1}
        samtools index {output.bam1}
        minimap2 -ax sr  --secondary-seq --MD --eqx --secondary=yes -N 2000 -t {threads} {input.mmi} {input.read2} | samtools view -b | samtools view -hF 4 | samtools sort > {output.bam2}
        samtools index {output.bam2}
        """       
rule pair_with_new_KIR4:
    input:
        bam1 = "{DATA_DIR}/{sample}/remapped1_new_kir4.bam",
        bam2 = "{DATA_DIR}/{sample}/remapped2_new_kir4.bam",
        bam  = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25.bam",
        fasta = "/gpfs/project/yogah100/cram_files/resources/kir_reference_2025/kir_gen_new_mod_with_utr_extended.fasta"
    output:
        bam1 = "{DATA_DIR}/{sample}/paired_new_kir_all4.bam",
        bam2 = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_all4.bam",
    threads: 72
    params:
          chr17="chr17"
    script:
        "../scripts/pairing.py"  
rule index_with_new_KIR4:
     input:
          "{DATA_DIR}/{sample}/paired_new_kir_all4.bam"
     output:
          sort="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam",
          index="{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam.bai"
     shell:
          """
          samtools sort {input} > {output.sort}
          samtools index {output.sort}
          """                 
rule sort_the_tag_and_index_with_new_KIR4:
     input:
         bam2 = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_all4.bam", 
     output:
         bam2 = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam",  
     shell:
         """
         samtools sort {input.bam2} > {output.bam2}
         samtools index {output.bam2}
         """ 
rule drop_out_with_new_KIR4:
     input:
         bam2 = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam"
     output:
         bam_half = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.5_all4.bam",
         bam_half_sort = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.5_sort_all4.bam",
         bam_0 = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.01_all4.bam",
         bam_0_sort = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.01_sort_all4.bam",
     shell:
         """
         samtools view -bs 42.5 {input.bam2} > {output.bam_half}
         samtools sort {output.bam_half} > {output.bam_half_sort}
         samtools view -bs 50.01 {output.bam_half_sort} > {output.bam_0}
         samtools sort {output.bam_0} > {output.bam_0_sort}
         samtools index {output.bam_0_sort}
         samtools index {output.bam_half_sort}
         """    
rule drop_out_0_with_new_KIR4:
     input:
         bam_half_sort = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.5_sort_all4.bam",
     output:
         bam_0 = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.0001_all4.bam",
         bam_0_sort = "{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_0.0001_sort_all4.bam",
     shell:
         """
         samtools view -bs 50.0001 {input.bam_half_sort} > {output.bam_0}
         samtools sort {output.bam_0} > {output.bam_0_sort}
         samtools index {output.bam_0_sort}
         """  
