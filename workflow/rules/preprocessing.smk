rule extract:
    input:
        expand("{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_tag_sort_all4.bam", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),
        expand("{DATA_DIR}/{sample}/paired_new_kir_sort_all4.bam.bai", DATA_DIR= config["Samples"]["samples_dir"],sample=config["Samples"]["sample_fastqs"]),


rule extract_reads:
    input:
        reads=lambda wildcards: config["Samples"]["sample_fastqs"][wildcards.sample],
        ref=config["Reference"]["fasta"],
        bed=config["Reference"]["KIR_regions_bed"]
    output:
        bam="{DATA_DIR}/{sample}/mapped.bam"
    threads: min(config["threads"], 30)
    params:
        tmpdir="{DATA_DIR}/{sample}/tmp",  
        background_region=config["background_region"] 
    run:
        import os

        reads = input.reads
        if reads[0].endswith(".cram"):
            shell("""
                mkdir -p {params.tmpdir}
                # 1. Reads overlapping KIR regions
                samtools index {reads}
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
                # 5. Extract background region
                samtools view -@ {threads} -b -T {input.ref} {reads} {params.background_region} \
                > {output.bam}.background.bam

                # 6. Merge + sort
                samtools merge -@ {threads} -f {output.bam}.merged.bam \
                {output.bam}.kir.named.bam  {output.bam}.background.bam {output.bam}.unmapped.bam

                samtools sort -@ {threads} -T {params.tmpdir} -o {output.bam} {output.bam}.merged.bam

                rm -f {output.bam}.kir.tmp* {output.bam}.names.tmp \
                    {output.bam}.kir.named.bam {output.bam}.unmapped.bam \
                    {output.bam}.merged.bam {output.bam}.background.bam 
            """)
        elif reads[0].endswith(".fastq") or reads[0].endswith(".fq") or reads[0].endswith(".fastq.gz"):
            bwa = config["bwa_path"]
            shell("""
                set -euo pipefail
                mkdir -p {params.tmpdir}
                # 0. Align reads
                {bwa} mem -t {threads} {input.ref} {reads} > {output.bam}.align.sam
                
                samtools view -@ {threads} -b {output.bam}.align.sam  | samtools sort -@ {threads} -o {output.bam}.align.tmp
                samtools index -@ {threads} {output.bam}.align.tmp

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

                # 5. Extract background region
                samtools view -@ {threads} -hb {output.bam}.align.tmp {params.background_region} \
                > {output.bam}.background.bam

                # 6. Merge + sort
                samtools merge -@ {threads} -f {output.bam}.merged.bam \
                {output.bam}.kir.named.bam {output.bam}.background.bam {output.bam}.unmapped.bam

                samtools sort -@ {threads} -T {params.tmpdir} -o {output.bam} {output.bam}.merged.bam
                samtools index {output.bam}

                # 7. Cleanup
                rm -f {output.bam}.align.sam \
                    {output.bam}.align.tmp \
                    {output.bam}.kir.tmp* \
                    {output.bam}.names.raw.tmp \
                    {output.bam}.names.tmp \
                    {output.bam}.kir.named.bam \
                    {output.bam}.unmapped.bam \
                    {output.bam}.merged.bam \
                    {output.bam}.sort.kir.tmp \
                    {output.bam}.sort.kir.tmp \
                    {output.bam}.background.bam
            """)
rule remove_dup_in_bam_files: 
    input:
       bam="{DATA_DIR}/{sample}/mapped.bam"
    output:
       bam="{DATA_DIR}/{sample}/dedup_mapped_noRG.bam",
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
     params:
         background_region=config["background_region"],
         background_region_chr=config["background_region"].split(":")[0]
     shell:
         """
         samtools index {input.bam}
         samtools view -b {input.bam} {params.background_region} > {output.bam}
         samtools index {output.bam}
         samtools view -b -e 'rname!="{params.background_region_chr}"' {input.bam} > {output.kir}
         samtools index {output.kir}
         """          
rule bam_to_fastq:
    input:
        bam="{DATA_DIR}/{sample}/mapped_filt_noChr17.bam"
    output:
        read1=temp("{DATA_DIR}/{sample}/mapped_filt.read1.fastq"),
        read2=temp("{DATA_DIR}/{sample}/mapped_filt.read2.fastq"),
        sing=temp("{DATA_DIR}/{sample}/singletons.fastq")
    threads: min(config["threads"], 20)
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
        fastq1="{DATA_DIR}/{sample}/mapped_filt.read1_mod.fastq",
        fastq2="{DATA_DIR}/{sample}/mapped_filt.read2_mod.fastq",  
        reference=config["Reference"]["fasta_2"],
    output:
        sam=temp("{DATA_DIR}/{sample}/mapped_the_second_time.sam")
    threads: min(config["threads"], 72)
    params: 
        bwa=config["bwa_path"]
    shell:
        """
        {params.bwa} mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} > {output.sam}
        """  
rule sam_to_bam2:
    input:
         "{data_dir}/{sample}/mapped_the_second_time.sam"
    output:
        temp("{data_dir}/{sample}/mapped_the_second_time.bam")
    threads: min(config["threads"], 72)
    shell:
        "samtools view -Sb {input} > {output}"

rule filter_unmapped_new_kir2_keep_utrs:
    input:
        bam = "{DATA_DIR}/{sample}/mapped_the_second_time.bam",
        bed = config["Reference"]["KIR_alleles_bed"]
    output:
        bam = "{DATA_DIR}/{sample}/mapped_filt_new_kir2.bam",
        roi = temp("{DATA_DIR}/{sample}/roi_tmp2.bam"),
    threads: min(config["threads"], 72)
    shell:
        """
        samtools view -b -L {input.bed} {input.bam} > {output.roi}
        samtools sort -@ {threads} -o {output.bam} {output.roi}
        samtools index -@ {threads} {output.bam}
        """ 
rule extract_mapped_reads_with_seqkit_new_kir2:
    input:
        bam="{DATA_DIR}/{sample}/mapped_filt_new_kir2.bam"
    output:
        read1=temp("{DATA_DIR}/{sample}/mapped_filt.read1_new_kir2.fastq"),
        read2=temp("{DATA_DIR}/{sample}/mapped_filt.read2_new_kir2.fastq"),
        sing=temp("{DATA_DIR}/{sample}/Sing_new_kir2.fastq"),
    threads: min(config["threads"], 20)
    shell:
        """
        picard SamToFastq I={input.bam} F={output.read1}  F2={output.read2} FU={output.sing}  VALIDATION_STRINGENCY=SILENT
        """                                  
rule proper_mapping_with_new_KIR_4:
    input:
        read1="{DATA_DIR}/{sample}/mapped_filt.read1_new_kir2.fastq",
        read2="{DATA_DIR}/{sample}/mapped_filt.read2_new_kir2.fastq",
        mmi= config["Reference"]["KIR_alleles"],
    output:
        bam1="{DATA_DIR}/{sample}/remapped1_new_kir4.bam",
        bam2="{DATA_DIR}/{sample}/remapped2_new_kir4.bam"
    threads: min(config["threads"], 72)
    shell:
        """
        minimap2 -ax sr --secondary-seq --MD --eqx --secondary=yes -N 2000 -t {threads} {input.mmi} {input.read1} | samtools view -b | samtools view -hF 4 | samtools sort > {output.bam1}
        samtools index {output.bam1}
        minimap2 -ax sr  --secondary-seq --MD --eqx --secondary=yes -N 2000 -t {threads} {input.mmi} {input.read2} | samtools view -b | samtools view -hF 4 | samtools sort > {output.bam2}
        samtools index {output.bam2}
        """       
rule tag_chr17_bam:
    """ZP-tag the background BAM and compute insert-size stats — runs once per sample."""
    input:
        bam="{DATA_DIR}/{sample}/mapped_filt_Chr17q25.bam"
    output:
        tagged_bam="{DATA_DIR}/{sample}/mapped_filt_Chr17q25_with_tag_new_kir_all4.bam",
        stats="{DATA_DIR}/{sample}/insert_stats.json"
    params:
        chr17=config["background_region"].split(":")[0]
    threads: min(config["threads"], 30)
    script:
        "../scripts/tag_chr17.py"

checkpoint split_kir_bams:
    """
    Split bam1 and bam2 into size-based contig chunks (target: ≤ max_chunk_gb each).
    Uses samtools view for extraction. Snakemake re-evaluates downstream rules after
    this checkpoint completes to discover the actual number of chunks produced.
    """
    input:
        bam1="{DATA_DIR}/{sample}/remapped1_new_kir4.bam",
        bam2="{DATA_DIR}/{sample}/remapped2_new_kir4.bam"
    output:
        directory("{DATA_DIR}/{sample}/kir_chunks")
    params:
        max_chunk_gb=config.get("max_chunk_gb", 0.5)
    run:
        import math, os
        os.makedirs(output[0], exist_ok=True)

        # Compute number of chunks so each is ≤ max_chunk_gb
        max_bytes  = max(os.path.getsize(input.bam1), os.path.getsize(input.bam2))
        n_chunks   = max(1, math.ceil(max_bytes / (params.max_chunk_gb * 1024 ** 3)))

        contigs_f  = os.path.join(output[0], "contigs.txt")
        chunks_pfx = os.path.join(output[0], "chunk_")

        # Extract KIR/IAG contig names from bam1 header using samtools
        shell(
            "samtools view -H {input.bam1} | grep '^@SQ' | "
            "awk '{{print $2}}' | sed 's/SN://' | grep -E '^(KIR|IAG)' > " + contigs_f
        )

        # Split contig list into n_chunks line-balanced groups
        shell(f"split -n l/{n_chunks} {contigs_f} {chunks_pfx}")

        # Build one shell command that launches all chunk extractions in parallel
        chunk_files = sorted(
            f for f in os.listdir(output[0]) if f.startswith("chunk_")
        )
        cmds = []
        for i, chunk_file in enumerate(chunk_files):
            chunk_path = os.path.join(output[0], chunk_file)
            bam1_out   = os.path.join(output[0], f"bam1_chunk{i}.bam")
            bam2_out   = os.path.join(output[0], f"bam2_chunk{i}.bam")
            # Each pair of samtools view calls is backgrounded with &
            # {input.bam1} / {input.bam2} are substituted by Snakemake's shell()
            cmds.append(
                "samtools view -b {input.bam1} "
                "$(cat " + chunk_path + " | tr '\\n' ' ') > " + bam1_out + " &"
            )
            cmds.append(
                "samtools view -b {input.bam2} "
                "$(cat " + chunk_path + " | tr '\\n' ' ') > " + bam2_out + " &"
            )
        # wait blocks until every background job finishes
        cmds.append("wait")
        shell("\n".join(cmds))

        # Remove temporary contig list and chunk files
        shell("rm -f " + contigs_f + " " + chunks_pfx + "*")

def _kir_chunk_ids(wildcards):
    """Return chunk indices discovered after the split_kir_bams checkpoint runs."""
    chunks_dir = checkpoints.split_kir_bams.get(**wildcards).output[0]
    ids, = glob_wildcards(os.path.join(chunks_dir, "bam1_chunk{chunk}.bam"))
    return ids

rule pair_kir_chunk:
    """Pair reads for one contig chunk — all chunks run in parallel."""
    input:
        bam1="{DATA_DIR}/{sample}/kir_chunks/bam1_chunk{chunk}.bam",
        bam2="{DATA_DIR}/{sample}/kir_chunks/bam2_chunk{chunk}.bam",
        stats="{DATA_DIR}/{sample}/insert_stats.json",
        fasta=config["Reference"]["KIR_alleles"]
    output:
        bam1=temp("{DATA_DIR}/{sample}/kir_chunks/paired_chunk{chunk}.bam")
    wildcard_constraints:
        chunk="\d+"
    threads:  min(config["threads"], 30)
    script:
        "../scripts/pairing.py"

rule merge_kir_pairs:
    """Merge all per-chunk paired BAMs into the final output."""
    input:
        bams=lambda wc: expand(
            "{DATA_DIR}/{sample}/kir_chunks/paired_chunk{chunk}.bam",
            DATA_DIR=wc.DATA_DIR,
            sample=wc.sample,
            chunk=_kir_chunk_ids(wc)
        )
    output:
        "{DATA_DIR}/{sample}/paired_new_kir_all4.bam"
    threads: min(config["threads"], 30)
    shell:
        "samtools merge -@ {threads} -f {output} {input.bams}"
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
