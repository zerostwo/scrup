rule align_starsolo:
    input:
        read1 = lambda wc: f"{OUTDIR}/fastqs/{wc.sample}_R1.fastq.gz",
        read2 = lambda wc: f"{OUTDIR}/fastqs/{wc.sample}_R2.fastq.gz",
        reference = f"{OUTDIR}/reference/{GENOME}/v{VERSION}/{TOOL}",
    output:
        summary = f"{OUTDIR}/alignments/{{sample}}/outs/GeneFull/Summary.csv"
    params:
        whitelist = "resources/whitelists",
        sample = lambda wc: wc.sample,
        outdir = lambda wc: f"{OUTDIR}/alignments",
        memory = ALIGN_MEMORY,
        create_bam = ALIGN_CREATE_BAM,
        keep_unmapped = ALIGN_KEEP_UNMAPPED,
    threads: ALIGN_THREADS
    conda: "../envs/alignment.yaml"
    log:
        stdout = f"{OUTDIR}/logs/align_starsolo/{{sample}}.stdout.log",
        stderr = f"{OUTDIR}/logs/align_starsolo/{{sample}}.stderr.log"
    benchmark: f"{OUTDIR}/benchmarks/align_starsolo/{{sample}}.benchmark.log"
    shell:
        """
        bash scripts/starsolo_10x.sh \
            --sample {params.sample} \
            --read1 {input.read1} \
            --read2 {input.read2} \
            --reference {input.reference} \
            --whitelists {params.whitelist} \
            --threads {threads} \
            --memory {params.memory} \
            --outdir {params.outdir} \
            {{"--bam" if params.create_bam else ""}} \
            {{"--unmapped" if params.keep_unmapped else ""}} \
            > {log.stdout} 2> {log.stderr}
        """
