rule sra2fastq:
    output:
        read1=f"{OUTDIR}/fastqs/{{sample}}_R1.fastq.gz",
        read2=f"{OUTDIR}/fastqs/{{sample}}_R2.fastq.gz",
        summary=f"{OUTDIR}/fastqs/{{sample}}.fastq_summary.csv"
    params:
        outdir=f"{OUTDIR}/fastqs"
    threads: DOWNLOAD_THREADS
    conda: "../envs/download.yaml"
    resources:
        download_slots=1
    log:
        stdout=f"{OUTDIR}/logs/sra2fastq/{{sample}}.log",
        stderr=f"{OUTDIR}/logs/sra2fastq/{{sample}}.err"
    shell:
        """
        bash scripts/sra2fastq.sh \
            --input {wildcards.sample} \
            --threads {threads} \
            --outdir {params.outdir} \
            --tmpdir {DOWNLOAD_TMPDIR} \
            --rename \
            > {log.stdout} 2> {log.stderr}
        """
