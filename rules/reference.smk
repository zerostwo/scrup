rule build_reference:
    output:
        directory(f"{OUTDIR}/reference/{GENOME}/v{VERSION}/{TOOL}")
    params:
        genome = GENOME,
        version = VERSION,
        tool = TOOL
    cache: True
    threads: REF_THREADS
    conda: "../envs/alignment.yaml"
    log:
        out = f"{OUTDIR}/logs/build_reference/{GENOME}_v{VERSION}_{TOOL}.out.log",
        err = f"{OUTDIR}/logs/build_reference/{GENOME}_v{VERSION}_{TOOL}.err.log"
    shell:
        """
        bash scripts/build_reference.sh \
            --genome {params.genome} \
            --version {params.version} \
            --tool {params.tool} \
            --outdir {output} \
            --filter \
            --threads {threads} \
            > {log.out} 2> {log.err}
        """
