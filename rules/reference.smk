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
        stdout = f"{OUTDIR}/logs/build_reference/{GENOME}_v{VERSION}_{TOOL}.stdout.log",
        stderr = f"{OUTDIR}/logs/build_reference/{GENOME}_v{VERSION}_{TOOL}.stderr.log"
    benchmark: f"{OUTDIR}/benchmarks/build_reference/{GENOME}_v{VERSION}_{TOOL}.benchmark.log"
    shell:
        """
        bash scripts/build_reference.sh \
            --genome {params.genome} \
            --version {params.version} \
            --tool {params.tool} \
            --outdir {output} \
            --filter \
            --threads {threads} \
            > {log.stdout} 2> {log.stderr}
        """
