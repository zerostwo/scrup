rule build_reference:
    output:
        directory(f"{OUTDIR}/reference/{GENOME}/v{VERSION}/{TOOL}")
    params:
        genome=GENOME,
        version=VERSION,
        tool=TOOL
    cache: True
    threads: REF_THREADS
    conda: "../envs/alignment.yaml"
    shell:
        """
        bash scripts/build_reference.sh \
            --genome {params.genome} \
            --version {params.version} \
            --tool {params.tool} \
            --outdir {output} \
            --filter \
            --threads {threads}
        """
