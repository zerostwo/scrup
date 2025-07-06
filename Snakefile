configfile: "config.yaml"

# Define global variables from config
SAMPLES_FILE = config["samples"]
OUTDIR = config["outdir"]

# Download
DOWNLOAD_THREADS = config["download"]["threads"]
DOWNLOAD_TMPDIR = config["download"]["tmpdir"]

# Reference genome build
GENOME = config["reference"]["genome"]
VERSION = config["reference"]["version"]
TOOL = config["reference"]["tool"]
REF_THREADS = config["reference"]["threads"]

# Alignment
ALIGN_MEMORY = config["alignment"]["memory"]
ALIGN_THREADS = config["alignment"]["threads"]

# Load sample list from file
with open(SAMPLES_FILE) as f:
    SAMPLE_IDS = [line.strip() for line in f if line.strip()]

include: "rules/reference.smk"
include: "rules/download.smk"
include: "rules/alignment.smk"

# Rule: all
rule all:
    input:
        # Reference genome build done
        f"{OUTDIR}/reference/{GENOME}/v{VERSION}/{TOOL}",
        # Fastq for each sample
        expand(f"{OUTDIR}/fastqs/{{sample}}.fastq_summary.csv", sample=SAMPLE_IDS),
        expand(f"{OUTDIR}/fastqs/{{sample}}_R1.fastq.gz", sample=SAMPLE_IDS),
        expand(f"{OUTDIR}/fastqs/{{sample}}_R2.fastq.gz", sample=SAMPLE_IDS),
        # Alignment for each sample
        expand(f"{OUTDIR}/alignments/{{sample}}/Aligned.sortedByCoord.out.bam", sample=SAMPLE_IDS)
