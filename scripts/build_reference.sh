#!/usr/bin/env bash
set -euo pipefail

# ---------------------- Default Parameters ----------------------
GENOME=""
VERSION=""
SPECIES=""
TOOL="STAR"
THREADS=12
OUTDIR="$PWD"
FILTER=false

# ---------------------- Help Message ----------------------
usage() {
  echo -e "\033[1mBuild human or mouse genome reference\033[0m"
  echo ""
  echo "Usage:"
  echo "  $0 [-g GENOME -v VERSION] | [-s SPECIES] [-T TOOL_PATH] [-t THREADS] [-o OUTDIR] [--filter]"
  echo ""
  echo "Required:"
  echo "  -g, --genome     Reference genome: GRCh38 | GRCh37 | GRCm39 | GRCm38"
  echo "  -v, --version    GENCODE version: e.g., 48 (human), M37 (mouse)"
  echo "     OR"
  echo "  -s, --species    Shortcut: human for GRCh38+48, mouse for GRCm39+M37"
  echo ""
  echo "Options:"
  echo "  -T, --tool       Tool path: STAR (default), or path to cellranger"
  echo "  -t, --threads    Threads to use (default: 12)"
  echo "  -o, --outdir     Output directory (default: current directory)"
  echo "  -f, --filter     Enable GTF filtering for scRNA-seq (default: off)"
  echo ""
  echo "Examples:"
  echo "  $0 --species human --filter"
  echo "  $0 --genome GRCh38 --version 48 --tool /path/to/cellranger -o ./ref_build"
  exit 1
}

# ---------------------- Logging ----------------------
log() { echo "[$(date +'%F %T')] $*" >&2; }

# ---------------------- Parse Arguments ----------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--genome) GENOME="$2"; shift 2 ;;
    -v|--version) VERSION="$2"; shift 2 ;;
    -s|--species) SPECIES="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -T|--tool) TOOL="$2"; shift 2 ;;
    -f|--filter) FILTER=true; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# ---------------------- Species Shortcut ----------------------
if [[ -n "$SPECIES" ]]; then
  case "$SPECIES" in
    human) GENOME="GRCh38"; VERSION="48" ;;
    mouse) GENOME="GRCm39"; VERSION="M37" ;;
    *) echo "[ERROR] Invalid species: $SPECIES (choose 'human' or 'mouse')"; exit 1 ;;
  esac
fi

# ---------------------- Validation ----------------------
if [[ -z "$GENOME" || -z "$VERSION" ]]; then
  echo "[ERROR] Must specify --genome and --version or use --species"
  exit 1
fi

# ---------------------- Paths ----------------------
SOURCE="$(dirname "$OUTDIR")/sources"
BUILD="$OUTDIR"
mkdir -p "$SOURCE" "$BUILD"

# ---------------------- Reference URLs ----------------------
if [[ "$GENOME" =~ ^GRCh ]]; then
  FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/${GENOME}.primary_assembly.genome.fa.gz"
  GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/gencode.v${VERSION}.primary_assembly.annotation.gtf.gz"
  FASTA_IN="${SOURCE}/${GENOME}.primary_assembly.genome.fa"
  GTF_IN="${SOURCE}/gencode.v${VERSION}.primary_assembly.annotation.gtf"
elif [[ "$GENOME" =~ ^GRCm ]]; then
  FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${VERSION}/${GENOME}.primary_assembly.genome.fa.gz"
  GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${VERSION}/gencode.v${VERSION}.primary_assembly.annotation.gtf.gz"
  FASTA_IN="${SOURCE}/${GENOME}.primary_assembly.genome.fa"
  GTF_IN="${SOURCE}/gencode.v${VERSION}.primary_assembly.annotation.gtf"
else
  log "Unsupported genome: $GENOME"
  exit 1
fi

# ---------------------- Download ----------------------
[[ -f "$FASTA_IN" ]] || { log "Downloading FASTA..."; curl -sSL "$FASTA_URL" | gunzip -c > "$FASTA_IN"; }
[[ -f "$GTF_IN" ]]   || { log "Downloading GTF..."; curl -sSL "$GTF_URL" | gunzip -c > "$GTF_IN"; }

# Use samtools to index the fasta file, if samtools is installed
if command -v samtools &> /dev/null; then
  log "Indexing FASTA file with samtools..."
  samtools faidx "$FASTA_IN"
fi

# ---------------------- GTF Filtering (Optional) ----------------------
TOOL_NAME=$(basename "$TOOL")
GTF_FILTER="$GTF_IN"

if [[ "$FILTER" == true ]]; then
  log "Filtering GTF for single-cell use..."

  BIOTYPES="(protein_coding|protein_coding_LoF|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene)"
  GENE_TAG="gene_type \"${BIOTYPES}\""
  TX_TAG="transcript_type \"${BIOTYPES}\""
  READTHROUGH="tag \"readthrough_transcript\""

  ALLOWLIST="${SOURCE}/gene_allowlist.txt"
  cat "$GTF_IN" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_TAG" \
    | grep -E "$TX_TAG" \
    | grep -Ev "$READTHROUGH" \
    | sed -n 's/.*gene_id "\([^"]*\)".*/\1/p' \
    | sort -u > "$ALLOWLIST"

  FILTER_CHRY=false
  if [[ "$GENOME" == "GRCh38" ]]; then
    VERSION_MAIN=$(echo "$VERSION" | sed 's/[^0-9]*//g')
    if [[ "$VERSION_MAIN" -ge 44 ]]; then
      log "chrY PAR filtering will be applied"
      FILTER_CHRY=true
    fi
  fi

  GTF_FILTER="${SOURCE}/filtered.gtf"
  grep -E "^#" "$GTF_IN" > "$GTF_FILTER"
  grep -Ff "$ALLOWLIST" "$GTF_IN" | {
    if [[ "$FILTER_CHRY" == true ]]; then
      awk -F "\t" '$1 != "chrY" || ($4 >= 2752083 && $4 < 56887903 && !/ENSG00000290840/)'
    else
      cat
    fi
  } >> "$GTF_FILTER"
fi

# ---------------------- Build Reference ----------------------
if [[ "$TOOL_NAME" == "STAR" ]]; then
  log "Building STAR index..."
  "$TOOL" \
    --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --outTmpDir "$BUILD/_STARtmp" \
    --genomeDir "$BUILD" \
    --genomeFastaFiles "$FASTA_IN" \
    --sjdbGTFfile "$GTF_FILTER"
  log "STAR reference built at $BUILD"

elif [[ "$TOOL_NAME" == "cellranger" ]]; then
  log "Building Cell Ranger reference..."
  "$TOOL" mkref \
    --genome="$GENOME" \
    --ref-version="$VERSION" \
    --fasta="$FASTA_IN" \
    --genes="$GTF_FILTER" \
    --nthreads="$THREADS"
  log "Cell Ranger reference built in $GENOME"

else
  log "Unsupported tool: $TOOL"
  exit 1
fi
