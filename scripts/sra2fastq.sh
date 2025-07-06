#!/usr/bin/env bash
set -euo pipefail

# ---------------------- Default Parameters ----------------------
THREADS=8
OUTDIR="./"
TMPROOT=""
USER_SET_TMPDIR=false
INPUT=""
RENAME=false

usage() {
  echo "Usage: $0 -i <ACCESSION | FILE> [-t THREADS] [-o OUTDIR] [-d TMPDIR]"
  echo "  -i, --input         Accession ID (SRP/SRX/SRS/SRR) or file with one per line"
  echo "  -t, --threads       Threads to use (default: 8)"
  echo "  -o, --outdir        Output directory (default: ./)"
  echo "  -d, --tmpdir        Temporary directory (default: OUTDIR)"
  echo "  -r, --rename        Rename output files to SRX_<role>.fastq.gz (default: false)"
  echo "  --prefetch          Prefetch command (default: prefetch)"
  echo "  --fasterq-dump      Fasterq-dump command (default: fasterq-dump)"
  echo "  --pigz              Pigz command (default: pigz)"
  echo "  --help              Show this help message"
  exit 1
}

# ---------------------- Parse Arguments ----------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -d|--tmpdir) TMPROOT="$2"; USER_SET_TMPDIR=true; shift 2 ;;
    -r|--rename) RENAME=true; shift ;;
    --prefetch) PREFETCH="$2"; shift 2 ;;
    --fasterq-dump) FASTERQ_DUMP="$2"; shift 2 ;;
    --pigz) PIGZ="$2"; shift 2 ;;
    --help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

[[ -z "$INPUT" ]] && { echo "ERROR: --input is required"; usage; }

# ---------------------- Prepare Environment ----------------------
mkdir -p "$OUTDIR"
SUCCESS_LOG="$OUTDIR/success.log"
FAIL_LOG="$OUTDIR/fail.log"
if [[ ! -f "$SUCCESS_LOG" ]]; then
  : > "$SUCCESS_LOG"
fi
if [[ ! -f "$FAIL_LOG" ]]; then
  : > "$FAIL_LOG"
fi

PREFETCH=${PREFETCH:-prefetch}
FASTERQ_DUMP=${FASTERQ_DUMP:-fasterq-dump}
PIGZ=${PIGZ:-pigz}

log() { echo "[$(date +'%F %T')] $*" >&2; }

# ---------------------- Resolver Functions ----------------------
resolve_to_srx() {
  local acc="$1"
  if [[ "$acc" =~ ^SRP ]]; then
    pysradb metadata "$acc" --assay | awk 'NR>1 {print $3}' | sort -u
  elif [[ "$acc" =~ ^SRS ]]; then
    pysradb srs-to-srx "$acc" | awk 'NR>1 {print $2}'
  elif [[ "$acc" =~ ^SRX ]]; then
    echo "$acc"
  elif [[ "$acc" =~ ^SRR ]]; then
    echo "$acc"
  else
    log "Unrecognized accession: $acc"
    exit 1
  fi
}

resolve_to_srr() {
  local acc="$1"
  if [[ "$acc" =~ ^SRX ]]; then
    pysradb srx-to-srr "$acc" | awk 'NR>1 {print $2}'
  elif [[ "$acc" =~ ^SRR ]]; then
    echo "$acc"
  else
    log "ERROR: Cannot resolve to SRR: $acc"
    exit 1
  fi
}


# ---------------------- Core Processing Functions ----------------------
process_one_srr() {
  local run=$1
  local tmp

  if [[ "$USER_SET_TMPDIR" == true ]]; then
    tmp="${TMPROOT}/${run}"
  else
    tmp="${OUTDIR}/${run}"
  fi

  if ls "$OUTDIR/${run}"_*.fastq.gz &>/dev/null; then
    log "$run outputs exist — skipping"
    return 0
  fi

  mkdir -p "$tmp"
  log "[$run] Downloading with prefetch..."
  "$PREFETCH" --max-size u --progress --output-file "$tmp/${run}.sra" "$run"

  local sra="$tmp/${run}.sra"
  [[ ! -f "$sra" ]] && { log "[$run] SRA file not found"; return 1; }

  log "[$run] Running fasterq-dump..."
  "$FASTERQ_DUMP" "$sra" -O "$tmp" -t "$tmp" -e "$THREADS" -p -S --include-technical -f

  log "[$run] Compressing..."
  find "$tmp" -name "*.fastq" -exec "$PIGZ" -p "$THREADS" {} +

  mv "$tmp"/*.fastq.gz "$OUTDIR/" || true
  rm -rf "$tmp"
  log "[$run] DONE"
}

merge_by_srx() {
  local srx="$1"
  local runs=($(resolve_to_srr "$srx"))
  [[ "${#runs[@]}" -eq 0 ]] && return

  if [[ "${#runs[@]}" -eq 1 ]]; then
    local r="${runs[0]}"
    for f in "$OUTDIR/${r}"_*.fastq.gz; do
      [[ -f "$f" ]] || continue
      mv "$f" "${f/${r}/${srx}}"
    done
  else
    log "[$srx] Merging ${#runs[@]} runs"
    local suffixes=$(ls "$OUTDIR/${runs[0]}"_*.fastq.gz | sed -E 's/.*_([0-9]+)\.fastq\.gz/\1/' | sort -u)
    for s in $suffixes; do
      cat $(for r in "${runs[@]}"; do echo "$OUTDIR/${r}_${s}.fastq.gz"; done) > "$OUTDIR/${srx}_${s}.fastq.gz"
    done
    for r in "${runs[@]}"; do
      rm -f "$OUTDIR/${r}"_*.fastq.gz
    done
  fi
  log "[$srx] Merge complete"
}

summarize_fastq() {
  local srx="$1"
  local summary="$OUTDIR/${srx}.fastq_summary.csv"

  echo "filename,read_length,role" > "$summary"

  local tmpfile
  tmpfile=$(mktemp)

  for f in "$OUTDIR/${srx}"_*.fastq.gz; do
    [[ -s "$f" ]] || { log "[$srx] Warning: file $f is missing or empty"; continue; }
    read_length=$(zcat "$f" | head -n 4 | sed -n '2p' | awk '{print length($0)}')
    echo "$(basename "$f") $read_length" >> "$tmpfile"
  done

  log "[$srx] Summary collected, inferring roles..."

  awk '
    BEGIN { OFS="," }
    {
      file = $1
      len = $2
      match(file, /_([0-9]+)\.fastq/, m)
      id = m[1]
      files[id] = file
      lengths[id] = len
      order[n++] = id
    }
    END {
      asort(order, sorted)

      # detect index reads
      idx = 0
      for (i = 1; i <= n; i++) {
        r = sorted[i]
        if (lengths[r] >= 7 && lengths[r] <= 10) {
          if (idx == 0) { roles[r] = "I1"; idx++ }
          else if (idx == 1) { roles[r] = "I2"; idx++ }
        }
      }

      # assign R1/R2
      rest = 0
      for (i = 1; i <= n; i++) {
        r = sorted[i]
        if (!(r in roles)) {
          nonidx[rest++] = r
        }
      }

      if (rest == 2) {
        a = nonidx[0]; b = nonidx[1]
        if (lengths[a] >= 25 && lengths[a] <= 35) {
          roles[a] = "R1"; roles[b] = "R2"
        } else if (lengths[b] >= 25 && lengths[b] <= 35) {
          roles[b] = "R1"; roles[a] = "R2"
        } else {
          roles[a] = "R1"; roles[b] = "R2"
        }
      } else if (rest == 1) {
        roles[nonidx[0]] = "R1"
      }

      # output
      for (i = 1; i <= n; i++) {
        r = sorted[i]
        print files[r], lengths[r], (r in roles ? roles[r] : "Unknown")
      }
    }
  ' "$tmpfile" >> "$summary"

  rm -f "$tmpfile"
  log "[$srx] Summary + role annotation complete"
}

rename_by_role() {
  local srx="$1"
  local summary="$OUTDIR/${srx}.fastq_summary.csv"
  [[ ! -f "$summary" ]] && { log "[$srx] No summary to rename from"; return; }

  log "[$srx] Renaming by role via symlink..."

  awk -F',' -v srx="$srx" '
    NR > 1 && $3 != "Unknown" {
      orig = $1              # e.g. SRX9207567_1.fastq.gz
      role = $3              # e.g. R1
      link = srx "_" role ".fastq.gz"
      printf("ln -sf %s %s\n", orig, link)
    }
  ' "$summary" | (cd "$OUTDIR" && bash)

  log "[$srx] Rename complete"
}

process_accession() {
  local acc="$1"
  if [[ "$acc" =~ ^SRR ]]; then
    process_one_srr "$acc"
    summarize_fastq "$acc"
  else
    local srxs=($(resolve_to_srx "$acc"))
    for srx in "${srxs[@]}"; do
      local srrs=($(resolve_to_srr "$srx"))
      for srr in "${srrs[@]}"; do
        process_one_srr "$srr"
      done
      merge_by_srx "$srx"
      summarize_fastq "$srx"
    done
  fi
}

process_accession_with_status() {
  local acc="$1"

  if grep -qxF "$acc" "$SUCCESS_LOG"; then
    log "[$acc] Already completed."
  else
    if process_accession "$acc"; then
      echo "$acc" >> "$SUCCESS_LOG"
    else
      echo "$acc" >> "$FAIL_LOG"
    fi
  fi

  if [[ "$RENAME" == true ]]; then
    rename_by_role "$acc"
  fi
}


# ---------------------- Main Loop ----------------------
if [[ -f "$INPUT" ]]; then
  while IFS= read -r line; do
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    process_accession_with_status "$line"
  done < "$INPUT"
else
  process_accession_with_status "$INPUT"
fi

log "All accessions processed."
