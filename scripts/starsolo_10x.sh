#!/usr/bin/env bash
set -euo pipefail

# ---------------------- Default Parameters ----------------------
SAMPLE=""
READ1=""
READ2=""
REFERENCE=""
WHITELISTS=""
THREADS=8
MEMORY="32G"
OUTDIR="./alignments"

# ---------------------- Help Message ----------------------
usage() {
  echo -e "\033[1mUsage: $0 [options]\033[0m"
  echo ""
  echo "Options:"
  echo "  -s, --sample       Sample tag/ID (required)"
  echo "  -1, --read1        Read 1 FASTQ file (required)"
  echo "  -2, --read2        Read 2 FASTQ file (required)"
  echo "  -r, --reference    STAR genome reference directory (required)"
  echo "  -l, --whitelists   Directory containing 10x whitelist files (required)"
  echo "  -t, --threads      Number of threads (default: 12)"
  echo "  -m, --memory       Max RAM, e.g., 64G (default: 64G)"
  echo "  -o, --outdir       Output directory (default: ./alignments)"
  echo "  -h, --help         Show this help message"
  exit 1
}

# ---------------------- Parse Arguments ----------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--sample) SAMPLE="$2"; shift 2 ;;
    -1|--read1) READ1="$2"; shift 2 ;;
    -2|--read2) READ2="$2"; shift 2 ;;
    -r|--reference) REFERENCE="$2"; shift 2 ;;
    -l|--whitelists) WHITELISTS="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -m|--memory) MEMORY="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Convert MEMORY to bytes (support G/M suffix)
MEMORY_IN_BYTES=""
if [[ "$MEMORY" =~ ^([0-9]+)[Gg]$ ]]; then
  MEMORY_IN_BYTES=$(( ${BASH_REMATCH[1]} * 1024 * 1024 * 1024 ))
elif [[ "$MEMORY" =~ ^([0-9]+)[Mm]$ ]]; then
  MEMORY_IN_BYTES=$(( ${BASH_REMATCH[1]} * 1024 * 1024 ))
else
  echo "Invalid memory format: $MEMORY (should be like 32G or 8000M)"
  exit 1
fi

# ---------------------- Check Required ----------------------
if [[ -z "$SAMPLE" || -z "$READ1" || -z "$READ2" || -z "$REFERENCE" || -z "$WHITELISTS" ]]; then
  echo "Error: Missing required arguments"
  usage
fi

# ---------------------- Load STAR + Tools ----------------------
STAR="STAR"
SEQTK="seqtk"
SAMTOOLS="samtools"
PIGZ="pigz"

# ---------------------- Output Setup ----------------------
FQDIR="$(dirname $READ1)"
ALIGNDIR="${OUTDIR}/${SAMPLE}"
mkdir -p "$ALIGNDIR"

# ---------------------- Main Script Logic ----------------------
log() { echo "[$(date +'%F %T')] $*" >&2; }

# Choose one of the two otions, depending on whether you need a BAM file 
BAM="--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM $MEMORY_IN_BYTES --outSAMunmapped Within --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"

# Define some key variables, in order to evaluate reads for being:
# 1) gzipped/bzipped/un-archived; 
# 2) having barcodes from the whitelist, and which; 
# 3) having consistent length; 
# 4) being single- or paired-end. 
GZIP="--readFilesCommand zcat"
ZCMD="zcat"
BC=""
NBC1=""
NBC2=""
NBC3=""
NBCA=""
R1LEN=""
R2LEN=""
R1DIS=""

# We need a small and random selection of reads. the solution below is a result 
# of much trial and error. In the end, we select 200k reads that represent all 
# of the files present in the FASTQ dir for this sample. Have to use numbers 
# because bamtofastq likes to make files with identical names in different folders..
log "Extracting reads for testing..."
COUNT=0
for i in `echo $READ1 | tr ',' ' '`
do
  $ZCMD $i | head -4000000 > ${ALIGNDIR}/$COUNT.R1_head &
  COUNT=$((COUNT+1))
done
wait 

COUNT=0
for i in `echo $READ2 | tr ',' ' ' `                                                
do 
  $ZCMD $i | head -4000000 > ${ALIGNDIR}/$COUNT.R2_head &
  COUNT=$((COUNT+1))
done
wait

# Same random seed makes sure you select same reads from R1 and R2
cat ${ALIGNDIR}/*.R1_head | $SEQTK sample -s100 - 200000 > ${ALIGNDIR}/test.R1.fastq &
cat ${ALIGNDIR}/*.R2_head | $SEQTK sample -s100 - 200000 > ${ALIGNDIR}/test.R2.fastq &
wait 
rm ${ALIGNDIR}/*.R1_head ${ALIGNDIR}/*.R2_head

# Elucidate the right barcode whitelist to use. Grepping out N saves us some trouble. 
# Note the special list for multiome experiments (737K-arc-v1.txt):
# 50k (out of 200,000) is a modified empirical number,
# matching only first 14-16 nt makes this more specific
log "Evaluating whitelists..."
NBC1=`cat ${ALIGNDIR}/test.R1.fastq | awk 'NR%4==2' | cut -c-14 | grep -F -f $WHITELISTS/737K-april-2014_rc.txt | wc -l`
NBC2=`cat ${ALIGNDIR}/test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f $WHITELISTS/737K-august-2016.txt | wc -l`
NBC3=`cat ${ALIGNDIR}/test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f $WHITELISTS/3M-february-2018_TRU.txt | wc -l`
NBC4=`cat ${ALIGNDIR}/test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f $WHITELISTS/3M-3pgex-may-2023_TRU.txt | wc -l`
NBC5=`cat ${ALIGNDIR}/test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f $WHITELISTS/3M-5pgex-jan-2023.txt | wc -l`
NBCA=`cat ${ALIGNDIR}/test.R1.fastq | awk 'NR%4==2' | cut -c-16 | grep -F -f $WHITELISTS/737K-arc-v1.txt | wc -l`
R1LEN=`cat ${ALIGNDIR}/test.R1.fastq | awk 'NR%4==2' | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}'`
R2LEN=`cat ${ALIGNDIR}/test.R2.fastq | awk 'NR%4==2' | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}'`
R1DIS=`cat ${ALIGNDIR}/test.R1.fastq | awk 'NR%4==2' | awk '{print length($0)}' | sort | uniq -c | wc -l`

if (( $NBC3 > 50000 )) 
then 
  BC=$WHITELISTS/3M-february-2018_TRU.txt
elif (( $NBC2 > 50000 ))
then
  BC=$WHITELISTS/737K-august-2016.txt
elif (( $NBCA > 50000 ))
then
  BC=$WHITELISTS/737K-arc-v1.txt
elif (( $NBC1 > 50000 )) 
then
  BC=$WHITELISTS/737K-april-2014_rc.txt
elif (( $NBC4 > 50000 )) 
then
  BC=$WHITELISTS/3M-3pgex-may-2023_TRU.txt
elif (( $NBC5 > 50000 )) 
then
  BC=$WHITELISTS/3M-5pgex-jan-2023.txt
else 
  >&2 echo "ERROR: No whitelist has matched a random selection of 200,000 barcodes! Match counts: $NBC1 (v1), $NBC2 (v2), $NBC3 (v3), $NBC4 (v4-3p), $NBC5 (v4-5p), $NBCA (multiome)."
  exit 1
fi 

# Check read lengths, fail if something funky is going on: 
PAIRED=False
UMILEN=""
CBLEN=""
if (( $R1DIS > 1 && $R1LEN <= 30 ))
then 
  >&2 echo "ERROR: Read 1 (barcode) has varying length; possibly someone thought it's a good idea to quality-trim it. Please check the fastq files."
  exit 1
elif (( $R1LEN < 24 )) 
then
  >&2 echo "ERROR: Read 1 (barcode) is less than 24 bp in length. Please check the fastq files."
  exit 1
elif (( $R2LEN < 40 )) 
then
  >&2 echo "ERROR: Read 2 (biological read) is less than 40 bp in length. Please check the fastq files."
  exit 1
fi

# Assign the necessary variables for barcode/UMI length/paired-end processing. 
# Script was changed to not rely on read length for the UMIs because of the epic Hassan case
# (v2 16bp barcodes + 10bp UMIs were sequenced to 28bp, effectively removing the effects of the UMIs)
if (( $R1LEN > 50 )) 
then
  PAIRED=True
fi

if [[ $BC == "$WHITELISTS/3M-february-2018_TRU.txt" || $BC == "$WHITELISTS/737K-arc-v1.txt" || $BC == "$WHITELISTS/3M-3pgex-may-2023_TRU.txt" || $BC == "$WHITELISTS/3M-5pgex-jan-2023.txt" ]] 
then 
  CBLEN=16
  UMILEN=12
elif [[ $BC == "$WHITELISTS/737K-august-2016.txt" ]] 
then
  CBLEN=16
  UMILEN=10
elif [[ $BC == "$WHITELISTS/737K-april-2014_rc.txt" ]] 
then
  CBLEN=14
  UMILEN=10
fi

# Yet another failsafe! Some geniuses managed to sequence v3 10x with a 26bp R1,
# which also causes STARsolo grief. This fixes it.
if (( $CBLEN + $UMILEN > $R1LEN ))
then
  NEWUMI=$((R1LEN-CBLEN))
  BCUMI=$((UMILEN+CBLEN))
  >&2 echo "WARNING: Read 1 length ($R1LEN) is less than the sum of appropriate barcode and UMI ($BCUMI). Changing UMI setting from $UMILEN to $NEWUMI!"
  UMILEN=$NEWUMI
elif (( $CBLEN + $UMILEN < $R1LEN ))
then
  BCUMI=$((UMILEN+CBLEN))
  >&2 echo "WARNING: Read 1 length ($R1LEN) is more than the sum of appropriate barcode and UMI ($BCUMI)."
fi

# It's hard to come up with a universal rule to correctly infer strand-specificity of the experiment. This is the best I could come up with: 
# 1) check if fraction of test reads (200k random ones) maps to GeneFull forward strand with higher than 50% probability; 
# 2) if not, run the same quantification with "--soloStand Reverse" and calculate the same stat; 
# 3) output a warning, and choose the strand with higher %;
# 4) if both percentages are below 10, 
log "Inferring strand-specificity..."
STRAND=Forward

mkdir -p ${ALIGNDIR}/test_forward/
$STAR --runThreadN $THREADS --genomeDir $REFERENCE \
  --readFilesIn ${ALIGNDIR}/test.R2.fastq ${ALIGNDIR}/test.R1.fastq \
  --runDirPerm All_RWX --outSAMtype None \
  --soloType CB_UMI_Simple --soloCBwhitelist $BC \
  --soloBarcodeReadLength 0 \
  --soloCBlen $CBLEN \
  --soloUMIstart $((CBLEN+1)) \
  --soloUMIlen $UMILEN --soloStrand Forward \
  --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloFeatures Gene GeneFull \
  --outTmpDir ${ALIGNDIR}/test_forward/_STARtmp \
  --outFileNamePrefix ${ALIGNDIR}/ \
  --soloOutFileNames test_forward/ features.tsv barcodes.tsv matrix.mtx &> /dev/null 

mkdir -p ${ALIGNDIR}/test_reverse/
$STAR --runThreadN $THREADS --genomeDir $REFERENCE \
  --readFilesIn ${ALIGNDIR}/test.R2.fastq ${ALIGNDIR}/test.R1.fastq \
  --runDirPerm All_RWX --outSAMtype None \
  --soloType CB_UMI_Simple --soloCBwhitelist $BC \
  --soloBarcodeReadLength 0 \
  --soloCBlen $CBLEN \
  --soloUMIstart $((CBLEN+1)) \
  --soloUMIlen $UMILEN --soloStrand Reverse \
  --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloFeatures Gene GeneFull \
  --outTmpDir ${ALIGNDIR}/test_reverse/_STARtmp \
  --outFileNamePrefix ${ALIGNDIR}/ \
  --soloOutFileNames test_reverse/ features.tsv barcodes.tsv matrix.mtx &> /dev/null

PCTFWD=`grep "Reads Mapped to GeneFull: Unique GeneFull" ${ALIGNDIR}/test_forward/GeneFull/Summary.csv | awk -F "," '{printf "%d\n",$2*100+0.5}'`
PCTREV=`grep "Reads Mapped to GeneFull: Unique GeneFull" ${ALIGNDIR}/test_reverse/GeneFull/Summary.csv | awk -F "," '{printf "%d\n",$2*100+0.5}'`

if (( $PCTREV > $PCTFWD )) 
then
  STRAND=Reverse
fi

if (( $PCTREV < 50 && $PCTFWD < 50)) 
then
  >&2 echo "WARNING: Low percentage of reads mapping to GeneFull: forward = $PCTFWD , reverse = $PCTREV"
fi 

# Finally, if paired-end experiment turned out to be 3' (yes, they do exist!), process it as single-end: 
if [[ $STRAND == "Forward" && $PAIRED == "True" ]]
then
  PAIRED=False
fi

# Write a file in the sample dir too, these metrics are not crucial but useful 
log "Done setting up the STARsolo run; here are final processing options:"
echo "============================================================================="
echo "Sample: $SAMPLE" | tee ${ALIGNDIR}/strand.txt
echo "Paired-end mode: $PAIRED" | tee -a ${ALIGNDIR}/strand.txt
echo "Strand (Forward = 3', Reverse = 5'): $STRAND, %reads mapped to GeneFull: forward = $PCTFWD , reverse = $PCTREV" | tee -a ${ALIGNDIR}/strand.txt
echo "CB whitelist: $BC, matches out of 200,000: $NBC3 (v3), $NBC2 (v2), $NBC1 (v1), $NBCA (multiome) " | tee -a ${ALIGNDIR}/strand.txt
echo "CB length: $CBLEN" | tee -a ${ALIGNDIR}/strand.txt
echo "UMI length: $UMILEN" | tee -a ${ALIGNDIR}/strand.txt
echo "GZIP: $GZIP" | tee -a ${ALIGNDIR}/strand.txt
echo "-----------------------------------------------------------------------------" | tee -a ${ALIGNDIR}/strand.txt
echo "Read 1 files: $READ1" | tee -a ${ALIGNDIR}/strand.txt
echo "-----------------------------------------------------------------------------" | tee -a ${ALIGNDIR}/strand.txt
echo "Read 2 files: $READ2" | tee -a ${ALIGNDIR}/strand.txt
echo "-----------------------------------------------------------------------------" | tee -a ${ALIGNDIR}/strand.txt

echo "sample,paired,strand,percent_forward,percent_reverse,cb_whitelist,cb_length,umi_length,gzip,read1_files,read2_files" > ${ALIGNDIR}/run_config.csv
echo "\"$SAMPLE\",\"$PAIRED\",\"$STRAND\",\"$PCTFWD\",\"$PCTREV\",\"$BC\",\"$CBLEN\",\"$UMILEN\",\"$GZIP\",\"$READ1\",\"$READ2\"" >> ${ALIGNDIR}/run_config.csv

log "Running STARsolo..."
if [[ $PAIRED == "True" ]]
then
  # Note the R1/R2 order of input fastq reads and --soloStrand Forward for 5' paired-end experiment
  $STAR \
    --runThreadN $THREADS \
    --genomeDir $REFERENCE \
    --readFilesIn $READ1 $READ2 \
    --runDirPerm All_RWX $GZIP \
    $BAM \
    --soloBarcodeMate 1 \
    --clip5pNbases 39 0 \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist $BC \
    --soloCBstart 1 \
    --soloCBlen $CBLEN \
    --soloUMIstart $((CBLEN+1)) \
    --soloUMIlen $UMILEN \
    --soloStrand Forward \
    --soloUMIdedup 1MM_CR \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter EmptyDrops_CR \
    --outFilterScoreMin 30 \
    --soloFeatures Gene GeneFull Velocyto \
    --outTmpDir ${ALIGNDIR}/_STARtmp \
    --outFileNamePrefix ${ALIGNDIR}/ \
    --soloOutFileNames outs/ features.tsv barcodes.tsv matrix.mtx \
    --soloMultiMappers EM
else 
  $STAR \
    --runThreadN $THREADS \
    --genomeDir $REFERENCE \
    --readFilesIn $READ2 $READ1 \
    --runDirPerm All_RWX \
    $GZIP \
    $BAM \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist $BC \
    --soloBarcodeReadLength 0 \
    --soloCBlen $CBLEN \
    --soloUMIstart $((CBLEN+1)) \
    --soloUMIlen $UMILEN \
    --soloStrand $STRAND \
    --soloUMIdedup 1MM_CR \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloCellFilter EmptyDrops_CR \
    --clipAdapterType CellRanger4 \
    --outFilterScoreMin 30 \
    --soloFeatures Gene GeneFull Velocyto \
    --outTmpDir ${ALIGNDIR}/_STARtmp \
    --outFileNamePrefix ${ALIGNDIR}/ \
    --soloOutFileNames outs/ features.tsv barcodes.tsv matrix.mtx \
    --soloMultiMappers EM
fi

log "Indexing BAM file..."
# Index the BAM file
if [[ -s ${ALIGNDIR}/Aligned.sortedByCoord.out.bam ]]
then
  $SAMTOOLS index -@8 ${ALIGNDIR}/Aligned.sortedByCoord.out.bam
fi

log "Compressing unmapped reads..."
# Max-CR bzip all unmapped reads with multicore pbzip2 
if [[ -s ${ALIGNDIR}/Unmapped.out.mate1 ]]
then
  $PIGZ -9 ${ALIGNDIR}/Unmapped.out.mate1 &
  $PIGZ -9 ${ALIGNDIR}/Unmapped.out.mate2 &
  wait
fi

log "Removing test files..."
rm -rf ${ALIGNDIR}/test.R?.fastq ${ALIGNDIR}/test_forward ${ALIGNDIR}/test_reverse

if [[ -d ${ALIGNDIR}/outs ]]; then
  log "Compressing STARsolo output..."

  for subdir in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered Velocyto/raw Velocyto/filtered; do
    TARGET_DIR="${ALIGNDIR}/outs/${subdir}"
    if [[ -d "$TARGET_DIR" ]]; then
      for file in "$TARGET_DIR"/*; do
        if [[ -f "$file" && ! "$file" =~ \.gz$ ]]; then
          # log "Compressing $file"
          gzip "$file" &
        fi
      done
    fi
  done

  wait
  log "Compression complete."
fi


wait
log "STARsolo completed for $SAMPLE"
