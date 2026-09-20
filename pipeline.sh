#!/usr/bin/env bash
# Single-sample paired-end NGS RNA-editing analysis pipeline
#
# Usage:
#   bash pipeline.sh <R1.fastq.gz> <R2.fastq.gz> <reference.fa> <target_position> <output_dir> [sample_name]
#
# Example:
#   bash pipeline.sh sample_1.fq.gz sample_2.fq.gz ref.fa 132 results sample
#
# Optional environment variables:
#   THREADS=6
#   FASTP_QUALITY=30
#   REDITOOLS_SPLITS=10
#   PYTHON_BIN=python
#   REDITOOLS_SCRIPT=/root/software/reditools2.0/src/cineca/reditools.py

set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash pipeline.sh <R1.fastq.gz> <R2.fastq.gz> <reference.fa> <target_position> <output_dir> [sample_name]

Required software:
  FastQC, fastp, Bowtie2, SAMtools, Python, and REDItools2
EOF
}

if [[ $# -lt 5 || $# -gt 6 ]]; then
    usage >&2
    exit 2
fi

R1=$1
R2=$2
REFERENCE=$3
TARGET_POSITION=$4
OUTPUT_DIR=$5
SAMPLE_NAME=${6:-}

THREADS=${THREADS:-6}
FASTP_QUALITY=${FASTP_QUALITY:-30}
REDITOOLS_SPLITS=${REDITOOLS_SPLITS:-10}
PYTHON_BIN=${PYTHON_BIN:-python}
REDITOOLS_SCRIPT=${REDITOOLS_SCRIPT:-/root/software/reditools2.0/src/cineca/reditools.py}

fail() {
    echo "[ERROR] $*" >&2
    exit 1
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

for file in "$R1" "$R2" "$REFERENCE"; do
    [[ -s "$file" ]] || fail "Input file does not exist or is empty: $file"
done

[[ "$TARGET_POSITION" =~ ^[1-9][0-9]*$ ]] || fail "target_position must be a positive integer"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || fail "THREADS must be a positive integer"
[[ "$FASTP_QUALITY" =~ ^[0-9]+$ ]] || fail "FASTP_QUALITY must be an integer"
[[ "$REDITOOLS_SPLITS" =~ ^[1-9][0-9]*$ ]] || fail "REDITOOLS_SPLITS must be a positive integer"

for command_name in fastqc fastp bowtie2 bowtie2-build samtools "$PYTHON_BIN"; do
    command -v "$command_name" >/dev/null 2>&1 || fail "Required command not found: $command_name"
done
[[ -f "$REDITOOLS_SCRIPT" ]] || fail "REDItools2 script not found: $REDITOOLS_SCRIPT"

if [[ -z "$SAMPLE_NAME" ]]; then
    SAMPLE_NAME=$(basename "$R1")
    SAMPLE_NAME=${SAMPLE_NAME%.fastq.gz}
    SAMPLE_NAME=${SAMPLE_NAME%.fq.gz}
    SAMPLE_NAME=${SAMPLE_NAME%.fastq}
    SAMPLE_NAME=${SAMPLE_NAME%.fq}
    SAMPLE_NAME=${SAMPLE_NAME%_R1_001}
    SAMPLE_NAME=${SAMPLE_NAME%_R1}
    SAMPLE_NAME=${SAMPLE_NAME%_1}
fi
[[ "$SAMPLE_NAME" =~ ^[A-Za-z0-9._-]+$ ]] || fail "sample_name may contain only letters, numbers, dots, underscores, and hyphens"

REFERENCE_COUNT=$(grep -c '^>' "$REFERENCE" || true)
[[ "$REFERENCE_COUNT" -eq 1 ]] || fail "reference.fa must contain exactly one reference sequence"

REFERENCE_ID=$(awk '/^>/{sub(/^>/, ""); split($0, fields, /[[:space:]]+/); print fields[1]; exit}' "$REFERENCE")
[[ -n "$REFERENCE_ID" ]] || fail "Could not read the reference ID from reference.fa"

REFERENCE_LENGTH=$(awk '
    /^>/ {if (seen) exit; seen=1; next}
    seen {gsub(/[[:space:]]/, ""); length_sum += length($0)}
    END {print length_sum + 0}
' "$REFERENCE")
[[ "$TARGET_POSITION" -le "$REFERENCE_LENGTH" ]] || fail "target_position exceeds the reference length (${REFERENCE_LENGTH} nt)"

mkdir -p \
    "$OUTPUT_DIR/01-fastqc" \
    "$OUTPUT_DIR/02-fastp" \
    "$OUTPUT_DIR/03-fastqc" \
    "$OUTPUT_DIR/04-bowtie2" \
    "$OUTPUT_DIR/05-reditools2" \
    "$OUTPUT_DIR/ref" \
    "$OUTPUT_DIR/logs"

REFERENCE_COPY="$OUTPUT_DIR/ref/reference.fa"
INDEX_PREFIX="$OUTPUT_DIR/ref/reference"
CLEAN_R1="$OUTPUT_DIR/02-fastp/${SAMPLE_NAME}_1.clean.fq.gz"
CLEAN_R2="$OUTPUT_DIR/02-fastp/${SAMPLE_NAME}_2.clean.fq.gz"
SORTED_BAM="$OUTPUT_DIR/04-bowtie2/${SAMPLE_NAME}.sorted.bam"
REDITOOLS_OUTPUT="$OUTPUT_DIR/05-reditools2/${SAMPLE_NAME}"
TARGET_OUTPUT="$OUTPUT_DIR/05-reditools2/${SAMPLE_NAME}_target.tsv"

cp "$REFERENCE" "$REFERENCE_COPY"

log "Step 1/5: FastQC on raw paired-end reads"
fastqc -t "$THREADS" -o "$OUTPUT_DIR/01-fastqc" "$R1" "$R2" \
    >"$OUTPUT_DIR/logs/01_fastqc.log" 2>&1

log "Step 2/5: Read preprocessing with fastp"
fastp \
    -q "$FASTP_QUALITY" \
    -w "$THREADS" \
    -i "$R1" \
    -I "$R2" \
    -o "$CLEAN_R1" \
    -O "$CLEAN_R2" \
    -h "$OUTPUT_DIR/02-fastp/${SAMPLE_NAME}.html" \
    -j "$OUTPUT_DIR/02-fastp/${SAMPLE_NAME}.json" \
    >"$OUTPUT_DIR/logs/02_fastp.log" 2>&1

log "Step 3/5: FastQC on processed reads"
fastqc -t "$THREADS" -o "$OUTPUT_DIR/03-fastqc" "$CLEAN_R1" "$CLEAN_R2" \
    >"$OUTPUT_DIR/logs/03_fastqc.log" 2>&1

log "Step 4/5: Paired-end alignment with Bowtie2"
bowtie2-build "$REFERENCE_COPY" "$INDEX_PREFIX" \
    >"$OUTPUT_DIR/logs/04_bowtie2_build.log" 2>&1

bowtie2 \
    -p "$THREADS" \
    -x "$INDEX_PREFIX" \
    -1 "$CLEAN_R1" \
    -2 "$CLEAN_R2" \
    2>"$OUTPUT_DIR/logs/04_bowtie2.log" \
    | samtools sort -@ "$THREADS" -O BAM -o "$SORTED_BAM" -

samtools index -@ "$THREADS" "$SORTED_BAM"

log "Step 5/5: Nucleotide counting with REDItools2"
"$PYTHON_BIN" "$REDITOOLS_SCRIPT" \
    -f "$SORTED_BAM" \
    -r "$REFERENCE_COPY" \
    -o "$REDITOOLS_OUTPUT" \
    --split-bam "$REDITOOLS_SPLITS" \
    >"$OUTPUT_DIR/logs/05_reditools2.log" 2>&1

[[ -s "$REDITOOLS_OUTPUT" ]] || fail "REDItools2 did not produce the expected output: $REDITOOLS_OUTPUT"

awk -F '\t' -v OFS='\t' -v ref="$REFERENCE_ID" -v pos="$TARGET_POSITION" '
    NR == 1 {print; next}
    $1 == ref && $2 == pos {print}
' "$REDITOOLS_OUTPUT" > "$TARGET_OUTPUT"

if [[ $(wc -l < "$TARGET_OUTPUT") -le 1 ]]; then
    log "Warning: no REDItools2 result was reported at ${REFERENCE_ID}:${TARGET_POSITION}"
else
    log "Target-site result written to: $TARGET_OUTPUT"
fi

log "Pipeline completed successfully"
log "Sorted BAM: $SORTED_BAM"
log "Complete REDItools2 table: $REDITOOLS_OUTPUT"
log "Target-site table: $TARGET_OUTPUT"
