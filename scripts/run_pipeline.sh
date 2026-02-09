# Step 1: Download FASTQs

#!/bin/bash

SECONDS=0

set -euo pipefail

# Path to existing FASTQ folder
FASTQ_DIR="../dataset/fastq_files"

# SRRs wanted (10 total)
SRR_LIST=(
  SRR36750771
  SRR36750772
  SRR36750775
  SRR36750776
  SRR36750778
  SRR36750779
  SRR36750780
  SRR36750781
  SRR36750783
  SRR36750784
)

for SRR in "${SRR_LIST[@]}"; do
  R1GZ="$FASTQ_DIR/${SRR}_1.fastq.gz"
  R2GZ="$FASTQ_DIR/${SRR}_2.fastq.gz"
  R1="$FASTQ_DIR/${SRR}_1.fastq"
  R2="$FASTQ_DIR/${SRR}_2.fastq"

  # Skip if already done
  if [[ -f "$R1GZ" && -f "$R2GZ" ]]; then
    echo "[OK] $SRR already exists — skipping"
    continue
  fi

  # Clean any partial leftovers from previous failures
  rm -f "$R1" "$R2" "$R1GZ" "$R2GZ"

  echo "[DL] Downloading $SRR"
  fasterq-dump --split-files "$SRR" -O "$FASTQ_DIR"

  echo "[ZIP] Compressing $SRR"
  gzip "$R1"
  gzip "$R2"

  # Safety check
  if [[ ! -s "$R1GZ" || ! -s "$R2GZ" ]]; then
    echo "[ERROR] $SRR gzip failed or empty — stopping"
    exit 1
  fi

  echo "[DONE] $SRR complete"
done

echo "[ALL DONE]"
# This beginning section was written with the help of AI, as there were MANY issues with downloads from previously written code. Some
# issues include: download failures due to lack of space, partial downloads, mostly disk issues. Eventually, I worked around this with 
# the help of google colab, google drive, and onedrive. I then bought a harddrive and offloaded the project, including fastqs, to 
# the said harddrive. 


# running this script from the project root ensures all relative paths work correctly, and that any intermediate files are created in the expected locations.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FASTQ_DIR="$PROJECT_ROOT/dataset/fastq_files"
ADAPTERS="$(brew --prefix trimmomatic)/share/trimmomatic/adapters/TruSeq3-PE.fa"

# ensure script runs from Project root so relative paths work
cd "$PROJECT_ROOT"

# Step 2: Pre-trim Quality Control

# Make directory structure for pre_trim QC outputs
mkdir -p qc/pre_trim/fastqc qc/pre_trim/multiqc

# FastQC on raw FastQs
if ls qc/pre_trim/fastqc/*_fastqc.html 1> /dev/null 2>&1; then
  echo "[OK] Pre-trim FastQC reports already exist — skipping"
else
  echo "[QC] Running FastQC on raw FASTQs"
  fastqc "$FASTQ_DIR"/*.fastq.gz -o qc/pre_trim/fastqc/
fi

# MultiQC aggregation
if [[ -f qc/pre_trim/multiqc/multiqc_report.html ]]; then
  echo "[OK] Pre-trim MultiQC report already exists — skipping"
else
  echo "[QC] Running MultiQC on pre-trim FastQC reports"
  multiqc qc/pre_trim/fastqc/ -o qc/pre_trim/multiqc/
fi


# Step 3: Trimming / Cleaning

mkdir -p trimmedReads
# ensure trimmomatic and adapter file are available (no hardcoded absolute paths)
command -v trimmomatic >/dev/null 2>&1 || { echo "[ERROR] trimmomatic not in PATH"; exit 1; }
if [[ ! -f "$ADAPTERS" ]]; then
  echo "[ERROR] Adapter file not found: $ADAPTERS"
  exit 1
fi

for SRR in "${SRR_LIST[@]}"; do
  PAIR1="trimmedReads/${SRR}_1_paired.fastq.gz"
  PAIR2="trimmedReads/${SRR}_2_paired.fastq.gz"

  if [[ -s "$PAIR1" && -s "$PAIR2" ]]; then
    echo "[OK] $SRR trimmed files exist — skipping"
    continue
  fi

  trimmomatic PE -threads 4 \
    "$FASTQ_DIR/${SRR}_1.fastq.gz" "$FASTQ_DIR/${SRR}_2.fastq.gz" \
    "$PAIR1" /dev/null \
    "$PAIR2" /dev/null \
    ILLUMINACLIP:"$ADAPTERS":2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 \
    2> trimmedReads/${SRR}_trimming.log

  if [[ ! -s "$PAIR1" || ! -s "$PAIR2" ]]; then
    echo "[ERROR] $SRR trimming failed — check trimmedReads/${SRR}_trimming.log"
    exit 1
  fi
done
# Visual Studio Text Editor was used to help debug the above code. 
# The main issue was that the trimmomatic command was not running correctly, and the output files were not being created. 
# After some debugging, I found that the issue was with the way the trimmomatic command was being called. 


# Step 4: Post-trim Quality Control

# Make directory structure for post_trim QC outputs
mkdir -p qc/post_trim/fastqc qc/post_trim/multiqc

# FastQC on trimmed FastQs
if ls qc/post_trim/fastqc/*_fastqc.html 1> /dev/null 2>&1; then
  echo "[OK] Post-trim FastQC reports already exist — skipping"
else
  echo "[QC] Running FastQC on trimmed FASTQs"
  fastqc trimmedReads/*_paired.fastq.gz -o qc/post_trim/fastqc/
fi

# MultiQC aggregation
if [[ -f qc/post_trim/multiqc/multiqc_report.html ]]; then
  echo "[OK] Post-trim MultiQC report already exists — skipping"
else
  echo "[QC] Running MultiQC on post-trim FastQC reports"
  multiqc qc/post_trim/fastqc/ -o qc/post_trim/multiqc/
fi


# Step 5: Quantification + Count Matrix




