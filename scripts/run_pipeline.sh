# Step 1: Download FASTQs

#!/bin/bash
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

# The following is now written FULLY by me with the help of class info and Bioinformagician on YouTube. 
# If more help outside of those two things happen, it will be noted by the code. Thanks!

# Step 2: Pre-trim Quality Control

# Step 3: Trimming / Cleaning

# Step 4: Post-trim Quality Control

# Step 5: Quantification + Count Matrix


