
# Intro to Bioinformatics — Pipeline Workspace


**Contents**
- **dataset/fastq_files/**: raw FASTQ files (paired-end) 
- **trimmed_fastq/**: output FASTQ files after trimming (empty here until pipeline run).
- **qc/pre_trim/** and **qc/post_trim/**: directories for quality control reports before and after trimming.
- **counts/**: downstream count files or outputs (empty here).
- **scripts/**: contains `run_pipeline.sh` to run the preprocessing pipeline.

Getting started
- Confirm you have the required tools installed (common tools: `fastqc`, `trim_galore` or `cutadapt`, `multiqc`, `bash`).
- From the repository root, review and run the pipeline script:

```bash
bash scripts/run_pipeline.sh
```

Notes
- `dataset/fastq_files/` currently contains paired FASTQ files (SRR*). 
- `qc/pre_trim/` and `qc/post_trim/` are placeholders for FastQC/MultiQC outputs; they may be empty until you run the pipeline.

Suggested next steps
- Open `scripts/run_pipeline.sh` and confirm tool paths and parameters match your system.
- Run the pipeline on a small subset first to validate configuration.
- Add a `requirements.txt` or documentation listing exact software versions if you plan to share or reproduce results.

Contact
- If you want, I can: (a) inspect `scripts/run_pipeline.sh` and suggest improvements, (b) add a `requirements.txt`, or (c) create a small example run command that processes a single sample.
