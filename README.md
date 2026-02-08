
# Intro to Bioinformatics — Pipeline Workspace
# following was written with Visual Studio AI 


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



- Add a `requirements.txt` or documentation listing exact software versions if you plan to share or reproduce results?


