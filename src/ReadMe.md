# batch_reads_type_classifier.py

## Overview

`batch_reads_type_classifier.py` is a script designed to classify sample types from an Excel file and run batch analysis accordingly. It supports two types of analysis workflows depending on the sample category and logs the results.

## Functionality

- Reads sample metadata from an Excel file.
- Automatically selects and executes the appropriate analysis module based on the sample type:
  - `change_KI_analysis_batch.py` for KI-type samples.
  - `change_T_analysis_batch.py` for T-type samples.
- Generates logs and outputs the results of the analysis.

## Dependencies

- `change_KI_analysis_batch.py`
- `change_T_analysis_batch.py`

## Usage

Run the script with a configuration file in YAML format:

```bash
python batch_reads_type_classifier.py --config ./config.yaml

# sgRNA-caused Off-target Extraction Pipeline

This script is designed to extract potential sgRNA off-target sequences from genome data for specific samples, perform window-based alignment, and output matching results and summary statistics.

## Input Files

### 1. `sgRNA.xlsx`

This Excel file contains:
- `sample`: Original sample name (e.g., AH2-KI)
- `Rename`: The renamed version used in folder names
- `sgRNA1`, `sgRNA2`, `sgRNA3`: sgRNA sequences for each sample

### 2. Per-sample `output.xlsx`

Located in:  
`/data5/wangxin/20241001_wcx/shuyu/20250307/{Rename}/output.xlsx`  
It must contain:
- Columns like `chr_seq`, `Start`, `Send`, etc.
- Only rows with the first column value `offtargetKI` will be used for analysis

### 3. Transgenic `.fa` files

Located in:  
`/data5/wangxin/20241001_wcx/shuyu/genome/{prefix}/{prefix}_transgene.fa`  
The `prefix` is extracted from the `sample` (e.g., AH2-KI → AH2)

## Output Files

For each sgRNA, the following will be generated:

- `extracted_sequences_{rna_number}.fasta`: Extracted target window sequences
- `filtered_results_{rna_number}.json`: Filtered matching results (mismatch <6, ends with GG or starts with CC)
- `window_seq_counts_{rna_number}.xlsx`: Sequence match statistics and normalized values (per 1000 UMIs)
- `summary.xlsx`: Summary count of sgRNA-caused vs random matches for each sample
- `new_output.xlsx`: Modified version of original `output.xlsx` with an additional `Category` column (`sgRNA_caused` or `random`)

## Script Entry Point

The main function used is:

```python
process_sample_blastn(sample_raw_name, sample, output_dir, target_seq, rna_number)
```

### Parameters:
- `sample_raw_name`: Original sample name (e.g., AH2-KI)
- `sample`: Renamed sample identifier (folder name)
- `output_dir`: Directory for storing results
- `target_seq`: sgRNA sequence
- `rna_number`: Index of the sgRNA (1, 2, or 3)

Called in a loop like this:

```python
for _, row in filtered_df.iterrows():
    ...
    for i, sgRNA in enumerate(sgRNAs, start=1):
        if pd.notna(sgRNA):
            process_sample_blastn(...)
```

## Notes

- `Csn-KI` is specially renamed as `CH3-KI`
- Both sense and antisense strand windows are supported (ends with GG or starts with CC)
- The PAM sequence (last 3 bp) is excluded from alignment
- Maximum mismatch allowed: <6

## Dependencies

Install required packages using:

```bash
pip install pandas biopython openpyxl
```
