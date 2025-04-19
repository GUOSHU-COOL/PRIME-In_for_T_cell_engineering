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
