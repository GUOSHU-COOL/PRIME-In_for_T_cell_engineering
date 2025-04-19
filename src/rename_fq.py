import os
import logging
import pandas as pd
import argparse
import yaml
import shutil
from common_analysis_utils import get_sample_info_from_excel, filter_fastq_by_length

def load_config(config_path: str) -> dict:
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def batch_process_from_excel(config: dict) -> None:
    df = pd.read_excel(config["excel_file"], dtype=str)

    merged_dir = "/data5/shuyu/data/wangchenxin/"
    os.makedirs(merged_dir, exist_ok=True)

    for _, row in df.iterrows():
        sample = row['sample']
        rename = row["Rename"]
        for read_type in ['R1', 'R2']:
            old_name = f"{sample}_{read_type}.fq.gz"
            new_name = f"{rename}_{read_type}.fq.gz"

            old_path = os.path.join(merged_dir, old_name)
            new_path = os.path.join(merged_dir, new_name)

            if os.path.exists(old_path):
                os.rename(old_path, new_path)
                print(f"Renamed: {old_name} → {new_name}")
            else:
                print(f"Missing file: {old_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch processor from Excel using config YAML.")
    parser.add_argument("--config", default='./config.yaml', help="Path to YAML config file")
    args = parser.parse_args()

    config = load_config(args.config)
    batch_process_from_excel(config)
