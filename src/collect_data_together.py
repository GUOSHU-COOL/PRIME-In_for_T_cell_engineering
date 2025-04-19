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
        clean_data_dir = row.iloc[7]  # Replace with column name if needed
        sample_data = get_sample_info_from_excel("/data5/wangxin/20241001_wcx/shuyu/20250307/need6colunm.xlsx", sample)

        cm_dir = config["cm_dir"]
        clean_data_relative = os.path.relpath(clean_data_dir, '/')
        fq_r1 = os.path.join(cm_dir, clean_data_relative, sample, f"{sample}_R1.fq.gz")
        fq_r2 = os.path.join(cm_dir, clean_data_relative, sample, f"{sample}_R2.fq.gz")

        # 拷贝到目标文件夹，改名防止覆盖
        dst_r1 = os.path.join(merged_dir, f"{sample}_R1.fq.gz")
        dst_r2 = os.path.join(merged_dir, f"{sample}_R2.fq.gz")

        try:
            shutil.copy2(fq_r1, dst_r1)
            shutil.copy2(fq_r2, dst_r2)
            print(f"Copied {sample} R1/R2 to {merged_dir}")
        except FileNotFoundError as e:
            print(f"Warning: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch processor from Excel using config YAML.")
    parser.add_argument("--config", default='./config.yaml', help="Path to YAML config file")
    args = parser.parse_args()

    config = load_config(args.config)
    batch_process_from_excel(config)
