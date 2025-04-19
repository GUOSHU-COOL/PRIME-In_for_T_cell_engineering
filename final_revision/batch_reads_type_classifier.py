import os
import logging
import pandas as pd
import argparse
import yaml
from datetime import datetime
from typing import Callable

import change_KI_analysis_batch
import change_T_analysis_batch


def setup_logging(output_dir: str, level: str = "INFO") -> None:
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"batch_process_{current_time}.log")

    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )


def load_config(config_path: str) -> dict:
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def validate_paths(config: dict) -> None:
    if not os.path.exists(config["excel_file"]):
        raise FileNotFoundError(f"Excel file not found: {config['excel_file']}")
    if not os.path.isdir(config["cm_dir"]):
        raise NotADirectoryError(f"CM directory not found: {config['cm_dir']}")
    os.makedirs(config["output_evaluate_dir"], exist_ok=True)
    if not os.access(config["output_evaluate_dir"], os.W_OK):
        raise PermissionError(f"Output directory not writable: {config['output_evaluate_dir']}")


def batch_process_from_excel(config: dict) -> None:
    df = pd.read_excel(config["excel_file"], dtype=str)
    logging.info(f"Loaded Excel with {len(df)} rows.")

    for _, row in df.iterrows():
        sample = row['sample']
        process_func: Callable = (
            change_KI_analysis_batch.process_sample if "KI" in sample
            else change_T_analysis_batch.process_sample
        )

        try:
            logging.info(f"Processing sample {sample} with {process_func.__module__}")
            process_func(sample, config["cm_dir"], config["output_evaluate_dir"],config["excel_file"])
        except Exception as e:
            logging.error(f"Error processing {sample}: {e}", exc_info=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch processor from Excel using config YAML.")
    parser.add_argument("--config",default='./config.yaml', help="Path to YAML config file")
    args = parser.parse_args()

    config = load_config(args.config)
    validate_paths(config)
    setup_logging(config["output_evaluate_dir"], config.get("log_level", "INFO"))
    batch_process_from_excel(config)
