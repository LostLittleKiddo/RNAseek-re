import os
import subprocess
import logging
from django.conf import settings
from rsa.models import Project, ProjectFiles

logger = logging.getLogger(__name__)

def check_fastqc_for_trimmomatic(data_txt_paths):
    """
    Parse fastqc_data.txt files to check 'Per base sequence quality' and 'Adapter Content'.
    Returns True if Trimmomatic should run (FAIL or WARN in either module), False otherwise.

    Args:
        data_txt_paths (list): List of paths to fastqc_data.txt files from FastQC.

    Returns:
        bool: True if Trimmomatic is needed, False if not.
    """
    modules_to_check = ["Per base sequence quality", "Adapter Content"]
    run_trimmomatic = False

    for data_txt_path in data_txt_paths:
        if not os.path.exists(data_txt_path):
            logger.error(f"fastqc_data.txt not found: {data_txt_path}")
            raise RuntimeError(f"fastqc_data.txt not found: {data_txt_path}")

        logger.debug(f"Parsing FastQC data file: {data_txt_path}")
        try:
            with open(data_txt_path, 'r') as f:
                lines = f.readlines()
                current_module = None
                for line in lines:
                    line = line.strip()
                    if line.startswith(">>"):
                        parts = line.split("\t")
                        print(parts)
                        if len(parts) >= 2 and parts[0] == ">>":
                            module_name = parts[1]
                            status = parts[2] if len(parts) > 2 else None
                            if module_name in modules_to_check and status in ["FAIL", "WARN"]:
                                logger.info(f"{module_name} status: {status} in {data_txt_path}, Trimmomatic required")
                                run_trimmomatic = True
                                break  # No need to check further for this file
                    if line == ">>END_MODULE":
                        current_module = None
                if run_trimmomatic:
                    break  # No need to check other files if one already requires Trimmomatic

        except Exception as e:
            logger.error(f"Error parsing {data_txt_path}: {str(e)}")
            raise RuntimeError(f"Failed to parse {data_txt_path}: {str(e)}")

    logger.info(f"Trimmomatic required: {run_trimmomatic}")
    return run_trimmomatic

def run_trimmomatic(project, data_txt_paths, output_dir):

    if check_fastqc_for_trimmomatic(data_txt_paths) == False:
        logger.info("Trimmomatic not required based on FastQC results.")
        return
    else:
        logger.info("Trimmomatic required based on FastQC results. Starting Trimmomatic...")