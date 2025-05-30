import os
import subprocess
import logging
from django.conf import settings
from rsa.models import Project, ProjectFiles

logger = logging.getLogger(__name__)

def run_fastqc(project, input_files, output_dir):
    """
    Run FastQC on input FASTQ files and register outputs.
    
    Args:
        project: Project instance.
        input_files: QuerySet of ProjectFiles (input FASTQ files).
        output_dir: Directory for FastQC output.
    
    Returns:
        list: Paths to fastqc_data.txt files.
    """
    os.makedirs(output_dir, exist_ok=True)
    data_txt_paths = []
    
    for input_file in input_files:
        fastq_path = input_file.path
        if not os.path.exists(fastq_path):
            logger.error(f"Input file not found: {fastq_path}")
            raise RuntimeError(f"Input file not found: {fastq_path}")
        
        fastqc_cmd = ['fastqc', fastq_path, '-o', output_dir, '--extract']
        logger.debug(f"FastQC command: {' '.join(fastqc_cmd)}")
        
        try:
            result = subprocess.run(fastqc_cmd, capture_output=True, text=True, check=True)
            logger.info(f"FastQC completed for {fastq_path}")
            
            base_name = os.path.splitext(os.path.basename(fastq_path))[0]
            if fastq_path.endswith('.gz'):
                base_name = os.path.splitext(base_name)[0]
            html_output = os.path.join(output_dir, f"{base_name}_fastqc.html")
            zip_output = os.path.join(output_dir, f"{base_name}_fastqc.zip")
            data_txt = os.path.join(output_dir, f"{base_name}_fastqc", "fastqc_data.txt")
            
            for output_path in [html_output, zip_output, data_txt]:
                if os.path.exists(output_path):
                    ProjectFiles.objects.create(
                        project=project,
                        type='fastqc_output',
                        path=output_path,
                        is_directory=False,
                        file_format=output_path.split('.')[-1]
                    )
                    logger.info(f"Registered FastQC output: {output_path}")
                    if output_path.endswith("fastqc_data.txt"):
                        data_txt_paths.append(output_path)
                else:
                    logger.warning(f"FastQC output not found: {output_path}")
        
        except subprocess.CalledProcessError as e:
            logger.error(f"FastQC failed for {fastq_path}: {e.stderr}")
            raise RuntimeError(f"FastQC failed: {e.stderr}")
    
    return data_txt_paths