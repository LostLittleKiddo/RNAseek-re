import os
import subprocess
import logging
from django.conf import settings
from rsa.models import Project, ProjectFiles
from weasyprint import HTML

logger = logging.getLogger(__name__)

def run_fastqc(project, input_files, output_dir):
    """
    Run FastQC on input FASTQ files, register outputs with file sizes, and convert HTML to PDF.
    
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
            pdf_output = os.path.join(output_dir, f"{base_name}_fastqc.pdf")
            zip_output = os.path.join(output_dir, f"{base_name}_fastqc.zip")
            data_txt = os.path.join(output_dir, f"{base_name}_fastqc", "fastqc_data.txt")
            
            # Register FastQC outputs (HTML, ZIP, data_txt)
            for output_path in [html_output, zip_output, data_txt]:
                if os.path.exists(output_path):
                    file_size = os.path.getsize(output_path) if os.path.isfile(output_path) else None
                    ProjectFiles.objects.create(
                        project=project,
                        type='fastqc_output',
                        path=output_path,
                        is_directory=False,
                        file_format=output_path.split('.')[-1],
                        size=file_size
                    )
                    logger.info(f"Registered FastQC output: {output_path} with size {file_size} bytes")
                    if output_path.endswith("fastqc_data.txt"):
                        data_txt_paths.append(output_path)
                else:
                    logger.warning(f"FastQC output not found: {output_path}")
            
            # Convert HTML to PDF using weasyprint
            if os.path.exists(html_output):
                try:
                    HTML(html_output).write_pdf(pdf_output)
                    if os.path.exists(pdf_output):
                        file_size = os.path.getsize(pdf_output)
                        ProjectFiles.objects.create(
                            project=project,
                            type='fastqc_output',
                            path=pdf_output,
                            is_directory=False,
                            file_format='pdf',
                            size=file_size
                        )
                        logger.info(f"Registered FastQC PDF output: {pdf_output} with size {file_size} bytes")
                    else:
                        logger.warning(f"PDF output not created: {pdf_output}")
                except Exception as e:
                    logger.error(f"Failed to convert HTML to PDF for {html_output}: {str(e)}")
                    # Continue pipeline even if PDF conversion fails
            else:
                logger.warning(f"HTML file not found for PDF conversion: {html_output}")
        
        except subprocess.CalledProcessError as e:
            logger.error(f"FastQC failed for {fastq_path}: {e.stderr}")
            raise RuntimeError(f"FastQC failed: {e.stderr}")
    
    return data_txt_paths