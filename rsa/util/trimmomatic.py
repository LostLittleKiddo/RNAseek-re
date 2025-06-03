import os
import logging
import subprocess
from django.conf import settings
from rsa.models import Project, ProjectFiles

logger = logging.getLogger(__name__)

def parse_fastqc_data(data_txt_path):
    """
    Parse FastQC data file to check for 'Per base sequence quality' or 'Adapter Content' status.
    
    Args:
        data_txt_path: Path to fastqc_data.txt file.
    
    Returns:
        bool: True if 'Per base sequence quality' or 'Adapter Content' is 'fail' or 'warn', else False.
    """
    try:
        with open(data_txt_path, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith(">>Per base sequence quality") or line.startswith(">>Adapter Content"):
                status = line.strip().split()[-1]
                if status in ["fail", "warn"]:
                    logger.info(f"Found {status} in {line.strip()} for {data_txt_path}")
                    return True
        return False
    except Exception as e:
        logger.error(f"Error parsing {data_txt_path}: {e}")
        raise RuntimeError(f"Failed to parse FastQC data file: {e}")

def find_paired_files(input_files):
    """
    Identify paired-end FASTQ files based on naming conventions (e.g., _R1 and _R2).
    
    Args:
        input_files: QuerySet of ProjectFiles (input FASTQ files).
    
    Returns:
        list: List of tuples (forward_path, reverse_path) for paired-end files, or None for single-end.
    """
    paired_files = []
    forward_files = []
    reverse_files = []

    for input_file in input_files:
        filename = os.path.basename(input_file.path)
        if '_R1' in filename:
            forward_files.append(input_file)
        elif '_R2' in filename:
            reverse_files.append(input_file)
        else:
            # Single-end file
            return None

    # Match forward and reverse files
    for f_file in forward_files:
        f_base = os.path.basename(f_file.path).replace('_R1', '')
        for r_file in reverse_files:
            r_base = os.path.basename(r_file.path).replace('_R2', '')
            if f_base == r_base:
                paired_files.append((f_file.path, r_file.path))
                break

    return paired_files if paired_files else None

def generate_trimmomatic_params(project, data_txt_paths, input_file_path, paired_file_path=None, quality_threshold=20, min_length=36):
    """
    Generate Trimmomatic command parameters as a list based on FastQC data.
    
    Args:
        project: Project instance.
        data_txt_paths: List of paths to FastQC data files.
        input_file_path: Path to the input FASTQ file (single-end) or forward read (paired-end).
        paired_file_path: Path to the reverse read FASTQ file (paired-end, optional).
        quality_threshold: Quality score for trimming (default: 20).
        min_length: Minimum read length after trimming (default: 36).
    
    Returns:
        tuple: (List of Trimmomatic parameters, input FASTQ path(s), output FASTQ file name(s)).
    """
    cmd = []
    adapters_path = os.path.join(os.path.dirname(__file__), "adapters.fa")
    
    # Validate input file path(s)
    if not os.path.exists(input_file_path):
        logger.error(f"Input file not found: {input_file_path}")
        raise RuntimeError(f"Input file not found: {input_file_path}")
    
    if paired_file_path and not os.path.exists(paired_file_path):
        logger.error(f"Paired input file not found: {paired_file_path}")
        raise RuntimeError(f"Paired input file not found: {paired_file_path}")
    
    # Generate output file names
    base_name = os.path.splitext(os.path.basename(input_file_path))[0]
    if input_file_path.endswith('.gz'):
        base_name = os.path.splitext(base_name)[0]
    
    if paired_file_path:
        # Paired-end mode
        forward_paired = f"{base_name}_tpaired_R1.fastq"
        forward_unpaired = f"{base_name}_tunpaired_R1.fastq"
        reverse_paired = f"{base_name.replace('_R1', '_R2')}_tpaired_R2.fastq"
        reverse_unpaired = f"{base_name.replace('_R1', '_R2')}_tunpaired_R2.fastq"
        output_files = (forward_paired, forward_unpaired, reverse_paired, reverse_unpaired)
        input_files = (input_file_path, paired_file_path)
    else:
        # Single-end mode
        output_files = f"{base_name}_trimmed.fastq"
        input_files = input_file_path
    
    # Validate FastQC data files
    for txt in data_txt_paths:
        if not os.path.exists(txt):
            logger.error(f"Required file not found: {txt}")
            raise RuntimeError(f"Required file not found: {txt}")
    
    # Check for adapter file existence
    if os.path.exists(adapters_path):
        cmd.append(f"ILLUMINACLIP:{adapters_path}:2:30:10")
    else:
        logger.warning(f"Adapter file {adapters_path} not found, skipping adapter trimming.")
    
    # Add fixed trimming parameters
    cmd.extend([
        "LEADING:3",
        "TRAILING:3",
        "SLIDINGWINDOW:4:15",
        f"MINLEN:{min_length}"
    ])
    
    return cmd, input_files, output_files

def run_trimmomatic(project, data_txt_paths, output_dir, input_files):
    """
    Run Trimmomatic if FastQC indicates issues with 'Per base sequence quality' or 'Adapter Content'.
    
    Args:
        project: Project instance.
        data_txt_paths: List of paths to fastqc_data.txt files.
        output_dir: Directory for Trimmomatic output.
        input_files: QuerySet of ProjectFiles (input FASTQ files).
    
    Returns:
        bool: True if Trimmomatic was executed, False if skipped.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Check FastQC results for any file
    needs_trimming = False
    for data_txt_path in data_txt_paths:
        if parse_fastqc_data(data_txt_path):
            needs_trimming = True
            logger.info(f"Trimmomatic required due to issues in {data_txt_path}")
            break
    
    if not needs_trimming:
        logger.info("No trimming required based on FastQC results.")
        return False
    
    # Check for paired-end files
    paired_files = find_paired_files(input_files)
    
    if paired_files:
        # Paired-end processing
        for forward_path, reverse_path in paired_files:
            cmd_params, input_files, output_files = generate_trimmomatic_params(project, data_txt_paths, forward_path, reverse_path)
            output_forward_paired = os.path.join(output_dir, output_files[0])
            output_forward_unpaired = os.path.join(output_dir, output_files[1])
            output_reverse_paired = os.path.join(output_dir, output_files[2])
            output_reverse_unpaired = os.path.join(output_dir, output_files[3])
            
            trimmomatic_cmd = [
                'trimmomatic', 'PE', '-phred33',
                input_files[0], input_files[1],
                output_forward_paired, output_forward_unpaired,
                output_reverse_paired, output_reverse_unpaired
            ] + cmd_params
            
            logger.debug(f"Trimmomatic command: {' '.join(trimmomatic_cmd)}")
            
            try:
                result = subprocess.run(trimmomatic_cmd, capture_output=True, text=True, check=True)
                logger.info(f"Trimmomatic completed for paired-end files {forward_path} and {reverse_path}")
                
                # Register paired output files
                for output_path, file_type in [
                    (output_forward_paired, 'trimmomatic_fastq_paired'),
                    (output_reverse_paired, 'trimmomatic_fastq_paired')
                ]:
                    if os.path.exists(output_path):
                        ProjectFiles.objects.create(
                            project=project,
                            type=file_type,
                            path=output_path,
                            is_directory=False,
                            file_format='fastq'
                        )
                        logger.info(f"Registered Trimmomatic output: {output_path}")
                    else:
                        logger.warning(f"Trimmomatic output not found: {output_path}")
                
                # Optionally register unpaired files
                for output_path in [output_forward_unpaired, output_reverse_unpaired]:
                    if os.path.exists(output_path):
                        ProjectFiles.objects.create(
                            project=project,
                            type='trimmomatic_fastq_unpaired',
                            path=output_path,
                            is_directory=False,
                            file_format='fastq'
                        )
                        logger.info(f"Registered Trimmomatic unpaired output: {output_path}")
            
            except subprocess.CalledProcessError as e:
                logger.error(f"Trimmomatic failed for paired-end files {forward_path} and {reverse_path}: {e.stderr}")
                raise RuntimeError(f"Trimmomatic failed: {e.stderr}")
    else:
        # Single-end processing
        for input_file in input_files:
            fastq_path = input_file.path
            if not os.path.exists(fastq_path):
                logger.error(f"Input file not found: {fastq_path}")
                raise RuntimeError(f"Input file not found: {fastq_path}")
            
            cmd_params, input_fastq, output_fastq_name = generate_trimmomatic_params(project, data_txt_paths, fastq_path)
            output_fastq = os.path.join(output_dir, output_fastq_name)
            
            trimmomatic_cmd = [
                'trimmomatic', 'SE', '-phred33',
                fastq_path,
                output_fastq,
            ] + cmd_params
            
            logger.debug(f"Trimmomatic command: {' '.join(trimmomatic_cmd)}")
            
            try:
                result = subprocess.run(trimmomatic_cmd, capture_output=True, text=True, check=True)
                logger.info(f"Trimmomatic completed for {fastq_path}")
                
                # Register output file
                if os.path.exists(output_fastq):
                    ProjectFiles.objects.create(
                        project=project,
                        type='trimmomatic_fastq',
                        path=output_fastq,
                        is_directory=False,
                        file_format='fastq'
                    )
                    logger.info(f"Registered Trimmomatic output: {output_fastq}")
                else:
                    logger.warning(f"Trimmomatic output not found: {output_fastq}")
            
            except subprocess.CalledProcessError as e:
                logger.error(f"Trimmomatic failed for {fastq_path}: {e.stderr}")
                raise RuntimeError(f"Trimmomatic failed: {e.stderr}")
    
    return True