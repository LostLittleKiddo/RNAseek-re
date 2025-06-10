import os
import subprocess
import logging
from django.conf import settings
from rsa.models import Project, ProjectFiles

logger = logging.getLogger(__name__)

def run_hisat2(project, input_files, output_dir, data_txt_paths):
    """
    Run HISAT2 alignment on trimmed or untrimmed FASTQ files for a project.
    
    Args:
        project: Project instance (contains species, genome_reference, sequencing_type).
        input_files: QuerySet of ProjectFiles (input FASTQ files, trimmed or untrimmed).
        output_dir: Directory for HISAT2 output (BAM files).
        data_txt_paths: List of paths to FastQC data files (for reference, not used directly).
    
    Returns:
        list: Paths to generated BAM files.
    """
    os.makedirs(output_dir, exist_ok=True)
    bam_files = []
    
    sequencing_type = project.sequencing_type.lower()
    # Map species to index prefix
    species_to_index_prefix = {
        'human': 'genome',
        'mouse': 'genome',
        'yeast': 'yeast_index'
    }
    index_prefix = species_to_index_prefix.get(project.species.lower(), 'genome')
    index_base = os.path.join(settings.BASE_DIR, 'rsa', 'references', 'index', project.species.lower(), f"{index_prefix}")
    index_file = f"{index_base}.1.ht2"
    logger.debug(f"Checking for HISAT2 index at: {index_file}")
    logger.debug(f"BASE_DIR resolved to: {settings.BASE_DIR}")
    if not os.path.exists(index_file):
        logger.error(f"HISAT2 index file does not exist at {index_file}")
        # Additional debug: List files in the directory
        index_dir = os.path.dirname(index_file)
        if os.path.exists(index_dir):
            logger.debug(f"Files in {index_dir}: {os.listdir(index_dir)}")
        else:
            logger.error(f"Index directory does not exist: {index_dir}")
        raise RuntimeError(f"HISAT2 index not found: {index_base}")
    
    if sequencing_type == 'paired':
        from .trimmomatic import find_paired_files
        paired_files = find_paired_files(input_files)
        if not paired_files:
            logger.error("No paired-end files found for paired-end project")
            raise RuntimeError("No paired-end files found for paired-end project")
        
        for forward_path, reverse_path in paired_files:
            # Extract base name by removing _R1 or _R2 (case-insensitive) from forward file
            base_name = os.path.splitext(os.path.basename(forward_path))[0]
            base_name = base_name.replace('_tpaired_R1', '').replace('_tpaired_r1', '')
            base_name = base_name.replace('_R1', '').replace('_r1', '')  # Additional removal for _R1/_r1
            output_bam = os.path.join(output_dir, f"{base_name}.bam")
            
            cmd = [
                'hisat2', '-x', index_base,
                '-1', forward_path, '-2', reverse_path
            ]
            samtools_cmd = ['samtools', 'view', '-b', '-o', output_bam]
            
            logger.debug(f"HISAT2 command: {' '.join(cmd)} | {' '.join(samtools_cmd)}")
            try:
                # Pipe HISAT2 output to samtools view to create BAM
                hisat2_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                samtools_process = subprocess.run(samtools_cmd, stdin=hisat2_process.stdout, capture_output=True, text=True, check=True)
                hisat2_stdout, hisat2_stderr = hisat2_process.communicate()
                
                if hisat2_process.returncode != 0:
                    logger.error(f"HISAT2 failed for pair {forward_path}, {reverse_path}: {hisat2_stderr}")
                    raise RuntimeError(f"HISAT2 failed: {hisat2_stderr}")
                
                logger.info(f"HISAT2 and samtools completed for pair {forward_path}, {reverse_path}: {output_bam}")
                file_size = os.path.getsize(output_bam) if os.path.isfile(output_bam) else None
                bam_files.append(output_bam)
                
                ProjectFiles.objects.create(
                    project=project,
                    type='hisat2_bam',
                    path=output_bam,
                    is_directory=False,
                    file_format='bam',
                    size=file_size
                )
                logger.info(f"Registered HISAT2 output: {output_bam} with size {file_size} bytes")
            
            except subprocess.CalledProcessError as e:
                logger.error(f"HISAT2 or samtools failed for pair {forward_path}, {reverse_path}: {e.stderr}")
                raise RuntimeError(f"HISAT2 or samtools failed: {e.stderr}")
    
    else:
        for input_file in input_files:
            fastq_path = input_file.path
            base_name = os.path.splitext(os.path.basename(fastq_path))[0].replace('_trimmed', '')
            output_bam = os.path.join(output_dir, f"{base_name}.bam")
            
            cmd = [
                'hisat2', '-x', index_base,
                '-U', fastq_path
            ]
            samtools_cmd = ['samtools', 'view', '-b', '-o', output_bam]
            
            logger.debug(f"HISAT2 command: {' '.join(cmd)} | {' '.join(samtools_cmd)}")
            try:
                # Pipe HISAT2 output to samtools view to create BAM
                hisat2_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                samtools_process = subprocess.run(samtools_cmd, stdin=hisat2_process.stdout, capture_output=True, text=True, check=True)
                hisat2_stdout, hisat2_stderr = hisat2_process.communicate()
                
                if hisat2_process.returncode != 0:
                    logger.error(f"HISAT2 failed for {fastq_path}: {hisat2_stderr}")
                    raise RuntimeError(f"HISAT2 failed: {hisat2_stderr}")
                
                logger.info(f"HISAT2 and samtools completed for {fastq_path}: {output_bam}")
                file_size = os.path.getsize(output_bam) if os.path.isfile(output_bam) else None
                bam_files.append(output_bam)
                
                ProjectFiles.objects.create(
                    project=project,
                    type='hisat2_bam',
                    path=output_bam,
                    is_directory=False,
                    file_format='bam',
                    size=file_size
                )
                logger.info(f"Registered HISAT2 output: {output_bam} with size {file_size} bytes")
            
            except subprocess.CalledProcessError as e:
                logger.error(f"HISAT2 or samtools failed for {fastq_path}: {e.stderr}")
                raise RuntimeError(f"HISAT2 or samtools failed: {e.stderr}")
    
    return bam_files