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
        output_dir: Directory for HISAT2 output (SAM files).
        data_txt_paths: List of paths to FastQC data files (for reference, not used directly).
    
    Returns:
        list: Paths to generated SAM files.
    """
    os.makedirs(output_dir, exist_ok=True)
    sam_files = []
    
    sequencing_type = project.sequencing_type.lower()
    index_base = os.path.join(settings.BASE_DIR, 'rsa', 'references', project.species, project.genome_reference)
    if not os.path.exists(f"{index_base}.1.ht2"):
        logger.error(f"HISAT2 index not found at {index_base}")
        raise RuntimeError(f"HISAT2 index not found: {index_base}")
    
    if sequencing_type == 'paired':
        from .trimmomatic import find_paired_files
        paired_files = find_paired_files(input_files)
        if not paired_files:
            logger.error("No paired-end files found for paired-end project")
            raise RuntimeError("No paired-end files found for paired-end project")
        
        for forward_path, reverse_path in paired_files:
            base_name = os.path.splitext(os.path.basename(forward_path))[0].replace('_tpaired_R1', '')
            output_sam = os.path.join(output_dir, f"{base_name}.sam")
            
            cmd = [
                'hisat2', '-x', index_base,
                '-1', forward_path, '-2', reverse_path,
                '-S', output_sam
            ]
            
            logger.debug(f"HISAT2 command: {' '.join(cmd)}")
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                logger.info(f"HISAT2 completed for pair {forward_path}, {reverse_path}: {output_sam}")
                sam_files.append(output_sam)
                
                ProjectFiles.objects.create(
                    project=project,
                    type='hisat2_sam',
                    path=output_sam,
                    is_directory=False,
                    file_format='sam'
                )
                logger.info(f"Registered HISAT2 output: {output_sam}")
            
            except subprocess.CalledProcessError as e:
                logger.error(f"HISAT2 failed for pair {forward_path}, {reverse_path}: {e.stderr}")
                raise RuntimeError(f"HISAT2 failed: {e.stderr}")
    
    else:
        for input_file in input_files:
            fastq_path = input_file.path
            base_name = os.path.splitext(os.path.basename(fastq_path))[0].replace('_trimmed', '')
            output_sam = os.path.join(output_dir, f"{base_name}.sam")
            
            cmd = [
                'hisat2', '-x', index_base,
                '-U', fastq_path,
                '-S', output_sam
            ]
            
            logger.debug(f"HISAT2 command: {' '.join(cmd)}")
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                logger.info(f"HISAT2 completed for {fastq_path}: {output_sam}")
                sam_files.append(output_sam)
                
                ProjectFiles.objects.create(
                    project=project,
                    type='hisat2_sam',
                    path=output_sam,
                    is_directory=False,
                    file_format='sam'
                )
                logger.info(f"Registered HISAT2 output: {output_sam}")
            
            except subprocess.CalledProcessError as e:
                logger.error(f"HISAT2 failed for {fastq_path}: {e.stderr}")
                raise RuntimeError(f"HISAT2 failed: {e.stderr}")
    
    return sam_files