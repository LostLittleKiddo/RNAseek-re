import os
import subprocess
import logging
from django.conf import settings
from rsa.models import Project, ProjectFiles

logger = logging.getLogger(__name__)

def run_samtools(project, input_files, output_dir):
    """
    Convert SAM files to sorted BAM files and generate index files using SAMtools.
    
    Args:
        project: Project instance.
        input_files: QuerySet of ProjectFiles (SAM files from HISAT2).
        output_dir: Directory for SAMtools output (BAM and BAI files).
    
    Returns:
        list: Paths to generated BAM files.
    """
    os.makedirs(output_dir, exist_ok=True)
    bam_files = []
    
    # Verify SAMtools is installed
    try:
        subprocess.run(['samtools', '--version'], capture_output=True, text=True, check=True)
        logger.debug("SAMtools is installed and accessible")
    except subprocess.CalledProcessError as e:
        logger.error("SAMtools is not installed or not found in PATH")
        raise RuntimeError("SAMtools is not installed or not found in PATH")
    
    for input_file in input_files:
        sam_path = input_file.path
        if not os.path.exists(sam_path):
            logger.error(f"SAM file not found: {sam_path}")
            raise RuntimeError(f"SAM file not found: {sam_path}")
        
        base_name = os.path.splitext(os.path.basename(sam_path))[0]
        bam_output = os.path.join(output_dir, f"{base_name}.bam")
        sorted_bam_output = os.path.join(output_dir, f"{base_name}.sorted.bam")
        bai_output = f"{sorted_bam_output}.bai"
        
        # Step 1: Convert SAM to BAM
        view_cmd = ['samtools', 'view', '-bS', sam_path, '-o', bam_output]
        logger.debug(f"SAMtools view command: {' '.join(view_cmd)}")
        try:
            result = subprocess.run(view_cmd, capture_output=True, text=True, check=True)
            logger.info(f"SAM to BAM conversion completed for {sam_path}: {bam_output}")
        except subprocess.CalledProcessError as e:
            logger.error(f"SAMtools view failed for {sam_path}: {e.stderr}")
            raise RuntimeError(f"SAMtools view failed: {e.stderr}")
        
        # Step 2: Sort BAM
        sort_cmd = ['samtools', 'sort', bam_output, '-o', sorted_bam_output]
        logger.debug(f"SAMtools sort command: {' '.join(sort_cmd)}")
        try:
            result = subprocess.run(sort_cmd, capture_output=True, text=True, check=True)
            logger.info(f"BAM sorting completed for {bam_output}: {sorted_bam_output}")
        except subprocess.CalledProcessError as e:
            logger.error(f"SAMtools sort failed for {bam_output}: {e.stderr}")
            raise RuntimeError(f"SAMtools sort failed: {e.stderr}")
        
        # Step 3: Index sorted BAM
        index_cmd = ['samtools', 'index', sorted_bam_output]
        logger.debug(f"SAMtools index command: {' '.join(index_cmd)}")
        try:
            result = subprocess.run(index_cmd, capture_output=True, text=True, check=True)
            logger.info(f"BAM indexing completed for {sorted_bam_output}: {bai_output}")
        except subprocess.CalledProcessError as e:
            logger.error(f"SAMtools index failed for {sorted_bam_output}: {e.stderr}")
            raise RuntimeError(f"SAMtools index failed: {e.stderr}")
        
        # Register sorted BAM and index files
        for output_path, file_type in [
            (sorted_bam_output, 'samtools_bam'),
            (bai_output, 'samtools_bai')
        ]:
            if os.path.exists(output_path):
                file_size = os.path.getsize(output_path) if os.path.isfile(output_path) else None
                ProjectFiles.objects.create(
                    project=project,
                    type=file_type,
                    path=output_path,
                    is_directory=False,
                    file_format=output_path.split('.')[-1],
                    size=file_size
                )
                logger.info(f"Registered SAMtools output: {output_path} with size {file_size} bytes")
                if output_path.endswith('.bam'):
                    bam_files.append(output_path)
            else:
                logger.warning(f"SAMtools output not found: {output_path}")
        
        # Clean up unsorted BAM file
        if os.path.exists(bam_output):
            os.remove(bam_output)
            logger.info(f"Removed unsorted BAM file: {bam_output}")
    
    return bam_files