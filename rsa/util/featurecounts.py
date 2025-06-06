# rsa/util/featurecounts.py
import os
import subprocess
import logging
from django.conf import settings
from rsa.models import Project, ProjectFiles
import re


logger = logging.getLogger(__name__)

def run_featurecounts(project, input_files, output_dir):
    """
    Run FeatureCounts to quantify reads from BAM files into a single counts.txt file.

    Args:
        project: Project instance (contains species, genome_reference, sequencing_type).
        input_files: QuerySet of ProjectFiles (sorted BAM files from SAMtools).
        output_dir: Directory for FeatureCounts output (counts.txt).

    Returns:
        list: Path to the generated counts.txt file (single file as a list for consistency).
    """
    os.makedirs(output_dir, exist_ok=True)
    counts_file = os.path.join(output_dir, "counts.txt")
    
    # Map species to GTF annotation file
    species_to_gtf = {
        'human': 'Homo_sapiens.GRCh38.gtf',
        'mouse': 'Mus_musculus.GRCm39.gtf',
        'yeast': 'Saccharomyces_cerevisiae.R64-1-1.gtf'
    }
    gtf_file = species_to_gtf.get(project.species.lower(), None)
    if not gtf_file:
        logger.error(f"No GTF annotation file defined for species: {project.species}")
        raise RuntimeError(f"No GTF annotation file defined for species: {project.species}")

    # Use the provided annotation directory
    gtf_path = os.path.join(settings.BASE_DIR, 'rsa', 'references', 'featurecounts', gtf_file)
    logger.debug(f"Checking for GTF annotation at: {gtf_path}")
    if not os.path.exists(gtf_path):
        logger.error(f"GTF annotation file not found: {gtf_path}")
        raise RuntimeError(f"GTF annotation file not found: {gtf_path}")

    # Verify FeatureCounts is installed
    try:
        subprocess.run(['featureCounts', '-v'], capture_output=True, text=True, check=True)
        logger.debug("FeatureCounts is installed and accessible")
    except subprocess.CalledProcessError as e:
        logger.error("FeatureCounts is not installed or not found in PATH")
        raise RuntimeError("FeatureCounts is not installed or not found in PATH")

    # Collect BAM file paths
    bam_files = [input_file.path for input_file in input_files if input_file.type == 'samtools_bam']
    if not bam_files:
        logger.error("No BAM files found for FeatureCounts")
        raise RuntimeError("No BAM files found for FeatureCounts")

    # Build FeatureCounts command
        # Build FeatureCounts command
    cmd = [
        'featureCounts',
        '-a', gtf_path,  # Annotation file
        '-o', counts_file,  # Output counts file
        '-g', 'gene_id',  # Count reads per gene ID
        '-t', 'exon',  # Feature type to count
    ]

    # Adjust parameters based on sequencing type
    sequencing_type = project.sequencing_type.lower()
    if sequencing_type == 'paired':
        cmd.append('-p')  # Paired-end mode
        cmd.append('--countReadPairs')  # Count read pairs instead of individual reads

    # Add BAM files to the command
    cmd.extend(bam_files)

    logger.debug(f"FeatureCounts command: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info(f"FeatureCounts completed: {counts_file}")
        
        # Post-process counts.txt to replace full paths with file names in header
                # Post-process counts.txt to replace full paths with file names in header
        if os.path.exists(counts_file):
            with open(counts_file, 'r') as f:
                lines = f.readlines()
            if lines:
                header = lines[1].strip().split('\t')  # Second line is the header
                for i, col in enumerate(header):
                    if i > 0:  # Skip first column (Geneid)
                        filename = os.path.basename(col)
                        # Remove .fastq and everything after it
                        filename = re.sub(r'\.fastq.*$', '', filename)
                        header[i] = filename
                lines[1] = '\t'.join(header) + '\n'
                with open(counts_file, 'w') as f:
                    f.writelines(lines)
                logger.info(f"Post-processed counts.txt to use file names in header")


        # Register counts.txt file
        if os.path.exists(counts_file):
            file_size = os.path.getsize(counts_file) if os.path.isfile(counts_file) else None
            ProjectFiles.objects.create(
                project=project,
                type='featurecounts_counts',
                path=counts_file,
                is_directory=False,
                file_format='txt',
                size=file_size
            )
            logger.info(f"Registered FeatureCounts output: {counts_file} with size {file_size} bytes")

        return [counts_file]
    
    except subprocess.CalledProcessError as e:
        logger.error(f"FeatureCounts failed: {e.stderr}")
        raise RuntimeError(f"FeatureCounts failed: {e.stderr}")