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
    Run FeatureCounts to quantify reads from BAM files into a single counts.csv file.

    Args:
        project: Project instance (contains species, genome_reference, sequencing_type).
        input_files: QuerySet of ProjectFiles (sorted BAM files from SAMtools).
        output_dir: Directory for FeatureCounts output (counts.csv).

    Returns:
        list: Path to the generated counts.csv file (single file as a list for consistency).
    """
    os.makedirs(output_dir, exist_ok=True)
    counts_file = os.path.join(output_dir, "counts.csv")
    
    # Map species to GFF3 annotation file
    species_to_gff3 = {
        'human': 'Homo_sapiens.GRCh38.114.gff3',
        'mouse': 'Mus_musculus.GRCm39.114.gff3',
        'yeast': 'Saccharomyces_cerevisiae.R64-1-1.114.gff3',
        'arabidopsis': 'Arabidopsis_thaliana.TAIR10.61.gff3',
        'worm': 'Caenorhabditis_elegans.WBcel235.114.gff3',
        'zebrafish': 'Danio_rerio.GRCz11.114.gff3',
        'fly': 'Drosophila_melanogaster.BDGP6.54.61.gff3',
        'rice': 'Oryza_sativa.IRGSP-1.0.61.gff3',
        'maize': 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.61.gff3'
    }
    gff3_file = species_to_gff3.get(project.species.lower(), None)
    if not gff3_file:
        logger.error(f"No GFF3 annotation file defined for species: {project.species}")
        raise RuntimeError(f"No GFF3 annotation file defined for species: {project.species}")

    # Use the provided annotation directory
    gff3_path = os.path.join(settings.BASE_DIR, 'rsa', 'references', 'gff3', gff3_file)
    logger.debug(f"Checking for GFF3 annotation at: {gff3_path}")
    if not os.path.exists(gff3_path):
        logger.error(f"GFF3 annotation file not found: {gff3_path}")
        raise RuntimeError(f"GFF3 annotation file not found: {gff3_path}")

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
    cmd = [
        'featureCounts',
        '-F', 'GFF',  # Specify GFF format for annotation
        '-a', gff3_path,  # Annotation file
        '-o', counts_file,  # Output counts file
        '-g', 'gene_id',  # Count reads per gene ID
        '-t', 'gene',  # Feature type to count
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
        
        # Post-process counts.csv to remove first row and use second row as header
        if os.path.exists(counts_file):
            with open(counts_file, 'r') as f:
                lines = f.readlines()
            if len(lines) > 1:  # Ensure there are at least two lines
                header = lines[1].strip().split('\t')  # Second line is the header
                for i, col in enumerate(header):
                    if i > 0:  # Skip first column (Geneid)
                        filename = os.path.basename(col)
                        # Remove .fastq and everything after it
                        filename = re.sub(r'\.fastq.*$', '', filename)
                        filename = re.sub(r'\.sorted.*$', '', filename)
                        header[i] = filename
                # Write header (second row) and data rows (third row onward)
                with open(counts_file, 'w') as f:
                    f.write('\t'.join(header) + '\n')  # Write modified header
                    f.writelines(lines[2:])  # Write data rows
                logger.info(f"Post-processed counts.csv to use second row as header and removed first row")

        # Register counts.csv file
        if os.path.exists(counts_file):
            file_size = os.path.getsize(counts_file) if os.path.isfile(counts_file) else None
            ProjectFiles.objects.create(
                project=project,
                type='featurecounts_counts',
                path=counts_file,
                is_directory=False,
                file_format='csv',
                size=file_size
            )
            logger.info(f"Registered FeatureCounts output: {counts_file} with size {file_size} bytes")

        return [counts_file]
    
    except subprocess.CalledProcessError as e:
        logger.error(f"FeatureCounts failed: {e.stderr}")
        raise RuntimeError(f"FeatureCounts failed: {e.stderr}")