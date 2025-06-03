import os
import logging
import subprocess
from django.conf import settings
from rsa.models import Project, ProjectFiles

logger = logging.getLogger(__name__)

def parse_fastqc_data(data_txt_path):
    """
    Parse FastQC data file to analyze 'Per base sequence quality' and 'Adapter Content' status.
    
    Args:
        data_txt_path: Path to fastqc_data.txt file.
    
    Returns:
        dict: Contains 'per_base_quality' and 'adapter_content' with status and details.
              - 'per_base_quality': {'status': 'pass'/'warn'/'fail', 'low_quality_positions': []}
              - 'adapter_content': {'status': 'pass'/'warn'/'fail', 'adapters': {name: max_percentage}}
    """
    try:
        result = {
            'per_base_quality': {'status': 'pass', 'low_quality_positions': []},
            'adapter_content': {'status': 'pass', 'adapters': {}}
        }
        with open(data_txt_path, 'r') as f:
            lines = f.readlines()
        
        in_base_quality_section = False
        for i, line in enumerate(lines):
            if line.startswith(">>Per base sequence quality"):
                result['per_base_quality']['status'] = line.strip().split()[-1]
                in_base_quality_section = True
            elif line.startswith(">>END_MODULE") and in_base_quality_section:
                in_base_quality_section = False
            elif in_base_quality_section and line.startswith("Base"):
                # Skip header
                continue
            elif in_base_quality_section and line.strip() and not line.startswith("#"):
                # Parse base quality data (e.g., "1 22.75470252723626 ...")
                fields = line.strip().split()
                if len(fields) > 1:
                    mean_quality = float(fields[1])
                    if mean_quality < 20:  # Threshold for low quality
                        result['per_base_quality']['low_quality_positions'].append(fields[0])
        
        in_adapter_section = False
        adapter_headers = []
        for i, line in enumerate(lines):
            if line.startswith(">>Adapter Content"):
                result['adapter_content']['status'] = line.strip().split()[-1]
                in_adapter_section = True
            elif line.startswith(">>END_MODULE") and in_adapter_section:
                in_adapter_section = False
            elif in_adapter_section and line.startswith("#Position"):
                adapter_headers = line.strip().split()[1:]  # e.g., ['Illumina Universal Adapter', ...]
            elif in_adapter_section and line.strip() and not line.startswith("#"):
                fields = line.strip().split()
                position = fields[0]
                for idx, value in enumerate(fields[1:], start=0):
                    if float(value) > 0.01:  # Threshold for significant adapter presence (1%)
                        adapter_name = adapter_headers[idx]
                        result['adapter_content']['adapters'][adapter_name] = max(
                            result['adapter_content']['adapters'].get(adapter_name, 0.0),
                            float(value)
                        )
        
        logger.info(f"Parsed {data_txt_path}: Per base quality {result['per_base_quality']['status']}, "
                    f"low quality positions {result['per_base_quality']['low_quality_positions']}, "
                    f"Adapter content {result['adapter_content']['status']}, adapters {result['adapter_content']['adapters']}")
        return result
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
    Generate tailored Trimmomatic command parameters based on FastQC data for a specific FASTQ file.
    
    Args:
        project: Project instance.
        data_txt_paths: List of paths to FastQC data files (one per FASTQ file).
        input_file_path: Path to the input FASTQ file (single-end) or forward read (paired-end).
        paired_file_path: Path to the reverse read FASTQ file (paired-end, optional).
        quality_threshold: Base quality score for trimming (default: 20).
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
    
    # Find the corresponding FastQC data file for this FASTQ file
    fastqc_data_path = None
    for txt in data_txt_paths:
        if base_name in txt:  # Match FASTQ file name to FastQC data file
            fastqc_data_path = txt
            break
    
    if not fastqc_data_path or not os.path.exists(fastqc_data_path):
        logger.error(f"FastQC data file not found for {input_file_path}")
        raise RuntimeError(f"FastQC data file not found for {input_file_path}")
    
    # Parse FastQC data to tailor parameters
    fastqc_results = parse_fastqc_data(fastqc_data_path)
    
    # Default parameters
    window_size = 4
    quality_threshold = quality_threshold
    leading_threshold = 3
    trailing_threshold = 3
    
    # Tailor parameters based on Per base sequence quality
    if fastqc_results['per_base_quality']['status'] in ['warn', 'fail']:
        low_quality_positions = fastqc_results['per_base_quality']['low_quality_positions']
        if low_quality_positions:
            # Increase leading threshold if early bases are low quality
            if any(int(pos.split('-')[0]) <= 5 for pos in low_quality_positions):
                leading_threshold = 10  # Stricter trimming for early bases
            window_size = 3  # Stricter window for quality scanning
            quality_threshold = max(quality_threshold, 25)  # Higher quality threshold
            logger.info(f"Stricter trimming for {input_file_path} due to low quality at positions: {low_quality_positions}")
    
    # Tailor adapter trimming based on Adapter Content
    if fastqc_results['adapter_content']['status'] in ['warn', 'fail'] and fastqc_results['adapter_content']['adapters']:
        if os.path.exists(adapters_path):
            # Use stricter ILLUMINACLIP parameters for significant adapter presence
            seed_mismatches = 2
            palindrome_clip_threshold = 15 if any(v > 0.1 for v in fastqc_results['adapter_content']['adapters'].values()) else 30
            simple_clip_threshold = 7
            cmd.append(f"ILLUMINACLIP:{adapters_path}:{seed_mismatches}:{palindrome_clip_threshold}:{simple_clip_threshold}:1:true")
            logger.info(f"Added adapter trimming for {input_file_path} targeting {fastqc_results['adapter_content']['adapters'].keys()}")
        else:
            logger.warning(f"Adapter file {adapters_path} not found, skipping adapter trimming.")
    
    # Add trimming parameters
    cmd.extend([
        f"LEADING:{leading_threshold}",
        f"TRAILING:{trailing_threshold}",
        f"SLIDINGWINDOW:{window_size}:{quality_threshold}",
        f"MINLEN:{min_length}"
    ])
    
    logger.debug(f"Generated Trimmomatic params for {input_file_path}: {' '.join(cmd)}")
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
    
    # Check sequencing type
    sequencing_type = getattr(project, 'sequencing_type', 'single').lower()
    if sequencing_type not in ['single', 'paired']:
        logger.error(f"Invalid sequencing type: {sequencing_type}")
        raise RuntimeError(f"Invalid sequencing type: {sequencing_type}")

    # Process based on sequencing type
    if sequencing_type == 'paired':
        # Paired-end processing
        paired_files = find_paired_files(input_files)
        if not paired_files:
            logger.error("No paired-end files found for paired-end project")
            raise RuntimeError("No paired-end files found for paired-end project")
        
        any_trimming_executed = False
        for forward_path, reverse_path in paired_files:
            # Find corresponding FastQC data files
            forward_base = os.path.splitext(os.path.basename(forward_path))[0]
            if forward_path.endswith('.gz'):
                forward_base = os.path.splitext(forward_base)[0]
            reverse_base = os.path.splitext(os.path.basename(reverse_path))[0]
            if reverse_path.endswith('.gz'):
                reverse_base = os.path.splitext(reverse_base)[0]
            
            forward_fastqc = next((txt for txt in data_txt_paths if forward_base in txt), None)
            reverse_fastqc = next((txt for txt in data_txt_paths if reverse_base in txt), None)
            
            if not forward_fastqc or not reverse_fastqc:
                logger.error(f"FastQC data file missing for pair: {forward_path}, {reverse_path}")
                raise RuntimeError("FastQC data file missing for paired-end files")
            
            # Check if trimming is needed for this pair
            forward_results = parse_fastqc_data(forward_fastqc)
            reverse_results = parse_fastqc_data(reverse_fastqc)
            if (forward_results['per_base_quality']['status'] == 'pass' and
                forward_results['adapter_content']['status'] == 'pass' and
                reverse_results['per_base_quality']['status'] == 'pass' and
                reverse_results['adapter_content']['status'] == 'pass'):
                logger.info(f"Skipping Trimmomatic for pair {forward_path}, {reverse_path}: "
                            "Both FastQC metrics pass")
                continue
            
            logger.info(f"Trimmomatic required for pair {forward_path}, {reverse_path}: "
                        f"Forward quality {forward_results['per_base_quality']['status']}, "
                        f"Forward adapters {forward_results['adapter_content']['status']}, "
                        f"Reverse quality {reverse_results['per_base_quality']['status']}, "
                        f"Reverse adapters {reverse_results['adapter_content']['status']}")

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
                any_trimming_executed = True
            
            except subprocess.CalledProcessError as e:
                logger.error(f"Trimmomatic failed for paired-end files {forward_path} and {reverse_path}: {e.stderr}")
                raise RuntimeError(f"Trimmomatic failed: {e.stderr}")
        
        return any_trimming_executed
    else:
        # Single-end processing
        any_trimming_executed = False
        for input_file in input_files:
            fastq_path = input_file.path
            if not os.path.exists(fastq_path):
                logger.error(f"Input file not found: {fastq_path}")
                raise RuntimeError(f"Input file not found: {fastq_path}")
            
            # Find corresponding FastQC data file
            base_name = os.path.splitext(os.path.basename(fastq_path))[0]
            if fastq_path.endswith('.gz'):
                base_name = os.path.splitext(base_name)[0]
            fastqc_data_path = next((txt for txt in data_txt_paths if base_name in txt), None)
            
            if not fastqc_data_path:
                logger.error(f"FastQC data file not found for {fastq_path}")
                raise RuntimeError(f"FastQC data file not found for {fastq_path}")
            
            # Check if trimming is needed
            fastqc_results = parse_fastqc_data(fastqc_data_path)
            if (fastqc_results['per_base_quality']['status'] == 'pass' and
                fastqc_results['adapter_content']['status'] == 'pass'):
                logger.info(f"Skipping Trimmomatic for {fastq_path}: Both FastQC metrics pass")
                continue
            
            logger.info(f"Trimmomatic required for {fastq_path}: "
                        f"Quality {fastqc_results['per_base_quality']['status']}, "
                        f"Adapters {fastqc_results['adapter_content']['status']}")

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
                any_trimming_executed = True
            
            except subprocess.CalledProcessError as e:
                logger.error(f"Trimmomatic failed for {fastq_path}: {e.stderr}")
                raise RuntimeError(f"Trimmomatic failed: {e.stderr}")
    
    return any_trimming_executed

def get_trimmomatic_file_ids(project, input_files, data_txt_paths):
    """
    Return ProjectFiles IDs for input files not trimmed and final trimmed output files.
    
    Args:
        project: Project instance.
        input_files: QuerySet of ProjectFiles (input FASTQ files).
        data_txt_paths: List of paths to FastQC data files.
    
    Returns:
        dict: {'untrimmed': [list of ProjectFiles.id], 'trimmed': [list of ProjectFiles.id]}
    """
    untrimmed_ids = []
    sequencing_type = getattr(project, 'sequencing_type', 'single').lower()
    
    if sequencing_type == 'paired':
        paired_files = find_paired_files(input_files)
        if not paired_files:
            logger.error("No paired-end files found for paired-end project")
            raise RuntimeError("No paired-end files found for paired-end project")
        
        for forward_path, reverse_path in paired_files:
            forward_base = os.path.splitext(os.path.basename(forward_path))[0]
            if forward_path.endswith('.gz'):
                forward_base = os.path.splitext(forward_base)[0]
            reverse_base = os.path.splitext(os.path.basename(reverse_path))[0]
            if reverse_path.endswith('.gz'):
                reverse_base = os.path.splitext(reverse_base)[0]
            
            forward_fastqc = next((txt for txt in data_txt_paths if forward_base in txt), None)
            reverse_fastqc = next((txt for txt in data_txt_paths if reverse_base in txt), None)
            
            if not forward_fastqc or not reverse_fastqc:
                logger.error(f"FastQC data file missing for pair: {forward_path}, {reverse_path}")
                raise RuntimeError("FastQC data file missing for paired-end files")
            
            forward_results = parse_fastqc_data(forward_fastqc)
            reverse_results = parse_fastqc_data(reverse_fastqc)
            
            forward_file = input_files.get(path=forward_path)
            reverse_file = input_files.get(path=reverse_path)
            
            if (forward_results['per_base_quality']['status'] == 'pass' and
                forward_results['adapter_content']['status'] == 'pass' and
                reverse_results['per_base_quality']['status'] == 'pass' and
                reverse_results['adapter_content']['status'] == 'pass'):
                untrimmed_ids.extend([forward_file.id, reverse_file.id])
                logger.info(f"Untrimmed paired-end files: {forward_path}, {reverse_path}")
    else:
        for input_file in input_files:
            fastq_path = input_file.path
            base_name = os.path.splitext(os.path.basename(fastq_path))[0]
            if fastq_path.endswith('.gz'):
                base_name = os.path.splitext(base_name)[0]
            fastqc_data_path = next((txt for txt in data_txt_paths if base_name in txt), None)
            
            if not fastqc_data_path:
                logger.error(f"FastQC data file not found for {fastq_path}")
                raise RuntimeError(f"FastQC data file not found for {fastq_path}")
            
            fastqc_results = parse_fastqc_data(fastqc_data_path)
            if (fastqc_results['per_base_quality']['status'] == 'pass' and
                fastqc_results['adapter_content']['status'] == 'pass'):
                untrimmed_ids.append(input_file.id)
                logger.info(f"Untrimmed single-end file: {fastq_path}")
    
    # Get trimmed file IDs
    trimmed_ids = list(ProjectFiles.objects.filter(
        project=project,
        type__in=['trimmomatic_fastq', 'trimmomatic_fastq_paired']
    ).values_list('id', flat=True))
    
    logger.info(f"Untrimmed file IDs: {untrimmed_ids}, Trimmed file IDs: {trimmed_ids}")
    return {'untrimmed': untrimmed_ids, 'trimmed': trimmed_ids}