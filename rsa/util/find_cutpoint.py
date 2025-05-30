import os
import logging

logger = logging.getLogger(__name__)

def find_cutpoint(data_txt_paths):
    """
    Process FastQC data.txt files to check 'Per base sequence quality' and 'Adapter Content' statuses.
    Generate Trimmomatic commands for files where either module has 'fail' or 'warn' status.
    
    Args:
        data_txt_paths (list): List of paths to fastqc_data.txt files.
    
    Returns:
        list: List of Trimmomatic commands for files needing trimming.
    """
    trimmomatic_commands = []
    
    for data_txt_path in data_txt_paths:
        if not os.path.exists(data_txt_path):
            logger.error(f"fastqc_data.txt not found: {data_txt_path}")
            continue
        
        per_base_quality_status = None
        adapter_content_status = None
        
        try:
            with open(data_txt_path, 'r') as f:
                for line in f:
                    if line.startswith('>>Per base sequence quality'):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            per_base_quality_status = parts[1].lower()
                    elif line.startswith('>>Adapter Content'):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            adapter_content_status = parts[1].lower()
                    if per_base_quality_status and adapter_content_status:
                        break
        except Exception as e:
            logger.error(f"Error reading {data_txt_path}: {e}")
            continue
        
        if per_base_quality_status in ['fail', 'warn'] or adapter_content_status in ['fail', 'warn']:
            logger.info(f"Trimming needed for {data_txt_path}: Per base quality={per_base_quality_status}, Adapter content={adapter_content_status}")
            
            fastq_dir = os.path.dirname(os.path.dirname(data_txt_path))
            base_name = os.path.basename(data_txt_path).replace('_fastqc/fastqc_data.txt', '')
            input_fastq = os.path.join(fastq_dir, f"{base_name}.fastq.gz")
            output_dir = os.path.join(fastq_dir, 'trimmomatic')
            os.makedirs(output_dir, exist_ok=True)
            output_fastq = os.path.join(output_dir, f"{base_name}_trimmed.fastq.gz")
            
            trimmomatic_cmd = [
                'trimmomatic', 'SE', '-phred33', input_fastq, output_fastq,
                'ILLUMINACLIP:adapters.fa:2:30:10', 'SLIDINGWINDOW:4:20', 'MINLEN:36'
            ]
            command_str = ' '.join(trimmomatic_cmd)
            trimmomatic_commands.append(command_str)
            logger.debug(f"Generated Trimmomatic command: {command_str}")
        else:
            logger.info(f"No trimming needed for {data_txt_path}")
    
    return trimmomatic_commands