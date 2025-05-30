import os
import subprocess
import logging
from rsa.models import Project, ProjectFiles

logger = logging.getLogger(__name__)

def run_trimmomatic(project, trimmomatic_commands, output_dir):
    """
    Execute Trimmomatic commands for single-end or paired-end sequencing.
    
    Args:
        project: Project instance.
        trimmomatic_commands: List of Trimmomatic command strings.
        output_dir: Directory for Trimmomatic output.
    
    Returns:
        bool: True if trimming was performed, None otherwise.
    """
    if not trimmomatic_commands:
        logger.info("No Trimmomatic commands to execute")
        return None
    
    os.makedirs(output_dir, exist_ok=True)
    
    if project.sequencing_type == 'single':
        for cmd in trimmomatic_commands:
            logger.debug(f"Executing Trimmomatic command: {cmd}")
            try:
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
                logger.info(f"Trimmomatic completed for command: {cmd}")
                
                output_path = cmd.split()[-3]
                if os.path.exists(output_path):
                    ProjectFiles.objects.create(
                        project=project,
                        type='trimmomatic_output',
                        path=output_path,
                        is_directory=False,
                        file_format='fastq.gz'
                    )
                    logger.info(f"Registered Trimmomatic output: {output_path}")
                else:
                    logger.warning(f"Trimmomatic output not found: {output_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Trimmomatic failed for command {cmd}: {e.stderr}")
                raise RuntimeError(f"Trimmomatic failed: {e.stderr}")
    
    else:  # Paired-end
        sample_commands = {}
        for cmd in trimmomatic_commands:
            input_path = cmd.split()[3]
            sample_name = os.path.basename(input_path).replace('_R1', '').replace('_R2', '')
            sample_name = os.path.splitext(sample_name)[0]
            if sample_name not in sample_commands:
                sample_commands[sample_name] = []
            sample_commands[sample_name].append(cmd)
        
        for sample_name, cmds in sample_commands.items():
            if len(cmds) != 2:
                logger.error(f"Paired-end expects R1 and R2 for {sample_name}, found {len(cmds)} commands")
                raise RuntimeError(f"Mismatched paired-end files for {sample_name}")
            
            r1_cmd = next((c for c in cmds if '_R1' in c), None)
            r2_cmd = next((c for c in cmds if '_R2' in c), None)
            if not r1_cmd or not r2_cmd:
                logger.error(f"Missing R1 or R2 command for {sample_name}")
                raise RuntimeError(f"Missing paired-end command for {sample_name}")
            
            r1_input = r1_cmd.split()[3]
            r1_output = r1_cmd.split()[4]
            r2_input = r2_cmd.split()[3]
            r2_output = r2_cmd.split()[4]
            unpaired_r1 = os.path.join(output_dir, f"{sample_name}_unpaired_R1.fastq.gz")
            unpaired_r2 = os.path.join(output_dir, f"{sample_name}_unpaired_R2.fastq.gz")
            
            paired_cmd = [
                'trimmomatic', 'PE', '-phred33', r1_input, r2_input,
                r1_output, unpaired_r1, r2_output, unpaired_r2,
                'ILLUMINACLIP:adapters.fa:2:30:10', 'SLIDINGWINDOW:4:20', 'MINLEN:36'
            ]
            command_str = ' '.join(paired_cmd)
            logger.debug(f"Executing paired-end Trimmomatic command: {command_str}")
            
            try:
                result = subprocess.run(command_str, shell=True, capture_output=True, text=True, check=True)
                logger.info(f"Trimmomatic completed for {sample_name}")
                
                for output_path in [r1_output, r2_output, unpaired_r1, unpaired_r2]:
                    if os.path.exists(output_path):
                        ProjectFiles.objects.create(
                            project=project,
                            type='trimmomatic_output',
                            path=output_path,
                            is_directory=False,
                            file_format='fastq.gz'
                        )
                        logger.info(f"Registered Trimmomatic output: {output_path}")
                    else:
                        logger.warning(f"Trimmomatic output not found: {output_path}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Trimmomatic failed for {sample_name}: {e.stderr}")
                raise RuntimeError(f"Trimmomatic failed: {e.stderr}")
    
    return True