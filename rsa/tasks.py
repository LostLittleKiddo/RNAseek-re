# rsa/tasks.py
from celery import shared_task
from channels.layers import get_channel_layer
from asgiref.sync import async_to_sync
from .models import Project, ProjectFiles
import time
import logging
import os
import subprocess
from django.conf import settings

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
    logger.info(f"find_cutpoint called with data.txt paths: {data_txt_paths}")
    trimmomatic_commands = []

    for data_txt_path in data_txt_paths:
        if not os.path.exists(data_txt_path):
            logger.error(f"fastqc_data.txt not found: {data_txt_path}")
            continue

        # Initialize flags for module statuses
        per_base_quality_status = None
        adapter_content_status = None

        # Parse fastqc_data.txt to find module statuses
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
                        break  # Exit once both statuses are found
        except Exception as e:
            logger.error(f"Error reading {data_txt_path}: {e}")
            continue

        # Check if trimming is needed
        if per_base_quality_status in ['fail', 'warn'] or adapter_content_status in ['fail', 'warn']:
            logger.info(f"Trimming needed for {data_txt_path}: Per base quality={per_base_quality_status}, Adapter content={adapter_content_status}")
            
            # Construct input and output paths
            fastq_dir = os.path.dirname(os.path.dirname(data_txt_path))  # Parent of fastqc_data.txt dir
            base_name = os.path.basename(data_txt_path).replace('_fastqc/fastqc_data.txt', '')
            input_fastq = os.path.join(fastq_dir, f"{base_name}.fastq.gz")
            output_dir = os.path.join(fastq_dir, 'trimmomatic')
            os.makedirs(output_dir, exist_ok=True)
            output_fastq = os.path.join(output_dir, f"{base_name}_trimmed.fastq.gz")

            # Construct Trimmomatic command (single-end, basic parameters)
            trimmomatic_cmd = [
                'trimmomatic',
                'SE',
                '-phred33',
                input_fastq,
                output_fastq,
                'ILLUMINACLIP:adapters.fa:2:30:10',
                'SLIDINGWINDOW:4:20',
                'MINLEN:36'
            ]
            command_str = ' '.join(trimmomatic_cmd)
            trimmomatic_commands.append(command_str)
            logger.debug(f"Generated Trimmomatic command: {command_str}")
        else:
            logger.info(f"No trimming needed for {data_txt_path}: Per base quality={per_base_quality_status}, Adapter content={adapter_content_status}")

    return trimmomatic_commands

@shared_task
def run_rnaseek_pipeline(project_id):
    """
    Celery task to run the RNAseek pipeline, starting with FastQC and conditionally Trimmomatic.
    """
    channel_layer = get_channel_layer()
    try:
        project = Project.objects.get(id=project_id)
        logger.info(f"Starting pipeline for project {project.name} (ID: {project_id})")

        # Helper function to update status and notify
        def update_status(new_status):
            project.status = new_status
            project.save()
            logger.debug(f"Project {project.name} status set to '{new_status}'")
            async_to_sync(channel_layer.group_send)(
                'project_status',
                {
                    'type': 'project_status_update',
                    'project_id': str(project.id),
                    'status': new_status,
                    'project_name': project.name
                }
            )

        # Set up output directory
        output_dir = os.path.join(settings.MEDIA_ROOT, 'output', str(project.session_id), str(project.id), 'fastqc')
        os.makedirs(output_dir, exist_ok=True)

        # Get input FASTQ files
        input_files = ProjectFiles.objects.filter(project=project, type='input_fastq')
        if not input_files:
            raise ValueError("No input FASTQ files found for the project")

        update_status('pending')
        time.sleep(1)  # Brief delay to ensure UI updates

        update_status('processing')
        data_txt_paths = []  # List to collect all fastqc_data.txt paths
        for input_file in input_files:
            fastq_path = input_file.path
            logger.info(f"Running FastQC on {fastq_path}")

            # Verify input file exists
            if not os.path.exists(fastq_path):
                logger.error(f"Input file not found: {fastq_path}")
                raise RuntimeError(f"Input file not found: {fastq_path}")

            # Construct FastQC command
            fastqc_cmd = [
                'fastqc',
                fastq_path,
                '-o', output_dir,
                '--extract'  # Keep ZIP file intact
            ]
            logger.debug(f"FastQC command: {' '.join(fastqc_cmd)}")

            # Execute FastQC
            try:
                result = subprocess.run(
                    fastqc_cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )
                logger.info(f"FastQC completed for {fastq_path}")
                logger.debug(f"FastQC stdout: {result.stdout}")
                logger.debug(f"FastQC stderr: {result.stderr}")

                # List output directory contents for debugging
                output_files = os.listdir(output_dir)
                logger.debug(f"Output directory contents: {output_files}")

                # Handle potential .gz extension in output naming
                base_name = os.path.splitext(os.path.basename(fastq_path))[0]
                if fastq_path.endswith('.gz'):
                    base_name = os.path.splitext(base_name)[0]  # Remove .gz if present
                html_output = os.path.join(output_dir, f"{base_name}_fastqc.html")
                zip_output = os.path.join(output_dir, f"{base_name}_fastqc.zip")
                fastqc_dir = os.path.join(output_dir, f"{base_name}_fastqc")
                data_txt = os.path.join(fastqc_dir, "fastqc_data.txt")

                # Register output files and collect fastqc_data.txt path
                for output_path in [html_output, zip_output, data_txt]:
                    if os.path.exists(output_path):
                        ProjectFiles.objects.create(
                            project=project,
                            type='fastqc_output',
                            path=output_path,
                            is_directory=False,
                            file_format=output_path.split('.')[-1]
                        )
                        logger.info(f"Registered FastQC output: {output_path}")
                        if output_path.endswith("fastqc_data.txt"):
                            data_txt_paths.append(output_path)
                    else:
                        logger.warning(f"Expected FastQC output not found: {output_path}")

            except subprocess.CalledProcessError as e:
                logger.error(f"FastQC failed for {fastq_path}: {e.stderr}")
                raise RuntimeError(f"FastQC failed: {e.stderr}")

        # Call find_cutpoint with all fastqc_data.txt paths
        trimmomatic_commands = find_cutpoint(data_txt_paths)
        logger.info(f"Generated {len(trimmomatic_commands)} Trimmomatic commands: {trimmomatic_commands}")

        # Handle Trimmomatic execution
        trimmed = None
        if trimmomatic_commands:
            trimmed = True
            update_status('trimming')
            trimmomatic_output_dir = os.path.join(settings.MEDIA_ROOT, 'output', str(project.session_id), str(project.id), 'trimmomatic')
            os.makedirs(trimmomatic_output_dir, exist_ok=True)

            if project.sequencing_type == 'single':
                # Single-end: Execute each command as-is
                for cmd in trimmomatic_commands:
                    logger.debug(f"Executing Trimmomatic command: {cmd}")
                    try:
                        result = subprocess.run(
                            cmd,
                            shell=True,
                            capture_output=True,
                            text=True,
                            check=True
                        )
                        logger.info(f"Trimmomatic completed for command: {cmd}")
                        logger.debug(f"Trimmomatic stdout: {result.stdout}")
                        logger.debug(f"Trimmomatic stderr: {result.stderr}")

                        # Register output file
                        output_path = cmd.split()[-3]  # Extract output_fastq from command
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
                # Group commands by sample name and find pairs
                sample_commands = {}
                for cmd in trimmomatic_commands:
                    input_path = cmd.split()[3]  # Extract input_fastq
                    sample_name = os.path.basename(input_path).replace('_R1', '').replace('_R2', '')
                    sample_name = os.path.splitext(sample_name)[0]  # Remove .fastq.gz
                    if sample_name not in sample_commands:
                        sample_commands[sample_name] = []
                    sample_commands[sample_name].append(cmd)

                for sample_name, cmds in sample_commands.items():
                    # Expect two commands (R1 and R2)
                    if len(cmds) != 2:
                        logger.error(f"Paired-end expects R1 and R2 for {sample_name}, found {len(cmds)} commands")
                        raise RuntimeError(f"Mismatched paired-end files for {sample_name}")

                    # Identify R1 and R2
                    r1_cmd = next((c for c in cmds if '_R1' in c), None)
                    r2_cmd = next((c for c in cmds if '_R2' in c), None)
                    if not r1_cmd or not r2_cmd:
                        logger.error(f"Missing R1 or R2 command for {sample_name}")
                        raise RuntimeError(f"Missing paired-end command for {sample_name}")

                    # Extract input and output paths
                    r1_input = r1_cmd.split()[3]
                    r1_output = r1_cmd.split()[4]
                    r2_input = r2_cmd.split()[3]
                    r2_output = r2_cmd.split()[4]
                    unpaired_r1 = os.path.join(trimmomatic_output_dir, f"{sample_name}_unpaired_R1.fastq.gz")
                    unpaired_r2 = os.path.join(trimmomatic_output_dir, f"{sample_name}_unpaired_R2.fastq.gz")

                    # Construct paired-end Trimmomatic command
                    paired_cmd = [
                        'trimmomatic',
                        'PE',
                        '-phred33',
                        r1_input,
                        r2_input,
                        r1_output,
                        unpaired_r1,
                        r2_output,
                        unpaired_r2,
                        'ILLUMINACLIP:adapters.fa:2:30:10',
                        'SLIDINGWINDOW:4:20',
                        'MINLEN:36'
                    ]
                    command_str = ' '.join(paired_cmd)
                    logger.debug(f"Executing paired-end Trimmomatic command: {command_str}")

                    try:
                        result = subprocess.run(
                            command_str,
                            shell=True,
                            capture_output=True,
                            text=True,
                            check=True
                        )
                        logger.info(f"Trimmomatic completed for {sample_name}")
                        logger.debug(f"Trimmomatic stdout: {result.stdout}")
                        logger.debug(f"Trimmomatic stderr: {result.stderr}")

                        # Register output files
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

        # Simulate remaining pipeline (replace with real steps later)
        time.sleep(1)

        update_status('completed')
        logger.info(f"Project {project.name} completed successfully")

    except Project.DoesNotExist:
        logger.error(f"Project with ID {project_id} not found")
        async_to_sync(channel_layer.group_send)(
            'project_status',
            {
                'type': 'project_status_update',
                'project_id': str(project_id),
                'status': 'failed',
                'project_name': 'Unknown'
            }
        )
    except Exception as e:
        logger.error(f"Error in pipeline for project {project_id}: {e}")
        project.status = 'failed'
        project.error_message = str(e)
        project.save()
        async_to_sync(channel_layer.group_send)(
            'project_status',
            {
                'type': 'project_status_update',
                'project_id': str(project.id),
                'status': 'failed',
                'project_name': project.name
            }
        )