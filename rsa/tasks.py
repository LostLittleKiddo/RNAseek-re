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

@shared_task
def run_rnaseek_pipeline(project_id):
    """
    Celery task to run the RNAseek pipeline, starting with FastQC on input files.
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

                for output_path in [html_output, zip_output]:
                    if os.path.exists(output_path):
                        ProjectFiles.objects.create(
                            project=project,
                            type='fastqc_output',
                            path=output_path,
                            is_directory=False,
                            file_format=output_path.split('.')[-1]
                        )
                        logger.info(f"Registered FastQC output: {output_path}")
                    else:
                        logger.warning(f"Expected FastQC output not found: {output_path}")

            except subprocess.CalledProcessError as e:
                logger.error(f"FastQC failed for {fastq_path}: {e.stderr}")
                raise RuntimeError(f"FastQC failed: {e.stderr}")

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