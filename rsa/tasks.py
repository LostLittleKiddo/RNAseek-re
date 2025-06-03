# rsa/tasks.py
from celery import shared_task
from channels.layers import get_channel_layer
from asgiref.sync import async_to_sync
from .models import Project, ProjectFiles
import time
import logging
from django.conf import settings
from .util.fastqc import run_fastqc
from .util.trimmomatic import run_trimmomatic
import os

logger = logging.getLogger(__name__)

@shared_task
def run_rnaseek_pipeline(project_id):
    """
    Celery task to run the RNAseek pipeline, starting with FastQC and conditionally Trimmomatic.
    """
    channel_layer = get_channel_layer()
    try:
        project = Project.objects.get(id=project_id)
        logger.info(f"Starting pipeline for project {project.name} (ID: {project_id})")
        
        def update_status(new_status, error_message=None):
            project.status = new_status
            if error_message:
                project.error_message = error_message
            project.save()
            logger.debug(f"Project {project.name} status set to '{new_status}'")
            event = {
                'type': 'project_status_update',
                'project_id': str(project.id),
                'status': new_status,
                'project_name': project.name
            }
            if error_message:
                event['error_message'] = error_message
            async_to_sync(channel_layer.group_send)(
                'project_status',
                event
            )
        
        output_dir = os.path.join(settings.MEDIA_ROOT, 'output', str(project.session_id), str(project.id), 'fastqc')
        input_files = ProjectFiles.objects.filter(project=project, type='input_fastq')
        if not input_files:
            raise ValueError("No input FASTQ files found for the project")
        
        update_status('pending')
        time.sleep(1) 
        
        update_status('processing')
        data_txt_paths = run_fastqc(project, input_files, output_dir)
        logger.info(f"FastQC data files generated: {data_txt_paths}")

        update_status('trimming')
        trimmomatic_output_dir = os.path.join(settings.MEDIA_ROOT, 'output', str(project.session_id), str(project.id), 'trimmomatic')
        run_trimmomatic(project, data_txt_paths, trimmomatic_output_dir, input_files)
        
        update_status('completed')
        logger.info(f"Project {project.name} completed successfully")
    
    except Project.DoesNotExist:
        error_msg = f"Project with ID {project_id} not found"
        logger.error(error_msg)
        async_to_sync(channel_layer.group_send)(
            'project_status',
            {
                'type': 'project_status_update',
                'project_id': str(project_id),
                'status': 'failed',
                'project_name': 'Unknown',
                'error_message': error_msg
            }
        )
        
    except Exception as e:
        error_msg = str(e)
        logger.error(f"Error in pipeline for project {project_id}: {error_msg}")
        project.status = 'failed'
        project.error_message = error_msg
        project.save()
        async_to_sync(channel_layer.group_send)(
            'project_status',
            {
                'type': 'project_status_update',
                'project_id': str(project.id),
                'status': 'failed',
                'project_name': project.name,
                'error_message': error_msg
            }
        )