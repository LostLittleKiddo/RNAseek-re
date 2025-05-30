from celery import shared_task
from channels.layers import get_channel_layer
from asgiref.sync import async_to_sync
from .models import Project, ProjectFiles
import time
import logging
from django.conf import settings
from .util.fastqc import run_fastqc
from .util.find_cutpoint import find_cutpoint
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
        
        output_dir = os.path.join(settings.MEDIA_ROOT, 'output', str(project.session_id), str(project.id), 'fastqc')
        input_files = ProjectFiles.objects.filter(project=project, type='input_fastq')
        if not input_files:
            raise ValueError("No input FASTQ files found for the project")
        
        update_status('pending')
        time.sleep(1)
        
        update_status('processing')
        data_txt_paths = run_fastqc(project, input_files, output_dir)
        
        trimmomatic_commands = find_cutpoint(data_txt_paths)
        logger.info(f"Generated {len(trimmomatic_commands)} Trimmomatic commands")
        
        if trimmomatic_commands:
            update_status('trimming')
            trimmomatic_output_dir = os.path.join(settings.MEDIA_ROOT, 'output', str(project.session_id), str(project.id), 'trimmomatic')
            run_trimmomatic(project, trimmomatic_commands, trimmomatic_output_dir)
        
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