# rsa/tasks.py
from celery import shared_task
from channels.layers import get_channel_layer
from asgiref.sync import async_to_sync
from .models import Project, ProjectFiles
import time
import logging
import os
from django.conf import settings

logger = logging.getLogger(__name__)

@shared_task
def run_rnaseek_pipeline(project_id):
    """
    Celery task to simulate the RNAseek pipeline with status updates and create a sample output file.
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

        # Simulate pipeline stages
        update_status('pending')
        time.sleep(2)

        update_status('processing')
        time.sleep(2)

        # Simulate creating an output file
        output_dir = os.path.join(settings.MEDIA_ROOT, 'output', str(project.session_id), str(project.id))
        os.makedirs(output_dir, exist_ok=True)
        output_file_path = os.path.join(output_dir, 'sample_output.txt')
        with open(output_file_path, 'w') as f:
            f.write(f"Sample output for project {project.name}")

        # Register the output file in ProjectFiles
        ProjectFiles.objects.create(
            project=project,
            type='output',
            path=output_file_path,
            is_directory=False,
            file_format='txt'
        )

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
        if project:
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