from celery import shared_task
from asgiref.sync import async_to_sync
from channels.layers import get_channel_layer
from .models import Project
import time
import logging

logger = logging.getLogger(__name__)

@shared_task
def run_rnaseek_pipeline(project_id):
    """
    Celery task to simulate the RNAseek pipeline with status updates over 5 seconds.
    """
    try:
        project = Project.objects.get(id=project_id)
        logger.info(f"Starting pipeline for project {project.name} (ID: {project_id})")

        # Simulate initial status
        project.status = 'pending'
        project.save()
        logger.debug(f"Project {project.name} status set to 'pending'")
        # Send status update via channels
        channel_layer = get_channel_layer()
        async_to_sync(channel_layer.group_send)(
            f'project_{project_id}',
            {
                'type': 'task_status_update',
                'status': project.status,
                'project_id': str(project.id)
            }
        )
        time.sleep(2)

        # Simulate processing
        project.status = 'processing'
        project.save()
        logger.debug(f"Project {project.name} status set to 'processing'")
        # Send status update via channels
        channel_layer = get_channel_layer()
        async_to_sync(channel_layer.group_send)(
            f'project_{project_id}',
            {
                'type': 'task_status_update',
                'status': project.status,
                'project_id': str(project.id)
            }
        )
        time.sleep(2)

        # Simulate completion
        project.status = 'completed'
        project.save()
        logger.info(f"Project {project.name} completed successfully")
        # Send status update via channels
        channel_layer = get_channel_layer()
        async_to_sync(channel_layer.group_send)(
            f'project_{project_id}',
            {
                'type': 'task_status_update',
                'status': project.status,
                'project_id': str(project.id)
            }
        )

    except Project.DoesNotExist:
        logger.error(f"Project with ID {project_id} not found")
    except Exception as e:
        logger.error(f"Error in pipeline for project {project_id}: {e}")
        project_to_update = None
        # Check if 'project' was defined in the try block and is not None
        if 'project' in locals() and project is not None:
            project_to_update = project
        else:
            try:
                project_to_update = Project.objects.get(id=project_id)
            except Project.DoesNotExist:
                logger.error(f"Project with ID {project_id} not found when attempting to set status to 'failed' in exception handler.")
            except Exception as fetch_exc: # Catching potential errors during the fetch
                logger.error(f"Further error fetching project {project_id} in exception handler: {fetch_exc}")

        if project_to_update:
            project_to_update.status = 'failed'
            project_to_update.save()
            logger.debug(f"Project {project_to_update.name} status set to 'failed' due to error: {e}")
            # Send status update via channels
            channel_layer = get_channel_layer()
            async_to_sync(channel_layer.group_send)(
                f'project_{project_id}', # Use project_id from task argument for group name
                {
                    'type': 'task_status_update',
                    'status': project_to_update.status,
                    'project_id': str(project_to_update.id)
                }
            )
        else:
            logger.error(f"Could not update project {project_id} to 'failed' as it could not be fetched.")
