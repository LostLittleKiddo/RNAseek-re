from celery import shared_task
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
        time.sleep(2)

        # Simulate processing
        project.status = 'processing'
        project.save()
        logger.debug(f"Project {project.name} status set to 'processing'")
        time.sleep(2)

        # Simulate completion
        project.status = 'completed'
        project.save()
        logger.info(f"Project {project.name} completed successfully")

    except Project.DoesNotExist:
        logger.error(f"Project with ID {project_id} not found")
    except Exception as e:
        logger.error(f"Error in pipeline for project {project_id}: {e}")
        if project:
            project.status = 'failed'
            project.save()