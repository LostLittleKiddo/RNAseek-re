# rsa/tasks.py
from celery import shared_task
from .models import Project
import time
import logging

logger = logging.getLogger(__name__)

@shared_task
def run_rnaseek_pipeline(project_id):
    """
    Celery task to simulate the RNAseek pipeline by updating project status.
    """
    try:
        # Fetch the project
        project = Project.objects.get(id=project_id)
        
        # Transition to running
        project.status = 'running'
        project.save()
        logger.info(f"Project {project.id} status updated to running")

        # Simulate pipeline processing with a 5-second delay
        time.sleep(5)

        # Transition to completed
        project.status = 'completed'
        project.save()
        logger.info(f"Project {project.id} status updated to completed")

    except Project.DoesNotExist:
        logger.error(f"Project {project_id} not found")
    except Exception as e:
        logger.error(f"Error in run_rnaseek_pipeline for project {project_id}: {e}")