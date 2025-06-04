from celery import shared_task
from channels.layers import get_channel_layer
from asgiref.sync import async_to_sync
from .models import Project, ProjectFiles
import time
import logging
from django.conf import settings
from .util.fastqc import run_fastqc
from .util.trimmomatic import run_trimmomatic, get_trimmomatic_file_ids
import os
from django.db import transaction
from django.core.exceptions import ValidationError

logger = logging.getLogger(__name__)

@shared_task
def run_rnaseek_pipeline(project_id):
    channel_layer = get_channel_layer()
    try:
        with transaction.atomic():
            project = Project.objects.select_for_update().get(id=project_id)
            if project.is_running:
                logger.warning(f"Pipeline already running for project {project.name} (ID: {project_id})")
                raise ValidationError(f"Pipeline already running for project {project_id}")
            project.is_running = True
            project.save()

        def update_status(new_status, error_message=None):
            with transaction.atomic():
                project.status = new_status
                if error_message:
                    project.error_message = error_message
                project.is_running = (new_status not in ['completed', 'failed'])
                project.save()
                logger.debug(f"Project {project.name} status set to '{new_status}'")
                event = {
                    'type': 'project_status_update',
                    'project_id': str(project.id),
                    'status': new_status,
                    'project_name': project.name,
                    'session_id': project.session_id,
                    'pvalue_cutoff': project.pvalue_cutoff,
                    'species': project.species,
                    'genome_reference': project.genome_reference,
                    'pipeline_version': project.pipeline_version,
                    'sequencing_type': project.sequencing_type
                }
                if error_message:
                    event['error_message'] = error_message
                logger.debug(f"Sending WebSocket update to group: project_status_{project.session_id}, event: {event}")
                async_to_sync(channel_layer.group_send)(
                    f'project_status_{project.session_id}',
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
        trimmomatic_results = run_trimmomatic(project, data_txt_paths, trimmomatic_output_dir, input_files)
        logger.info(f"Trimmomatic results: {trimmomatic_results}")

        if trimmomatic_results['trimmed']:
            update_status('post_trimmomatic_check')
            post_trimmomatic_fastqc_dir = os.path.join(settings.MEDIA_ROOT, 'output', str(project.session_id), str(project.id), 'post_trimmomatic_fastqc')
            post_trimmomatic_files = ProjectFiles.objects.filter(
                project=project,
                path__in=trimmomatic_results['trimmed']
            )
            post_trimmomatic_data_txt_paths = run_fastqc(project, post_trimmomatic_files, post_trimmomatic_fastqc_dir)
            logger.info(f"Post-Trimmomatic FastQC data files generated: {post_trimmomatic_data_txt_paths}")
        else:
            logger.info("No trimmed files to run post-Trimmomatic FastQC on")

        update_status('completed')
        logger.info(f"Project {project.name} completed successfully")

    except Exception as e:
        error_msg = str(e)
        logger.error(f"Error in pipeline for project {project_id}: {error_msg}")
        with transaction.atomic():
            project.status = 'failed'
            project.error_message = error_msg
            project.is_running = False
            project.save()
        async_to_sync(channel_layer.group_send)(
            f'project_status_{project.session_id}',
            {
                'type': 'project_status_update',
                'project_id': str(project.id),
                'status': 'failed',
                'project_name': project.name,
                'session_id': project.session_id,
                'pvalue_cutoff': project.pvalue_cutoff,
                'species': project.species,
                'genome_reference': project.genome_reference,
                'pipeline_version': project.pipeline_version,
                'sequencing_type': project.sequencing_type,
                'error_message': error_msg
            }
        )