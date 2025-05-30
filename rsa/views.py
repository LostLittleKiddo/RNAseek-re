# rsa/views.py
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from django.http import FileResponse, JsonResponse
from .models import User, Project, ProjectFiles
from .forms import RNAseekForm, DeseqMetadataForm
from .tasks import run_rnaseek_pipeline
import uuid
import logging
import os
import csv
import shutil
from django.conf import settings
from django.core.exceptions import PermissionDenied

# Set up logging for debugging
logger = logging.getLogger(__name__)

# Mock file class to simulate file objects for DeseqMetadataForm
class MockFile:
    def __init__(self, name):
        self.name = name

# Existing home and results views remain unchanged
def home(request):
    session_id = request.COOKIES.get('session_id')
    user = None
    is_new_user = False

    if session_id:
        try:
            user = User.objects.get(session_id=session_id)
            logger.debug(f"Found user: {user.username} with session_id: {session_id}")
        except User.DoesNotExist:
            logger.warning(f"Invalid session_id: {session_id}")
            pass
        except Exception as e:
            logger.error(f"Error querying user: {e}")
            pass

    if not user:
        session_id = str(uuid.uuid4())
        user = User.objects.create(
            username=f"guest_{session_id[:8]}",
            session_id=session_id
        )
        is_new_user = True
        logger.info(f"Created new user: {user.username} with session_id: {session_id}")
    else:
        User.objects.filter(id=user.id).exclude(session_id=session_id).update(session_id=None)

    if str(user.session_id) != session_id:
        user.session_id = session_id
        user.save()
        logger.debug(f"Updated session_id for user: {user.username}")

    form = RNAseekForm()
    deseq_form = None

    if request.method == 'POST':
        form = RNAseekForm(request.POST, request.FILES)
        if form.is_valid():
            try:
                uploaded_files = [f.name for f in form.cleaned_data['files']]
                logger.debug(f"Uploaded files: {uploaded_files}")
                logger.debug(f"Sequencing type: {form.cleaned_data['sequencing_type']}")

                deseq_form = DeseqMetadataForm(
                    request.POST,
                    files=form.cleaned_data['files'],
                    sequencing_type=form.cleaned_data['sequencing_type']
                )
                if deseq_form.is_valid():
                    project = Project.objects.create(
                        user=user,
                        session_id=session_id,
                        name=form.cleaned_data['project_name'],
                        status='pending',
                        species=form.cleaned_data['genome_of_interest'],
                        genome_reference={
                            'yeast': 'Saccharomyces cerevisiae (R64-1-1)',
                            'human': 'Homo sapiens (GRCh38)',
                            'mouse': 'Mus musculus (GRCm39)'
                        }.get(form.cleaned_data['genome_of_interest'], 'Unknown'),
                        pipeline_version='1.0.0',
                        sequencing_type=form.cleaned_data['sequencing_type'],
                        pvalue_cutoff=form.cleaned_data['pvalue_cutoff']
                    )

                    project_dir = os.path.join(settings.MEDIA_ROOT, 'r_fastq', str(session_id), str(project.id))
                    os.makedirs(project_dir, exist_ok=True)

                    for file in form.cleaned_data['files']:
                        file_extension = os.path.splitext(file.name)[1].lower()
                        file_format = 'fastq.gz' if file_extension == '.gz' else 'fastq'
                        file_path = os.path.join(project_dir, file.name)

                        with open(file_path, 'wb+') as destination:
                            for chunk in file.chunks():
                                destination.write(chunk)

                        ProjectFiles.objects.create(
                            project=project,
                            type='input_fastq',
                            path=file_path,
                            is_directory=False,
                            file_format=file_format
                        )

                    deseq_dir = os.path.join(settings.MEDIA_ROOT, 'deseq', str(session_id), str(project.id))
                    os.makedirs(deseq_dir, exist_ok=True)
                    metadata_path = os.path.join(deseq_dir, 'metadata.csv')

                    with open(metadata_path, 'w', newline='') as csvfile:
                        writer = csv.writer(csvfile)
                        writer.writerow(['sample', 'condition'])
                        for sample_name in deseq_form.sample_names:
                            condition_field = deseq_form.cleaned_data[f'condition_{sample_name}']
                            condition = deseq_form.cleaned_data['condition1'] if condition_field == 'condition1' else deseq_form.cleaned_data['condition2']
                            writer.writerow([sample_name, condition])
                            logger.debug(f"Metadata entry: sample={sample_name}, condition={condition}")

                    run_rnaseek_pipeline.delay(project.id)
                    logger.info(f"Triggered Celery task for project {project.name} (ID: {project.id})")

                    return JsonResponse({'project_id': str(project.id), 'message': f"Project '{project.name}' created successfully! Analysis is pending."})

                else:
                    logger.warning(f"DESeq2 metadata form validation failed: {deseq_form.errors}")
                    messages.error(request, "Please correct the errors in the DESeq2 metadata form.")
                    return JsonResponse({'error': 'Invalid DESeq2 metadata'}, status=400)

            except Exception as e:
                logger.error(f"Error processing form: {e}")
                messages.error(request, "An error occurred while processing your submission. Please try again.")
                return JsonResponse({'error': 'Form processing failed'}, status=400)
        else:
            logger.warning(f"RNAseek form validation failed: {form.errors}")
            messages.error(request, "Please correct the errors in the form.")
            return JsonResponse({'error': 'Invalid form data'}, status=400)

    response = render(request, 'home.html', {
        'is_new_user': is_new_user,
        'form': form,
        'deseq_form': deseq_form
    })

    response['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    response['Pragma'] = 'no-cache'
    response['Expires'] = '0'

    response.set_cookie(
        'session_id',
        session_id,
        max_age=30 * 24 * 60 * 60,  # 30 days
        httponly=True,
        secure=request.is_secure(),
        samesite='Lax'
    )

    return response

def example_analysis(request):
    if request.method != 'POST':
        logger.warning("Invalid method for example-analysis endpoint")
        return JsonResponse({'error': 'Method not allowed'}, status=405)

    session_id = request.COOKIES.get('session_id')
    if not session_id:
        logger.error("No session_id provided for example analysis")
        return JsonResponse({'error': 'Session expired. Please start a new session.'}, status=401)

    try:
        user = User.objects.get(session_id=session_id)
        logger.debug(f"Found user: {user.username} for example analysis")
    except User.DoesNotExist:
        logger.error("Invalid session_id for example analysis")
        return JsonResponse({'error': 'Invalid session. Please start a new session.'}, status=401)

    try:
        # Predefined parameters
        project_data = {
            'project_name': 'Example Analysis',
            'genome_of_interest': 'yeast',
            'sequencing_type': 'single',
            'pvalue_cutoff': 0.05,
            'example_analysis': 'true',  # Added to bypass files validation
        }
        deseq_data = {
            'condition1': 'control',
            'condition2': 'treatment',
            'condition_sample1': 'condition1',
            'condition_sample2': 'condition2',
        }

        # Validate forms
        form = RNAseekForm(project_data)
        # Create mock file objects for DeseqMetadataForm
        mock_files = [MockFile('sample1.fastq.gz'), MockFile('sample2.fastq.gz')]
        deseq_form = DeseqMetadataForm(deseq_data, files=mock_files, sequencing_type='single')
        if not form.is_valid():
            logger.warning(f"RNAseek form validation failed for example analysis: {form.errors}")
            return JsonResponse({'error': f"Invalid form data: {form.errors.as_json()}"}, status=400)
        if not deseq_form.is_valid():
            logger.warning(f"DESeq2 metadata form validation failed for example analysis: {deseq_form.errors}")
            return JsonResponse({'error': f"Invalid DESeq2 metadata: {deseq_form.errors.as_json()}"}, status=400)

        # Verify sample files exist
        sample_files = ['sample1.fastq.gz', 'sample2.fastq.gz']
        for file_name in sample_files:
            source_path = os.path.join(settings.STATIC_ROOT, 'example', file_name)
            if not os.path.exists(source_path):
                logger.error(f"Sample file not found: {source_path}")
                return JsonResponse({'error': f"Sample file {file_name} not found"}, status=400)

        # Create project
        project = Project.objects.create(
            user=user,
            session_id=session_id,
            name=form.cleaned_data['project_name'],
            status='pending',
            species=form.cleaned_data['genome_of_interest'],
            genome_reference='Saccharomyces cerevisiae (R64-1-1)',
            pipeline_version='1.0.0',
            sequencing_type=form.cleaned_data['sequencing_type'],
            pvalue_cutoff=form.cleaned_data['pvalue_cutoff']
        )

        # Copy sample files from static/example to project directory
        project_dir = os.path.join(settings.MEDIA_ROOT, 'r_fastq', str(session_id), str(project.id))
        os.makedirs(project_dir, exist_ok=True)
        for file_name in sample_files:
            source_path = os.path.join(settings.STATIC_ROOT, 'example', file_name)
            dest_path = os.path.join(project_dir, file_name)
            shutil.copy2(source_path, dest_path)
            logger.debug(f"Copied {file_name} to {dest_path}")

            ProjectFiles.objects.create(
                project=project,
                type='input_fastq',
                path=dest_path,
                is_directory=False,
                file_format='fastq.gz'
            )

        # Create DESeq2 metadata
        deseq_dir = os.path.join(settings.MEDIA_ROOT, 'deseq', str(session_id), str(project.id))
        os.makedirs(deseq_dir, exist_ok=True)
        metadata_path = os.path.join(deseq_dir, 'metadata.csv')
        with open(metadata_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['sample', 'condition'])
            writer.writerow(['sample1', deseq_data['condition1']])
            writer.writerow(['sample2', deseq_data['condition2']])
            logger.debug(f"Created metadata CSV at {metadata_path}")

        # Trigger Celery task
        run_rnaseek_pipeline.delay(project.id)
        logger.info(f"Triggered Celery task for example project {project.name} (ID: {project.id})")

        return JsonResponse({
            'project_id': str(project.id),
            'message': f"Example Analysis '{project.name}' started successfully!"
        })

    except Exception as e:
        logger.error(f"Error in example analysis: {str(e)}")
        return JsonResponse({'error': f"Error starting example analysis: {str(e)}"}, status=500)

def results(request):
    session_id = request.COOKIES.get('session_id')
    if not session_id:
        messages.error(request, "Session expired. Please start a new session.")
        return redirect('home')

    try:
        user = User.objects.get(session_id=session_id)
        projects = Project.objects.filter(user=user).order_by('-created_at')
        return render(request, 'results.html', {'projects': projects})
    except User.DoesNotExist:
        messages.error(request, "Invalid session. Please start a new session.")
        return redirect('home')

def project_detail(request, project_id):
    session_id = request.COOKIES.get('session_id')
    if not session_id:
        messages.error(request, "Session expired. Please start a new session.")
        return redirect('home')

    try:
        user = User.objects.get(session_id=session_id)
        project = get_object_or_404(Project, id=project_id, user=user)
        if project.status != 'completed':
            messages.warning(request, "Project analysis is not yet completed.")
            return redirect('results')

        files = ProjectFiles.objects.filter(project=project).order_by('created_at')
        return render(request, 'project_detail.html', {
            'project': project,
            'files': files
        })
    except User.DoesNotExist:
        messages.error(request, "Invalid session. Please start a new session.")
        return redirect('home')

def download_file(request, file_id):
    session_id = request.COOKIES.get('session_id')
    if not session_id:
        logger.error("No session_id provided for file download")
        raise PermissionDenied("Session expired. Please start a new session.")

    try:
        user = User.objects.get(session_id=session_id)
        project_file = get_object_or_404(ProjectFiles, id=file_id, project__user=user)
        file_path = project_file.path
        logger.debug(f"Attempting to serve file: {file_path}")
        logger.debug(f"Does file exist? {os.path.exists(file_path)}")
        logger.debug(f"File path details: absolute={os.path.abspath(file_path)}, is_file={os.path.isfile(file_path)}")

        if not os.path.exists(file_path):
            logger.error(f"File not found: {file_path}")
            messages.error(request, "File not found.")
            return redirect('project_detail', project_id=project_file.project.id)

        file = open(file_path, 'rb')
        file_name = os.path.basename(file_path)
        response = FileResponse(file, as_attachment=True, filename=file_name)
        logger.info(f"Serving file: {file_name} for project {project_file.project.id}")
        return response
    except User.DoesNotExist:
        logger.error("Invalid session_id for file download")
        raise PermissionDenied("Invalid session. Please start a new session.")
    except Exception as e:
        logger.error(f"Error serving file {file_id}: {str(e)}")
        messages.error(request, "An error occurred while downloading the file.")
        return redirect('project_detail', project_id=project_file.project.id)