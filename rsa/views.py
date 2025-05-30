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
from django.conf import settings
from django.core.exceptions import PermissionDenied

# Set up logging for debugging
logger = logging.getLogger(__name__)

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