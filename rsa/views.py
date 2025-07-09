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
from django.db.models import Q
import csv

# Set up logging for debugging
logger = logging.getLogger(__name__)

# Mock file class to simulate file objects for DeseqMetadataForm
class MockFile:
    def __init__(self, name):
        self.name = name

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
                            'arabidopsis': 'Arabidopsis thaliana (TAIR10)',
                            'worm': 'Caenorhabditis elegans (WBcel235)',
                            'zebrafish': 'Danio rerio (GRCz11)',
                            'fly': 'Drosophila melanogaster (BDGP6)',
                            'human': 'Homo sapiens (GRCh38)',
                            'mouse': 'Mus musculus (GRCm39)',
                            'rice': 'Oryza sativa (IRGSP-1.0)',
                            'yeast': 'Saccharomyces cerevisiae (R64-1-1)',
                            'maize': 'Zea mays (Zm-B73-REFERENCE-NAM-5.0)'

                        }.get(form.cleaned_data['genome_of_interest'], 'Unknown'),
                        pipeline_version='1.0.0',
                        sequencing_type=form.cleaned_data['sequencing_type'],
                        pvalue_cutoff=form.cleaned_data['pvalue_cutoff']
                    )

                    project_dir = os.path.join(settings.MEDIA_ROOT, 'r_fastq', str(session_id), str(project.id))
                    os.makedirs(project_dir, exist_ok=True)

                    total_size = 0
                    for file in form.cleaned_data['files']:
                        file_extension = os.path.splitext(file.name)[1].lower()
                        file_format = 'fastq.gz' if file_extension == '.gz' else 'fastq'
                        file_path = os.path.join(project_dir, file.name)

                        with open(file_path, 'wb+') as destination:
                            for chunk in file.chunks():
                                destination.write(chunk)

                        file_size = os.path.getsize(file_path) if os.path.isfile(file_path) else 0
                        total_size += file_size
                        ProjectFiles.objects.create(
                            project=project,
                            type='input_fastq',
                            path=file_path,
                            is_directory=False,
                            file_format=file_format,
                            size=file_size
                        )
                        logger.info(f"Registered input FASTQ file: {file_path} with size {file_size} bytes")

                        deseq_dir = os.path.join(settings.MEDIA_ROOT, 'deseq', str(session_id), str(project.id))
                        os.makedirs(deseq_dir, exist_ok=True)
                        metadata_path = os.path.join(deseq_dir, 'metadata.csv')
                        with open(metadata_path, 'w', newline='') as csvfile:
                            writer = csv.writer(csvfile)
                            writer.writerow(['sample', 'condition'])
                            for condition in sorted(grouped_samples.keys()):
                                for sample_name in sorted(grouped_samples[condition]):
                                    writer.writerow([sample_name, condition])
                                    logger.debug(f"Writing grouped metadata: sample={sample_name}, condition={condition}")
                            csvfile.flush()  # Ensure file is written
                            os.fsync(csvfile.fileno())  # Force to disk
                        if not os.path.isfile(metadata_path):
                            logger.error(f"Metadata file not created: {metadata_path}")
                            raise RuntimeError(f"Metadata file not created: {metadata_path}")
                        file_size = os.path.getsize(metadata_path)
                        total_size += file_size
                        ProjectFiles.objects.create(
                            project=project,
                            type='deseq_metadata',
                            path=metadata_path,
                            is_directory=False,
                            file_format='csv',
                            size=file_size
                        )
                        logger.info(f"Registered DESeq2 metadata file: {metadata_path} with size {file_size} bytes")

                    project.project_size = total_size
                    project.save()
                    logger.info(f"Updated project {project.name} (ID: {project.id}) with total size {total_size} bytes")

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
            'condition_sample1_control': 'condition1',
            'condition_sample2_control': 'condition1',
            'condition_sample3_control': 'condition1',
            'condition_sample4_treatment': 'condition2',
            'condition_sample5_treatment': 'condition2',
            'condition_sample6_treatment': 'condition2',
            'condition_sample7_treatment': 'condition2',
        }

        # Validate forms
        form = RNAseekForm(project_data)
        # Create mock file objects for DeseqMetadataForm
        mock_files = [MockFile('sample1_control.fastq.gz'), MockFile('sample2_control.fastq.gz'), MockFile('sample3_control.fastq.gz'), MockFile('sample4_treatment.fastq.gz'), MockFile('sample5_treatment.fastq.gz'), MockFile('sample6_treatment.fastq.gz'), MockFile('sample7_treatment.fastq.gz')]
        deseq_form = DeseqMetadataForm(deseq_data, files=mock_files, sequencing_type='single')
        if not form.is_valid():
            logger.warning(f"RNAseek form validation failed for example analysis: {form.errors}")
            return JsonResponse({'error': f"Invalid form data: {form.errors.as_json()}"}, status=400)
        if not deseq_form.is_valid():
            logger.warning(f"DESeq2 metadata form validation failed for example analysis: {deseq_form.errors}")
            return JsonResponse({'error': f"Invalid DESeq2 metadata: {deseq_form.errors.as_json()}"}, status=400)

        # Verify sample files exist
        sample_files = ['sample1_control.fastq.gz', 'sample2_control.fastq.gz', 'sample3_control.fastq.gz', 'sample4_treatment.fastq.gz', 'sample5_treatment.fastq.gz', 'sample6_treatment.fastq.gz', 'sample7_treatment.fastq.gz']
        for file_name in sample_files:
            source_path = os.path.join(settings.BASE_DIR, 'rsa', 'references', 'example' , file_name)
            if not os.path.exists(source_path):
                logger.error(f"Sample file not found: {source_path}")
                return JsonResponse({'error': f"Sample file {file_name} not found"}, status=400)

        if Project.objects.filter(user=user, name=form.cleaned_data['project_name'], is_running=True).exists():
            logger.warning(f"Project {form.cleaned_data['project_name']} is already running")
            return JsonResponse({'error': 'Project is already running'}, status=400)
        
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

        project_dir = os.path.join(settings.MEDIA_ROOT, 'r_fastq', str(session_id), str(project.id))
        os.makedirs(project_dir, exist_ok=True)
        total_size = 0
        for file_name in sample_files:
            source_path = os.path.join(settings.BASE_DIR, 'rsa', 'references', 'example' , file_name)
            dest_path = os.path.join(project_dir, file_name)
            shutil.copy2(source_path, dest_path)
            logger.debug(f"Copied {file_name} to {dest_path}")

            file_size = os.path.getsize(dest_path) if os.path.isfile(dest_path) else 0
            total_size += file_size
            ProjectFiles.objects.create(
                project=project,
                type='input_fastq',
                path=dest_path,
                is_directory=False,
                file_format='fastq.gz',
                size=file_size
            )
            logger.info(f"Registered input FASTQ file: {dest_path} with size {file_size} bytes")

        # Create DESeq2 metadata
        deseq_dir = os.path.join(settings.MEDIA_ROOT, 'deseq', str(session_id), str(project.id))
        os.makedirs(deseq_dir, exist_ok=True)
        metadata_path = os.path.join(deseq_dir, 'metadata.csv')
        with open( metadata_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['sample', 'condition'])
            writer.writerow(['sample1_control', deseq_data['condition1']])
            writer.writerow(['sample2_control', deseq_data['condition1']])
            writer.writerow(['sample3_control', deseq_data['condition1']])
            writer.writerow(['sample4_treatment', deseq_data['condition2']])
            writer.writerow(['sample5_treatment', deseq_data['condition2']])
            writer.writerow(['sample6_treatment', deseq_data['condition2']])
            writer.writerow(['sample7_treatment', deseq_data['condition2']])
            logger.debug(f"Created metadata CSV at {metadata_path}")
            csvfile.flush()  # Ensure file is written
            os.fsync(csvfile.fileno())  # Force to disk
        if not os.path.isfile(metadata_path):
            logger.error(f"Metadata file not created: {metadata_path}")
            raise RuntimeError(f"Metadata file not created: {metadata_path}")
        file_size = os.path.getsize(metadata_path)
        total_size += file_size
        ProjectFiles.objects.create(
            project=project,
            type='deseq_metadata',
            path=metadata_path,
            is_directory=False,
            file_format='csv',
            size=file_size
        )
        logger.info(f"Registered DESeq2 metadata file: {metadata_path} with size {file_size} bytes")
        # Update project_size
        project.project_size = total_size
        project.save()
        logger.info(f"Updated project {project.name} (ID: {project.id}) with total size {total_size} bytes")

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
    logger.debug(f"Results view: session_id from cookie: {session_id}")
    if not session_id:
        logger.error("No session_id in cookies for results view")
        messages.error(request, "Session expired. Please start a new session.")
        return redirect('home')

    try:
        user = User.objects.get(session_id=session_id)
        projects = Project.objects.filter(user=user).order_by('-created_at')
        logger.debug(f"Results view: Found user {user.username} with {projects.count()} projects")
        response = render(request, 'results.html', {
            'projects': projects,
            'session_id': session_id  # Pass session_id to template
        })
        
        # Add cache-control headers to prevent caching
        response['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
        response['Pragma'] = 'no-cache'
        response['Expires'] = '0'
        
        return response
    except User.DoesNotExist:
        logger.error("Invalid session_id for results view: {session_id}")
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
        files = ProjectFiles.objects.filter(project=project).exclude(
            Q(type__in=['samtools_bam', 'samtools_bai']) |
            Q(type='fastqc_output', file_format__in=['txt', 'html'])
        ).order_by('created_at')
        
        # Read metadata.csv
        metadata_content = None
        try:
            metadata_file = ProjectFiles.objects.get(project=project, type='deseq_metadata')
            with open(metadata_file.path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                metadata_content = list(reader)  # List of rows, each row is a list of values
                logger.debug(f"Read metadata.csv for project {project.id}: {metadata_content}")
        except ProjectFiles.DoesNotExist:
            logger.warning(f"No metadata.csv found for project {project.id}")
            metadata_content = []
        except Exception as e:
            logger.error(f"Error reading metadata.csv for project {project.id}: {str(e)}")
            metadata_content = []

        # Read deseq_output.csv
        deseq_output_content = None
        try:
            deseq_output_file = ProjectFiles.objects.get(project=project, type='deseq_output')
            with open(deseq_output_file.path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                deseq_output_content = list(reader)[:6]  # Limit to first 5 rows + header for preview
                logger.debug(f"Read deseq_output.csv for project {project.id}: {deseq_output_content}")
        except ProjectFiles.DoesNotExist:
            logger.warning(f"No deseq_output.csv found for project {project.id}")
            deseq_output_content = []
        except Exception as e:
            logger.error(f"Error reading deseq_output.csv for project {project.id}: {str(e)}")
            deseq_output_content = []

        # Read go_gsea_results.csv
        go_gsea_output_content = None
        try:
            go_gsea_output_file = ProjectFiles.objects.get(project=project, type='go_gsea_output')
            with open(go_gsea_output_file.path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                go_gsea_output_content = list(reader)[:6]  # Limit to first 5 rows + header for preview
                logger.debug(f"Read go_gsea_results.csv for project {project.id}: {go_gsea_output_content}")
        except ProjectFiles.DoesNotExist:
            logger.warning(f"No go_gsea_results.csv found for project {project.id}")
            go_gsea_output_content = []
        except Exception as e:
            logger.error(f"Error reading go_gsea_results.csv for project {project.id}: {str(e)}")
            go_gsea_output_content = []

        # Read kegg_gsea_results.csv
        kegg_gsea_output_content = None
        try:
            kegg_gsea_output_file = ProjectFiles.objects.get(project=project, type='kegg_gsea_output')
            with open(kegg_gsea_output_file.path, 'r') as csvfile:
                reader = csv.reader(csvfile)
                kegg_gsea_output_content = list(reader)[:6]  # Limit to first 5 rows + header for preview
                logger.debug(f"Read kegg_gsea_results.csv for project {project.id}: {kegg_gsea_output_content}")
        except ProjectFiles.DoesNotExist:
            logger.warning(f"No kegg_gsea_results.csv found for project {project.id}")
            kegg_gsea_output_content = []
        except Exception as e:
            logger.error(f"Error reading kegg_gsea_results.csv for project {project.id}: {str(e)}")
            kegg_gsea_output_content = []

        return render(request, 'project_detail.html', {
            'project': project,
            'files': files,
            'metadata_content': metadata_content,
            'deseq_output_content': deseq_output_content,
            'go_gsea_output_content': go_gsea_output_content,
            'kegg_gsea_output_content': kegg_gsea_output_content
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