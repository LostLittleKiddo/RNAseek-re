# rsa/views.py
from django.shortcuts import render, redirect
from django.contrib import messages
from .models import User, Project, ProjectFiles
from .forms import RNAseekForm, DeseqMetadataForm
import uuid
import logging
import os
import csv
from django.conf import settings

# Set up logging for debugging
logger = logging.getLogger(__name__)

def home(request):
    # Get session_id from cookie
    session_id = request.COOKIES.get('session_id')
    user = None
    is_new_user = False

    if session_id:
        try:
            # Find user by session_id
            user = User.objects.get(session_id=session_id)
            logger.debug(f"Found user: {user.username} with session_id: {session_id}")
        except User.DoesNotExist:
            logger.warning(f"Invalid session_id: {session_id}")
            pass
        except Exception as e:
            logger.error(f"Error querying user: {e}")
            pass

    if not user:
        # Generate a new session_id
        session_id = str(uuid.uuid4())
        # Create a new user
        user = User.objects.create(
            username=f"guest_{session_id[:8]}",
            session_id=session_id
        )
        is_new_user = True
        logger.info(f"Created new user: {user.username} with session_id: {session_id}")
    else:
        # Invalidate other sessions for this user
        User.objects.filter(id=user.id).exclude(session_id=session_id).update(session_id=None)

    # Update the user's session_id if necessary
    if str(user.session_id) != session_id:
        user.session_id = session_id
        user.save()
        logger.debug(f"Updated session_id for user: {user.username}")

    # Initialize the forms
    form = RNAseekForm()
    deseq_form = None

    if request.method == 'POST':
        form = RNAseekForm(request.POST, request.FILES)
        if form.is_valid():
            try:
                # Log uploaded files for debugging
                uploaded_files = [f.name for f in form.cleaned_data['files']]
                logger.debug(f"Uploaded files: {uploaded_files}")
                logger.debug(f"Sequencing type: {form.cleaned_data['sequencing_type']}")

                # Create DESeq2 metadata form with uploaded files and sequencing type
                deseq_form = DeseqMetadataForm(
                    request.POST,
                    files=form.cleaned_data['files'],
                    sequencing_type=form.cleaned_data['sequencing_type']
                )
                if deseq_form.is_valid():
                    # Create a new Project instance
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
                        pipeline_version='1.0.0',  # Adjust as needed
                        sequencing_type=form.cleaned_data['sequencing_type'],
                        pvalue_cutoff=form.cleaned_data['pvalue_cutoff']
                    )

                    # Create directory structure: r_fastq/sessionid/projectid
                    project_dir = os.path.join(settings.MEDIA_ROOT, 'r_fastq', str(session_id), str(project.id))
                    os.makedirs(project_dir, exist_ok=True)

                    # Save uploaded files
                    for file in form.cleaned_data['files']:
                        file_extension = os.path.splitext(file.name)[1].lower()
                        file_format = 'fastq.gz' if file_extension == '.gz' else 'fastq'
                        file_path = os.path.join(project_dir, file.name)

                        # Write file to disk
                        with open(file_path, 'wb+') as destination:
                            for chunk in file.chunks():
                                destination.write(chunk)

                        # Save file metadata to ProjectFiles
                        ProjectFiles.objects.create(
                            project=project,
                            type='input_fastq',
                            path=file_path,
                            is_directory=False,
                            file_format=file_format
                        )

                    # Create DESeq2 metadata.csv
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

                    logger.info(f"Project {project.name} created for user {user.username} with {len(form.cleaned_data['files'])} files")
                    logger.info(f"DESeq2 metadata saved to {metadata_path}")
                    messages.success(request, f"Project '{project.name}' created successfully! Analysis is pending.")
                    return redirect('home')

                else:
                    logger.warning(f"DESeq2 metadata form validation failed: {deseq_form.errors}")
                    messages.error(request, "Please correct the errors in the DESeq2 metadata form.")

            except Exception as e:
                logger.error(f"Error processing form: {e}")
                messages.error(request, "An error occurred while processing your submission. Please try again.")
        else:
            logger.warning(f"RNAseek form validation failed: {form.errors}")
            messages.error(request, "Please correct the errors in the form.")

    # Render the home page
    response = render(request, 'home.html', {
        'is_new_user': is_new_user,
        'form': form,
        'deseq_form': deseq_form
    })

    # Set cache-control headers to prevent browser from caching form data
    response['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    response['Pragma'] = 'no-cache'
    response['Expires'] = '0'

    # Set the session_id cookie
    response.set_cookie(
        'session_id',
        session_id,
        max_age=30 * 24 * 60 * 60,  # 30 days
        httponly=True,
        secure=request.is_secure(),
        samesite='Lax'
    )

    return response