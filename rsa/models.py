# rsa/models.py
import uuid
from django.db import models

class User(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    username = models.CharField(max_length=150, unique=True)
    session_id = models.UUIDField(unique=True, editable=False, null=True)  # Allow null for initial creation
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'user'

    def __str__(self):
        return self.username

class Project(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey('User', on_delete=models.CASCADE, related_name='projects', null=True, blank=True)
    session_id = models.UUIDField(editable=False)
    name = models.CharField(max_length=255)
    status = models.CharField(max_length=50)
    species = models.CharField(max_length=100)
    genome_reference = models.CharField(max_length=255)
    pipeline_version = models.CharField(max_length=50)
    sequencing_type = models.CharField(max_length=20, choices=[('single', 'Single-End'), ('paired', 'Paired-End')], default='single')
    pvalue_cutoff = models.FloatField(default=0.05, help_text="P-value cutoff for differential expression analysis")
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'project'

    def __str__(self):
        return self.name

class ProjectFiles(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    project = models.ForeignKey(Project, on_delete=models.CASCADE, related_name='inputs')
    type = models.CharField(max_length=50)
    path = models.CharField(max_length=255)
    is_directory = models.BooleanField(default=False)
    file_format = models.CharField(max_length=50)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'project_input'

    def __str__(self):
        return f"{self.type} - {self.project.id}"