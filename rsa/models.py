from django.db import models
import uuid

class User(models.Model):
    username = models.CharField(max_length=100, unique=True)
    session_id = models.CharField(max_length=36, null=True, blank=True, unique=True)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.username

class Project(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    session_id = models.CharField(max_length=36)
    name = models.CharField(max_length=200)
    status = models.CharField(max_length=20, default='pending')
    species = models.CharField(max_length=100)
    genome_reference = models.CharField(max_length=200)
    pipeline_version = models.CharField(max_length=20)
    sequencing_type = models.CharField(max_length=20)
    pvalue_cutoff = models.FloatField(default=0.05)
    created_at = models.DateTimeField(auto_now_add=True)
    error_message = models.TextField(null=True, blank=True)  # New field for error messages

    def __str__(self):
        return f"{self.name} ({self.status})"

class ProjectFiles(models.Model):
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    type = models.CharField(max_length=50)
    path = models.CharField(max_length=500)
    is_directory = models.BooleanField(default=False)
    file_format = models.CharField(max_length=50)

    def __str__(self):
        return f"{self.project.name} - {self.type}"