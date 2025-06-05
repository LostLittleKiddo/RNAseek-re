# rsa/admin.py
# from django.contrib import admin
# from .models import User, Project, ProjectFiles

# @admin.register(User)
# class UserAdmin(admin.ModelAdmin):
#     list_display = ('username', 'session_id', 'created_at')
#     search_fields = ('username', 'session_id')

# @admin.register(Project)
# class ProjectAdmin(admin.ModelAdmin):
#     list_display = ('name', 'user', 'status', 'species', 'sequencing_type', 'created_at')
#     list_filter = ('status', 'species', 'sequencing_type')
#     search_fields = ('name', 'user__username', 'session_id')

# @admin.register(ProjectFiles)
# class ProjectFilesAdmin(admin.ModelAdmin):
#     list_display = ('project', 'type', 'path', 'file_format', 'size', 'created_at')
#     list_filter = ('type', 'file_format')
#     search_fields = ('project__name', 'path')