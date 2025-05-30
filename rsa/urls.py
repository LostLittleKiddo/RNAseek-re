# rsa/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('results/', views.results, name='results'),
    path('result/<int:project_id>/', views.project_detail, name='project_detail'),
    path('download/<int:file_id>/', views.download_file, name='download_file'),
    path('example-analysis/', views.example_analysis, name='example_analysis'),
]