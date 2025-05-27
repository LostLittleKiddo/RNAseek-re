# rsa/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('results/', views.results, name='results'),
    path('results/<uuid:project_id>/', views.results, name='result'),
]