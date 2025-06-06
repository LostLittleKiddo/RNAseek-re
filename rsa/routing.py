# rsa/routing.py
from django.urls import re_path
from . import consumers

websocket_urlpatterns = [
    re_path(r'ws/projects/(?P<session_id>[^/]+)/$', consumers.ProjectStatusConsumer.as_asgi()),
]