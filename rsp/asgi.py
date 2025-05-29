# import os

# from django.core.asgi import get_asgi_application
# from channels.routing import ProtocolTypeRouter

# os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'rsp.settings')

# application = ProtocolTypeRouter({
    
#     'http':get_asgi_application()

# })

# rsp/asgi.py
import os
from django.core.asgi import get_asgi_application
from channels.routing import ProtocolTypeRouter, URLRouter
from channels.auth import AuthMiddlewareStack
import importlib

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'rsp.settings')

# Initialize Django ASGI application
django_asgi_app = get_asgi_application()

# Lazily import rsa.routing to avoid AppRegistryNot Potter
application = ProtocolTypeRouter({
    'http': django_asgi_app,
    'websocket': AuthMiddlewareStack(
        URLRouter(
            importlib.import_module('rsa.routing').websocket_urlpatterns
        )
    ),
})