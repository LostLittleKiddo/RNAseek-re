# rsa/consumers.py
import json
from channels.generic.websocket import AsyncWebsocketConsumer
from channels.db import database_sync_to_async
from .models import Project
import logging

logger = logging.getLogger(__name__)

class ProjectStatusConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        self.group_name = 'project_status'
        await self.channel_layer.group_add(self.group_name, self.channel_name)
        await self.accept()
        logger.debug(f"WebSocket connected for channel: {self.channel_name}")

    async def disconnect(self, close_code):
        await self.channel_layer.group_discard(self.group_name, self.channel_name)
        logger.debug(f"WebSocket disconnected for channel: {self.channel_name}")

    async def receive(self, text_data):
        pass

    async def project_status_update(self, event):
        project_id = event['project_id']
        status = event['status']
        project_name = event.get('project_name', '')
        error_message = event.get('error_message', '')

        message = {
            'project_id': project_id,
            'status': status,
            'project_name': project_name
        }
        if error_message:
            message['error_message'] = error_message

        await self.send(text_data=json.dumps(message))
        logger.debug(f"Sent status update: project_id={project_id}, status={status}, error_message={error_message}")