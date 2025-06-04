# rsa/consumers.py
import json
from channels.generic.websocket import AsyncWebsocketConsumer
from channels.db import database_sync_to_async
from .models import Project
import logging

logger = logging.getLogger(__name__)

class ProjectStatusConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        # Get session_id from query string
        self.session_id = self.scope['url_route']['kwargs'].get('session_id')
        if not self.session_id:
            logger.error("No session_id provided in WebSocket connection")
            await self.close()
            return
        
        self.group_name = f'project_status_{self.session_id}'
        await self.channel_layer.group_add(self.group_name, self.channel_name)
        await self.accept()
        logger.info(f"WebSocket connected for session: {self.session_id}, group: {self.group_name}, channel: {self.channel_name}")

    async def disconnect(self, close_code):
        if hasattr(self, 'group_name'):
            await self.channel_layer.group_discard(self.group_name, self.channel_name)
            logger.info(f"WebSocket disconnected for session: {self.session_id}, group: {self.group_name}, channel: {self.channel_name}")

    async def receive(self, text_data):
        logger.debug(f"Received WebSocket message for session: {self.session_id}: {text_data}")

    async def project_status_update(self, event):
        project_id = event['project_id']
        status = event['status']
        project_name = event.get('project_name', '')
        error_message = event.get('error_message', '')
        session_id = event.get('session_id', '')

        # Log the entire event for debugging
        logger.debug(f"Received project_status_update for session: {self.session_id}, event: {event}")

        # Only process updates matching the session_id
        if session_id != self.session_id:
            logger.debug(f"Ignoring update for session {session_id}; expected {self.session_id}")
            return

        message = {
            'project_id': project_id,
            'status': status,
            'project_name': project_name,
            'session_id': session_id
        }
        if error_message:
            message['error_message'] = error_message

        await self.send(text_data=json.dumps(message))
        logger.info(f"Sent status update: project_id={project_id}, status={status}, session_id={session_id}")