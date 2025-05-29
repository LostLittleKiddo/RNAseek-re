import json
from channels.generic.websocket import AsyncWebsocketConsumer
from channels.db import database_sync_to_async # Import for database operations
from .models import Project # Import your Project model

class StatusConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        self.project_id = self.scope['url_route']['kwargs']['project_id']
        self.project_group_name = f'project_{self.project_id}'

        # Join room group
        await self.channel_layer.group_add(
            self.project_group_name,
            self.channel_name
        )

        await self.accept()

        # Optional: Send current status on connect
        current_status = await self.get_project_status(self.project_id)
        if current_status:
            await self.send(text_data=json.dumps({
                'status': current_status,
                'project_id': self.project_id
            }))

    async def disconnect(self, close_code):
        # Leave room group
        await self.channel_layer.group_discard(
            self.project_group_name,
            self.channel_name
        )

    # Receive message from WebSocket (if needed, not strictly necessary for status updates from server)
    # async def receive(self, text_data):
    #     text_data_json = json.loads(text_data)
    #     message = text_data_json['message']
    #
    #     # Send message to room group (example of broadcasting)
    #     # await self.channel_layer.group_send(
    #     #     self.project_group_name,
    #     #     {
    #     #         'type': 'chat_message', # This would call a method named chat_message
    #     #         'message': message
    #     #     }
    #     # )

    # Receive status update from Celery task (via channel layer)
    async def task_status_update(self, event):
        status = event['status']
        project_id = event['project_id']

        # Send message to WebSocket
        await self.send(text_data=json.dumps({
            'status': status,
            'project_id': project_id
        }))

    @database_sync_to_async
    def get_project_status(self, project_id):
        try:
            project = Project.objects.get(id=project_id)
            return project.status
        except Project.DoesNotExist:
            return None
