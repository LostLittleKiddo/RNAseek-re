import pytest
import json
from channels.testing import WebsocketCommunicator
from django.test import TestCase
from rsa.consumers import StatusConsumer
from rsp.asgi import application
from rsa.models import Project, User  # Import the Project and User models


class StatusConsumerTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        # Create a User instance for the test
        test_user = User.objects.create(username="testuser")
        # Create a Project instance for the test, associated with the user
        Project.objects.create(id=1, user=test_user, name="Test Project", status="initial_status", session_id="test_session") # Added session_id

    @pytest.mark.asyncio
    @pytest.mark.django_db(transaction=True)  # Ensure DB transaction handling for async tests
    async def test_websocket_connection_and_initial_message(self):
        communicator = WebsocketCommunicator(application, "/ws/status/1/")
        connected, subprotocol = await communicator.connect()
        assert connected

        # Check for the initial status message
        response = await communicator.receive_from() # Use receive_from for text data
        expected_response = {
            "status": "initial_status",
            "project_id": "1" # project_id in the message is a string
        }
        assert json.loads(response) == expected_response

        await communicator.disconnect()

# Create your tests here.
