from app.models import Torque

class Config:
    DB = 'msdserver.sqlite'
    QUEUE = Torque
    USER = 'msdserver'
