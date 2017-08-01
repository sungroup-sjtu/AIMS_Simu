from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from config import Config
from mstools.jobmanager import *

app = Flask(__name__)
app.config.from_object(Config)
app.jinja_env.auto_reload = True

db = SQLAlchemy(app)

if Config.JOB_MANAGER == 'local':
    jobmanager = Local(nprocs=Config.NPROC_PER_JOB, env_cmd=Config.ENV_CMD)
elif Config.JOB_MANAGER == 'torque':
    jobmanager = Torque(queue_dict=Config.QUEUE_DICT, env_cmd=Config.ENV_CMD)
elif Config.JOB_MANAGER == 'slurm':
    jobmanager = Slurm(queue_dict=Config.QUEUE_DICT, env_cmd=Config.ENV_CMD)
else:
    raise Exception('Job manager not supported')

from .main import main as main_blueprint
from .api import api as api_blueprint

app.register_blueprint(main_blueprint)
app.register_blueprint(api_blueprint, url_prefix='/api')
