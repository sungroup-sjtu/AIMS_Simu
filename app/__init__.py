from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from config import Config
from mstools.jobmanager import *
from mstools.simulation import *

app = Flask(__name__)
app.config.from_object(Config)
app.jinja_env.auto_reload = True

db = SQLAlchemy(app)

if Config.JOB_MANAGER == 'local':
    jobmanager = Local()
elif Config.JOB_MANAGER == 'torque':
    jobmanager = Torque(queue=Config.JOB_QUEUE, nprocs=Config.NPROC_PER_JOB)
else:
    raise Exception('Job manager not supported')

if Config.SIMULATION_ENGINE == 'gmx':
    simulation = GmxSimulation(packmol_bin=Config.PACKMOL_BIN, dff_root=Config.DFF_ROOT,
                               gmx_bin=Config.GMX_BIN, jobmanager=jobmanager)
elif Config.SIMULATION_ENGINE == 'lammps':
    simulation = LammpsSimulation(packmol_bin=Config.PACKMOL_BIN, dff_root=Config.DFF_ROOT,
                                  lmp_bin=Config.LMP_BIN, jobmanager=jobmanager)
else:
    raise Exception('Simulation engine not supported')

from .main import main as main_blueprint
from .api import api as api_blueprint

app.register_blueprint(main_blueprint)
app.register_blueprint(api_blueprint, url_prefix='/api')
