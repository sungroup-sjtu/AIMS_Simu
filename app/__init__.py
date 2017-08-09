import logging
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from config import Config

app = Flask(__name__)
app.config.from_object(Config)
app.jinja_env.auto_reload = True

db = SQLAlchemy(app)

log = logging.getLogger('app')
log.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s    %(levelname)-10s %(message)s')
fh = logging.FileHandler(Config.LOG)
ch = logging.StreamHandler()
fh.setFormatter(formatter)
ch.setFormatter(formatter)
log.addHandler(fh)
log.addHandler(ch)

from .main import main as main_blueprint
from .api import api as api_blueprint

app.register_blueprint(main_blueprint)
app.register_blueprint(api_blueprint, url_prefix='/api')
