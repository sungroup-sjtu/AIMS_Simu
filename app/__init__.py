from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from config import Config

app=Flask(__name__)
app.config.from_object(Config)
app.jinja_env.auto_reload = True

db=SQLAlchemy(app)

from .main import main as main_blueprint
from .api import api as api_blueprint

app.register_blueprint(main_blueprint)
app.register_blueprint(api_blueprint, url_prefix='/api')
