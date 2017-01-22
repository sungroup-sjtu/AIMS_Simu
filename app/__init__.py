from flask import Flask
from config import Config

app=Flask(__name__)
app.config.from_object(Config)
app.jinja_env.auto_reload = True

from .main import main as main_blueprint
from .api import api as api_blueprint

app.register_blueprint(main_blueprint)
app.register_blueprint(api_blueprint)
