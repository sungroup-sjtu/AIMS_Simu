import os
import logging
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

from config import configs

db = SQLAlchemy()


def create_app(config_name):
    app = Flask(__name__)
    conf = configs[config_name]
    app.config.from_object(conf)
    conf.init_app(app)

    db.init_app(app)

    from .main import main
    from .api import api

    app.register_blueprint(main)
    app.register_blueprint(api, url_prefix='/api')
    app.jinja_env.auto_reload = True

    ### Default app.logger works weried. Replace it with a new logger
    logger = logging.getLogger(config_name)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s    %(levelname)-10s %(message)s')
    fh = logging.FileHandler(conf.LOG)
    ch = logging.StreamHandler()
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    app.logger = logger

    if not os.path.exists(conf.DB):
        with app.app_context():
            db.create_all()

    return app
