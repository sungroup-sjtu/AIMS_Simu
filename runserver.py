#!/usr/bin/env python3
# coding=utf-8

import os, sys
from config import Config
from app import app
from app import db

if not os.path.exists(Config.DB_PATH):
    db.create_all()

app.run(host='0.0.0.0', port=5050, debug=True, threaded=True)
