#!/usr/bin/env python3
# coding=utf-8

import os
from app import app, db, Config

if not os.path.exists(Config.DB):
    db.create_all()

app.run(host='0.0.0.0', port=5050, debug=True, threaded=True)
