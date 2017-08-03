#!/usr/bin/env python3
# coding=utf-8

import os, sys

cmd = sys.argv[1]

if cmd == 'init':
    from app import db

    db.create_all()
elif cmd == 'run':
    from app import app

    app.run(host='0.0.0.0', port=5050, debug=True, threaded=True)
