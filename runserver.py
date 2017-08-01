#!/usr/bin/env python3
# coding=utf-8

from app import app

app.run(host='0.0.0.0', port=5050, debug=True, threaded=True)

