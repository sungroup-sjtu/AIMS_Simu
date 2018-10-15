#!/usr/bin/env python3
# coding=utf-8

import sys
from app import create_app

app = create_app(sys.argv[1])

app.run(host='0.0.0.0', port=5050, debug=True, threaded=True)
