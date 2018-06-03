#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')

from app.models import *


tasks = Task.query
n_to_build = tasks.filter(Task.stage==0).filter(Task.status==9).count()
n_to_run   = tasks.filter(Task.stage==1).filter(Task.status==9).count()
n_running  = tasks.filter(Task.stage==2).filter(Task.status==1).count()
n_done     = tasks.filter(Task.stage==2).filter(Task.status==9).count()
print(n_to_build, 'to build\n',
      n_to_run, 'to run\n',
      n_running, 'running\n',
      n_done, 'done\n',
      )
