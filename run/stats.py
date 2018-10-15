#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')

from app import create_app
from app.models import *

app = create_app(sys.argv[1])
app.app_context().push()

tasks = Task.query

n_tasks___ = tasks.count()
n_to_build = tasks.filter(Task.stage == 0).filter(Task.status == 9).count()
n_to_run__ = tasks.filter(Task.stage == 1).filter(Task.status == 9).count()
n_running_ = tasks.filter(Task.stage == 2).filter(Task.status == 1).count()
n_done____ = tasks.filter(Task.stage == 2).filter(Task.status == 9).count()
n_analyzed = tasks.filter(Task.stage == 2).filter(Task.status == 10).count()
n_failed__ = tasks.filter(Task.status == -1).count()

print('%6i  tasks' % n_tasks___)
print('%6i  to build' % n_to_build)
print('%6i  to run' % n_to_run__)
print('%6i  running' % n_running_)
print('%6i  done' % n_done____)
print('%6i  analyzed' % n_analyzed)
print('%6i  failed' % n_failed__)
