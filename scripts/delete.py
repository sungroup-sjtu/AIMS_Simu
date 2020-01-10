#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models import *

import argparse
parser = argparse.ArgumentParser(description='This is a code to delete all the information related to a specific compute')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
parser.add_argument('-i', '--id', type=int, help='The specific compute id to delete', default=0)
parser.add_argument('-ti', '--taskid', type=int, help='The specific task id to delete', default=0)
parser.add_argument('--unstarted', type=bool, help='delete all unstarted tasks', default=False)
parser.add_argument('--unfinished', type=bool, help='delete all unfinished tasks', default=False)
opt = parser.parse_args()

procedure = opt.procedure
app = create_app(procedure)
app.app_context().push()

if opt.id != 0:
    compute = Compute.query.get(opt.id)
    tasks = compute.tasks
    for task in tasks:
        if task.procedure == procedure:
            task.remove()
        else:
            raise Exception('Compute %i is not a %s (procedure) compute' % (opt.id, procedure))
    db.session.delete(compute)
elif opt.taskid != 0:
    task = Task.query.filter(Task.id == opt.idtask).first()
    task.remove()
elif opt.unstarted:
    tasks = Task.query.filter(Task.procedure == procedure)
    for task in tasks:
        if task.stage != Compute.Stage.RUNNING:
            task.remove()
elif opt.unfinished:
    tasks = Task.query.filter(Task.procedure == procedure)
    for task in tasks:
        if task.stage != Compute.Stage.RUNNING or task.status not in [Compute.Status.DONE, Compute.Status.ANALYZED]:
            task.remove()
db.session.commit()
