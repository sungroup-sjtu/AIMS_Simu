#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models import *


import argparse
parser = argparse.ArgumentParser(description='This is a code to reset tasks')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
parser.add_argument('-i', '--id', type=int, help='The specific compute id to reset', default=0)
parser.add_argument('-ti', '--taskid', type=int, help='The specific task id to reset', default=0)
parser.add_argument('--all', type=bool, help='Reset all task with specific procedure', default=False)
parser.add_argument('--failed_tasks', type=bool, help='Reset all failed tasks', default=False)
opt = parser.parse_args()

procedure = opt.procedure
app = create_app(procedure)
app.app_context().push()

if opt.id != 0:
    compute = Compute.query.get(opt.id)
    tasks = compute.tasks
    print('Total task = %i' % (tasks.count()))
    count = 1
    for task in tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.reset()
        else:
            raise Exception('Compute %i is not a %s (procedure) compute' % (opt.id, procedure))
    db.session.commit()
if opt.taskid != 0:
    tasks = Task.query.filter(Task.procedure == procedure).filter(Task.id == opt.taskid)
    for task in tasks:
        if task.procedure == procedure:
            task.reset()
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (opt.taskid, procedure))
    db.session.commit()
if opt.failed_tasks:
    failed_tasks = Task.query.filter(Task.status == Compute.Status.FAILED).filter(Task.procedure==procedure)
    print('\nTotal task = %i' % (failed_tasks.count()))
    count = 1
    for task in failed_tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.reset()
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (task.id, procedure))
    db.session.commit()
if opt.all:
    tasks = Task.query.filter(Task.procedure==procedure)
    print('Total task = %i' % (tasks.count()))
    count = 1
    for task in tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.reset()
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (task.id, procedure))
    db.session.commit()
