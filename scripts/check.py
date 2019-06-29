#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models import *

import argparse
parser = argparse.ArgumentParser(description='This is a code to reanalyze a specific compute')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
opt = parser.parse_args()

procedure = opt.procedure
app = create_app(procedure)
app.app_context().push()

jobs = Job.query.filter()
print(jobs.count())
for job in jobs:
    t1 = job.t
    t2 = int(job.name.split('-')[-2])
    p1 = job.p
    p2 = int(job.name.split('-')[-1])
    if t1 != t2:
        print(t1, t2 ,job.task.t_min, job.task.t_max, p1, p2)

'''
if opt.id != 0:
    compute = Compute.query.get(opt.id)
    tasks = compute.tasks
    print('Total task = %i' % (tasks.count()))
    count = 1
    for task in tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.status = Compute.Status.STARTED
            for job in task.jobs:
                job.converged = False
                job.status = Compute.Status.STARTED
                job.result = None
        else:
            raise Exception('Compute %i is not a %s (procedure) compute' % (opt.id, procedure))
    db.session.commit()
elif opt.taskid != 0:
    tasks = Task.query.filter(Task.procedure == procedure).filter(Task.id == opt.taskid)
    for task in tasks:
        if task.procedure == procedure:
            task.status = Compute.Status.STARTED
            for job in task.jobs:
                job.converged = False
                job.status = Compute.Status.STARTED
                job.result = None
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (opt.taskid, procedure))
    db.session.commit()
elif opt.failed_jobs:
    failed_jobs = Job.query.filter(Job.status == Compute.Status.FAILED)
    print('Total job = %i' % (failed_jobs.count()))
    count = 1
    for job in failed_jobs:
        sys.stdout.write('\rpresent job = %i' % (count))
        count += 1
        if job.task.procedure == procedure:
            job.converged = False
            job.status = Compute.Status.STARTED
            job.result = None
            job.task.status = Compute.Status.STARTED
    failed_tasks = Task.query.filter(Task.status == Compute.Status.FAILED).filter(Task.procedure==procedure)
    print('\nTotal task = %i' % (failed_tasks.count()))
    count = 1
    for task in failed_tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.status = Compute.Status.STARTED
            for job in task.jobs:
                job.converged = False
                job.status = Compute.Status.STARTED
                job.result = None
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (task.id, procedure))
    db.session.commit()
else:
    tasks = Task.query.filter(Task.procedure==procedure)
    print('Total task = %i' % (tasks.count()))
    count = 1
    for task in tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.status = Compute.Status.STARTED
            for job in task.jobs:
                job.converged = False
                job.status = Compute.Status.STARTED
                job.result = None
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (task.id, procedure))
    db.session.commit()
    '''


