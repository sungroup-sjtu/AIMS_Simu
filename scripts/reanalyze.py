#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models import *

import argparse
parser = argparse.ArgumentParser(description='This is a code to reanalyze a specific compute')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
parser.add_argument('-i', '--id', type=int, help='The specific compute id to reanalyze', default=0)
parser.add_argument('-ti', '--taskid', type=int, help='The specific task id to reanalyze', default=0)
parser.add_argument('-ji', '--jobid', type=int, help='The specific job id to reanalyze', default=0)
parser.add_argument('--failed_jobs', type=bool, help='Reanalyze all failed jobs', default=False)
parser.add_argument('--failed_tasks', type=bool, help='Reanalyze all failed tasks', default=False)
parser.add_argument('--unfinished_tasks', type=bool, help='Reanalyze all unfinished tasks', default=False)
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
            task.stage = Compute.Stage.RUNNING
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
            task.stage = Compute.Stage.RUNNING
            task.status = Compute.Status.STARTED
            for job in task.jobs:
                job.converged = False
                job.status = Compute.Status.STARTED
                job.result = None
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (opt.taskid, procedure))
    db.session.commit()
elif opt.jobid != 0:
    job = Job.query.filter(Job.id == opt.jobid).first()
    if job.task.procedure == procedure:
        job.converged = False
        job.status = Compute.Status.STARTED
        job.result = None
        job.task.stage = Compute.Stage.RUNNING
        job.task.status = Compute.Status.STARTED
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
            job.task.stage = Compute.Stage.RUNNING
            job.task.status = Compute.Status.STARTED
    db.session.commit()
elif opt.failed_tasks:
    failed_tasks = Task.query.filter(Task.status == Compute.Status.FAILED).filter(Task.procedure == procedure)
    print('\nTotal task = %i' % (failed_tasks.count()))
    count = 1
    for task in failed_tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.stage = Compute.Stage.RUNNING
            task.status = Compute.Status.STARTED
            for job in task.jobs:
                job.converged = False
                job.status = Compute.Status.STARTED
                job.result = None
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (task.id, procedure))
    db.session.commit()
elif opt.unfinished_tasks:
    unfinished_tasks = Task.query.filter(Task.procedure == procedure).filter(Task.stage == Compute.Stage.RUNNING).filter(Task.status == Compute.Status.STARTED)
    print('\nTotal task = %i' % (unfinished_tasks.count()))
    count = 1
    for task in unfinished_tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.stage = Compute.Stage.RUNNING
            task.status = Compute.Status.STARTED
            for job in task.jobs:
                job.converged = False
                job.status = Compute.Status.STARTED
                job.result = None
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (task.id, procedure))
    db.session.commit()
else:
    tasks = Task.query.filter(Task.procedure == procedure)
    print('Total task = %i' % (tasks.count()))
    count = 1
    for task in tasks:
        sys.stdout.write('\rpresent task = %i' % (count))
        count += 1
        if task.procedure == procedure:
            task.stage = Compute.Stage.RUNNING
            task.status = Compute.Status.STARTED
            for job in task.jobs:
                job.converged = False
                job.status = Compute.Status.STARTED
                job.result = None
        else:
            raise Exception('Task %i is not a %s (procedure) compute' % (task.id, procedure))
    db.session.commit()


