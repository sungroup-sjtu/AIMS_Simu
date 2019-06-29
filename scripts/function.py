#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models import *

def version_gz2xy():
    procedure = 'npt'
    app = create_app(procedure)
    app.app_context().push()

    jobs = Job.query
    for job in jobs:
        job.pbs_jobs_id = json.dumps([job.pbs_job_id])
        job.repeat_id = 1
    db.session.commit()

    # npt
    tasks = Task.query.filter(Task.procedure == procedure)
    for task in tasks:
        t_list = []
        p_list = []
        for job in task.jobs:
            if job.t not in t_list:
                t_list.append(job.t)
            if job.p not in p_list:
                p_list.append(job.p)
        task.n_mol_ratio = json.dumps([1])
        if len(t_list) * len(p_list) == task.jobs.count():
            task.t_list = json.dumps(t_list)
            task.p_list = json.dumps(p_list)
        else:
            print(task.id)
        db.session.commit()
    #nvt-slab
    tasks = Task.query.filter(Task.procedure == procedure)
    for task in tasks:
        t_list = []
        for job in task.jobs:
            if job.t not in t_list:
                t_list.append(job.t)
        task.n_mol_ratio = json.dumps([1])
        task.t_list = json.dumps(t_list)
        db.session.commit()

procedure = 'npt'
app = create_app(procedure)
app.app_context().push()

tasks = Task.query.filter(Task.procedure==procedure)
for task in tasks:
    if task.jobs.count() != 56:
        print(task.id, task.jobs.count())