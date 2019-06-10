#!/usr/bin/env python3
# coding=utf-8

import os, sys, time
from sqlalchemy import func

sys.path.append('..')
from app import create_app
from app.models import Task, Compute, PbsJob

app = create_app(sys.argv[1])
app.app_context().push()

CWD = os.getcwd()


def detect_exit():
    os.chdir(CWD)
    if os.path.exists('EXIT-' + sys.argv[1]):
        print('EXIT file detected')
        os.remove('EXIT-' + sys.argv[1])
        sys.exit()


def process_pbs_job(n_pbs=20):
    for pbs_job in PbsJob.query.filter(PbsJob.submitted == False).limit(n_pbs):
        detect_exit()
        pbs_job.submit()


def process_task_build(n_task=20, random=False):
    tasks = Task.query.filter(Task.stage == Compute.Stage.SUBMITTED).filter(Task.status == Compute.Status.DONE)
    if random:
        tasks = tasks.order_by(func.random())
    for task in tasks.limit(n_task):
        detect_exit()
        task.build()
        if task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
            task.run()


def process_task_run(n_task=20, random=False):
    tasks = Task.query.filter(Task.stage == Compute.Stage.BUILDING).filter(Task.status == Compute.Status.DONE)
    if random:
        tasks = tasks.order_by(func.random())
    for task in tasks.limit(n_task):
        detect_exit()
        if task.run() == -1:
            break


def process_task_check(n_task=20):
    app.jobmanager.update_stored_jobs()
    app.jm_extend.update_stored_jobs()
    tasks = Task.query.filter(Task.stage == Compute.Stage.RUNNING).filter(Task.status == Compute.Status.STARTED)
    for task in tasks.limit(n_task):
        detect_exit()
        task.check_finished_multiprocessing()


def process_task_extend():
    tasks = Task.query.filter(Task.stage == Compute.Stage.RUNNING).filter(Task.status == Compute.Status.STARTED)
    for task in tasks:
        detect_exit()
        if task.extend() == -1:
            break


if __name__ == '__main__':
    while True:
        process_pbs_job(n_pbs=50)
        process_task_run(n_task=20, random=True)
        process_task_build(n_task=20, random=True)
        process_task_check(n_task=500)
        process_task_extend()

        app.logger.info('Sleep 1800 seconds ...')
        time.sleep(1800)
