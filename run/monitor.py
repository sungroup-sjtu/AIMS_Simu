#!/usr/bin/env python3
# coding=utf-8

import sys, time

sys.path.append('..')
from app.models import Task, Compute, PbsJob
from app import log


def process_pbs_job():
    for pbs_job in PbsJob.query.filter(PbsJob.submitted == False):
        pbs_job.submit()


def process_task_build(procedure=None, n_task=20):
    tasks = Task.query.filter(Task.stage == Compute.Stage.SUBMITTED).filter(Task.status == Compute.Status.DONE)
    if procedure != None:
        tasks = tasks.filter(Task.procedure == procedure)
    for task in tasks.limit(n_task):
        task.build()
        # task.run()


def process_task_run(procedure=None, n_task=20):
    tasks = Task.query.filter(Task.stage == Compute.Stage.BUILDING).filter(Task.status == Compute.Status.DONE)
    if procedure != None:
        tasks = tasks.filter(Task.procedure == procedure)
    for task in tasks.limit(n_task):
        if task.run() == -1:
            break


def process_task_check(procedure=None, n_task=20):
    tasks = Task.query.filter(Task.stage == Compute.Stage.RUNNING).filter(Task.status == Compute.Status.STARTED)
    if procedure != None:
        tasks = tasks.filter(Task.procedure == procedure)
    for task in tasks.limit(n_task):
        # task.check_finished()
        task.check_finished_multiprocessing()


def process_task_extend(procedure=None):
    tasks = Task.query
    if procedure != None:
        tasks = tasks.filter(Task.procedure == procedure)
    for task in tasks:
        if task.ready_to_extend:
            task.extend()


if __name__ == '__main__':
    while True:
        process_pbs_job()
        process_task_extend(procedure='nvt-slab')
        process_task_run(procedure='nvt-slab', n_task=20)
        process_task_build(procedure='nvt-slab', n_task=20)
        process_task_check(procedure='nvt-slab', n_task=1000)
        log.info('Sleep 1800 seconds ...')
        time.sleep(1800)
