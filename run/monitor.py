#!/usr/bin/env python3
# coding=utf-8

import sys, time

sys.path.append('..')
from app.models import Task, Compute, PbsJob
from app import log


def process_pbs_jobs():
    for pbs_job in PbsJob.query.filter(PbsJob.submitted == False):
        pbs_job.submit()


def process_tasks():
    n_build_tasks = 0
    for task in Task.query.all():
        if task.stage == Compute.Stage.SUBMITTED and task.status == Compute.Status.DONE:
            # do not build too much tasks once
            if n_build_tasks <= 25:
                task.build()
                n_build_tasks += 1

        if task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
            task.run()


        elif task.stage == Compute.Stage.RUNNING and task.status == Compute.Status.STARTED:
            task.check_finished()

        if task.ready_to_extend:
            task.extend()


if __name__ == '__main__':
    while True:
        process_pbs_jobs()
        process_tasks()
        log.info('Sleep 300 seconds ...')
        time.sleep(300)
