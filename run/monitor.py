#!/usr/bin/env python3
# coding=utf-8

import sys, time

sys.path.append('..')
from app.models import Task, Compute
from app import log


def process_tasks():
    for task in Task.query.all():
        if task.stage == Compute.Stage.SUBMITTED and task.status == Compute.Status.DONE:
            task.build()

        if task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
            task.run()

        elif task.stage == Compute.Stage.RUNNING and task.status == Compute.Status.STARTED:
            task.check_finished()

        if task.ready_to_extend:
            task.extend()


if __name__ == '__main__':
    while True:
        process_tasks()
        log.info('Sleep 300 seconds ...')
        time.sleep(300)
