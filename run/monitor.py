#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')

from app.models import *


def monitor_tasks():
    # TODO optimize logic
    tasks = Task.query.all()
    for task in tasks:
        if task.stage == Compute.Stage.SUBMITTED and task.status == Compute.Status.DONE:
            try:
                task.build()
            except Exception as e:
                raise
                print('Error when build task %s: %s' % (repr(task), str(e)))
            try:
                task.run()
            except Exception as e:
                raise
                print('Error when run task %s: %s' % (repr(task), str(e)))

        elif task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
            try:
                task.run()
            except Exception as e:
                raise
                print('Error when run task %s: %s' % (repr(task), str(e)))

        elif task.stage == Compute.Stage.RUNNING and task.status == Compute.Status.STARTED:
            task.check_finished()

        for job in Job.query:
            if job.status == Compute.Status.DONE and not job.converged:
                job.analyze()

            elif job.status == Compute.Status.ANALYZED and job.converged == False:
                job.extend()


if __name__ == '__main__':
    while True:
        monitor_tasks()
        print('Sleeping......')
        time.sleep(300)
