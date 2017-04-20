#!/usr/bin/env python3

from app.models import *


def monitor_tasks():
    tasks = Task.query.all()
    for task in tasks:
        if task.stage == Compute.Stage.SUBMITTED and task.status == Compute.Status.DONE:
            try:
                task.build()
            except Exception as e:
                print('Error when build task %s: %s' % (repr(task), str(e)))

        elif task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
            try:
                task.run()
            except Exception as e:
                print('Error when run task %s: %s' % (repr(task), str(e)))

        elif task.stage == Compute.Stage.RUNNING and task.status == Compute.Status.STARTED:
            for job in task.jobs:
                if job.status == Compute.Status.STARTED:
                    job.check_finished()

                elif job.status == Compute.Status.DONE and job.converged == None:
                    job.analyze()

                elif job.status == Compute.Status.DONE and job.converged == False:
                    if not job.next_cycle_started:
                        job.start_next_cycle()


if __name__ == '__main__':
    while True:
        monitor_tasks()
        time.sleep(300)
