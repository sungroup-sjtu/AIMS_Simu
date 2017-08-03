#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
import warnings
import logging
from app.models import *

log = logging.getLogger()
log.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s    %(levelname)-10s %(message)s')
fh = logging.FileHandler('log.txt')
ch = logging.StreamHandler()
fh.setFormatter(formatter)
ch.setFormatter(formatter)
log.addHandler(fh)
log.addHandler(ch)


def process_tasks():
    for task in Task.query.all():
        if task.stage == Compute.Stage.SUBMITTED and task.status == Compute.Status.DONE:
            try:
                log.info('Build task %s' % task)
                task.build()
            except Exception as e:
                log.error('Build task failed %s: %s' % (task, repr(e)))
            else:
                try:
                    log.info('Run task %s' % task)
                    task.run()
                except Exception as e:
                    log.error('Run task failed %s: %s' % (task, repr(e)))

        elif task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
            try:
                log.info('Run task %s' % task)
                task.run()
            except Exception as e:
                log.error('Run task failed %s: %s' % (task, repr(e)))

        elif task.stage == Compute.Stage.RUNNING and task.status == Compute.Status.STARTED:
            try:
                log.info('Check task status %s' % task)
                task.check_finished()
            except Exception as e:
                log.error('Check task status failed %s: %s' % (task, repr(e)))


def process_jobs():
    for job in Job.query.all():
        if job.status == Compute.Status.DONE and not job.converged:
            try:
                log.info('Analyze job %s' % job)
                job.analyze()
            except Exception as e:
                log.error('Analyze job failed %s: %s' % (job, repr(e)))

        elif job.status == Compute.Status.ANALYZED and job.converged == False:
            try:
                log.info('Extend job %s' % job)
                job.extend()
            except Exception as e:
                log.error('Extend job failed %s: %s' % (job, repr(e)))


if __name__ == '__main__':
    while True:
        with warnings.catch_warnings(record=True) as ws:
            process_tasks()
            process_jobs()
            for w in ws:
                log.warning(w.message)
        print('Sleep 300 seconds ...')
        time.sleep(300)