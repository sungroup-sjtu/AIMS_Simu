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
            log.info('Build task %s' % task)
            try:
                task.build()
            except Exception as e:
                log.error('Build task failed %s: %s' % (task, repr(e)))

        if task.stage == Compute.Stage.BUILDING and task.status == Compute.Status.DONE:
            log.info('Run task %s' % task)
            try:
                task.run()
            except Exception as e:
                log.error('Run task failed %s: %s' % (task, repr(e)))

        elif task.stage == Compute.Stage.RUNNING and task.status == Compute.Status.STARTED:
            log.info('Check task status %s' % task)
            try:
                task.check_finished()
            except Exception as e:
                log.error('Check task status failed %s: %s' % (task, repr(e)))

        if task.ready_to_extend:
            log.info('Extend task %s' % task)
            try:
                task.extend()
            except Exception as e:
                log.error('Extend task failed %s: %s' % (task, repr(e)))


if __name__ == '__main__':
    while True:
        with warnings.catch_warnings(record=True) as ws:
            process_tasks()
            for w in ws:
                log.warning(w.message)
        print('Sleep 300 seconds ...')
        time.sleep(300)
