#!/usr/bin/env python3

import sys

sys.path.append('..')
from app.models import *
from sqlalchemy import func

import pandas as pd
import tarfile
import pathlib


def dump_tasks(task_ids):
    tasks = Task.query.filter(Task.id.in_(task_ids))
    df = pd.read_sql(tasks.statement, tasks.session.bind)
    df.to_csv('tasks.csv', index=False)


def dump_jobs(task_ids):
    jobs = Job.query.filter(Job.task_id.in_(task_ids))
    df = pd.read_sql(jobs.statement, jobs.session.bind)
    df.to_csv('jobs.csv', index=False)


def tar_tasks(task_ids):
    WORK_DIR = pathlib.Path(Config.WORK_DIR).absolute()
    tasks = Task.query.filter(Task.id.in_(task_ids))
    for task in tasks:
        print('Tar task: %s' % task)
        tar = tarfile.open('archieve/%s.tar' % task.name, 'w')
        abs_dir = pathlib.Path(task.dir).absolute()
        rel_dir = abs_dir.relative_to(WORK_DIR)
        tar.add(abs_dir, rel_dir)
        tar.close()


def create_compute(remark) -> int:
    compute = Compute()
    compute.web_id = 0
    compute.web_user_id = 0
    compute.web_ip = '127.0.0.1'
    compute.json_detail = json.dumps(None)
    compute.remark = remark

    db.session.add(compute)
    db.session.commit()
    return compute.id


def create_pbs() -> int:
    pbs = PbsJob()
    pbs.name = 'import_' + random_string(4)
    pbs.sh_file = pbs.name
    pbs.submitted = True

    db.session.add(pbs)
    db.session.commit()
    return pbs.id


def load_tasks(compute_id, task_offset):
    df = pd.read_csv('tasks.csv')
    df.compute_id = compute_id
    df.id += task_offset

    df.to_sql('task', db.session)


def load_jobs(pbs_job_id, task_offset, job_offset):
    df = pd.read_csv('jobs.csv')
    df.pbs_job_id = pbs_job_id
    df.task_id += task_offset
    df.id += job_offset

    df.to_sql('job', db.session)


def untar_tasks():
    files = os.listdir('archieve')
    for i, f in enumerate(files):
        cmd = 'tar xvf archieve/%s -C %s' % (f, Config.WORK_DIR)
        print('%i / %i' % (i + 1, len(files)), cmd)
        os.system(cmd)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'dump':
        compute_ids = map(int, sys.argv[2:])
        tasks = Task.query.filter(Task.compute_id.in_(compute_ids))
        task_ids = [task.id for task in tasks]
        dump_tasks(task_ids)
        dump_jobs(task_ids)
        # tar_tasks(task_ids)

        task_offset = db.session.query(func.max(Task.id)).scalar()
        job_offset = db.session.query(func.max(Job.id)).scalar()
        print(task_offset, job_offset)

    elif cmd == 'load':
        remark = sys.argv[2]
        compute_id = create_compute(remark)
        pbs_id = create_pbs()
        task_offset = db.session.query(func.max(Task.id)).scalar()
        job_offset = db.session.query(func.max(Job.id)).scalar()
        load_tasks(compute_id, task_offset)
        load_jobs(pbs_id, task_offset, job_offset)
        # untar_tasks()

    else:
        raise Exception('ERROR: Cmd unknown')
