#!/usr/bin/env python3

'''
./transfer.py dump config compute_id [compute_id ...]
./transfer.py load config remark
'''

import sys
import pathlib
import tarfile
import pandas as pd

sys.path.append('..')
from app import create_app
from app.models import *
from sqlalchemy import func

app = create_app(sys.argv[1])
app.app_context().push()


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
    df.procedure += '-import'
    df.compute_id = compute_id
    df.id += (task_offset + 1 - df.id.min())

    df.to_sql('task', db.engine, if_exists='append', index=False)
    df.to_csv('tasks-load.csv', index=False)

    print(df.shape[0], 'tasks imported')


def load_jobs(pbs_job_id, task_offset, job_offset):
    df = pd.read_csv('jobs.csv')
    df.pbs_job_id = pbs_job_id
    df.task_id += (task_offset + 1 - df.task_id.min())
    df.id += (job_offset - df.id.min() + 1)

    df.to_sql('job', db.engine, if_exists='append', index=False)
    df.to_csv('jobs-load.csv', index=False)

    print(df.shape[0], 'jobs imported')


def untar_tasks():
    files = os.listdir('archieve')
    for i, f in enumerate(files):
        cmd = 'tar xvf archieve/%s -C %s' % (f, Config.WORK_DIR)
        print('%i / %i' % (i + 1, len(files)), cmd)
        os.system(cmd)


if __name__ == '__main__':
    cmd = sys.argv[2]
    if cmd == 'dump':
        compute_ids = map(int, sys.argv[3:])
        tasks = Task.query.filter(Task.compute_id.in_(compute_ids))
        task_ids = [task.id for task in tasks]
        dump_tasks(task_ids)
        dump_jobs(task_ids)
        # tar_tasks(task_ids)

        task_offset = db.session.query(func.max(Task.id)).scalar()
        job_offset = db.session.query(func.max(Job.id)).scalar()
        print(task_offset, job_offset)

    elif cmd == 'load':
        remark = sys.argv[3]
        compute_id = create_compute(remark)
        pbs_id = create_pbs()
        task_offset = db.session.query(func.max(Task.id)).scalar() or 0
        job_offset = db.session.query(func.max(Job.id)).scalar() or 0
        load_tasks(compute_id, task_offset)
        load_jobs(pbs_id, task_offset, job_offset)
        # untar_tasks()

    else:
        raise Exception('ERROR: Cmd unknown')
