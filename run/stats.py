#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')

from app import create_app
from app.models import *

procedure = sys.argv[1]
app = create_app(procedure)
app.app_context().push()

tasks = Task.query

n_tasks___ = tasks.filter(Task.procedure == sys.argv[1]).count()
n_to_build = tasks.filter(Task.procedure == sys.argv[1]).filter(Task.stage == 0).filter(Task.status == 9).count()
n_to_run__ = tasks.filter(Task.procedure == sys.argv[1]).filter(Task.stage == 1).filter(Task.status == 9).count()
n_running_ = tasks.filter(Task.procedure == sys.argv[1]).filter(Task.stage == 2).filter(Task.status == 1).count()
n_done____ = tasks.filter(Task.procedure == sys.argv[1]).filter(Task.stage == 2).filter(Task.status == 9).count()
n_analyzed = tasks.filter(Task.procedure == sys.argv[1]).filter(Task.stage == 2).filter(Task.status == 10).count()
n_failed__ = tasks.filter(Task.procedure == sys.argv[1]).filter(Task.status == -1).count()

def get_extend_info():
    jobs_extend_tmp = Job.query.filter(Job.status == Compute.Status.ANALYZED).filter(Job.converged == False)\
        .filter(Job.cycle < Config.EXTEND_CYCLE_LIMIT)
    jobs_extend = []
    for job in jobs_extend_tmp:
        if job.task.procedure == procedure and job.result is not None:
            jobs_extend.append(job)
    if len(jobs_extend) == 0:
        return
    name_list = json.loads(jobs_extend[0].result).get('name')
    extend_jobs_dict = {}
    for name in name_list:
        extend_jobs_dict.update({name: {}})
    for i, job in enumerate(jobs_extend):
        sys.stdout.write('\r%i / %i' % (i, len(jobs_extend)))
        result = json.loads(job.result)
        continue_list = result.get('continue')
        continue_n = result.get('continue_n')
        for i, name in enumerate(name_list):
            if continue_list[i] == True:
                simulation_part_n = job.get_simulation_part(name)
                if extend_jobs_dict.get(name).get(str(simulation_part_n)) is None:
                    extend_jobs_dict.get(name)[str(simulation_part_n)] = []
                extend_jobs_dict.get(name).get(str(simulation_part_n)).append(ExtendJob(job, continue_n[i]))
    for name in name_list:
        for simulation_part_n in extend_jobs_dict.get(name).keys():
            print('%6i  %s-%s jobs need to extend' % (len(extend_jobs_dict.get(name).get(simulation_part_n)), name, simulation_part_n))

print('%6i  tasks' % n_tasks___)
print('%6i  to build' % n_to_build)
print('%6i  to run' % n_to_run__)
print('%6i  running' % n_running_)
print('%6i  done' % n_done____)
print('%6i  analyzed' % n_analyzed)
print('%6i  failed' % n_failed__)
get_extend_info()
